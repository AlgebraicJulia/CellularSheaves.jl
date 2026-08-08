############################################################################################
# Compensated (as-if-doubled-precision) matrix-vector products for BlockSparseMatrix.
#
# Include this file after utils.jl and block_sparse_matrix.jl (it uses unwrapadj,
# unwrapsym, intriangle, promote_eltype, and the wrapper type aliases), e.g. add
#
#     include("comp/comp.jl")
#
# after include("blas/blas.jl") in BlockSparseArrays.jl.
#
# Naming convention, applied uniformly:
#
#   two[verb]   returns (or maintains) a value/error PAIR.
#   ext[verb]   returns a SINGLE result, accurate to extended (doubled) precision.
#
# For any vector x, its shadow cx satisfies: x + cx is the accurate value.
#
# Design: two-phase.
#
#   phase 1:  (s, cs) ← op(A) * b     accumulated with error-free transformations
#                                     (twogemv_impl!, twosymv_impl!); the block
#                                     kernels never see α or β.
#   phase 2:  either fold:            c ← α (s + cs) + β c            (extaxpby!)
#             or keep the pair:  (c, cc) ← α (s + cs) + β c           (twoaxpby!)
#
# Accuracy guarantees, from strongest to weakest (details in the docstrings):
#
#   twosum, twomul        exact — the pair reconstructs the true result with no
#                         error at all.
#   extaxpby!, extmul!    the single result is as accurate as if the whole
#                         computation were done in doubled working precision and
#                         rounded once (Ogita–Rump–Oishi, SISC 2005).
#   twoaxpby!, twomul!    the pair carries the doubled-precision value to first
#                         order; O(u²) terms are unaccounted. The pair is
#                         normalized, and its hi lane equals the extmul! result
#                         bitwise.
#
# IMPORTANT:
#   - No @fastmath anywhere in these files, ever. Reassociation legally reduces
#     twosum's correction term to zero. Julia's default semantics are safe.
#   - Real element types only. twosum is valid componentwise for Complex, but
#     twomul is not (fma does not distribute over complex products); complex
#     support needs the complex EFTs of Graillat et al. The tA arguments are kept
#     in the signatures so that adding complex later is a local change.
############################################################################################

############################################################################################
# Error-free transformations
############################################################################################

@inline function twosum(x, y)
    σ  = x + y
    x′ = σ - y
    y′ = σ - x′
    return σ, (x - x′) + (y - y′)
end

@inline function twomul(x, y)
    π = x * y
    return π, fma(x, y, -π)
end

############################################################################################
# Compensated accumulation helpers (phase-1 machinery)
############################################################################################

# Register pair: (Δ, δ) += a * x. Returns the updated pair.
@inline function twoacc(Δ, δ, a, x)
    π, ep = twomul(a, x)
    Δ, es = twosum(Δ, π)
    return Δ, δ + (ep + es)
end

# Memory pair: (s[i], cs[i]) += a * x.
@inline function twoaxpy!(s::AbstractVector, cs::AbstractVector, i::Integer, a, x)
    @inbounds begin
        π, ep = twomul(a, x)
        σ, es = twosum(s[i], π)
        s[i]  = σ
        cs[i] += ep + es
    end
    return
end

# Fold a register pair into a memory pair: (s[j], cs[j]) += (Δ, δ).
@inline function twofold!(s::AbstractVector, cs::AbstractVector, j::Integer, Δ, δ)
    @inbounds begin
        σ, e  = twosum(s[j], Δ)
        s[j]  = σ
        cs[j] += δ + e
    end
    return
end

############################################################################################
# Epilogues (phase 2)
############################################################################################

"""
    extaxpby!(c, α, s, cs, β) -> c

Folded epilogue: `c ← α (s + cs) + β c`, committing a single working-precision
rounding per entry. Given a phase-1 pair `(s, cs)`, the result is as accurate
as if `α * op(A) * b + β * c` were computed in doubled working precision and
rounded once. All α/β fastpaths live here.
"""
function extaxpby!(c::AbstractVector, α::Number, s::AbstractVector, cs::AbstractVector, β::Number)
    if iszero(β)
        if isone(α)
            @inbounds for i in eachindex(c)
                c[i] = s[i] + cs[i]
            end
        else
            @inbounds for i in eachindex(c)
                π, ep = twomul(α, s[i])
                c[i]  = π + fma(α, cs[i], ep)
            end
        end
    elseif isone(β)
        if isone(α)
            @inbounds for i in eachindex(c)
                σ, es = twosum(c[i], s[i])
                c[i]  = σ + (cs[i] + es)
            end
        elseif isone(-α)
            @inbounds for i in eachindex(c)
                σ, es = twosum(c[i], -s[i])
                c[i]  = σ + (es - cs[i])
            end
        else
            @inbounds for i in eachindex(c)
                π, ep = twomul(α, s[i])
                σ, es = twosum(c[i], π)
                c[i]  = σ + (fma(α, cs[i], ep) + es)
            end
        end
    elseif isone(-β)
        if isone(α)
            @inbounds for i in eachindex(c)
                σ, es = twosum(-c[i], s[i])
                c[i]  = σ + (cs[i] + es)
            end
        elseif isone(-α)
            @inbounds for i in eachindex(c)
                σ, es = twosum(-c[i], -s[i])
                c[i]  = σ + (es - cs[i])
            end
        else
            @inbounds for i in eachindex(c)
                π, ep = twomul(α, s[i])
                σ, es = twosum(-c[i], π)
                c[i]  = σ + (fma(α, cs[i], ep) + es)
            end
        end
    else
        @inbounds for i in eachindex(c)
            π₁, e₁ = twomul(α, s[i])
            π₂, e₂ = twomul(β, c[i])
            σ,  e₃ = twosum(π₁, π₂)
            c[i]   = σ + (fma(α, cs[i], e₁) + e₂ + e₃)
        end
    end

    return c
end

"""
    twoaxpby!(c, cc, α, s, cs, β) -> (c, cc)

Pair epilogue: `(c, cc) ← α (s + cs) + β c`, skipping the final rounding that
[`extaxpby!`](@ref) commits — analogous to how [`twosum`](@ref) returns the
error a bare `+` would discard. Reads the incoming `c` before overwriting it.

Guarantees: `c + cc` carries the doubled-precision value to first order (O(u²)
terms are unaccounted, unlike the exact scalar EFTs). The pair is normalized
via a fused final `twosum`, so `(c, cc)` is a proper hi/lo pair and `c` equals
the output of [`extaxpby!`](@ref) (hence of [`extmul!`](@ref)) **bitwise** —
the pair form is a strict refinement of the folded form.

Consumers that merge pairs downstream and need to re-normalize can broadcast
`twosum` over the merged lanes.
"""
function twoaxpby!(
        c::AbstractVector, cc::AbstractVector,
        α::Number, s::AbstractVector, cs::AbstractVector, β::Number,
    )
    if iszero(β)
        if isone(α)
            @inbounds for i in eachindex(c)
                c[i], cc[i] = twosum(s[i], cs[i])
            end
        else
            @inbounds for i in eachindex(c)
                π, ep = twomul(α, s[i])
                c[i], cc[i] = twosum(π, fma(α, cs[i], ep))
            end
        end
    elseif isone(β)
        if isone(α)
            @inbounds for i in eachindex(c)
                σ, es = twosum(c[i], s[i])
                c[i], cc[i] = twosum(σ, cs[i] + es)
            end
        elseif isone(-α)
            @inbounds for i in eachindex(c)
                σ, es = twosum(c[i], -s[i])
                c[i], cc[i] = twosum(σ, es - cs[i])
            end
        else
            @inbounds for i in eachindex(c)
                π, ep = twomul(α, s[i])
                σ, es = twosum(c[i], π)
                c[i], cc[i] = twosum(σ, fma(α, cs[i], ep) + es)
            end
        end
    elseif isone(-β)
        if isone(α)
            @inbounds for i in eachindex(c)
                σ, es = twosum(-c[i], s[i])
                c[i], cc[i] = twosum(σ, cs[i] + es)
            end
        elseif isone(-α)
            @inbounds for i in eachindex(c)
                σ, es = twosum(-c[i], -s[i])
                c[i], cc[i] = twosum(σ, es - cs[i])
            end
        else
            @inbounds for i in eachindex(c)
                π, ep = twomul(α, s[i])
                σ, es = twosum(-c[i], π)
                c[i], cc[i] = twosum(σ, fma(α, cs[i], ep) + es)
            end
        end
    else
        @inbounds for i in eachindex(c)
            π₁, e₁ = twomul(α, s[i])
            π₂, e₂ = twomul(β, c[i])
            σ,  e₃ = twosum(π₁, π₂)
            c[i], cc[i] = twosum(σ, fma(α, cs[i], e₁) + e₂ + e₃)
        end
    end

    return c, cc
end

############################################################################################
# Kernels (phase 1)
############################################################################################

include("twogemv.jl")
include("twosymv.jl")

############################################################################################
# Entry points — mirror the plain mul! dispatch, restricted to real elements
############################################################################################

# Dispatch helper: route a (possibly wrapped) matrix to the right phase-1 driver.
# Returns the filled pair (s, cs).
function twomv_impl!(
        c::AbstractVector, A::MaybeAdjOrTransBlockSparseMatrix{T, I}, b::AbstractVector,
        s::AbstractVector, cs::AbstractVector,
    ) where {T <: Real, I}
    @assert length(c) == size(A, 1)
    @assert length(b) == size(A, 2)
    uA, tA = unwrapadj(A)
    return twogemv_impl!(tA, uA, b, s, cs)
end

function twomv_impl!(
        c::AbstractVector, A::HermOrSymBlockSparseMatrix{T, I}, b::AbstractVector,
        s::AbstractVector, cs::AbstractVector,
    ) where {T <: Real, I}
    @assert length(c) == size(A, 1)
    @assert length(b) == size(A, 2)
    uA, tA = unwrapsym(A)

    if A.uplo == 'L'
        uplo = :L
    else
        uplo = :U
    end

    return twosymv_impl!(tA, uplo, uA, b, s, cs)
end

"""
    extmul!(c, A, b, α, β[, s, cs]) -> c

Extended-precision `c = α * A * b + β * c` for real-element `BlockSparseMatrix`
`A`, optionally wrapped in `Transpose`, `Adjoint`, `Symmetric`, or `Hermitian`.

Guarantee: the result is as accurate as if the entire computation were carried
out in doubled working precision and rounded once. Workspace vectors `s`, `cs`
must have the same length as `c`; if omitted they are allocated. After the
call, `s + cs` holds the unscaled product `op(A) * b` to doubled precision
(first order).
"""
function extmul!(
        c::AbstractVector, A::Union{MaybeAdjOrTransBlockSparseMatrix, HermOrSymBlockSparseMatrix},
        b::AbstractVector, α::Number, β::Number,
        s::AbstractVector = similar(c, promote_eltype(A, b)),
        cs::AbstractVector = similar(c, promote_eltype(A, b)),
    )
    twomv_impl!(c, A, b, s, cs)
    return extaxpby!(c, α, s, cs, β)
end

"""
    twomul!(c, cc, A, b, α, β[, s, cs]) -> (c, cc)

Pair-returning variant of [`extmul!`](@ref): writes a normalized hi/lo pair
`(c, cc)` such that `c + cc == α * A * b + β * c₀` to doubled working
precision, without the final rounding.

Guarantees: first-order accurate (O(u²) terms are unaccounted — unlike the
scalar EFT [`twomul`](@ref), whose pair is exact); normalized; and `c` equals
the output of [`extmul!`](@ref) bitwise, so the pair form strictly refines the
folded form.
"""
function twomul!(
        c::AbstractVector, cc::AbstractVector,
        A::Union{MaybeAdjOrTransBlockSparseMatrix, HermOrSymBlockSparseMatrix},
        b::AbstractVector, α::Number, β::Number,
        s::AbstractVector = similar(c, promote_eltype(A, b)),
        cs::AbstractVector = similar(c, promote_eltype(A, b)),
    )
    twomv_impl!(c, A, b, s, cs)
    return twoaxpby!(c, cc, α, s, cs, β)
end
