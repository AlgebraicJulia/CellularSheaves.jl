# Augmentation KKT backend (rebuilt — rebuild_spec.md §1). There is no settings object: α is a
# per-call argument (a trajectory-level parameter owned by the IPM, held as solver state, fixed at construction), and the ρ-shift ladder bounds are baked into the workspace at construction. The
# augmented system is F = (1/α)·H + BᵀB. There is no CRAIG iteration budget of our own: CRAIG terminates
# in exact arithmetic by iteration n+m, and Krylov.jl caps at its dimension-based default when itmax is
# unset — the hang guard belongs to the library that owns the iteration, at its termination bound.
#
# A different KKT backend is a different workspace type behind the same α/atol argument interface —
# there is no settings-level polymorphism (the KKTWorkspace abstract type is the only extension point).

struct UzawaWorkspace{UPLO, T, I <: Integer} <: KKTWorkspace{T}
    F::FChordalTriangular{:N, UPLO, T, I}
    L::BlockSparseMatrix{T, I}
    facwrk::FactorizationWorkspace{T, I}
    divwrk::DivisionWorkspace{T, I}
    itrwrk::CraigWorkspace{T, T, Vector{T}}
    r::Vector{T}
    α::Scalar{T}
    rgmin::T                    # baked: ρ-shift ladder lower bound
    rgmax::T                    # baked: ρ-shift ladder upper bound
end

function UzawaWorkspace(F::FChordalTriangular{:N, UPLO, T, I}, L::BlockSparseMatrix{T, I}, B::BlockSparseMatrix{T, I};
                        rgmin::T, rgmax::T) where {UPLO, T, I <: Integer}
    m, n = size(B)
    facwrk = FactorizationWorkspace(F)
    divwrk = DivisionWorkspace(F, 1)
    itrwrk = CraigWorkspace(m, n, Vector{T})
    r = zeros(T, m)
    α = ones(T)
    return UzawaWorkspace(F, L, facwrk, divwrk, itrwrk, r, α, rgmin, rgmax)
end

function make_kkt(B::BlockSparseMatrix{T, I}; elim::EliminationAlgorithm = DEFAULT_ELIMINATION_ALGORITHM,
                  rgmin::T, rgmax::T) where {T, I}
    weights, graph = weightedgraph(B)
    R, P, S = symbolic(weights, graph; alg=elim)
    B = selectvtxs(B, R.perm)
    F = FChordalTriangular{:N, :L, T, I}(S)
    L = B' * B
    wrk = UzawaWorkspace(F, L, B; rgmin, rgmax)
    return R, P, B, wrk
end

# build & factor F = (1/α)·A + BᵀB (BᵀB precomputed as wrk.L), with the ρ-shift ladder on failure.
function initkkt!(wrk::UzawaWorkspace{UPLO, T}, A::BlockSparseMatrix; α::T) where {UPLO, T}
    wrk.α[] = α
    return init_uzw!(wrk.facwrk, wrk.F, wrk.L, A, α, wrk.rgmin, wrk.rgmax)
end

function init_uzw!(
        facwrk::FactorizationWorkspace{T},
        F::ChordalTriangular{:N, UPLO, T},
        L::BlockSparseMatrix{T},
        A::BlockSparseMatrix{T},
        α::T,
        rgmin::T,
        rgmax::T
    ) where {UPLO, T}
    @assert size(F, 1) == size(L, 1) == size(A, 1)

    β = inv(α)
    ρ = rgmin
    #
    # assemble and factor the augmented matrix
    #
    #   F = β A + Bᵀ B
    #
    copyto!(F, L)
    axpy!(β, A, F)
    info = cholesky!(facwrk, F; check=false)
    #
    # on failure, add a diagonal shift ρ I and retry:
    #
    #   F ← β A + Bᵀ B + ρ I
    #
    if !iszero(info)
        copyto!(F, L)
        axpy!(β, A, F)
        axpy!(ρ, I, F)
        info = cholesky!(facwrk, F; check=false)
    end

    while !iszero(info) && 8ρ ≤ rgmax
        ρ *= 8
        copyto!(F, L)
        axpy!(β, A, F)
        axpy!(ρ, I, F)
        info = cholesky!(facwrk, F; check=false)
    end

    return iszero(info), ρ
end

# base KKT solve to `atol`; returns the base CRAIG iteration count (Krylov's dimension-based cap on hang).
function solve_kkt!(
    wrk::UzawaWorkspace{UPLO, T},
    x::AbstractVector{T},
    y::AbstractVector{T},
    A::BlockSparseMatrix{T},
    B::BlockSparseMatrix{T},
    f::AbstractVector{T},
    g::AbstractVector{T},
    y0 = nothing;
    atol::T,
) where {UPLO, T}
    return solve_uzw!(wrk.divwrk, wrk.itrwrk, x, y, wrk.r, wrk.F, A, B,
                      f, g, wrk.α[], atol, y0)
end

#
# Solve the KKT system
#
#   [ A -Bᵀ ] [ x ] = [ f ]
#   [ B  0  ] [ y ]   [ g ]
#
function solve_uzw!(
        divwrk::DivisionWorkspace{T},
        itrwrk::CraigWorkspace{T},
        x::AbstractVector{T},
        y::AbstractVector{T},
        r::AbstractVector{T},
        F::ChordalTriangular{:N, UPLO, T},
        A::BlockSparseMatrix{T},
        B::BlockSparseMatrix{T},
        f::AbstractVector{T},
        g::AbstractVector{T},
        α::T,
        atol::T,
        y0 = nothing,
    ) where {UPLO, T}
    m, n = size(B)

    @assert length(x) == n
    @assert length(y) == m
    @assert length(r) == m
    @assert length(f) == n
    @assert length(g) == m
    @assert size(F, 1) == n
    @assert size(A, 1) == n
    niter = 1; β = inv(α)
    #
    # solve for x:
    #
    #   (β A + Bᵀ B) x =  β f + Bᵀ (g + β y₀)
    #
    copyto!(r, g)

    if !isnothing(y0)
        axpy!(β, y0, r)
    end

    copyto!(x, f)
    mul!(x, B', r, one(T), β)
    ldiv!(divwrk, F,  x)
    ldiv!(divwrk, F', x)
    #
    # compute the residual
    #
    #   r = g - B x
    #
    copyto!(r, g)
    mul!(r, B, x, -one(T), one(T))
    #
    # solve for δx and δy:
    #
    #   [ β A + Bᵀ B  -Bᵀ ] [ δx ] = [ 0 ]
    #   [ B            0  ] [ δy ]   [ r ]
    #
    function prec!(u, v)
        #
        # solve for u:
        #
        #   (β A + Bᵀ B) u = v
        #
        copyto!(u, v)
        ldiv!(divwrk, F,  u)
        ldiv!(divwrk, F', u)
    end

    N = LinearOperator(T, n, n, true, true, prec!)
    craig!(itrwrk, B, r; ldiv = false, btol = zero(T), N, atol)
    niter += itrwrk.stats.niter
    #
    # update x:
    #
    #   x = x + δx
    #
    axpy!(one(T), itrwrk.x, x)
    #
    # recover y:
    #
    #   y = y₀ + 1/β (δy + (g - B x))
    #
    copyto!(r, g)
    mul!(r, B, x, -one(T), one(T))
    copyto!(y, itrwrk.y)
    axpy!(one(T), r, y)
    lmul!(α, y)

    if !isnothing(y0)
        axpy!(one(T), y0, y)
    end

    return niter
end
