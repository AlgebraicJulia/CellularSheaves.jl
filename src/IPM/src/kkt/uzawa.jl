# Augmentation KKT backend (rebuilt — rebuild_spec.md §1). There is no settings object: α is a
# per-call argument (a trajectory-level parameter owned by the IPM, held as solver state and mutated
# only by updateaug!), and the ρ-shift ladder bounds are baked into the workspace at construction. The
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
end

function UzawaWorkspace(F::FChordalTriangular{:N, UPLO, T, I}, L::BlockSparseMatrix{T, I}, B::BlockSparseMatrix{T, I}) where {UPLO, T, I <: Integer}
    m, n = size(B)
    facwrk = FactorizationWorkspace(F)
    divwrk = DivisionWorkspace(F, 1)
    itrwrk = CraigWorkspace(m, n, Vector{T})
    r = zeros(T, m)
    α = ones(T)
    return UzawaWorkspace(F, L, facwrk, divwrk, itrwrk, r, α)
end

function make_kkt(B::BlockSparseMatrix{T, I}; elim::EliminationAlgorithm = DEFAULT_ELIMINATION_ALGORITHM) where {T, I}
    weights, graph = weightedgraph(B)
    R, P, S = symbolic(weights, graph; alg=elim)
    B = selectvtxs(B, R.perm)
    F = FChordalTriangular{:N, :L, T, I}(S)
    L = B' * B
    wrk = UzawaWorkspace(F, L, B)
    return R, P, B, wrk
end

# build & factor F = (1/α)·A + BᵀB (BᵀB precomputed as wrk.L), with the ρ-shift ladder on failure.
# `rgmin` is the floor the ladder starts from (the running s.ρ[]).
function initkkt!(wrk::UzawaWorkspace{UPLO, T}, A::BlockSparseMatrix; α::T, rgmin::T) where {UPLO, T}
    wrk.α[] = α
    return init_uzw!(wrk.facwrk, wrk.F, wrk.L, A, α, rgmin)
end

function init_uzw!(
        facwrk::FactorizationWorkspace{T},
        F::ChordalTriangular{:N, UPLO, T},
        L::BlockSparseMatrix{T},
        A::BlockSparseMatrix{T},
        α::T,
        rgmin::T,
    ) where {UPLO, T}
    @assert size(F, 1) == size(L, 1) == size(A, 1)

    β = inv(α)
    #
    # assemble and factor the augmented matrix
    #
    #   F = β A + Bᵀ B
    #
    copyto!(F, L)
    axpy!(β, A, F)
    info = cholesky!(facwrk, F; check=false)
    #
    # unshifted factorization succeeded — no shift applied, so the APPLIED shift is 0. Returning the
    # actual applied shift (not the floor `rgmin`) is what the history records; the caller keeps the
    # monotonic floor separately via s.ρ[] = max(s.ρ[], applied). The ladder below starts from the
    # floor `rgmin` (= s.ρ[]), so it never restarts below a shift a previous solve already needed.
    #
    iszero(info) && return true, zero(T)
    #
    # on failure, add a diagonal shift ρ I and retry, climbing the ρ-ladder from the floor `rgmin`:
    #
    #   F ← β A + Bᵀ B + ρ I,   ρ ← rgmin, 8·rgmin, 64·rgmin, …
    #
    # The ×10 rung count is a panic cap against an infinite loop, not a tuned bound: with the floor at
    # u·‖B‖² a single rung fixes the roundoff-tipped pivot in practice, and needing a shift 8^10 above
    # the floor means the matrix is genuinely singular, so failing there is correct.
    #
    ρ = rgmin
    for _ in 1:10
        copyto!(F, L)
        axpy!(β, A, F)
        axpy!(ρ, I, F)
        iszero(cholesky!(facwrk, F; check=false)) && return true, ρ
        ρ *= 8
    end

    return false, ρ
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
    # compute the residual norm
    #
    #   nr₀ = ‖ g - B x ‖
    #
    nr0 = norm(r)
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

    return niter, nr0
end
