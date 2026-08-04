# Augmentation KKT backend (rebuilt — rebuild_spec.md §1). There is no settings object: α is a
# per-call argument (a trajectory-level parameter owned by the IPM, held as solver state and mutated
# only by updateaug!), and the ρ-shift ladder bounds are baked into the workspace at construction. The
# augmented system is F = (1/α)·H + BᵀB. There is no CRAIG iteration budget of our own: CRAIG terminates
# in exact arithmetic by iteration n+m, and Krylov.jl caps at its dimension-based default when itmax is
# unset — the hang guard belongs to the library that owns the iteration, at its termination bound.
#
# A different KKT backend is a different workspace type behind the same α/atol argument interface —
# there is no settings-level polymorphism (the KKTSolver abstract type is the only extension point).

# Inexact-Uzawa interior discipline: an interior refinement pass drives its CRAIG correction only to
# θ·(its entry residual), not to atol — each pass cuts the residual ~10× and the outer refinement loop
# converges it. Base solves and the final exit run exact (θ = 0). See pricing_augmentation §1.
const INTERIOR_THETA = 0.1

struct UzawaSolver{UPLO, T, I <: Integer} <: KKTSolver{T}
    F::FChordalTriangular{:N, UPLO, T, I}
    L::BlockSparseMatrix{T, I}
    facwrk::FactorizationWorkspace{T, I}
    divwrk::DivisionWorkspace{T, I}
    itrwrk::CraigWorkspace{T}
    r::Vector{T}
    α::Scalar{T}
    sd::Vector{T}   # refinement scratch: dual residual (n)
    sp::Vector{T}   # refinement scratch: primal residual (m)
    dp::Vector{T}   # refinement scratch: primal correction (n)
    dy::Vector{T}   # refinement scratch: dual correction (m)
end

function UzawaSolver(F::FChordalTriangular{:N, UPLO, T, I}, L::BlockSparseMatrix{T, I}, B::BlockSparseMatrix{T, I}) where {UPLO, T, I <: Integer}
    m, n = size(B)
    facwrk = FactorizationWorkspace(F)
    divwrk = DivisionWorkspace(F, 1)
    itrwrk = CraigWorkspace{T}(m, n)
    r = zeros(T, m)
    α = ones(T)
    return UzawaSolver(F, L, facwrk, divwrk, itrwrk, r, α, zeros(T, n), zeros(T, m), zeros(T, n), zeros(T, m))
end

function makekkt(B::BlockSparseMatrix{T, I}; elim::EliminationAlgorithm = DEFAULT_ELIMINATION_ALGORITHM) where {T, I}
    weights, graph = weightedgraph(B)
    R, P, S = symbolic(weights, graph; alg=elim)
    B = selectvtxs(B, R.perm)
    F = FChordalTriangular{:N, :L, T, I}(S)
    L = B' * B
    wrk = UzawaSolver(F, L, B)
    return R, P, B, wrk
end

# build & factor F = (1/α)·A + BᵀB (BᵀB precomputed as wrk.L), with the ρ-shift ladder on failure.
# `rgmin` is the floor the ladder starts from (the running s.ρ[]).
function initkkt!(wrk::UzawaSolver{UPLO, T}, A::BlockSparseMatrix; α::T, rgmin::T) where {UPLO, T}
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

# Solve the 2-row KKT system to atol = max(force_tol, floor_tol): a base solve plus the solver-owned
# iterative-refinement loop, returning a typed exit (Govaerts–Pryce BE+1). This is the guarantee — the
# residual is driven to atol or the exit says why not (floor / stalled / itmax). Refinement scratch is
# workspace-resident. Returns (nbase, nrefine, npass, status, fres, entry_pres, entry_dres).
function solvekkt!(
    wrk::UzawaSolver{UPLO, T},
    Δp::AbstractVector{T},
    Δy::AbstractVector{T},
    A::BlockSparseMatrix{T},
    B::BlockSparseMatrix{T},
    f::AbstractVector{T},
    rp::AbstractVector{T},
    nc::T,
    ng::T,
    y0 = nothing;
    force_tol::T,
    floor_tol::T,
    stall::T,
    itmax::Int,
) where {UPLO, T}
    atol = max(force_tol, floor_tol)
    #
    # base solve
    #
    nbase, fres = solveuzw!(wrk.divwrk, wrk.itrwrk, Δp, Δy, wrk.r, wrk.F, A, B,
                             f, rp, wrk.α[], atol, y0)
    #
    # iterative refinement of the 2-row residual
    #
    sd = wrk.sd; sp = wrk.sp; dp = wrk.dp; dy = wrk.dy
    nrefine = 0
    npass = 0
    status = REFINE_ITMAX
    prv = typemax(T)
    pres1 = T(NaN)
    dres1 = T(NaN)

    for i in 1:itmax
        #
        #   [ sd ] = [ f  ] - [ A -Bᵀ ] [ Δp ]
        #   [ sp ]   [ rp ]   [ B  0  ] [ Δy ]
        #
        mulkkt!(sd, sp, A, B, Δp, Δy)
        axpby!(one(T), f,  -one(T), sd)
        axpby!(one(T), rp, -one(T), sp)

        dres = norm(sd, Inf) / (one(T) + nc)
        pres = norm(sp, Inf) / (one(T) + ng)
        res = max(dres, pres)

        if isone(i)
            pres1 = pres
            dres1 = dres
        end

        if res ≤ force_tol
            status = REACHED_FORCE
            break
        end
        if res ≤ floor_tol
            status = REACHED_FLOOR
            break
        end
        if res > stall * prv
            status = REFINE_STALLED
            break
        end

        prv = res
        #
        #   [ A -Bᵀ ] [ dp ] = [ sd ]
        #   [ B  0  ] [ dy ]   [ sp ]
        #
        n, _ = solveuzw!(wrk.divwrk, wrk.itrwrk, dp, dy, wrk.r, wrk.F, A, B,
                          sd, sp, wrk.α[], atol; θ = T(INTERIOR_THETA))
        nrefine += n
        npass += 1

        axpy!(one(T), dp, Δp)
        axpy!(one(T), dy, Δy)
    end

    return nbase, nrefine, npass, status, fres, pres1, dres1
end

#
# Solve the KKT system
#
#   [ A -Bᵀ ] [ x ] = [ f ]
#   [ B  0  ] [ y ]   [ g ]
#
function solveuzw!(
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
        y0 = nothing;
        θ::T = zero(T),
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
    # craig applies the preconditioner (F Fᵀ)⁻¹ = (β A + BᵀB)⁻¹ inline via F + divwrk. Absolute stop
    # only (rtol = 0): base solves (θ = 0) drive to exactly atol; interior inexact-Uzawa passes (θ > 0)
    # stop at max(atol, θ·nr0), each reducing its own entry residual by ~1/θ and no further. max(atol, 0)
    # = atol, so the one formula covers both.
    nc, _ = craig!(itrwrk, B, F, divwrk, r; btol = zero(T), atol = max(atol, θ * nr0), rtol = zero(T))
    niter += nc
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
