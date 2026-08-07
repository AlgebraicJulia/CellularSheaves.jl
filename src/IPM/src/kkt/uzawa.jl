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
    r::FVector{T}
    α::FScalar{T}
    sd::FVector{T}   # refinement scratch: dual residual (n)
    sp::FVector{T}   # refinement scratch: primal residual (m)
    dp::FVector{T}   # refinement scratch: primal correction (n)
    dy::FVector{T}   # refinement scratch: dual correction (m)
end

function UzawaSolver(F::FChordalTriangular{:N, UPLO, T, I}, L::BlockSparseMatrix{T, I}, B::BlockSparseMatrix{T, I}) where {UPLO, T, I <: Integer}
    m, n = size(B)
    facwrk = FactorizationWorkspace(F)
    divwrk = DivisionWorkspace(F, 1)
    itrwrk = CraigWorkspace{T}(m, n)
    r = FVector{T}(undef, m)
    α = FScalar{T}(undef)
    α[] = one(T)
    return UzawaSolver(F, L, facwrk, divwrk, itrwrk, r, α,
                       FVector{T}(undef, n), FVector{T}(undef, m), FVector{T}(undef, n), FVector{T}(undef, m))
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
    return inituzw!(wrk.facwrk, wrk.F, wrk.L, A, α, rgmin)
end

function inituzw!(
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

#
# One iteration of the 2-row iterative-
# refinement loop. Forms the residual of
# (Δp, Δy)
#
#   [ sd ] = [ f  ] - [ A -Bᵀ ] [ Δp ]
#   [ sp ]   [ rp ]   [ B  0  ] [ Δy ]
#
# and, unless it already meets its per-row floor
# (ptol/ytol raised to the terms' roundoff) or
# has stalled against the previous residual
# norms (pprv, dprv), solves the correction
#
#   [ A -Bᵀ ] [ dp ] = [ sd ]
#   [ B  0  ] [ dy ]   [ sp ]
#
# to θ·(its entry residual) and folds it in.
# Returns (status, pres, dres, n): KKT_SOLVED or
# KKT_STAGNATED (no correction, n = 0), or
# KKT_CONTINUE (corrected, n = craig count).
# The stall test compares res = max(pres/pfloor,
# dres/dfloor) against its previous value rebuilt
# from (pprv, dprv).
#
function refinekkt!(
        wrk::UzawaSolver{UPLO, T},
        Δp::AbstractVector{T},
        Δy::AbstractVector{T},
        A::BlockSparseMatrix{T},
        B::BlockSparseMatrix{T},
        f::AbstractVector{T},
        rp::AbstractVector{T},
        pprv::T,
        dprv::T;
        ptol::T,
        ytol::T,
        stall::T,
    ) where {UPLO, T}
    sd = wrk.sd; sp = wrk.sp; dp = wrk.dp; dy = wrk.dy

    #
    # form the residual with the three matvecs kept SEPARATE so each row can be floored at its own
    # evaluation limit — u × the sum of the magnitudes differenced to form it. ‖A Δp‖ is H-amplified but
    # cancels against ‖Bᵀ Δy‖ in A Δp - Bᵀ Δy, so a fused ‖A Δp - Bᵀ Δy‖ floor would undercount exactly when
    # it matters (see refinebdr!). dp holds Bᵀ Δy just long enough to floor and reassemble; the correction
    # below overwrites it. The floors are the linear solve's own — the IPM passes force_tol, not a guessed
    # floor_tol.
    #
    mul!(sp, B, Δp)                  # sp = B Δp
    mul!(sd, Symmetric(A, :L), Δp)   # sd = A Δp
    mul!(dp, B', Δy)                 # dp = Bᵀ Δy
    pfloor = max(ptol, 100 * eps(T) * (norm(rp) + norm(sp)))
    dfloor = max(ytol, 100 * eps(T) * (norm(f) + norm(sd) + norm(dp)))
    axpby!(one(T), rp, -one(T), sp)                        # sp = rp - B Δp
    axpby!(one(T), f, -one(T), sd); axpy!(one(T), dp, sd)  # sd = f - A Δp + Bᵀ Δy

    dres = norm(sd)
    pres = norm(sp)
    res  = max(pres / pfloor, dres / dfloor)

    if res ≤ one(T)
        # met the floor: SOLVED iff that floor was the requested tol, else STAGNATED (bottomed out at
        # roundoff short of tol — the direction is as accurate as arithmetic allows, but not to tol).
        st = (pres ≤ ptol && dres ≤ ytol) ? KKT_SOLVED : KKT_STAGNATED
        return st, pres, dres, 0
    end
    res > stall * max(pprv / pfloor, dprv / dfloor) && return KKT_STAGNATED, pres, dres, 0

    n, _ = solveuzw!(wrk.divwrk, wrk.itrwrk, dp, dy, wrk.r, wrk.F, A, B,
                     sd, sp, wrk.α[]; atol = ptol, rtol = T(INTERIOR_THETA))
    axpy!(one(T), dp, Δp)
    axpy!(one(T), dy, Δy)

    return KKT_CONTINUE, pres, dres, n
end

#
# Solve the 2-row KKT system
#
#   [ A  -Bᵀ ] [ Δp ]   [ f  ]
#   [ B   0  ] [ Δy ] = [ rp ]
#
# to the caller's residual targets, via one
# base solve (solveuzw!) plus the solver-
# owned iterative-refinement loop (Govaerts–
# Pryce BE+1). Drives
#
#   ‖ rp - B Δp ‖        ≤ ptol
#   ‖ f - A Δp + Bᵀ Δy ‖ ≤ ytol
#
# and returns a typed exit: KKT_SOLVED, or
# KKT_STAGNATED / KKT_ITMAX saying why
# not. Refinement scratch is workspace-
# resident. Returns (nbase, nrefine, npass,
# status, fres, entry_pres, entry_dres).
#
function solvekkt!(
    wrk::UzawaSolver{UPLO, T},
    Δp::AbstractVector{T},
    Δy::AbstractVector{T},
    A::BlockSparseMatrix{T},
    B::BlockSparseMatrix{T},
    f::AbstractVector{T},
    rp::AbstractVector{T};
    pwarm::Bool = false,
    ywarm::Bool = false,
    ptol::T,
    ytol::T,
    stall::T,
    itmax::Int,
) where {UPLO, T}
    #
    # craig's exit residual is the next refinement pass's primal residual sp (rp − B(Δp+dp) = sp − B·dp),
    # and the loop's convergence test is ‖sp‖₂ ≤ ptol. So craig's target IS ptol — no over-solving to a
    # bare tolerance. Only ptol appears: craig owns the primal row; the dual row is craig-independent.
    #
    atol = ptol
    sd = wrk.sd; sp = wrk.sp; dp = wrk.dp; dy = wrk.dy
    #
    # base solve. pwarm/ywarm declare whether Δp/Δy already hold a live seed. A cold buffer is zeroed and
    # the base runs directly on (f, rp); a warm buffer is corrected — its residual drives one base-quality
    # (rtol = 0) solveuzw!, folded back in. Both are the same solveuzw!: with a zeroed buffer the residual
    # is (f, rp) bitwise, so the direct form is identical and just skips a matvec on zero.
    #
    copyto!(sd, f)
    copyto!(sp, rp)

    if pwarm
        mul!(sd, Symmetric(A, :L), Δp, -one(T), one(T))
        mul!(sp,           B,      Δp, -one(T), one(T))
    end

    if ywarm
        mul!(sd, B', Δy, one(T), one(T))
    end

    nbase, fres = solveuzw!(wrk.divwrk, wrk.itrwrk, dp, dy, wrk.r, wrk.F, A, B, sd, sp, wrk.α[]; atol)

    if pwarm
        axpy!(one(T), dp, Δp)
    else
        copyto!(Δp, dp)
    end

    if ywarm
        axpy!(one(T), dy, Δy)
    else
        copyto!(Δy, dy)
    end

    #
    # iterative refinement of the 2-row residual
    #
    nrefine = 0
    npass = 0
    status = KKT_ITMAX
    pprv = typemax(T)
    dprv = typemax(T)
    pres1 = T(NaN)
    dres1 = T(NaN)

    for i in 1:itmax
        st, pres, dres, n = refinekkt!(wrk, Δp, Δy, A, B, f, rp, pprv, dprv; ptol, ytol, stall)

        if isone(i)
            pres1 = pres
            dres1 = dres
        end

        if st != KKT_CONTINUE
            status = st
            break
        end

        pprv = pres
        dprv = dres
        nrefine += n
        npass += 1
    end

    return nbase, nrefine, npass, status, fres, pres1, dres1
end

#
# Solve the KKT system
#
#   [ A -Bᵀ ] [ x ] = [ f ]
#   [ B  0  ] [ y ]   [ g ]
#
# where A is a n x n positive-definite
# matrix and B is a m × n matrix, α is a
# positive number, and L is the n × n
# Cholesky factor of the augmented matrix
#
#   1/α A + Bᵀ B = L Lᵀ.
#
# The solution (x, y) meets the
# following tolerances
#
#   ‖ f - A x + Bᵀ y ‖ ≤ c ( α ‖ |L| |Lᵀ| |x| ‖ + ‖ |A| |x| + |Bᵀ| |y| + |f| ‖ )
#   ‖ Bx - g ‖         ≤ ϵ
#
# where ϵ is a specified tolerance (`atol`) and c
# is a small constant depending on m and n.
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
        α::T;
        atol::T = √eps(T),
        rtol::T = zero(T),
        itmax::Int = sum(size(B)),
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
    #   (1/α A + Bᵀ B) x =  1/α f + Bᵀ g
    #
    copyto!(r, g)
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
    #   [ 1/α A + Bᵀ B  -Bᵀ ] [ δx ] = [ 0 ]
    #   [            B   0  ] [ δy ]   [ r ]
    #
    # and update
    #
    #   x ← x +   δx
    #   r ← r - B δx
    #
    nc, _ = craig!(itrwrk, B, F, divwrk, x, y, r; atol = max(atol, rtol * nr0), itmax)
    niter += nc
    #
    # recover y:
    #
    #   y = α (δy + r)
    #
    axpy!(one(T), r, y)
    lmul!(α, y)

    return niter, nr0
end
