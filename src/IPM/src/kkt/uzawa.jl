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
# Min-cost α-window for the 2-row solve, read off the base solve's own residual decay — the α-
# controller's level structure. Two thresholds bound the zero-refinement operating range:
#
#   α₀ = α · fres / ptol         primal one-shot floor: base primal residual scales ∝ 1/α, so α ≥ α₀
#                                is where the bare base already meets ptol (level 0).
#   α_c = α · ytol / s0          dual ceiling: base dual residual s0 scales ∝ α, so α ≤ α_c is where
#                                the dual meets ytol with no refinement.
#
# When α₀ > α_c no bare α clears both rows; CRAIG buys primal accuracy, letting a smaller α still meet
# ptol. Level k folds k CRAIG steps of the observed decay ‖r_k‖/‖r_0‖ into the primal floor:
#
#   α_k = ( α₀ · αᵏ · ‖r_k‖/‖r_0‖ )^{1/(k+1)}          (α_0 = α₀)
#
# decreasing in k. αmin is the α at the cheapest level that still clears the dual ceiling (smallest k
# with α_k ≤ α_c); αmax = α_c. hist[k+1] = ‖r_k‖ is the base CRAIG decay, still resident from the base
# solve — compute the window BEFORE any refinement pass overwrites it. Computed in log space so αᵏ does
# not overflow for a large α.
#
function kktwindow!(
        wrk::UzawaSolver{UPLO, T},
        A::BlockSparseMatrix{T},
        B::BlockSparseMatrix{T},
        f::AbstractVector{T},
        Δp::AbstractVector{T},
        Δy::AbstractVector{T},
        fres::T,
        nbase::Int;
        ptol::T,
        ytol::T,
    ) where {UPLO, T}
    #
    # base dual residual s0 = ‖ f - A Δp + Bᵀ Δy ‖ (unfused, as refinement forms it), accumulated into sd
    # ALONE — no separate Bᵀ Δy buffer. sd is dead scratch after the base solve (the bordered caller's
    # refinehsd! overwrites it via mulkkt! before reading), but wrk.dp is NOT free downstream: the HSD
    # column path reuses the inner solver's dp buffer, and clobbering it regresses the solve (verified: a
    # dp-clobber costs 07 87→56 CRAIG, an sd-clobber costs nothing). Hence the β=1 accumulate.
    #
    sd = wrk.sd
    mul!(sd, Symmetric(A, :L), Δp)   # sd = A Δp
    axpby!(one(T), f, -one(T), sd)   # sd = f - A Δp
    mul!(sd, B', Δy, one(T), one(T)) # sd = f - A Δp + Bᵀ Δy
    s0 = norm(sd)

    α  = wrk.α[]
    αc = s0 > zero(T) ? α * ytol / s0 : T(Inf)
    r0 = fres                        # ‖r_0‖ = base primal residual = hist[1]
    (r0 > zero(T) && α > zero(T)) || return zero(T), αc

    hist  = wrk.itrwrk.hist
    n0    = nbase - 1                # CRAIG steps the base took
    logα  = log(α)
    logr0 = log(r0)
    logα0 = logα + logr0 - log(ptol)
    logαc = αc == T(Inf) ? T(Inf) : log(αc)

    αmin = exp(logα0)                # level 0 fallback
    for k in 0:n0
        k + 1 ≤ length(hist) || break
        rk    = hist[k + 1]
        logrk = rk > zero(T) ? log(rk) : logr0 + log(eps(T))   # clamp an exact-zero residual
        logαk = (logα0 + k * logα + (logrk - logr0)) / (k + 1)
        αmin  = exp(logαk)
        logαk ≤ logαc && break
    end

    return αmin, αc
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
# status, fres, entry_pres, entry_dres, αmin,
# αmax) — the last two the min-cost α-window
# (kktwindow!) the controller geo-means.
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
    # min-cost α-window, read off the base solve's CRAIG decay (hist, still resident) and base dual
    # residual — BEFORE the refinement loop below overwrites hist with its correction solves.
    #
    αmin, αmax = kktwindow!(wrk, A, B, f, Δp, Δy, fres, nbase; ptol, ytol)

    #
    # iterative refinement of the 2-row residual (refinekkt! inlined). Each pass forms the residual with
    # the three matvecs kept SEPARATE so each row is floored at its own evaluation limit — u × the sum of
    # the magnitudes differenced. ‖A Δp‖ is H-amplified but cancels against ‖Bᵀ Δy‖ in A Δp - Bᵀ Δy, so a
    # fused floor would undercount exactly when it matters. Unless a row meets its floor (met tol →
    # SOLVED, roundoff-limited short of tol → STAGNATED) or the residual has stalled against the previous
    # pass, solve the correction to θ·(its entry residual) and fold it in.
    #
    nrefine = 0
    npass = 0
    status = KKT_ITMAX
    pprv = typemax(T)
    dprv = typemax(T)
    pres1 = T(NaN)
    dres1 = T(NaN)

    for i in 1:itmax
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

        if isone(i)
            pres1 = pres
            dres1 = dres
        end

        if res ≤ one(T)
            status = (pres ≤ ptol && dres ≤ ytol) ? KKT_SOLVED : KKT_STAGNATED
            break
        end
        if res > stall * max(pprv / pfloor, dprv / dfloor)
            status = KKT_STAGNATED
            break
        end

        n, _ = solveuzw!(wrk.divwrk, wrk.itrwrk, dp, dy, wrk.r, wrk.F, A, B,
                         sd, sp, wrk.α[]; atol = ptol, rtol = T(INTERIOR_THETA))
        axpy!(one(T), dp, Δp)
        axpy!(one(T), dy, Δy)

        pprv = pres
        dprv = dres
        nrefine += n
        npass += 1
    end

    return nbase, nrefine, npass, status, fres, pres1, dres1, αmin, αmax
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
