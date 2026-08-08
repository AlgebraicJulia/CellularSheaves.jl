# Bordered KKT solver — wraps a 2-row KKTSolver and adds the HSD τ-border via a Woodbury/Schur update.
# The Woodbury column w₂ = KKT⁻¹[c; g] and the capacitance scalar S are cached for a given factorization
# and reused by the predictor and corrector at that factorization. The column starts as a bare backsolve
# (1 application) and is driven per consumer to targets derived from dt = num/S (Change 2), so S moves as
# the column tightens — cache (Δp2, Δy2, S) together and recompute S whenever the column moves. The border
# data (c, g, Q, Hc) is passed per call, not stored. initkkt! factors + bare-backsolves w₂ + caches
# capacitance; solvekkt! does the 2-row base solve + column loop + Schur lift + 3-row refinement net.
struct BorderedSolver{T, W <: KKTSolver{T}} <: KKTSolver{T}
    inner::W          # the 2-row solver (factor + base solve primitive)
    Δp2::FVector{T}    # cached Woodbury column, primal part — carried across iterations as the seed (n)
    Δy2::FVector{T}    # cached Woodbury column, dual part — re-recovered per consumer, carried as seed (m)
    rcol::FVector{T}   # MAINTAINED column residual (the correction's r, never re-formed as g - B Δp2) (m)
    dya::FVector{T}    # accumulated RAW CRAIG dual Σδyᵢ across the column's continued solves (m)
    yseed::FVector{T}  # cross-iteration multiplier seed (= carried Δy2), used by the dual recovery (m)
    Δp0::FVector{T}    # direction base solve primal Δp₀, kept for the Schur lift (n)
    Δy0::FVector{T}    # direction base solve dual   Δy₀ (m)
    aτ::FVector{T}     # cached border row c - 2Qp/τ (n)
    S::FScalar{T}      # cached capacitance scalar
    δ::FScalar{T}      # cached τ-row diagonal pᵀQp/τ² + κ/τ
end

function BorderedSolver(inner::W, B::BlockSparseMatrix{T, I}) where {T, I, W <: KKTSolver{T}}
    m, n = size(B)
    S = FScalar{T}(undef); S[] = one(T)
    δ = FScalar{T}(undef); δ[] = zero(T)
    Δp2 = FVector{T}(undef, n); fill!(Δp2, zero(T))   # zero seed on the first factorization (cold bootstrap)
    Δy2 = FVector{T}(undef, m); fill!(Δy2, zero(T))
    return BorderedSolver(inner, Δp2, Δy2,
                          FVector{T}(undef, m),                        # rcol
                          FVector{T}(undef, m), FVector{T}(undef, m),  # dya, yseed
                          FVector{T}(undef, n), FVector{T}(undef, m),  # Δp0, Δy0
                          FVector{T}(undef, n), S, δ)                  # aτ, S, δ
end

############################################################################################
# initkkt! (factor + bare-backsolve column + cache S) / solvekkt! (base + column loop + lift)
############################################################################################

#
# Factor the border once and cache what the
# per-solve Schur lift reuses. Factors the
# inner F, bare-backsolves the Woodbury
# column
#
#   [ H  -Bᵀ ] [ Δp2 ]   [ c ]
#   [ B   0  ] [ Δy2 ] = [ g ]
#
# (one F-apply, no CRAIG — its accuracy is
# driven per consumer, not here), then caches
# the capacitance S and border row aτ. Runs
# once per factorization, shared by the
# predictor and corrector. Returns
# (ok, ρ, wtuple), wtuple = (wbase, wrefn,
# wpass, wstat) — the column base's counts.
#
function initkkt!(
        bw::BorderedSolver{T},
        H::BlockSparseMatrix{T},
        B::BlockSparseMatrix{T},
        c::AbstractVector{T},
        g::AbstractVector{T};
        α::T,
        rgmin::T,
        p::AbstractVector{T},
        τ::T,
        κ::T,
        Qp::AbstractVector{T},
        y0,
    ) where {T}
    ok, ρ = initkkt!(bw.inner, H; α, rgmin)
    ok || return ok, ρ, (0, 0, 0, KKT_ITMAX)
    #
    # column base as a WARM CORRECTION seeded by the carried column (x0 = Δp2, y0 = Δy2 from the previous
    # factorization; both 0 on the first). The primal seed cancels — Δp2 = G⁻¹(βc + Bᵀ(g + βy0)) is identical
    # to a cold solve (§2 of the p-warm note) — but the residual we KEEP is the correction's residual
    # rcol = sp - B·xb, NOT g - B·Δp2. THAT is the quantity the dual recovery consumes (α-amplified), and
    # solvekkt!'s CRAIG maintains rcol in place; it is never re-formed against the accumulated column — that
    # would re-cancel O(1) quantities to u‖g‖ and throw the gain away. ONE F-apply, no CRAIG here. y0 (the
    # carried Δy2) is the cross-iteration multiplier seed; stash it as yseed and reset the δy accumulator.
    #
    #   sp = g - B·x0 ,   sd = c - H·x0 + Bᵀ·y0
    #   xb = G⁻¹(β·sd + Bᵀ·sp) ,   Δp2 = x0 + xb
    #   rcol = sp - B·xb ,   Δy2 = yseed + α·rcol
    #
    inr = bw.inner
    β = inv(inr.α[])
    copyto!(bw.yseed, bw.Δy2)                              # carried multiplier seed (before the recovery overwrites Δy2)
    mulkkt!(inr.sd, inr.sp, H, B, bw.Δp2, bw.Δy2)          # inr.sd = H·x0 - Bᵀ·y0 ,  inr.sp = B·x0
    axpby!(one(T), c, -one(T), inr.sd)                     # inr.sd = c - H·x0 + Bᵀ·y0 = sd
    copyto!(bw.rcol, g); axpy!(-one(T), inr.sp, bw.rcol)   # bw.rcol = g - B·x0 = sp
    copyto!(inr.dp, inr.sd); rmul!(inr.dp, β)              # inr.dp = β·sd
    mul!(inr.dp, B', bw.rcol, one(T), one(T))              # inr.dp = β·sd + Bᵀ·sp
    ldiv!(inr.divwrk, inr.F,  inr.dp)
    ldiv!(inr.divwrk, inr.F', inr.dp)                      # inr.dp = xb = G⁻¹(β·sd + Bᵀ·sp)
    axpy!(one(T), inr.dp, bw.Δp2)                          # Δp2 = x0 + xb
    mul!(bw.rcol, B, inr.dp, -one(T), one(T))              # bw.rcol = sp - B·xb  (MAINTAINED residual)
    copyto!(bw.Δy2, bw.rcol); lmul!(inr.α[], bw.Δy2); axpy!(one(T), bw.yseed, bw.Δy2)   # Δy2 = yseed + α·rcol
    fill!(bw.dya, zero(T))
    #
    # cache the border row aτ = c - 2Qp/τ and the capacitance — the Schur complement of the τ row,
    #
    #   S = aτᵀ Δp2 + gᵀ Δy2 + δ,   δ = pᵀQp/τ² + κ/τ
    #
    # by inner product rather than the algebraically-equal quadratic form: sharing the aτᵀΔp2, gᵀΔy2
    # rounding with the lift's numerator keeps the τ row exact however loosely the column was solved.
    #
    bw.δ[] = dot(p, Qp) / τ^2 + κ / τ
    copyto!(bw.aτ, c)
    axpby!(-2 / τ, Qp, one(T), bw.aτ)
    bw.S[] = dot(bw.aτ, bw.Δp2) + dot(g, bw.Δy2) + bw.δ[]
    return ok, ρ, (1, 0, 0, KKT_SOLVED)
end

#
# Solve the 3-row HSD bordered system
#
#   [ H    -Bᵀ  -c ] [ Δp ]   [ f  ]
#   [ B     0   -g ] [ Δy ] = [ rp ]
#   [ aτᵀ   gᵀ   δ ] [ Δτ ]   [ fτ ]
#
# by Schur complement on the τ row: base-
# solve the 2-row part, drive the cached
# Woodbury column to the accuracy the lift
# needs, apply the lift Δτ = num/S, then run
# the 3-row refinement net. Refinement scratch
# is the inner solver's. Returns (niter, ncol,
# npass, status, Δτ, αmin, αmax, wamin): niter
# the direction CRAIG (base + refinement), αmin
# from the base CRAIG levels, αmax from the
# 3-row dual, wamin the Woodbury column floor
# (kktwindow! on the column).
#
function solvekkt!(
        bw::BorderedSolver{T},
        Δp::AbstractVector{T},
        Δy::AbstractVector{T},
        H::BlockSparseMatrix{T},
        B::BlockSparseMatrix{T},
        c::AbstractVector{T},
        g::AbstractVector{T},
        ng::T,                 # cached ‖g‖ (solver-lifetime constant)
        f::AbstractVector{T},
        rp::AbstractVector{T},
        fτ::T;
        pwarm::Bool = false,
        ywarm::Bool = false,
        ptol::T,
        ytol::T,
        τtol::T,
        stall::T,
        itmax::Int,
    ) where {T}
    inner = bw.inner
    sp = inner.sp; sd = inner.sd; dp = inner.dp; dy = inner.dy
    α = inner.α[]
    #
    # direction base for the column-free part Δp₀, Δy₀ — the bare 2-row base solve (itmax=0, refinehsd! is
    # the sole refinement), inlined. pwarm/ywarm: a warm caller has pre-seeded Δp/Δy (corrected against the
    # residual), a cold one gets them copied straight from the base correction.
    #
    #   [ H  -Bᵀ ] [ Δp₀ ]   [ f  ]
    #   [ B   0  ] [ Δy₀ ] = [ rp ]
    #
    copyto!(sd, f); copyto!(sp, rp)
    if pwarm
        mul!(sd, Symmetric(H, :L), Δp, -one(T), one(T))
        mul!(sp, B,               Δp, -one(T), one(T))
    end
    if ywarm
        mul!(sd, B', Δy, one(T), one(T))
    end
    nbase, fres = solveuzw!(inner.divwrk, inner.itrwrk, dp, dy, inner.r, inner.F, H, B, sd, sp, α; atol = ptol)
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
    # primal window floor αmin from the base CRAIG decay (hist still resident — read it BEFORE the column
    # CRAIG below overwrites it). kktwindow!'s αmax is the 2-row dual ceiling; the bordered ceiling is the
    # 3-row dual, computed in the refinement net below, so take only αmin here.
    #
    αmin, _ = kktwindow!(inner, H, B, f, Δp, Δy, fres, nbase; ptol, ytol)
    #
    # preserve the base solve Δp₀, Δy₀: the column loop re-lifts Δp = Δp₀ + Δτ Δp2 from them every pass, so
    # they must survive the in-place lift into (Δp, Δy).
    #
    copyto!(bw.Δp0, Δp)
    copyto!(bw.Δy0, Δy)
    #
    # the base solve's primal budget rb and the τ-row numerator num (both column-free):
    #
    #   rb  = ‖ rp - B Δp₀ ‖
    #   num = fτ - aτᵀ Δp₀ - gᵀ Δy₀
    #
    mul!(sp, B, Δp)
    axpby!(one(T), rp, -one(T), sp)
    rb  = norm(sp)
    num = fτ - dot(bw.aτ, Δp) - dot(g, Δy)
    #
    # column tighten — continue the column's CRAIG (base done once in initkkt!) to the accuracy THIS
    # consumer's lift needs, then re-recover the dual and refresh S. CRAIG continues on the MAINTAINED
    # residual bw.rcol — it is never re-formed as g - B Δp2 (that re-cancels to u‖g‖ and loses the dual
    # gain, p-warm note §3). Only the primal is driven; the dual is refinehsd!'s job — the recovery keeps
    # Δy2 feasible. δy accumulates in dya (craig resets Δy2 each call), so a corrector continuing from the
    # predictor carries δy_pred + δy_corr. tgt is first-order in S; the floor uses ‖B Δp2‖ = ‖g - rcol‖.
    #
    #   dt  = num / S ,   tgt = max( (ptol - rb)/|dt| , 100u(‖g‖ + ‖B Δp2‖) )
    #   craig!(Δp2, Δy2, rcol; atol = tgt)           Δp2 += δx ; rcol maintained ; Δy2 = δy_incr
    #   dya += Δy2 ;  Δy2 = yseed + α (dya + rcol) ;  S = aτᵀ Δp2 + gᵀ Δy2 + δ
    #
    dt  = num / bw.S[]
    copyto!(sp, g); axpy!(-one(T), bw.rcol, sp)          # sp = g - rcol = B Δp2  (floor only; no matvec)
    tgt = max((ptol - rb) / abs(dt), 100 * eps(T) * (ng + norm(sp)))
    ngcol = norm(bw.rcol)                                # column entry residual (= hist[1]) for the wamin floor
    ncol, _ = craig!(inner.itrwrk, B, inner.F, inner.divwrk, bw.Δp2, bw.Δy2, bw.rcol; atol = tgt)   # Woodbury/column CRAIG — reported separately, not folded into the direction cost
    axpy!(one(T), bw.Δy2, bw.dya)                        # dya += δy_incr
    copyto!(bw.Δy2, bw.dya); axpy!(one(T), bw.rcol, bw.Δy2)   # Δy2 = dya + rcol
    lmul!(inner.α[], bw.Δy2)                             # Δy2 = α (dya + rcol)
    axpy!(one(T), bw.yseed, bw.Δy2)                      # Δy2 = yseed + α (dya + rcol)
    bw.S[] = dot(bw.aτ, bw.Δp2) + dot(g, bw.Δy2) + bw.δ[]
    #
    # Woodbury floor wamin — the α below which the column base no longer meets ptol (needs CRAIG), computed
    # the SAME way as the predictor/corrector floors: kktwindow! on the column solve's own CRAIG decay (now
    # resident in the inner workspace). Only the floor is meaningful for the column, so αmax is discarded —
    # the column has no refinement/dual ceiling. Runs before the refinement net overwrites the decay.
    #
    wamin, _ = kktwindow!(inner, H, B, c, bw.Δp2, bw.Δy2, ngcol, ncol + 1; ptol, ytol)
    #
    # final lift with the converged S — Δp, Δy already carry it on a satisfied/stalled exit; this makes them
    # consistent with S on an itmax exit too. Δτ pairs with the same S (τ-row exactness):
    #
    #   Δτ ← num / S ,   Δp ← Δp₀ + Δτ Δp2 ,   Δy ← Δy₀ + Δτ Δy2
    #
    Δτ = num / bw.S[]
    copyto!(Δp, bw.Δp0); axpy!(Δτ, bw.Δp2, Δp)
    copyto!(Δy, bw.Δy0); axpy!(Δτ, bw.Δy2, Δy)
    #
    # 3-row refinement net (refinehsd! inlined — Govaerts–Pryce BE+1). Each pass forms the residual
    #
    #   [ sd ]   [ f  ]   [ H    -Bᵀ  -c ] [ Δp ]
    #   [ sp ] = [ rp ] - [ B     0   -g ] [ Δy ]
    #   [ sτ ]   [ fτ ]   [ aτᵀ   gᵀ   δ ] [ Δτ ]
    #
    # solves the 2-row correction (solveuzw!) to θ·(entry residual), Schur-lifts dτ = (sτ - aτᵀdp -
    # gᵀdy)/S, and updates — until ‖sp‖ ≤ ptol ∧ ‖sd‖ ≤ ytol ∧ |sτ| ≤ τtol (KKT_SOLVED), or it stalls /
    # runs out of passes. craig works in unscaled 2-norms so its target is ptol (first-order: the applied
    # dp is Schur-lifted, so sp' = craig_resid + dτ·r₂). Scratch (sp, sd, dp, dy) is the inner solver's.
    #
    atol = ptol
    status = KKT_ITMAX
    npass = 0; nrefine = 0
    prv = typemax(T)
    dres1 = T(NaN)   # entry 3-row dual residual, kept for the αmax ceiling

    for i in 1:itmax
        mulkkt!(sd, sp, H, B, Δp, Δy)
        axpy!(-Δτ, c, sd)
        axpy!(-Δτ, g, sp)
        axpby!(one(T), f,  -one(T), sd)
        axpby!(one(T), rp, -one(T), sp)
        sτ = fτ - dot(bw.aτ, Δp) - dot(g, Δy) - bw.δ[] * Δτ

        pres = norm(sp)
        dres = norm(sd)
        τres = abs(sτ)

        isone(i) && (dres1 = dres)

        # converged iff all three rows are within their absolute targets
        res = max(pres / ptol, dres / ytol, τres / τtol)

        if res ≤ one(T)
            status = KKT_SOLVED
            break
        end
        if res > stall * prv
            status = KKT_STAGNATED
            break
        end
        prv = res

        n, _ = solveuzw!(inner.divwrk, inner.itrwrk, dp, dy, inner.r, inner.F, H, B, sd, sp, α; atol, rtol = T(INTERIOR_THETA))
        nrefine += n
        npass += 1
        #
        # Schur lift + border correction:  dτ = (sτ - aτᵀdp - gᵀdy)/S ; dp += dτ Δp2 ; dy += dτ Δy2
        #
        dτ = (sτ - dot(bw.aτ, dp) - dot(g, dy)) / bw.S[]
        axpy!(dτ, bw.Δp2, dp)
        axpy!(dτ, bw.Δy2, dy)
        axpy!(one(T), dp, Δp)
        axpy!(one(T), dy, Δy)
        Δτ += dτ
    end
    #
    # 3-row dual ceiling: αmax = α·ytol / (base 3-row dual residual). dres1 is the entry (base-direction)
    # dual residual ‖f - H Δp + Bᵀ Δy + Δτ c‖ — the bordered ceiling, distinct from the 2-row dual.
    #
    αmax = (isfinite(dres1) && dres1 > zero(T)) ? α * ytol / dres1 : T(Inf)
    return nbase + nrefine, ncol, npass, status, Δτ, αmin, αmax, wamin
end

