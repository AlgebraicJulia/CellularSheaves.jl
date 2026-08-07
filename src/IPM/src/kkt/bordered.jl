# Bordered KKT solver — wraps a 2-row KKTSolver and adds the HSD τ-border via a Woodbury/Schur update.
# The Woodbury column w₂ = KKT⁻¹[c; g] and the capacitance scalar S are cached for a given factorization
# and reused by the predictor and corrector at that factorization. The column starts as a bare backsolve
# (1 application) and is driven per consumer to targets derived from dt = num/S (Change 2), so S moves as
# the column tightens — cache (Δp2, Δy2, S) together and recompute S whenever the column moves. The border
# data (c, g, Q, Hc) is passed per call, not stored. initkkt! factors + bare-backsolves w₂ + caches
# capacitance; solvekkt! does the 2-row base solve + column loop + Schur lift + 3-row refinement net.
struct BorderedSolver{T, W <: KKTSolver{T}} <: KKTSolver{T}
    inner::W          # the 2-row solver (factor + base solve primitive)
    Δp2::FVector{T}    # cached Woodbury column, primal part — CRAIG accumulates δx into it (n)
    Δy2::FVector{T}    # cached Woodbury column, dual part — re-recovered per consumer (m)
    dya::FVector{T}    # accumulated RAW CRAIG dual Σδyᵢ across the column's continued solves (m)
    yseed::FVector{T}  # cross-iteration multiplier seed y0 (= s.Δy0), used by the dual recovery (m)
    Δp0::FVector{T}    # base solve primal Δp₀ (n)
    Δy0::FVector{T}    # base solve dual   Δy₀ (m)
    aτ::FVector{T}     # cached border row c - 2Qp/τ (n)
    S::FScalar{T}      # cached capacitance scalar
    δ::FScalar{T}      # cached τ-row diagonal pᵀQp/τ² + κ/τ
end

function BorderedSolver(inner::W, B::BlockSparseMatrix{T, I}) where {T, I, W <: KKTSolver{T}}
    m, n = size(B)
    S = FScalar{T}(undef); S[] = one(T)
    δ = FScalar{T}(undef); δ[] = zero(T)
    return BorderedSolver(inner,
                          FVector{T}(undef, n), FVector{T}(undef, m),   # Δp2, Δy2
                          FVector{T}(undef, m), FVector{T}(undef, m),   # dya, yseed
                          FVector{T}(undef, n), FVector{T}(undef, m),   # Δp0, Δy0
                          FVector{T}(undef, n), S, δ)                   # aτ, S, δ
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
# wpass, wstat, wfres, wpres, wdres) — the
# column solve's history row.
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
    ok || return ok, ρ, (0, 0, 0, KKT_ITMAX, T(NaN), T(NaN), T(NaN))
    #
    # bare backsolve for the Woodbury column [ H -Bᵀ; B 0 ][Δp2; Δy2] = [c; g] — ONE application (itmax=0
    # skips CRAIG), no refinement (Change 2). The column is not solved to a tolerance here: its required
    # accuracy is set per consumer by dt = num/S (computable only after S, which this backsolve yields).
    # Solving it up front to the Newton tolerances over-solves fixed data (c,g) by up to 10 decades late.
    #
    inr = bw.inner
    β = inv(inr.α[]); r = inr.r
    #
    # base solve of the column — the inlined `solveuzw!` base (no CRAIG here; itmax=0). solvekkt! continues
    # this same CRAIG per consumer, to the accuracy the lift needs. The cross-iteration seed y0 enters both
    # the RHS and the dual recovery; stash it (yseed) for solvekkt!, and reset the δy accumulator.
    #
    #   G·Δp2 = β·c + Bᵀ(g + β·y0) ,   G = βH + BᵀB = F Fᵀ
    #   r     = g - B·Δp2
    #   Δy2   = y0 + α·r            (base: δy_total = 0)
    #
    copyto!(r, g)
    isnothing(y0) || axpy!(β, y0, r)
    copyto!(bw.Δp2, c)
    mul!(bw.Δp2, B', r, one(T), β)
    ldiv!(inr.divwrk, inr.F,  bw.Δp2)
    ldiv!(inr.divwrk, inr.F', bw.Δp2)
    copyto!(r, g)
    mul!(r, B, bw.Δp2, -one(T), one(T))
    r2 = norm(r)
    copyto!(bw.Δy2, r); lmul!(inr.α[], bw.Δy2)
    isnothing(y0) || axpy!(one(T), y0, bw.Δy2)
    fill!(bw.dya, zero(T))
    isnothing(y0) ? fill!(bw.yseed, zero(T)) : copyto!(bw.yseed, y0)
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
    return ok, ρ, (1, 0, 0, KKT_SOLVED, r2, r2, T(NaN))
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
# the 3-row refinement net. The column re-
# solve CRAIG folds into nbase; refinement
# scratch is the inner solver's. Returns
# (nbase, nrefine, npass, status, fres,
# entry_pres, entry_dres, Δτ).
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
    sp = inner.sp; sd = inner.sd
    #
    # base solve for the column-free part Δp₀, Δy₀ — the 2-row system to (ptol, ytol). pwarm/ywarm just
    # forward to the inner solver: a warm caller has pre-seeded Δp/Δy, a cold one gets them zeroed.
    #
    #   [ H  -Bᵀ ] [ Δp₀ ]   [ f  ]
    #   [ B   0  ] [ Δy₀ ] = [ rp ]
    #
    nb, nr, _, _, fres, _, _ = solvekkt!(inner, Δp, Δy, H, B, f, rp; pwarm, ywarm, ptol, ytol, stall, itmax)
    nbase = nb + nr
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
    # consumer's lift needs, then re-recover the dual and refresh S. ONE craig call, not a loop: tgt is
    # first-order in S (the pre-tighten value) and craig iterates internally to it. Only the primal is
    # driven; the dual is refinehsd!'s job (§ pass-economy) — the recovery just keeps Δy2 feasible for the
    # tightened Δp2. δy accumulates in dya (craig resets Δy2 each call), so a corrector continuing from the
    # predictor carries δy_pred + δy_corr.
    #
    #   dt  = num / S ,   tgt = max( (ptol - rb)/|dt| , 100u(‖g‖ + ‖B Δp2‖) )
    #   r   = g - B Δp2
    #   craig!(Δp2, Δy2, r; atol = tgt)              Δp2 += δx ; r ← g - B Δp2 ; Δy2 = δy_incr
    #   dya += Δy2 ;  Δy2 = yseed + α (dya + r) ;  S = aτᵀ Δp2 + gᵀ Δy2 + δ
    #
    dt  = num / bw.S[]
    mul!(sp, B, bw.Δp2)                                   # sp = B Δp2
    tgt = max((ptol - rb) / abs(dt), 100 * eps(T) * (ng + norm(sp)))
    axpby!(one(T), g, -one(T), sp)                       # sp = g - B Δp2  (column residual)
    nc, _ = craig!(inner.itrwrk, B, inner.F, inner.divwrk, bw.Δp2, bw.Δy2, sp; atol = tgt)
    nbase += nc
    axpy!(one(T), bw.Δy2, bw.dya)                        # dya += δy_incr
    copyto!(bw.Δy2, bw.dya); axpy!(one(T), sp, bw.Δy2)   # Δy2 = dya + r
    lmul!(inner.α[], bw.Δy2)                             # Δy2 = α (dya + r)
    axpy!(one(T), bw.yseed, bw.Δy2)                      # Δy2 = yseed + α (dya + r)
    bw.S[] = dot(bw.aτ, bw.Δp2) + dot(g, bw.Δy2) + bw.δ[]
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
    # 3-row refinement net (a single verifying pass when both column targets were met):
    #
    npass, nrefine, status, Δτ, pres1, dres1 = refinehsd!(
        Δp, Δy, Δτ, inner, H, B, c, g, bw.aτ, bw.δ[], rp, f, fτ, bw.Δp2, bw.Δy2, bw.S[],
        sp, sd, inner.dp, inner.dy;
        itmax, ptol, ytol, τtol, stall,
    )
    return nbase, nrefine, npass, status, fres, pres1, dres1, Δτ
end

############################################################################################
# refinehsd! — Govaerts–Pryce BE+1 refinement on the 3-row bordered system
############################################################################################

#
# Iterative refinement of the 3-row bordered
# system. Each pass forms the residual
#
#   [ sd ]   [ f  ]   [ H   -Bᵀ  -c ] [ Δp ]
#   [ sp ] = [ rp ] - [ B    0   -g ] [ Δy ]
#   [ sτ ]   [ fτ ]   [ aτᵀ  gᵀ   δ ] [ Δτ ]
#
# solves the 2-row correction (solveuzw!),
# Schur-lifts dτ = (sτ - aτᵀdp - gᵀdy)/S,
# and updates — until ‖sp‖ ≤ ptol ∧
# ‖sd‖ ≤ ytol ∧ |sτ| ≤ τtol (KKT_SOLVED), or it
# stalls / runs out of passes. Scratch
# (sp, sd, dp, dy) is the caller's. Returns
# (npass, nrefine, status, Δτ, entry_pres,
# entry_dres).
#
function refinehsd!(
    Δp::AbstractVector{T},
    Δy::AbstractVector{T},
    Δτ::T,
    wrk::UzawaSolver{UPLO, T},
    H::BlockSparseMatrix{T},
    B::BlockSparseMatrix{T},
    c::AbstractVector{T},
    g::AbstractVector{T},
    aτ::AbstractVector{T}, # border row: c - 2Qp/τ
    δ::T,                  # τ-row diagonal pᵀQp/τ² + κ/τ
    rp::AbstractVector{T},  # RHS for row P
    f::AbstractVector{T},   # RHS for row D
    fτ::T,                  # RHS for row T
    Δp2::AbstractVector{T}, # Woodbury direction
    Δy2::AbstractVector{T},
    S::T,                   # Woodbury capacitance scalar
    sp::AbstractVector{T},  # scratch for primal residual
    sd::AbstractVector{T},  # scratch for dual residual
    dp::AbstractVector{T},  # scratch for correction
    dy::AbstractVector{T};
    itmax::Int,
    ptol::T,
    ytol::T,
    τtol::T,
    stall::T,
) where {UPLO, T}
    niter = 0
    npass = 0
    #
    # craig works in unscaled 2-norms; the loop measures pres = ‖sp‖₂ against ptol, so its craig tolerance
    # IS ptol. First-order in the bordered path: craig's exit residual is sp − B·dp_uzw, but the applied dp
    # is Schur-lifted (dp_uzw + dτ·Δp2), so sp' = craig_resid + dτ·r₂ (r₂ = g − B·Δp2, §3). So ptol is the
    # correct first-order translation, not an equality. Only ptol — the dual row is border-free.
    #
    atol = ptol
    status = KKT_ITMAX
    prv = typemax(T)
    pres1 = T(NaN)
    dres1 = T(NaN)

    for i in 1:itmax
        #
        # compute the residuals
        #
        #   [ sd ]   [ f  ]   [ H    -Bᵀ  -c ] [ Δp ]
        #   [ sp ] = [ rp ] - [ B     0   -g ] [ Δy ]
        #   [ sτ ]   [ fτ ]   [ aτᵀ   gᵀ   δ ] [ Δτ ]
        #
        mulkkt!(sd, sp, H, B, Δp, Δy)
        axpy!(-Δτ, c, sd)
        axpy!(-Δτ, g, sp)
        axpby!(one(T), f,  -one(T), sd)
        axpby!(one(T), rp, -one(T), sp)
        sτ = fτ - dot(aτ, Δp) - dot(g, Δy) - δ * Δτ

        pres = norm(sp)
        dres = norm(sd)
        τres = abs(sτ)

        if isone(i)
            pres1 = pres
            dres1 = dres
        end

        # τ (border) row: τtol is the caller's absolute target for the border residual (= tol·(1+nc+ng),
        # matching the old abs(sτ)/(1+nc+ng) ≤ tol test exactly). Converged iff all three rows are within
        # their targets: ‖sp‖ ≤ ptol ∧ ‖sd‖ ≤ ytol ∧ abs(sτ) ≤ τtol.
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
        #
        # solve for dp and dy:
        #
        #   [ H -Bᵀ ] [ dp ] = [ sd ]
        #   [ B  0  ] [ dy ]   [ sp ]
        #
        n, _ = solveuzw!(wrk.divwrk, wrk.itrwrk, dp, dy, wrk.r, wrk.F, H, B, sd, sp, wrk.α[]; atol, rtol = T(INTERIOR_THETA))
        niter += n
        npass += 1
        #
        # apply the Schur lift (aτ = c - 2Qp/τ):
        #
        #   dτ = (sτ - aτᵀdp - gᵀdy) / S
        #
        dτ = (sτ - dot(aτ, dp) - dot(g, dy)) / S
        #
        # apply border correction:
        #
        #   dp ← dp + dτ Δp2
        #   dy ← dy + dτ Δy2
        #
        axpy!(dτ, Δp2, dp)
        axpy!(dτ, Δy2, dy)
        #
        # update directions:
        #
        #   Δp ← Δp + dp
        #   Δy ← Δy + dy
        #   Δτ ← Δτ + dτ
        #
        axpy!(one(T), dp, Δp)
        axpy!(one(T), dy, Δy)
        Δτ += dτ
    end

    return npass, niter, status, Δτ, pres1, dres1
end
