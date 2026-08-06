# Bordered KKT solver — wraps a 2-row KKTSolver and adds the HSD τ-border via a Woodbury/Schur update.
# The Woodbury column w₂ = KKT⁻¹[c; g] and the capacitance scalar S are cached for a given factorization
# and reused by the predictor and corrector at that factorization. The column starts as a bare backsolve
# (1 application) and is driven per consumer to targets derived from dt = num/S (Change 2), so S moves as
# the column tightens — cache (Δp2, Δy2, S) together and recompute S whenever the column moves. The border
# data (c, g, Q, Hc) is passed per call, not stored. initkkt! factors + bare-backsolves w₂ + caches
# capacitance; solvekkt! does the 2-row base solve + column loop + Schur lift + 3-row refinement net.
struct BorderedSolver{T, W <: KKTSolver{T}} <: KKTSolver{T}
    inner::W          # the 2-row solver (factor + base solve primitive)
    Δp2::FVector{T}    # cached Woodbury column, primal part (n)
    Δy2::FVector{T}    # cached Woodbury column, dual part (m)
    aτ::FVector{T}     # cached border row c - 2Qp/τ (n)
    QΔp2::FVector{T}   # capacitance scratch (n)
    S::FScalar{T}      # cached capacitance scalar
end

function BorderedSolver(inner::W, B::BlockSparseMatrix{T, I}) where {T, I, W <: KKTSolver{T}}
    m, n = size(B)
    S = FScalar{T}(undef)
    S[] = one(T)
    return BorderedSolver(inner, FVector{T}(undef, n), FVector{T}(undef, m),
                          FVector{T}(undef, n), FVector{T}(undef, n), S)
end

############################################################################################
# initkkt! (bare backsolve) / solvekkt! (base + column loop + lift) / capacitance!
############################################################################################

# solve for the Woodbury auxiliary directions
#
#   [ H  -Bᵀ ] [ Δp2 ]   [ c ]
#   [ B   0  ] [ Δy2 ] = [ g ]
#
# BorderedSolver do-once: factor F (inner, combined H), solve & cache the Woodbury column
# w₂ = KKT⁻¹[c; g] to atol, and compute & cache the capacitance scalar S and border row aτ. Returns
# (ok, ρ, (wbase, wrefn, wpass, wstat, wfres, wpres, wdres)) — the border solve's counts for the history.
function initkkt!(
        bw::BorderedSolver{T},
        H::BlockSparseMatrix{T},
        Hc::BlockSparseMatrix{T},
        B::BlockSparseMatrix{T},
        c::AbstractVector{T},
        g::AbstractVector{T},
        Q::BlockSparseMatrix{T};
        α::T,
        rgmin::T,
        p::AbstractVector{T},
        τ::T,
        κ::T,
        Qp::AbstractVector{T},
        y0,
        ptol::T,
        ytol::T,
        stall::T,
        itmax::Int,
    ) where {T}
    ok, ρ = initkkt!(bw.inner, H; α, rgmin)
    ok || return ok, ρ, (0, 0, 0, REFINE_ITMAX, T(NaN), T(NaN), T(NaN))
    #
    # bare backsolve for the Woodbury column [ H -Bᵀ; B 0 ][Δp2; Δy2] = [c; g] — ONE application, no CRAIG,
    # no refinement (Change 2). The column is not solved to a tolerance here: its required accuracy is set
    # per consumer by dt = num/S (computable only after S, which this backsolve yields). Solving it up
    # front to the Newton tolerances over-solves fixed data (c,g) by up to 10 decades late in a run.
    #
    napp, r2 = backsolve!(bw.inner, bw.Δp2, bw.Δy2, H, B, c, g, y0)
    #
    # cache the capacitance scalar S and border row aτ = c - 2Qp/τ (stable split form)
    #
    bw.S[] = capacitance!(bw.QΔp2, bw.aτ, τ, κ, bw.Δp2, bw.Δy2, c, g, Qp, p, Hc, Q)
    return ok, ρ, (napp, 0, 0, REACHED, r2, r2, T(NaN))
end

# BorderedSolver solve: 2-row base solve to (ptol, ytol), then the column loop (tighten the Woodbury
# column to the consumer's dt-derived targets, recomputing S), then the Schur lift, then refinehsd!
# (3-row refinement net). Refinement scratch is the inner solver's. Column re-solve CRAIG folds into
# nbase. Returns (nbase, nrefine, npass, status, fres, entry_pres, entry_dres, Δτ).
function solvekkt!(
        bw::BorderedSolver{T},
        Δp::AbstractVector{T},
        Δy::AbstractVector{T},
        H::BlockSparseMatrix{T},
        B::BlockSparseMatrix{T},
        c::AbstractVector{T},
        g::AbstractVector{T},
        Qp::AbstractVector{T},
        p::AbstractVector{T},
        τ::T,
        κ::T,
        f::AbstractVector{T},
        rp::AbstractVector{T},
        fτ::T,
        y0 = nothing;
        ptol::T,
        ytol::T,
        τtol::T,
        stall::T,
        itmax::Int,
    ) where {T}
    inner = bw.inner
    #
    # base solve: the 2-row Newton system to (ptol, ytol) — backsolve + CRAIG (primal) + IR (dual). num
    # and rb come from this and need no column; this is the one reordering — base solve BEFORE the column.
    #
    nb, nr, _, _, fres, _, _ = solvekkt!(inner, Δp, Δy, H, B, f, rp, y0; ptol, ytol, stall, itmax)
    nbase = nb + nr
    #
    # rb = ‖rp - B·Δp0‖ (primal budget already spent); num = fτ - aτᵀΔp0 - gᵀΔy0 (needs no column).
    #
    sp = inner.sp; sd = inner.sd
    mul!(sp, B, Δp)
    axpby!(one(T), rp, -one(T), sp)
    rb = norm(sp)
    num = fτ - dot(bw.aτ, Δp) - dot(g, Δy)
    #
    # column loop: drive the Woodbury column to targets set by the consumer's dt = num/S. The lift injects
    # dt·r₂ into the primal row (must fit ptol - rb) and dt·r_d2 into the dual row (must fit ytol), where
    # r₂ = g - B·Δp2 and r_d2 = H·Δp2 - Bᵀ·Δy2 - c. Tightening the column moves Δp2,Δy2 → S → dt →
    # targets, so recompute S and re-check (cap 6, observed depth ≤ 2). δ and the CRAIG floor are fixed.
    #
    δ = dot(p, Qp) / τ^2 + κ / τ
    flr = eps(T) * (norm(c) + norm(g)) * κ
    for _ in 1:6
        dt = num / bw.S[]
        adt = abs(dt)
        ptgt = (ptol - rb) / adt
        ytgt = ytol / adt
        #
        # a target under the CRAIG floor cannot be reached by tightening the column — proceed and let the
        # 3-row refinement net see the (small) injected residual rather than thrashing CRAIG on it.
        #
        (ptgt < flr || ytgt < flr) && break
        #
        # r₂ = ‖g - B·Δp2‖ (primal) and r_d2 = ‖H·Δp2 - Bᵀ·Δy2 - c‖ (dual); mulkkt! yields both pieces.
        #
        mulkkt!(sd, sp, H, B, bw.Δp2, bw.Δy2)
        axpby!(one(T), g, -one(T), sp)
        r2 = norm(sp)
        axpy!(-one(T), c, sd)
        r_d2 = norm(sd)
        (r2 ≤ ptgt && r_d2 ≤ ytgt) && break
        #
        # tighten the column to (ptgt, ytgt) and recompute S from the moved column (bare inner product —
        # a tightened column barely cancels; also pre-figures the capacitance guard removal).
        #
        nb2, nr2, _, _, _, _, _ = solvekkt!(inner, bw.Δp2, bw.Δy2, H, B, c, g; ptol=ptgt, ytol=ytgt, stall, itmax)
        nbase += nb2 + nr2
        bw.S[] = dot(bw.aτ, bw.Δp2) + dot(g, bw.Δy2) + δ
    end
    #
    # lift (dt paired with the FINAL S — the τ-row exactness of Change 1 needs S formed from the same
    # Δp2,Δy2 the lift adds).
    #
    dt = num / bw.S[]
    Δτ = dt
    axpy!(dt, bw.Δp2, Δp)
    axpy!(dt, bw.Δy2, Δy)
    #
    # 3-row refinement safety net (should pass at pass 0 when both targets were met).
    #
    npass, nrefine, status, Δτ, pres1, dres1 = refinehsd!(
        Δp, Δy, Δτ, inner, H, B, c, g, Qp, p, bw.aτ, τ, κ, rp, f, fτ, bw.Δp2, bw.Δy2, bw.S[],
        inner.sp, inner.sd, inner.dp, inner.dy;
        itmax, ptol, ytol, τtol, stall,
    )
    return nbase, nrefine, npass, status, fres, pres1, dres1, Δτ
end

function capacitance!(
        QΔp2::AbstractVector{T},
        aτ::AbstractVector{T},
        τ::T,
        κ::T,
        Δp2::AbstractVector{T},
        Δy2::AbstractVector{T},
        c::AbstractVector{T},
        g::AbstractVector{T},
        Qp::AbstractVector{T},
        p::AbstractVector{T},
        W::BlockSparseMatrix{T},
        Q::BlockSparseMatrix{T},
    ) where {T}
    #
    # compute aτ = c - 2Qp/τ  (the border row; used by the Schur lift in newton!/refinehsd!)
    #
    copyto!(aτ, c)
    axpby!(-2 / τ, Qp, one(T), aτ)
    #
    # Woodbury capacitance by inner product (units spec, HSD bordered §1):
    #
    #   S = aτᵀΔp2 + gᵀΔy2 + δ,   δ = pᵀQp/τ² + κ/τ
    #
    # This is the Schur complement of the τ row and is ALGEBRAICALLY identical to the quadratic form
    # S = Δp2ᵀWΔp2 + (Δp2-p/τ)ᵀQ(Δp2-p/τ) + κ/τ (cross terms cancel via HΔp2-BᵀΔy2=c, BΔp2=g, W=H-Q).
    # Numerically it is BETTER: the same aτᵀΔp2/gᵀΔy2 rounding appears here and in the lift's numerator
    # num = fτ - aτᵀΔp0 - gᵀΔy0, so they cancel and the τ row is satisfied to machine precision regardless
    # of how loosely the Woodbury column was solved. It is also cheaper (two dot products, no matvecs).
    #
    δ = dot(p, Qp) / τ^2 + κ / τ
    S = dot(aτ, Δp2) + dot(g, Δy2) + δ
    #
    # guard: S is divided by (the Schur lift dτ = num/S). The true capacitance is ≥ κ/τ — the quadratic
    # form is Δp2ᵀWΔp2 + (Δp2-p/τ)ᵀQ(Δp2-p/τ) + κ/τ, a sum of PSD terms plus κ/τ (W, Q ⪰ 0). So an inner
    # product landing below that floor has been corrupted by cancellation; fall back to the manifestly
    # positive quadratic form. Above the floor the inner product is used, keeping the τ-row exact.
    #
    if S > κ / τ
        return S
    end

    mul!(QΔp2, Symmetric(Q, :L), Δp2)
    Sq = dot(Δp2, Symmetric(W, :L), Δp2) + κ / τ
    @inbounds for i in eachindex(Δp2)
        Sq += (Δp2[i] - p[i] / τ) * (QΔp2[i] - Qp[i] / τ)
    end

    return Sq
end

############################################################################################
# refinehsd! — Govaerts–Pryce BE+1 refinement on the 3-row bordered system
############################################################################################

function refinehsd!(
    Δp::AbstractVector{T},
    Δy::AbstractVector{T},
    Δτ::T,
    wrk::UzawaSolver{UPLO, T},
    H::BlockSparseMatrix{T},
    B::BlockSparseMatrix{T},
    c::AbstractVector{T},
    g::AbstractVector{T},
    Qp::AbstractVector{T},
    p::AbstractVector{T},
    aτ::AbstractVector{T}, # border row: c - 2Qp/τ
    τ::T,
    κ::T,
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
    status = REFINE_ITMAX
    prv = typemax(T)
    #
    # compute the sum
    #
    #   1/τ pᵀQp + κ
    #
    pQpτκ = dot(p, Qp) / τ + κ
    pres1 = T(NaN)
    dres1 = T(NaN)

    for i in 1:itmax
        #
        # compute the residuals
        #
        #   [ sd ]   [ f  ]   [ H  -Bᵀ  -c ] [ Δp ]
        #   [ sp ] = [ rp ] - [ B   0   -g ] [ Δy ]
        #   [ sτ ]   [ fτ ]   [ cᵀ  gᵀ   0 ] [ Δτ ]
        #
        sτ = fτ - mulhsd!(sd, sp, H, B, c, g, Δp, Δy, Δτ)
        axpby!(one(T), f,  -one(T), sd)
        axpby!(one(T), rp, -one(T), sp)
        #
        # correct sτ:
        #
        #   sτ ← sτ + 2pᵀQΔp/τ - (pᵀQp/τ² + κ/τ) Δτ
        #
        sτ += (2dot(Qp, Δp) - Δτ * pQpτκ) / τ

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
            status = REACHED
            break
        end

        if res > stall * prv
            status = REFINE_STALLED
            break
        end

        prv = res
        #
        # solve for dp and dy:
        #
        #   [ H -Bᵀ ] [ dp ] = [ sd ]
        #   [ B  0  ] [ dy ]   [ sp ]
        #
        n, _ = solveuzw!(wrk.divwrk, wrk.itrwrk, dp, dy, wrk.r, wrk.F, H, B, sd, sp, wrk.α[], atol; θ = T(INTERIOR_THETA))
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
