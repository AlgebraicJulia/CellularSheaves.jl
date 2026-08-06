# Bordered KKT solver — wraps a 2-row KKTSolver and adds the HSD τ-border via a Woodbury/Schur update.
# The Woodbury column w₂ = KKT⁻¹[c; g] and the capacitance scalar S are cached: they are fixed for a
# given factorization and reused by the predictor and corrector solves at that factorization. The border
# data (c, g, Q, Hc) is passed per call, not stored. initkkt! factors + solves/caches w₂ + capacitance;
# solvekkt! does the bordered base solve (newton!) + the 3-row refinement (refinehsd!).
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
# woodbury! / capacitance! / newton!
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
    # solve & cache the Woodbury column   [ H -Bᵀ; B 0 ] [Δp2; Δy2] = [c; g]  to the working tolerance.
    # (The units spec's §3 solved it to floor_tol — the theory being that its residual r₂ = g − B·Δp2 is
    # re-injected as dτ·r₂ by every Schur lift, and the capacitance bw.S[] is computed from Δp2. But §2
    # showed dτ·r₂ is NOT the refinement bottleneck here — npass didn't rise — and forcing the column to
    # ~100eps burned 22–28 CRAIG on it, drove Δp2/S to the noise floor, and broke e04-hsd to
    # NUMERICAL_FAILURE. So it stays at working tol.) NB the RHS here IS (c, g) with nc/ng = ‖c‖/‖g‖, so
    # (1+nc)/(1+ng) is a genuine relative-residual test, not fixed problem-data constants — don't "fix" it.
    #
    wtuple = solvekkt!(bw.inner, bw.Δp2, bw.Δy2, H, B, c, g, y0;
                       ptol, ytol, stall, itmax)
    #
    # cache the capacitance scalar S and border row aτ = c - 2Qp/τ (stable split form)
    #
    bw.S[] = capacitance!(bw.QΔp2, bw.aτ, τ, κ, bw.Δp2, bw.Δy2, c, g, Qp, p, Hc, Q)
    return ok, ρ, wtuple
end

# BorderedSolver solve: the bordered guarantee — newton! (2-row base + Schur lift, reusing the cached
# Woodbury column and capacitance) then refinehsd! (3-row iterative refinement). Refinement scratch is
# the inner solver's. Returns (nbase, nrefine, npass, status, fres, entry_pres, entry_dres, Δτ).
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
    # craig units: the base solve's primal tolerance IS ptol — one place. newton! forwards this straight
    # to solveuzw!; do not re-scale inside it. First-order in the bordered path (the Schur lift injects
    # dτ·r₂ each pass — see refinehsd! and §3).
    atol = ptol
    nbase, Δτ, fres = newton!(Δp, Δy, bw.inner, H, B, g, f, rp, fτ, bw.Δp2, bw.Δy2, bw.aτ, bw.S[], y0; atol)
    npass, nrefine, status, Δτ, pres1, dres1 = refinehsd!(
        Δp, Δy, Δτ, bw.inner, H, B, c, g, Qp, p, bw.aτ, τ, κ, rp, f, fτ, bw.Δp2, bw.Δy2, bw.S[],
        bw.inner.sp, bw.inner.sd, bw.inner.dp, bw.inner.dy;
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

#
# solve for the directions Δp, Δy, and Δτ
#
#   [  H          -Bᵀ             -c ] [ Δp ]   [ f  ]
#   [  B           0              -g ] [ Δy ] = [ rp ]
#   [ cᵀ - 2pᵀQ/τ  gᵀ  pᵀQp/τ² + κ/τ ] [ Δτ ]   [ fτ ]
#
function newton!(
        Δp::AbstractVector{T},
        Δy::AbstractVector{T},
        wrk::UzawaSolver{UPLO, T},
        H::BlockSparseMatrix{T},
        B::BlockSparseMatrix{T},
        g::AbstractVector{T},
        f::AbstractVector{T},
        rp::AbstractVector{T},
        fτ::T,
        Δp2::AbstractVector{T},
        Δy2::AbstractVector{T},
        aτ::AbstractVector{T},
        S::T,
        y0 = nothing;
        atol::T,
    ) where {UPLO, T}
    #
    # base solve for Δp and Δy (the 2-row primitive; the 3-row refinement is refinehsd!):
    #
    #   [ H -Bᵀ ] [ Δp ] = [ f  ]
    #   [ B  0  ] [ Δy ]   [ rp ]
    #
    niter, nr0 = solveuzw!(wrk.divwrk, wrk.itrwrk, Δp, Δy, wrk.r, wrk.F, H, B, f, rp, wrk.α[], atol, y0)
    #
    # apply the Schur lift:
    #
    #   Δτ = (fτ - aτᵀΔp - gᵀΔy) / S
    #
    Δτ = (fτ - dot(aτ, Δp) - dot(g, Δy)) / S
    #
    # apply Woodbury update:
    #
    #   Δp ← Δp + Δτ Δp2
    #   Δy ← Δy + Δτ Δy2
    #
    axpy!(Δτ, Δp2, Δp)
    axpy!(Δτ, Δy2, Δy)

    return niter, Δτ, nr0
end
