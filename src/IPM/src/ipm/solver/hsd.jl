struct HSDSolver{T, I, W, C} <: AbstractSolver{T}
    Q::BlockSparseMatrix{T, I}
    H::BlockSparseMatrix{T, I}
    Hc::BlockSparseMatrix{T, I}  # cone-only Hessian (before Q added)
    B::BlockSparseMatrix{T, I}
    c::FVector{T}
    g::FVector{T}
    p::FVector{T}
    d::FVector{T}
    y::FVector{T}
    K::FVector{C}
    scaling::HSDScaling{T}
    P::FPermutation{I}
    wrk::HSDWorkspace{T}
    caches::Caches{T, I}
    conewrk::ConeWorkspace{T}
    kkt::W
    hist::HSDHistory{T}
    ν::Int
    settings::HSDSettings{T}
    ρ::FScalar{T}
    τ::FScalar{T}
    κ::FScalar{T}
    Δy0::FVector{T}
    nc::FScalar{T}
    ng::FScalar{T}
    α::FScalar{T}       # effective augmentation penalty; fixed at construction (raug); no controller
    timers::TimerOutput
end

function result(s::HSDSolver{T}, status::IPMStatus) where {T}
    p = Vector{T}(undef, length(s.p))
    d = Vector{T}(undef, length(s.d))
    y = Vector{T}(undef, length(s.y))

    ldiv!(p, s.P, s.p)
    ldiv!(d, s.P, s.d)
    copyto!(y, s.y)

    τ, κ = unscale!(p, d, y, s.τ[], s.κ[], s.scaling)   # interior unscale + user-frame (τ, κ)

    if status in (OPTIMAL, NEAR_OPTIMAL, STALLED, ITERATION_LIMIT, NUMERICAL_FAILURE)
        ldiv!(τ, p)
        ldiv!(τ, d)
        ldiv!(τ, y)
    elseif status in (PRIMAL_INFEASIBLE, NEAR_PRIMAL_INFEASIBLE)
        ldiv!(norm(y), y)
    elseif status in (DUAL_INFEASIBLE, NEAR_DUAL_INFEASIBLE)
        np = norm(p)
        ldiv!(np, p)
        ldiv!(np, d)
    end

    niter = 0
    nsolve = 0

    for row in s.hist
        niter += 1
        nsolve += row.pbase + row.prefn + row.cbase + row.crefn + row.wbase + row.wrefn
    end
    #
    # C5 (addendum): the answer at the returned point, in user frame. pres/dres are recomputed
    # from a fresh residual pass (the last step's w.rp/w.rd describe the pre-update iterate), then
    # the residual VECTOR is un-scaled per-element (÷rscl for rows, ÷cscl for cols) BEFORE norming
    # — the E/D equilibration reweighting destroys the mass distribution, so un-scaling the stored
    # embedding norm would be wrong. Objectives are bilinear (equilibration-invariant); mu is
    # evaluated directly (the τκ / ν+1→ν wedges are full-size at off-path floors).
    #
    w = s.wrk; scl = s.scaling
    if status in (PRIMAL_INFEASIBLE, NEAR_PRIMAL_INFEASIBLE)
        # Certificate value gᵀy at the returned ‖y‖ = 1, in the user frame. The HSD embedding
        # scales g (the border) by d₂ = scaling.bscl, so undo it here (as unscale! does for τ).
        mu = pres = dres = dobj = T(NaN); pobj = dot(s.g, y) / scl.bscl[]
    elseif status in (DUAL_INFEASIBLE, NEAR_DUAL_INFEASIBLE)
        # Certificate value cᵀp at the returned ‖p‖ = 1 (symmetric border-scale undo; no
        # dual-infeasible example exercises this branch yet — see e14 for the primal analogue).
        mu = pres = dres = pobj = T(NaN); dobj = dot(s.c, p) / scl.bscl[]
    elseif status == ILL_POSED || status == NEAR_ILL_POSED || isempty(s.hist)
        mu = pres = dres = pobj = dobj = T(NaN)
    else
        residuals!(s)   # refresh w.rp (m), w.rd (n), w.Qp = Q·p at the terminal iterate
        pres = norm(w.rp ./ scl.rscl) / τ / (1 + norm(s.g ./ scl.rscl))
        dres = norm(w.rd ./ scl.cscl) / τ / (1 + norm(s.c ./ scl.cscl))
        mu   = dot(s.p, s.d) / (s.ν * τ^2)
        pQp  = dot(s.p, w.Qp) / τ^2
        pobj = pQp / 2 - dot(s.c, s.p) / τ
        dobj = dot(s.g, s.y) / τ - pQp / 2
    end

    return HSDResult{T}(p, d, y, status, niter, nsolve, τ, κ, s.hist, s.timers,
                        mu, pres, dres, pobj, dobj)
end

############################################################################################
# mulhsd!
############################################################################################

#
# compute the matrix-vector product
#
#   [ u ]   [ H   -Bᵀ  -c ] [ x ]
#   [ v ] = [ B    0   -g ] [ y ]
#   [ w ]   [ cᵀ   gᵀ   0 ] [ z ]
#
function mulhsd!(
        u::AbstractVector{T},
        v::AbstractVector{T},
        H::BlockSparseMatrix{T},
        B::BlockSparseMatrix{T},
        c::AbstractVector{T},
        g::AbstractVector{T},
        x::AbstractVector{T},
        y::AbstractVector{T},
        z::T,
    ) where {T}
    #
    # compute the matrix-vector product
    #
    #   [ u ] ← [ H  -Bᵀ ] [ x ]
    #   [ v ]   [ B   0  ] [ y ]
    #
    mulkkt!(u, v, H, B, x, y)
    #
    # subtract cz and gz:
    #
    #   u ← u - cz
    #   v ← v - gz
    #
    axpy!(-z, c, u)
    axpy!(-z, g, v)
    #
    # compute w:
    #
    #   w ← cᵀx + gᵀy
    #
    return dot(c, x) + dot(g, y)
end

############################################################################################
# residuals!
############################################################################################

#
# compute the negated residuals
#
#   [ rd ]   [ d  ]   [ Q   -Bᵀ  -c ] [ p ]
#   [ rp ] = [ 0  ] - [ B    0   -g ] [ y ]
#   [ rτ ]   [ κ  ]   [ cᵀ   gᵀ   0 ] [ τ ]
#
# then correct rτ:
#
#   rτ ← rτ + pᵀQp/τ
#
function residuals!(
        rd::AbstractVector{T},
        rp::AbstractVector{T},
        Qp::AbstractVector{T},
        τ::T,
        κ::T,
        B::BlockSparseMatrix,
        p::AbstractVector,
        d::AbstractVector,
        y::AbstractVector,
        c::AbstractVector,
        g::AbstractVector,
        Q::BlockSparseMatrix,
    ) where {T}
    #
    # compute the matrix-vector product
    #
    #   [ rd ]   [ Q   -Bᵀ  -c ] [ p ]
    #   [ rp ] = [ B    0   -g ] [ y ]
    #   [ rτ ]   [ cᵀ   gᵀ   0 ] [ τ ]
    #
    rτ = mulhsd!(rd, rp, Q, B, c, g, p, y, τ)
    #
    # correct rd and rp
    #
    #   rd ← -rd + d
    #   rp ← -rp  
    #
    axpby!(one(T), d, -one(T), rd)
    lmul!(-one(T), rp)
    #
    # correct rτ:
    #
    #   gap ← pᵀQp/τ - rτ
    #
    return dot(p, Qp) / τ - rτ
end

function residuals!(s::HSDSolver)
    w = s.wrk
    mul!(w.Qp, Symmetric(s.Q, :L), s.p)
    return residuals!(w.rd, w.rp, w.Qp, s.τ[], s.κ[], s.B, s.p, s.d, s.y, s.c, s.g, s.Q)
end

############################################################################################
# woodbury! / capacitance! / newton!
############################################################################################

# solve for the Woodbury auxiliary directions
#
#   [ H  -Bᵀ ] [ Δp2 ]   [ c ]
#   [ B   0  ] [ Δy2 ] = [ g ]
#
function woodbury!(s::HSDSolver{T}; force_tol::T, floor_tol::T, y0 = nothing) where {T}
    return woodbury!(
        s.wrk, s.kkt, s.settings, s.H, s.B, s.c, s.g, s.nc[], s.ng[];
        force_tol, floor_tol, y0,
    )
end

function woodbury!(
        w::HSDWorkspace{T},
        kkt::KKTWorkspace{T},
        set::HSDSettings{T},
        H::BlockSparseMatrix{T},
        B::BlockSparseMatrix{T},
        c::AbstractVector{T},
        g::AbstractVector{T},
        nc::T,
        ng::T;
        force_tol::T,
        floor_tol::T,
        y0 = nothing,
    ) where {T}
    atol = max(force_tol, floor_tol)
    #
    # solve for the Woodbury auxiliary directions
    #
    #   [ H  -Bᵀ ] [ Δp2 ]   [ c ]
    #   [ B   0  ] [ Δy2 ] = [ g ]
    #
    wbase = solve_kkt!(kkt, w.Δp2, w.Δy2, H, B, c, g, y0; atol)
    #
    # use iterative refinement to improve tha
    # accuracy of the solutions Δp2, Δy2
    #
    wpass, wrefn, wstat = refinekkt!(
        w.Δp2, w.Δy2, kkt, H, B,
        c, g, w.sy, w.sp, w.dp, w.dy, nc, ng;
        itmax=set.refine_itmax, force_tol, floor_tol, stall=set.refine_stall,
    )

    return wbase, wrefn, wpass, wstat
end

function capacitance!(
        QΔp2::AbstractVector{T},
        aτ::AbstractVector{T},
        τ::T,
        κ::T,
        Δp2::AbstractVector{T},
        c::AbstractVector{T},
        Qp::AbstractVector{T},
        p::AbstractVector{T},
        W::BlockSparseMatrix{T},
        Q::BlockSparseMatrix{T},
    ) where {T}
    #
    # compute aτ = c - 2Qp/τ
    #
    copyto!(aτ, c)
    axpby!(-2 / τ, Qp, one(T), aτ)
    #
    # compute Q Δp2
    #
    mul!(QΔp2, Symmetric(Q, :L), Δp2)
    #
    # compute the Woodbury capacitance scalar
    #
    #   S = Δp2ᵀ W Δp2 + (Δp2 - p/τ)ᵀ Q (Δp2 - p/τ) + κ/τ
    #
    S = dot(Δp2, Symmetric(W, :L), Δp2) + κ / τ

    @inbounds for i in eachindex(Δp2)
        S += (Δp2[i] - p[i] / τ) * (QΔp2[i] - Qp[i] / τ)
    end

    return S
end

function capacitance!(s::HSDSolver)
    w = s.wrk
    return capacitance!(w.QΔp2, w.aτ, s.τ[], s.κ[], w.Δp2, s.c, w.Qp, s.p, s.Hc, s.Q)
end

############################################################################################
# refinehsd! — Govaerts–Pryce BE+1 refinement on the 3-row bordered system
############################################################################################

function refinehsd!(
    Δp::AbstractVector{T},
    Δy::AbstractVector{T},
    Δτ::T,
    wrk::KKTWorkspace{T},
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
    dy::AbstractVector{T},
    nc::T,
    ng::T;
    itmax::Int,
    force_tol::T,
    floor_tol::T,
    stall::T,
) where {T}
    niter = 0
    npass = 0
    status = REFINE_ITMAX
    prv = typemax(T)
    #
    # compute the sum
    #
    #   1/τ pᵀQp + κ
    #
    pQpτκ = dot(p, Qp) / τ + κ

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

        pres = norm(sp, Inf) / (1 + ng)
        dres = norm(sd, Inf) / (1 + nc)
        τres = abs(sτ) / (1 + nc + ng)

        res = max(pres, dres, τres)

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
        # solve for dp and dy:
        #
        #   [ H -Bᵀ ] [ dp ] = [ sd ]
        #   [ B  0  ] [ dy ]   [ sp ]
        #
        niter += solve_kkt!(wrk, dp, dy, H, B, sd, sp; atol=max(force_tol, floor_tol))
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

    return npass, niter, status, Δτ
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
        wrk::KKTWorkspace{T},
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
        atol::T = T(0.1),
    ) where {T}
    #
    # solve for Δp and Δy:
    #
    #   [ H -Bᵀ ] [ Δp ] = [ f  ]
    #   [ B  0  ] [ Δy ]   [ rp ]
    #
    niter = solve_kkt!(wrk, Δp, Δy, H, B, f, rp, y0; atol)
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

    return niter, Δτ
end

############################################################################################
# solvepredictor! / solvecorrector!
############################################################################################

#
# solve for the Mehrotra predictor direction
#
#   [  H          -Bᵀ             -c ] [ Δpa ]   [ rd - d ]
#   [  B           0              -g ] [ Δya ] = [ rp     ]
#   [ cᵀ - 2pᵀQ/τ  gᵀ  pᵀQp/τ² + κ/τ ] [ Δτa ]   [ rτ - κ ]
#
# and recover
#
#   Δκa = -κ (1 + Δτa/τ)
#
function solvepredictor!(s::HSDSolver{T}, gap::T, aτ::AbstractVector{T}, S::T; force_tol::T, floor_tol::T) where {T}
    return solvepredictor!(
        s.wrk, s.kkt, s.settings, s.H, s.B, s.Q, s.c, s.g, s.p, s.d,
        s.τ[], s.κ[], s.nc[], s.ng[], gap, aτ, S;
        force_tol, floor_tol,
    )
end

function solvepredictor!(
        w::HSDWorkspace{T},
        kkt::KKTWorkspace{T},
        set::HSDSettings{T},
        H::BlockSparseMatrix{T},
        B::BlockSparseMatrix{T},
        Q::BlockSparseMatrix{T},
        c::AbstractVector{T},
        g::AbstractVector{T},
        p::AbstractVector{T},
        d::AbstractVector{T},
        τ::T,
        κ::T,
        nc::T,
        ng::T,
        gap::T,
        aτ::AbstractVector{T},
        S::T;
        force_tol::T,
        floor_tol::T,
    ) where {T}
    atol = max(force_tol, floor_tol)

    axpby!(-1,  d, 0, w.f)
    axpby!( 1, w.rd, 1, w.f)
    #
    # solve for the directions Δpa, Δya, and Δτa
    #
    #   [  H          -Bᵀ             -c ] [ Δpa ]   [ rd - d ]
    #   [  B           0              -g ] [ Δya ] = [ rp     ]
    #   [ cᵀ - 2pᵀQ/τ  gᵀ  pᵀQp/τ² + κ/τ ] [ Δτa ]   [ rτ - κ ]
    #
    pbase, Δτa = newton!(
        w.Δpa, w.Δya,
        kkt, H, B, g,
        w.f, w.rp, gap, w.Δp2, w.Δy2, aτ, S;
        atol
    )
    #
    # use iterative refinement to improve
    # the solutions Δpa, Δya, and Δτa
    #
    ppass, prefn, pstat, Δτa = refinehsd!(
        w.Δpa, w.Δya, Δτa,
        kkt, H, B, c, g, w.Qp, p, aτ,
        τ, κ, w.rp, w.f, gap, w.Δp2, w.Δy2, S,
        w.sy, w.sp, w.dp, w.dy, nc, ng;
        itmax=set.refine_itmax, force_tol, floor_tol, stall=set.refine_stall
    )
    #
    # recover Δda:
    #
    #   Δda ← Q Δpa - Δτa c - Bᵀ Δya - rd
    #
    mul!(w.Δda, Symmetric(Q, :L), w.Δpa)
    axpy!(-Δτa, c, w.Δda)
    mul!(w.Δda, B', w.Δya, -one(T), one(T))
    axpy!(-one(T), w.rd, w.Δda)
    #
    # recover Δκa:
    #
    #   Δκa = -κ (1 + 1/τ Δτa)
    #
    Δκa = -κ * (τ + Δτa) / τ
    return pbase, prefn, ppass, pstat, Δτa, Δκa
end

#
# solve for the Mehrotra combined direction
#
#   [ H           -Bᵀ               -c ] [ Δp ]   [ rd*                         ]
#   [ B            0                -g ] [ Δy ] = [ rp                          ]
#   [ cᵀ - 2pᵀQ/τ  gᵀ    pᵀQp/τ² + κ/τ ] [ Δτ ]   [ rτ - κ + (σμ - Δτa·Δκa) / τ ]
#
# where rd* is the corrected dual residual, and recover
#
#   Δκ = (σμ - τκ - Δτa·Δκa - κ·Δτ) / τ
#
function solvecorrector!(s::HSDSolver{T}, μ::T, gap::T, Δτa::T, Δκa::T, aτ::AbstractVector{T}, S::T; force_tol::T, floor_tol::T) where {T}
    return solvecorrector!(
        s.wrk, s.kkt, s.settings, s.H, s.B, s.Q, s.c, s.g, s.K, s.p, s.d,
        s.caches, s.conewrk, s.ν, s.τ[], s.κ[], s.nc[], s.ng[],
        μ, gap, Δτa, Δκa, aτ, S;
        force_tol, floor_tol,
    )
end

function solvecorrector!(
        w::HSDWorkspace{T},
        kkt::KKTWorkspace{T},
        set::HSDSettings{T},
        H::BlockSparseMatrix{T},
        B::BlockSparseMatrix{T},
        Q::BlockSparseMatrix{T},
        c::AbstractVector{T},
        g::AbstractVector{T},
        K::AbstractVector,
        p::AbstractVector{T},
        d::AbstractVector{T},
        caches::Caches{T},
        conewrk::ConeWorkspace{T},
        ν::Integer,
        τ::T,
        κ::T,
        nc::T,
        ng::T,
        μ::T,
        gap::T,
        Δτa::T,
        Δκa::T,
        aτ::AbstractVector{T},
        S::T;
        force_tol::T,
        floor_tol::T,
    ) where {T}
    atol = max(force_tol, floor_tol)
    #
    # compute the largest step length αa ∈ (0, 1]
    # such that the perturbed iterates
    #
    #   p + αa Δpa ∈ K
    #   d + αa Δda ∈ K*
    #
    #   τ + αa Δτa ≥ 0
    #   κ + αa Δκa ≥ 0
    #
    # lie within their respective cones
    #
    αa = one(T)

    for v in vtxs(B)
        τpv, τdv = maxsteps(K[v], v, p, d, w.Δpa, w.Δda, caches, B, conewrk)
        αa = min(αa, τpv, τdv)
    end

    if Δτa < 0
        αa = min(αa, -τ / Δτa)
    end

    if Δκa < 0
        αa = min(αa, -κ / Δκa)
    end

    #
    # compute the affine-step complementarity
    #
    #   μa = (⟨p + αa Δpa, d + αa Δda⟩ + (τ + αa Δτa)(κ + αa Δκa)) / (ν + 1)
    #
    μa = (τ + αa * Δτa) * (κ + αa * Δκa)

    for j in cols(B)
        μa += (p[j] + αa * w.Δpa[j]) * (d[j] + αa * w.Δda[j])
    end

    μa /= (ν + 1)
    #
    # compute the centering parameter
    #
    #   σμ ← clamp(μa (μa / μ)², 0, μ)
    #
    σμ = clamp(μa * (μa / μ)^2, zero(T), μ)
    #
    # compute fκ and fτ:
    #
    #   fκ = σμ - τκ -       Δτa Δκa
    #   fτ = rτ -  κ + (σμ - Δτa Δκa) / τ
    #
    fκ = σμ - τ * κ - Δτa * Δκa
    fτ = gap + (σμ - Δτa * Δκa) / τ
    #
    # set f to the Mehrotra corrector term:
    #
    #   f ← -d + (σμ e - Δpa ∘ Δda) / p
    #
    # where e is the Jordan identity element e ∈ K.
    #
    for v in vtxs(B)
        initcorrector!(K[v], v, w.f, caches, p, d, w.Δpa, w.Δda, σμ, B, conewrk)
    end

    axpy!(1, w.rd, w.f)

    axpy!(-Δτa, w.Δy2, w.Δya)
    #
    # solve for the directions Δp, Δy, and Δτ
    #
    #   [ H           -Bᵀ               -c ] [ Δp ]   [ rd*                         ]
    #   [ B            0                -g ] [ Δy ] = [ rp                          ]
    #   [ cᵀ - 2pᵀQ/τ  gᵀ    pᵀQp/τ² + κ/τ ] [ Δτ ]   [ rτ - κ + (σμ - Δτa·Δκa) / τ ]
    #
    cbase, Δτ = newton!(
        w.Δp, w.Δy,
        kkt, H, B, g,
        w.f, w.rp, fτ, w.Δp2, w.Δy2, aτ, S, w.Δya;
        atol
    )
    #
    # use iterative refinement to improve
    # the solutions Δp, Δy, and Δτ
    #
    cpass, crefn, cstat, Δτ = refinehsd!(
        w.Δp, w.Δy, Δτ,
        kkt, H, B, c, g, w.Qp, p, aτ,
        τ, κ, w.rp, w.f, fτ, w.Δp2, w.Δy2, S,
        w.sy, w.sp, w.dp, w.dy, nc, ng;
        itmax=set.refine_itmax, force_tol, floor_tol, stall=set.refine_stall
    )
    #
    # recover Δd:
    #
    #   Δd ← Q Δp - Δτ c - Bᵀ Δy - rd
    #
    mul!(w.Δd, Symmetric(Q, :L), w.Δp)
    axpy!(-Δτ, c, w.Δd)
    mul!(w.Δd, B', w.Δy, -one(T), one(T))
    axpy!(-one(T), w.rd, w.Δd)
    #
    # recover Δκ:
    #
    #   Δκ ← (fκ - κ Δτ) / τ
    #
    Δκ = (fκ - κ * Δτ) / τ
    return cbase, crefn, cpass, cstat, Δτ, Δκ
end

############################################################################################
# startingpoint!
############################################################################################

function identitypoint(B::BlockSparseMatrix{T}, K) where {T}
    n = size(B, 2)
    m = size(B, 1)

    p = FVector{T}(undef, n)
    d = FVector{T}(undef, n)
    y = FVector{T}(undef, m)

    for v in vtxs(B)
        r = colrange(B, v)
        identity!(view(p, r), K[v])
        identity!(view(d, r), K[v])
    end

    fill!(y, zero(T))

    return p, d, y
end

############################################################################################
# init
############################################################################################

function CommonSolve.init(prob::IPMProblem{T, I}, settings::HSDSettings{T}) where {T, I}
    n = size(prob.B, 2)
    m = size(prob.B, 1)
    ν = conedegree(prob.K, prob.B)

    scaling = HSDScaling{T}(n, m)

    if settings.scale_itmax > 0
        B = copy(prob.B)
        Q = copy(prob.Q)
        c = copy(prob.c)
        g = copy(prob.g)

        equilibrate!(scaling, B, Q, c, g; itmax=settings.scale_itmax)
        update!(scaling, c, g)   # Stage 2: equilibrate the c/g embedding border → τ ≈ 1
    else
        B = prob.B
        Q = prob.Q
        c = prob.c
        g = prob.g
    end

    R, P, B, kkt = make_kkt(B; elim=settings.elim, rgmin=settings.rgmin, rgmax=settings.rgmax)

    c = P * c
    Q = halfselectvtxs(halfselectvtxs(Q, R.perm), R.perm)
    cones = tounion(prob.K, R.perm)

    p, d, y = identitypoint(B, cones)

    caches = Caches(cones, B)

    for v in vtxs(B)
        initcache!(cache(caches, v, cones[v]))
    end

    H = allocblockdiag(B)
    Hc = allocblockdiag(B)
    conewrk = ConeWorkspace{T}(cones, B)
    hsdwrk = HSDWorkspace{T}(m, n)
    hist = HSDHistory{T}()
    ρ = FScalar{T}(undef)
    τ = FScalar{T}(undef)
    κ = FScalar{T}(undef)
    nc = FScalar{T}(undef)
    ng = FScalar{T}(undef)
    α = FScalar{T}(undef)

    ρ[] = settings.rgmin
    nc[] = norm(c)
    ng[] = norm(g)
    Δy0 = FVector{T}(undef, m); fill!(Δy0, false)

    τ[] = one(T)
    κ[] = one(T)

    nH = sqrt(n)
    nQ = norm(Symmetric(Q, :L))
    nB = norm(B)
    α[] = settings.aaug + settings.raug * (nH + nQ) / nB^2

    solver = HSDSolver(Q, H, Hc, B, c, g, p, d, y, cones,
        scaling, P, hsdwrk, caches, conewrk, kkt,
        hist, ν, settings, ρ, τ, κ, Δy0, nc, ng, α, TimerOutput()
    )

    return solver
end

############################################################################################
# mu / isoptimal / infeasibility / ill-posed
############################################################################################
#
# compute the centrality parameter
#
#   μ = (pᵀd + τκ) / (ν + 1)
#
function mu(s::HSDSolver{T}) where {T}
    pd = dot(s.p, s.d)
    return (pd + s.τ[] * s.κ[]) / (s.ν + 1)
end

function isoptimal(s::HSDSolver{T}, μ::T, pres::T, dres::T) where {T}
    τ = s.τ[]
    w = s.wrk

    pQp = dot(s.p, w.Qp) / τ^2
    pobj = pQp / 2 - dot(s.c, s.p) / τ
    dobj = dot(s.g, s.y) / τ - pQp / 2

    max(pres, dres) < s.settings.feas_tol && (μ < s.settings.gap_tol * τ^2 || pobj - dobj < s.settings.gap_tol * (1 + abs(pobj) + abs(dobj)))
end

function isprimalinfeasible(w, τ, κ, B, d, g, y, ng, rtol, atol)
    flag = τ / κ < rtol

    if flag
        gy = dot(g, y)
        ny = norm(y)
        flag = gy > atol * ny * (1 + ng)

        if flag
            nQp = norm(w.Qp)
            copyto!(w.f, w.Qp)
            axpy!(-1, d, w.f)
            mul!(w.f, B', y, -1, 1)
            flag = max(nQp, norm(w.f)) < rtol * gy * (1 + norm(d) / ny)
        end
    end

    return flag
end

function isprimalinfeasible(s::HSDSolver, rtol, atol)
    return isprimalinfeasible(s.wrk, s.τ[], s.κ[], s.B, s.d, s.g, s.y, s.ng[], rtol, atol)
end

function isprimalinfeasible(s::HSDSolver)
    return isprimalinfeasible(s, s.settings.infeas_rel, s.settings.infeas_abs)
end

function isnearprimalinfeasible(s::HSDSolver)
    f = s.settings.near_factor
    return isprimalinfeasible(s, f * s.settings.infeas_rel, s.settings.infeas_abs)
end

function isdualinfeasible(w, τ, κ, B, p, c, nc, rtol, atol)
    flag = τ / κ < rtol

    if flag
        cp = dot(c, p)
        np = norm(p)
        flag = cp > atol * np * (1 + nc)

        if flag
            mul!(w.sy, B, p)
            flag = max(norm(w.sy), norm(w.Qp)) < rtol * abs(cp)
        end
    end

    return flag
end

function isdualinfeasible(s::HSDSolver, rtol, atol)
    return isdualinfeasible(s.wrk, s.τ[], s.κ[], s.B, s.p, s.c, s.nc[], rtol, atol)
end

function isdualinfeasible(s::HSDSolver)
    return isdualinfeasible(s, s.settings.infeas_rel, s.settings.infeas_abs)
end

function isneardualinfeasible(s::HSDSolver)
    f = s.settings.near_factor
    return isdualinfeasible(s, f * s.settings.infeas_rel, s.settings.infeas_abs)
end

function isillposed(s::HSDSolver)
    return isillposed(s.hist; tol=s.settings.illposed_tol)
end

function isnearillposed(s::HSDSolver)
    tol = s.settings.near_factor * s.settings.illposed_tol
    return max(s.τ[], s.κ[]) <= tol
end

function nearstatus(s::HSDSolver, status::IPMStatus)
    residuals!(s)

    if isnearoptimal(s)
        status = NEAR_OPTIMAL
    elseif isnearprimalinfeasible(s)
        status = NEAR_PRIMAL_INFEASIBLE
    elseif isneardualinfeasible(s)
        status = NEAR_DUAL_INFEASIBLE
    elseif isnearillposed(s)
        status = NEAR_ILL_POSED
    end

    return status
end

############################################################################################
# step!
############################################################################################

function step!(s::HSDSolver{T}) where {T}
    status = CONTINUE

    pbase = cbase = wbase = 0   # base-solve CRAIG per role (woodbury counted into craig1, archive-wins)
    prefn = crefn = wrefn = 0   # refinement CRAIG per role
    ppass = cpass = wpass = 0   # refinement passes per role
    pstat = cstat = wstat = REACHED_FORCE   # refinement exit status per role

    step = zero(T)

    w = s.wrk
    τ = s.τ[]
    κ = s.κ[]
    #
    # compute negated residuals
    #
    #   [ rd ]   [ d ]   [ Q   -Bᵀ  -c       ] [ p ]
    #   [ rp ] = [ 0 ] - [ B    0   -g       ] [ y ]
    #   [ rτ ]   [ κ ]   [ cᵀ   gᵀ  -pᵀQp/τ² ] [ τ ]
    #
    gap = residuals!(s)
    #
    # compute the centrality parameter
    #
    #   μ = (pᵀd + τκ) / (ν + 1)
    #
    μ = mu(s)
    pres = norm(w.rp) / τ / (1 + s.ng[])
    dres = norm(w.rd) / τ / (1 + s.nc[])

    if isoptimal(s, μ, pres, dres)
        status = OPTIMAL
    elseif isprimalinfeasible(s)
        status = PRIMAL_INFEASIBLE
    elseif isdualinfeasible(s)
        status = DUAL_INFEASIBLE
    elseif isillposed(s)
        status = ILL_POSED
    else
        #
        # compute the Hessian
        #
        #   f''(w)
        #
        # of the primal barrier function f
        # at the Nestorov-Todd scaling point w
        #
        # for non-symmetric cones, no such point
        # exists, so the Hessian is replaced
        # by a Tuncel scaling matrix
        #
        @timeit s.timers "scale" flag = scale!(s)
        #
        # cache H and add the quadratic term
        #
        #   H ← H + Q
        #
        copyto!(s.Hc, s.H)

        for v in vtxs(s.B)
            axpy!(true, block(s.Q, v, v, v), block(s.H, v, v, v))
        end

        if !flag
            if s.settings.verbose > 1
                @warn "Scaling failed."
            end

            status = nearstatus(s, NUMERICAL_FAILURE)
        else
            if !@timeit s.timers "initkkt" initkkt!(s)
                if s.settings.verbose > 1
                    @warn "Failed to initialize KKT solver."
                end

                status = nearstatus(s, NUMERICAL_FAILURE)
            else
                #
                # compute tolerances for predictor and corrector solves
                #
                #   force: min(θ μ/μ₁, ceil)
                #   floor: 100ϵ (1 + max(‖rp‖, ‖rd‖))
                #
                if isempty(s.hist.μ)
                    μ1 = μ
                else
                    μ1 = first(s.hist.μ)
                end

                force_tol = min(s.settings.forcing_frac * μ / μ1, s.settings.forcing_ceil)
                floor_tol = 100eps(T) * (1 + max(norm(w.rp, Inf), norm(w.rd, Inf)))
                #
                # solve for the Woodbury auxiliary directions
                #
                #   [ H  -Bᵀ ] [ Δp2 ]   [ c ]
                #   [ B   0  ] [ Δy2 ] = [ g ]
                #
                wbase, wrefn, wpass, wstat = @timeit s.timers "woodbury" woodbury!(s; force_tol, floor_tol, y0 = s.Δy0)
                copyto!(s.Δy0, w.Δy2)
                #
                # compute the Woodbury capacitance scalar
                #
                #   S = Δp2ᵀ W Δp2 + (Δp2 - p/τ)ᵀ Q (Δp2 - p/τ) + κ/τ
                #
                S = capacitance!(s)
                #
                # solve for the Mehrota predictor direction
                #
                #   [  H          -Bᵀ             -c ] [ Δpa ]   [ rd - d ]
                #   [  B           0              -g ] [ Δya ] = [ rp     ]
                #   [ cᵀ - 2pᵀQ/τ  gᵀ  pᵀQp/τ² + κ/τ ] [ Δτa ]   [ rτ - κ ]
                #
                # and recover
                #
                #   Δκa = (-τκ - κ Δτa) / τ
                #
                pbase, prefn, ppass, pstat, Δτa, Δκa = @timeit s.timers "predictor" solvepredictor!(s, gap, w.aτ, S; force_tol, floor_tol)

                for v in vtxs(s.B)
                    if s.K[v] isa CofreeCone
                        fill!(view(w.Δda, colrange(s.B, v)), zero(T))
                    end
                end
                #
                # solve for the Mehrotra combined direction
                #
                #   [  H          -Bᵀ             -c ] [ Δp ]   [ rd*                         ]
                #   [  B           0              -g ] [ Δy ] = [ rp                          ]
                #   [ cᵀ - 2pᵀQ/τ  gᵀ  pᵀQp/τ² + κ/τ ] [ Δτ ]   [ rτ - κ + (σμ - Δτa·Δκa) / τ ]
                #
                # where rd* is the corrected dual residual, and recover
                #
                #   Δκ = (σμ - τκ - Δτa·Δκa - κ·Δτ) / τ
                #
                cbase, crefn, cpass, cstat, Δτ, Δκ = @timeit s.timers "corrector" solvecorrector!(s, μ, gap, Δτa, Δκa, w.aτ, S; force_tol, floor_tol)

                for v in vtxs(s.B)
                    if s.K[v] isa CofreeCone
                        fill!(view(w.Δd, colrange(s.B, v)), zero(T))
                    end
                end

                if s.settings.verbose > 1 && (pstat != REACHED_FORCE || cstat != REACHED_FORCE || wstat != REACHED_FORCE)
                    @info "KKT solve above target tolerance" pstat cstat wstat
                end
                #
                # find the largest step ∈ (0, 1] such that
                #
                #   p + step Δp ∈ K
                #   d + step Δd ∈ K*
                #
                #   τ + step Δτ ≥ 0
                #   κ + step Δκ ≥ 0
                #
                @timeit s.timers "maxsteps" begin
                    τp, τd = maxsteps(s, w.Δp, w.Δd; step_frac=s.settings.step_frac)
                    step = min(τp, τd)
                end

                if Δτ < 0
                    step = min(step, s.settings.step_frac * (-τ / Δτ))
                end

                if Δκ < 0
                    step = min(step, s.settings.step_frac * (-κ / Δκ))
                end
                #
                # compute the updated iterates
                #
                #   p ← p + step Δp ∈ K
                #   d ← d + step Δd ∈ K*
                #
                #   τ ← step Δτ ≥ 0
                #   κ ← step Δκ ≥ 0
                #
                axpy!(step, w.Δp, s.p)
                axpy!(step, w.Δd, s.d)
                axpy!(step, w.Δy, s.y)

                s.τ[] = τ + step * Δτ
                s.κ[] = κ + step * Δκ

                if isstalled(s)
                    if s.settings.verbose > 1
                        @warn "Stalling detected."
                    end

                    status = nearstatus(s, STALLED)
                elseif isnumfail(s)
                    if s.settings.verbose > 1
                        @warn "Step collapse detected."
                    end

                    status = nearstatus(s, NUMERICAL_FAILURE)
                end
            end
        end
    end

    push!(s.hist, (; μ, step, pres, dres, gap, α=s.α[], ρ=s.ρ[], τ=s.τ[], κ=s.κ[],
        pbase, prefn, ppass, pstat, cbase, crefn, cpass, cstat, wbase, wrefn, wpass, wstat))
    if status == CONTINUE && atfloor(s.hist; patience=s.settings.floor_patience)
        if s.settings.verbose > 1
            @warn "Refinement floor reached $(s.settings.floor_patience) consecutive times"
        end

        status = nearstatus(s, NUMERICAL_FAILURE)
    end

    return status
end

function reinit!(solver::HSDSolver{T}, prob::IPMProblem{T}; frac::Real=0.1) where {T}
    @assert ncols(prob.B) == ncols(solver.B)
    @assert nrows(prob.B) == nrows(solver.B)
    @assert nvtxs(prob.B) == nvtxs(solver.B)
    @assert nouts(prob.B) == nouts(solver.B)

    n = ncols(solver.B)

    c = copy(prob.c)
    g = copy(prob.g)

    for j in cols(solver.B)
        c[j] *= solver.scaling.cscl[j]
    end

    for i in rows(solver.B)
        g[i] *= solver.scaling.rscl[i]
    end

    update!(solver.scaling, c, g)   # Stage 2: refresh the border scale for the new c, g

    mul!(solver.c, solver.P, c)
    copyto!(solver.g, g)

    solver.nc[] = norm(solver.c)
    solver.ng[] = norm(solver.g)

    τ = solver.τ[]
    κ = solver.κ[]

    if τ > sqrt(eps(T)) * κ
        f = convert(T, frac)
    else
        f = one(T)
    end

    if f < 1
        ldiv!(τ, solver.p)
        ldiv!(τ, solver.d)
        ldiv!(τ, solver.y)
    end

    p, d, y = identitypoint(solver.B, solver.K)
    axpby!(f, p, 1 - f, solver.p)
    axpby!(f, d, 1 - f, solver.d)
    axpby!(f, y, 1 - f, solver.y)

    solver.τ[] = one(T)
    solver.κ[] = dot(solver.p, solver.d) / solver.ν

    empty!(solver.hist)

    for v in vtxs(solver.B)
        initcache!(cache(solver.caches, v, solver.K[v]))
    end

    solver.ρ[] = solver.settings.rgmin

    nH = sqrt(T(n))
    nQ = norm(Symmetric(solver.Q, :L))
    nB = norm(solver.B)
    solver.α[] = solver.settings.aaug + solver.settings.raug * (nH + nQ) / nB^2

    return solver
end
