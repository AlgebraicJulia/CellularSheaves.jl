struct HSDSolver{T, I, V} <: AbstractSolver{T}
    Q::BlockSparseMatrix{T, I}
    H::BlockSparseMatrix{T, I}
    B::BlockSparseMatrix{T, I}
    f::FVector{T}
    g::FVector{T}
    p::FVector{T}
    d::FVector{T}
    y::FVector{T}
    K::FVector{V}
    scaling::HSDScaling{T}
    C::FPermutation{I}
    R::FPermutation{I}
    wrk::HSDWorkspace{T}
    caches::Caches{T, I}
    sched::ConeSchedule{T, I}
    kkt::BorderedSolver{T, :L, I}
    hist::HSDHistory{T}
    ν::Int
    settings::HSDSettings{T}
    ρ::FScalar{T}
    τ::FScalar{T}
    κ::FScalar{T}
    Δp2::FVector{T}
    Δy2::FVector{T}
    Δe2::FVector{T}
    nf::FScalar{T}
    ng::FScalar{T}
    sg::FScalar{T}     # ‖g‖ in original (unscaled) units — stopping-test primal denominator
    sf::FScalar{T}     # ‖f‖ in original (unscaled) units — stopping-test dual denominator
    nB::FScalar{T}      # ‖B‖ — fixed for the solver's lifetime; the cold-start augmentation anchor
    δ::FScalar{T}       # reciprocal augmentation 1/α; owned by setaug!
    timers::TimerOutput
end

function result(s::HSDSolver{T}, status::IPMStatus) where {T}
    pu = copy(s.p)
    du = copy(s.d)
    yu = copy(s.y)

    τ, κ = unscale!(pu, du, yu, s.τ[], s.κ[], s.scaling)   # interior unscale + user-frame (τ, κ)

    if status in (OPTIMAL, NEAR_OPTIMAL, STALLED, ITERATION_LIMIT, NUMERICAL_FAILURE)
        ldiv!(τ, pu)
        ldiv!(τ, du)
        ldiv!(τ, yu)
    elseif status in (PRIMAL_INFEASIBLE, NEAR_PRIMAL_INFEASIBLE)
        ldiv!(norm(yu), yu)
    elseif status in (DUAL_INFEASIBLE, NEAR_DUAL_INFEASIBLE)
        np = norm(pu)
        ldiv!(np, pu)
        ldiv!(np, du)
    end

    niter = 0
    nsolve = 0

    for row in s.hist
        niter += 1
        nsolve += row.piter + row.citer + row.witer
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
        mu = pres = dres = dobj = T(NaN); pobj = dot(s.g, yu) / scl.bscl[]
    elseif status in (DUAL_INFEASIBLE, NEAR_DUAL_INFEASIBLE)
        # Certificate value fᵀp at the returned ‖p‖ = 1 (symmetric border-scale undo; no
        # dual-infeasible example exercises this branch yet — see e14 for the primal analogue).
        mu = pres = dres = pobj = T(NaN); dobj = dot(s.f, pu) / scl.bscl[]
    elseif status == ILL_POSED || status == NEAR_ILL_POSED || isempty(s.hist)
        mu = pres = dres = pobj = dobj = T(NaN)
    else
        residuals!(s)   # refresh w.rp (m), w.rd (n), w.Qp = Q·p at the terminal iterate
        pres = norm(w.rp ./ scl.rscl) / τ / (1 + norm(s.g ./ scl.rscl))
        dres = norm(w.rd ./ scl.cscl) / τ / (1 + norm(s.f ./ scl.cscl))
        mu   = dot(s.p, s.d) / (s.ν * τ^2)
        pQp  = dot(s.p, w.Qp) / τ^2
        # f, g carry the embedding border scale bscl (update!); undo it on the linear terms,
        # as the infeasibility-certificate branches above already do.
        pobj = pQp / 2 - dot(s.f, s.p) / τ / scl.bscl[]
        dobj = dot(s.g, s.y) / τ / scl.bscl[] - pQp / 2
    end

    p = Vector{T}(undef, length(pu))
    d = Vector{T}(undef, length(du))
    y = Vector{T}(undef, length(yu))

    ldiv!(p, s.C, pu)
    ldiv!(d, s.C, du)
    ldiv!(y, s.R, yu)

    return HSDResult{T}(p, d, y, status, niter, nsolve, τ, κ, s.hist, s.timers,
                        mu, pres, dres, pobj, dobj)
end

############################################################################################
# residuals!
############################################################################################

#
# compute the negated residuals
#
#   [ rd ]   [ d  ]   [ Q   -Bᵀ  -f ] [ p ]
#   [ rp ] = [ 0  ] - [ B    0   -g ] [ y ]
#   [ rτ ]   [ κ  ]   [ fᵀ   gᵀ   0 ] [ τ ]
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
        f::AbstractVector,
        g::AbstractVector,
        Q::BlockSparseMatrix,
    ) where {T}
    #
    # compute the matrix-vector product
    #
    #   [ rd ]   [ Q   -Bᵀ  -f ] [ p ]
    #   [ rp ] = [ B    0   -g ] [ y ]
    #   [ rτ ]   [ fᵀ   gᵀ   0 ] [ τ ]
    #
    mulkkt!(rd, rp, Q, B, p, y)
    axpy!(-τ, f, rd)
    axpy!(-τ, g, rp)
    rτ = dot(f, p) + dot(g, y)
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
    mul!(w.Qp, s.Q, s.p)
    return residuals!(w.rd, w.rp, w.Qp, s.τ[], s.κ[], s.B, s.p, s.d, s.y, s.f, s.g, s.Q)
end

############################################################################################
# solvepredictor! / solvecorrector!
############################################################################################

#
# solve for the Mehrotra predictor direction
#
#   [  H          -Bᵀ             -f ] [ Δpa ]   [ rd - d ]
#   [  B           0              -g ] [ Δya ] = [ rp     ]
#   [ fᵀ - 2pᵀQ/τ  gᵀ  pᵀQp/τ² + κ/τ ] [ Δτa ]   [ rτ - κ ]
#
# and recover
#
#   Δκa = -κ (1 + Δτa/τ)
#
function solvepredictor!(s::HSDSolver{T}, gap::T; ptol::T, ytol::T, τtol::T) where {T}
    return solvepredictor!(
        s.wrk, s.kkt, s.Δp2, s.Δy2, s.Δe2, s.settings, s.H, s.B, s.Q, s.f, s.g, s.d,
        s.τ[], s.κ[], gap;
        ptol, ytol, τtol,
    )
end

function solvepredictor!(
        w::HSDWorkspace{T},
        kkt::BorderedSolver{T},
        Δp2::AbstractVector{T},  # Woodbury column primal (HSD-owned; threaded to solvekkt!) (n)
        Δy2::AbstractVector{T},  # Woodbury column dual, accreted                             (m)
        Δe2::AbstractVector{T},  # Woodbury column residual e - B Δp2                         (m)
        set::HSDSettings{T},
        H::BlockSparseMatrix{T},
        B::BlockSparseMatrix{T},
        Q::BlockSparseMatrix{T},
        f::AbstractVector{T},
        g::AbstractVector{T},
        d::AbstractVector{T},
        τ::T,
        κ::T,
        gap::T;
        ptol::T,
        ytol::T,
        τtol::T,
    ) where {T}
    axpby!(-1,  d, 0, w.f)
    axpby!( 1, w.rd, 1, w.f)
    #
    # solve for the directions Δpa, Δya, and Δτa (bordered base + 3-row refinement)
    #
    #   [  H          -Bᵀ             -f ] [ Δpa ]   [ rd - d ]
    #   [  B           0              -g ] [ Δya ] = [ rp     ]
    #   [ fᵀ - 2pᵀQ/τ  gᵀ  pᵀQp/τ² + κ/τ ] [ Δτa ]   [ rτ - κ ]
    #
    piter, pncol, ppass, pstat, Δτa, dmin, dmax = solvekkt!(
        kkt, w.Δpa, w.Δya, Δp2, Δy2, Δe2, H, B, f, g, w.f, w.rp, gap;
        gtol=ptol, ftol=ytol, ηtol=τtol, stall=set.refine_stall_tol, rfmax=set.refine_max_iter, cgmax=set.newton_max_iter,
    )
    #
    # recover Δda:
    #
    #   Δda ← Q Δpa - Δτa f - Bᵀ Δya - rd
    #
    mul!(w.Δda, Q, w.Δpa)
    axpy!(-Δτa, f, w.Δda)
    mul!(w.Δda, B', w.Δya, -one(T), one(T))
    axpy!(-one(T), w.rd, w.Δda)
    #
    # recover Δκa:
    #
    #   Δκa = -κ (1 + 1/τ Δτa)
    #
    Δκa = -κ * (τ + Δτa) / τ
    return piter, pncol, ppass, pstat, Δτa, Δκa, dmin, dmax
end

#
# solve for the Mehrotra combined direction
#
#   [ H           -Bᵀ               -f ] [ Δp ]   [ rd*                         ]
#   [ B            0                -g ] [ Δy ] = [ rp                          ]
#   [ fᵀ - 2pᵀQ/τ  gᵀ    pᵀQp/τ² + κ/τ ] [ Δτ ]   [ rτ - κ + (σμ - Δτa·Δκa) / τ ]
#
# where rd* is the corrected dual residual, and recover
#
#   Δκ = (σμ - τκ - Δτa·Δκa - κ·Δτ) / τ
#
function solvecorrector!(s::HSDSolver{T}, μ::T, gap::T, Δτa::T, Δκa::T; ptol::T, ytol::T, τtol::T) where {T}
    return solvecorrector!(
        s.wrk, s.kkt, s.Δp2, s.Δy2, s.Δe2, s.settings, s.H, s.B, s.Q, s.f, s.g, s.K, s.p, s.d,
        s.caches, s.sched, s.ν, s.τ[], s.κ[], s.δ[],
        μ, gap, Δτa, Δκa;
        ptol, ytol, τtol,
    )
end

function solvecorrector!(
        w::HSDWorkspace{T},
        kkt::BorderedSolver{T},
        Δp2::AbstractVector{T},  # Woodbury column primal (HSD-owned; threaded to solvekkt!) (n)
        Δy2::AbstractVector{T},  # Woodbury column dual, accreted                             (m)
        Δe2::AbstractVector{T},  # Woodbury column residual e - B Δp2                         (m)
        set::HSDSettings{T},
        H::BlockSparseMatrix{T},
        B::BlockSparseMatrix{T},
        Q::BlockSparseMatrix{T},
        f::AbstractVector{T},
        g::AbstractVector{T},
        K::AbstractVector,
        p::AbstractVector{T},
        d::AbstractVector{T},
        caches::Caches{T},
        sched::ConeSchedule{T},
        ν::Integer,
        τ::T,
        κ::T,
        δ::T,
        μ::T,
        gap::T,
        Δτa::T,
        Δκa::T;
        ptol::T,
        ytol::T,
        τtol::T,
    ) where {T}
    α = inv(δ)
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
    αa = maxsteps(sched, K, p, d, w.Δpa, w.Δda, caches, B, w.step; step_frac=one(T))

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
    initcorrector!(sched, K, w.f, caches, p, d, w.Δpa, w.Δda, σμ, B)

    axpy!(1, w.rd, w.f)
    #
    # warm-start the corrector from the predictor's column-free directions (the column dual uses the
    # physical value Δy2 + α Δe2, folding the transient tail):
    #
    #   Δp₀a = Δpa - Δτa Δp2,   Δy₀a = Δya - Δτa (Δy2 + α Δe2)
    #
    # (the 2-row base's answer differs from these only by the corrector's second-order RHS) — seed both
    # into the output buffers and let the base refine from there.
    #
    copyto!(w.Δp, w.Δpa); axpy!(-Δτa, Δp2, w.Δp)
    copyto!(w.Δy, w.Δya); axpy!(-Δτa, Δy2, w.Δy); axpy!(-Δτa * α, Δe2, w.Δy)
    #
    # solve for the directions Δp, Δy, and Δτ
    #
    #   [ H           -Bᵀ               -f ] [ Δp ]   [ rd*                         ]
    #   [ B            0                -g ] [ Δy ] = [ rp                          ]
    #   [ fᵀ - 2pᵀQ/τ  gᵀ    pᵀQp/τ² + κ/τ ] [ Δτ ]   [ rτ - κ + (σμ - Δτa·Δκa) / τ ]
    #
    citer, cncol, cpass, cstat, Δτ, _, _ = solvekkt!(
        kkt, w.Δp, w.Δy, Δp2, Δy2, Δe2, H, B, f, g, w.f, w.rp, fτ;
        xwarm=true, ywarm=true, gtol=ptol, ftol=ytol, ηtol=τtol, stall=set.refine_stall_tol, rfmax=set.refine_max_iter, cgmax=set.newton_max_iter,
    )
    #
    # recover Δd:
    #
    #   Δd ← Q Δp - Δτ f - Bᵀ Δy - rd
    #
    mul!(w.Δd, Q, w.Δp)
    axpy!(-Δτ, f, w.Δd)
    mul!(w.Δd, B', w.Δy, -one(T), one(T))
    axpy!(-one(T), w.rd, w.Δd)
    #
    # recover Δκ:
    #
    #   Δκ ← (fκ - κ Δτ) / τ
    #
    Δκ = (fκ - κ * Δτ) / τ
    return citer, cncol, cpass, cstat, Δτ, Δκ
end

############################################################################################
# identitypoint!
############################################################################################

function identitypoint!(p::AbstractVector{T}, d::AbstractVector{T}, B::BlockSparseMatrix{T}, K) where {T}
    for v in vtxs(B)
        r = colrange(B, v)
        identity!(view(p, r), K[v])
        identity!(view(d, r), K[v])
    end

    return p, d
end

############################################################################################
# constructor / reinit! / init
############################################################################################

function reinit!(s::HSDSolver{T}) where {T}
    identitypoint!(s.p, s.d, s.B, s.K)
    fill!(s.y, zero(T))

    for v in vtxs(s.B)
        initcache!(cache(s.caches, v, s.K[v]))
    end

    empty!(s.hist)

    s.τ[] = one(T)
    s.κ[] = one(T)
    fill!(s.Δp2, false)
    fill!(s.Δy2, false)
    fill!(s.Δe2, false)

    return s
end

function HSDSolver(prob::IPMProblem{T, I}, settings::HSDSettings{T}) where {T, I}
    n = size(prob.B, 2)
    m = size(prob.B, 1)
    ν = conedegree(prob.K, prob.B)

    weights, graph = weightedgraph(prob.B, prob.Q)
    D, C, S = symbolic(weights, graph; alg=settings.elim_alg)

    B = selectvtxs(prob.B, D.perm)
    Q = halfselectvtxs(halfselectvtxs(prob.Q, D.perm), D.perm)
    f = C * prob.f
    g = copy(prob.g)
    cones = tounion(prob.K, D.perm)

    scaling = HSDScaling{T}(n, m)

    if settings.scale_max_iter > 0
        equilibrate!(scaling, B, Q, f, g; itmax=settings.scale_max_iter)
        update!(scaling, f, g)   # Stage 2: equilibrate the f/g embedding border → τ ≈ 1
    end

    kkt = BorderedSolver(S, B)   # HSD solves the 3-row bordered system

    C = C * prob.C
    R = prob.R

    p = FVector{T}(undef, n)
    d = FVector{T}(undef, n)
    y = FVector{T}(undef, m)

    caches = Caches(cones, B)

    H = copy(Q)
    sched = ConeSchedule{T}(cones, B, nthreads())
    hsdwrk = HSDWorkspace{T}(m, n, nvtxs(B))
    hist = HSDHistory{T}()
    ρ = FScalar{T}(undef)
    τ = FScalar{T}(undef)
    κ = FScalar{T}(undef)
    nf = FScalar{T}(undef)
    ng = FScalar{T}(undef)
    sg = FScalar{T}(undef)
    sf = FScalar{T}(undef)
    nB = FScalar{T}(undef)
    δ = FScalar{T}(undef)

    nB[] = norm(B)
    nf[] = norm(f)
    ng[] = norm(g)
    sg[] = scalenorm(g, scaling.rscl)
    sf[] = scalenorm(f, scaling.cscl)
    ρ[] = initreg(nB[])
    Δp2 = FVector{T}(undef, n)
    Δy2 = FVector{T}(undef, m)
    Δe2 = FVector{T}(undef, m)

    solver = HSDSolver(Q, H, B, f, g, p, d, y, cones,
        scaling, C, R, hsdwrk, caches, sched, kkt,
        hist, ν, settings, ρ, τ, κ, Δp2, Δy2, Δe2, nf, ng, sg, sf, nB, δ, TimerOutput()
    )

    return reinit!(solver)
end

function HSDSolver(prob::IPMProblem{T}; kw...) where {T}
    return HSDSolver(prob, HSDSettings{T}(; kw...))
end

function CommonSolve.init(prob::IPMProblem{T}, settings::HSDSettings{T}) where {T}
    return HSDSolver(prob, settings)
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

function isdualinfeasible(w, τ, κ, B, p, f, nf, rtol, atol)
    flag = τ / κ < rtol

    if flag
        cp = dot(f, p)
        np = norm(p)
        flag = cp > atol * np * (1 + nf)

        if flag
            mul!(w.sy, B, p)
            flag = max(norm(w.sy), norm(w.Qp)) < rtol * abs(cp)
        end
    end

    return flag
end

function isdualinfeasible(s::HSDSolver, rtol, atol)
    return isdualinfeasible(s.wrk, s.τ[], s.κ[], s.B, s.p, s.f, s.nf[], rtol, atol)
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
# initkkt! (HSD) / step!
############################################################################################

function initkkt!(s::HSDSolver{T}) where {T}
    τ  = s.τ[]
    Qp = s.wrk.Qp
    d  = s.wrk.aτ
    copyto!(d, s.f); axpby!(-2 / τ, Qp, one(T), d)
    φ  = dot(s.p, Qp) / τ^2 + s.κ[] / τ
    flag, ρ, wbase = initkkt!(
        s.kkt, s.Δp2, s.Δy2, s.Δe2, s.H, s.B, s.f, s.g, d, φ;
        δ=s.δ[], rgmin=s.ρ[],
    )
    s.ρ[] = max(s.ρ[], ρ)
    return flag, ρ, wbase
end

function step!(s::HSDSolver{T}) where {T}
    status = CONTINUE

    piter = citer = witer = 0   # total CG per role (base + refinement; woodbury = column base + tightening)
    ppass = cpass = 0   # refinement passes per role
    pstat = cstat = KKT_SOLVED   # refinement exit status per role

    step = zero(T)
    dmin = dmax = T(NaN)   # per-solve min-cost δ-window (bordered solvekkt!)
    ρ = zero(T)      # ρ-shift actually applied this step (0 = none / no factorization); recorded below

    w = s.wrk
    τ = s.τ[]
    κ = s.κ[]
    #
    # compute negated residuals
    #
    #   [ rd ]   [ d ]   [ Q   -Bᵀ  -f       ] [ p ]
    #   [ rp ] = [ 0 ] - [ B    0   -g       ] [ y ]
    #   [ rτ ]   [ κ ]   [ fᵀ   gᵀ  -pᵀQp/τ² ] [ τ ]
    #
    gap = residuals!(s)
    #
    # compute the centrality parameter
    #
    #   μ = (pᵀd + τκ) / (ν + 1)
    #
    μ = mu(s)
    τu = s.scaling.bscl[] * τ
    pres = scalenorm(w.rp, s.scaling.rscl) / τu / (1 + s.sg[])
    dres = scalenorm(w.rd, s.scaling.cscl) / τu / (1 + s.sf[])

    pQp = dot(s.p, w.Qp)
    pobj = (pQp / (2 * τ^2) - dot(s.f, s.p) / τ) / s.scaling.bscl[]^2
    dobj = (dot(s.g, s.y) / τ - pQp / (2 * τ^2)) / s.scaling.bscl[]^2

    if isoptimal(s, pobj, dobj, pres, dres)
        status = OPTIMAL
    elseif isprimalinfeasible(s)
        status = PRIMAL_INFEASIBLE
    elseif isdualinfeasible(s)
        status = DUAL_INFEASIBLE
    elseif isillposed(s)
        status = ILL_POSED
    elseif !(μ > 0)
        if s.settings.verbose > 1
            @warn "Nonpositive μ."
        end

        status = nearstatus(s, NUMERICAL_FAILURE)
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

        if !flag
            if s.settings.verbose > 1
                @warn "Scaling failed."
            end

            status = nearstatus(s, NUMERICAL_FAILURE)
        else
            #
            # choose augmentation parameter δ
            #
            setaug!(s)
            #
            # compute the KKT-solve tolerances (the border w₂ solve inside initkkt! needs them)
            #
            #   force: min(θ μ/μ₁, ceil)
            #
            μ1 = isempty(s.hist.μ) ? μ : first(s.hist.μ)

            tol = min(FORCING_FRAC * μ / μ1, FORCING_CEIL)
            ptol = tol * (1 + s.ng[])
            ytol = tol * (1 + s.nf[])
            τtol = tol * (1 + s.ng[] + s.nf[])   # border-row target (bordered predictor/corrector only)
            #
            # factor F and solve/cache the Woodbury column w₂ + capacitance (the bordered do-once)
            #
            #   [ H  -Bᵀ ] [ Δp2 ]   [ f ]
            #   [ B   0  ] [ Δy2 ] = [ g ] ,   S = Δp2ᵀ W Δp2 + (Δp2 - p/τ)ᵀ Q (Δp2 - p/τ) + κ/τ
            #
            initok, ρ, wbase = @timeit s.timers "initkkt" initkkt!(s)
            if !initok
                if s.settings.verbose > 1
                    @warn "Failed to initialize KKT solver."
                end

                status = nearstatus(s, NUMERICAL_FAILURE)
            else
                #
                # solve for the Mehrota predictor direction
                #
                #   [  H          -Bᵀ             -f ] [ Δpa ]   [ rd - d ]
                #   [  B           0              -g ] [ Δya ] = [ rp     ]
                #   [ fᵀ - 2pᵀQ/τ  gᵀ  pᵀQp/τ² + κ/τ ] [ Δτa ]   [ rτ - κ ]
                #
                piter, pncol, ppass, pstat, Δτa, Δκa, dmin, dmax = @timeit s.timers "predictor" solvepredictor!(s, gap; ptol, ytol, τtol)

                for v in vtxs(s.B)
                    if s.K[v] isa CofreeCone
                        fill!(view(w.Δda, colrange(s.B, v)), zero(T))
                    end
                end
                #
                # solve for the Mehrotra combined direction
                #
                #   [  H          -Bᵀ             -f ] [ Δp ]   [ rd*                         ]
                #   [  B           0              -g ] [ Δy ] = [ rp                          ]
                #   [ fᵀ - 2pᵀQ/τ  gᵀ  pᵀQp/τ² + κ/τ ] [ Δτ ]   [ rτ - κ + (σμ - Δτa·Δκa) / τ ]
                #
                citer, cncol, cpass, cstat, Δτ, Δκ = @timeit s.timers "corrector" solvecorrector!(s, μ, gap, Δτa, Δκa; ptol, ytol, τtol)
                #
                # materialize the physical column dual into s.Δy2 (fold the transient tail once, now that
                # the column is final) so it seeds next iteration's base. s.Δp2/s.Δe2 are already physical.
                #
                #   Δy2 ← Δy2 + α Δe2
                #
                axpy!(inv(s.δ[]), s.Δe2, s.Δy2)
                #
                # The Woodbury column (s.Δp2/s.Δy2) is threaded in/out of initkkt!/predictor/corrector, so
                # the final column already lives in it — no copy-back. It seeds next iteration's column base
                # (the KKT solver is memoryless; initkkt! reseeds/wipes from it).
                #
                # Woodbury role: wbase (1, the initkkt! column base) + wrefn (predictor + corrector column
                # tightening) → witer. piter/citer report the DIRECTION solves only, deconfounded from it.
                #
                wrefn = pncol + cncol
                witer = wbase + wrefn

                for v in vtxs(s.B)
                    if s.K[v] isa CofreeCone
                        fill!(view(w.Δd, colrange(s.B, v)), zero(T))
                    end
                end

                if s.settings.verbose > 1 && (pstat !== KKT_SOLVED || cstat !== KKT_SOLVED)
                    @info "KKT solve above target tolerance" pstat cstat
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
                step = @timeit s.timers "maxsteps" maxsteps(s, w.Δp, w.Δd; step_frac=s.settings.step_frac)

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
                end
            end
        end
    end

    push!(s.hist, (; μ, step, pres, dres, pobj, dobj, δ=s.δ[], ρ, τ=s.τ[], κ=s.κ[],
        piter, ppass, pstat, citer, cpass, cstat, witer,
        dmin, dmax))

    return status
end

