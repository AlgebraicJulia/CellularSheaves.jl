# α-window controller constants (alpha_window_laws.md §7): floor budget, ceiling gap, controller cap.
# budget_hsd=0 (universal minus the hidden ~0.3 woodbury workload); cap=1.5 = 1 + 0.5 near-edge bias.
const CTRL_BDG_HSD =  0.0
const CTRL_GAP_HSD = -0.4
const CTRL_CAP_HSD =  1.5

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
    nB::FScalar{T}      # ‖B‖ — fixed for the solver's lifetime; the cold-start augmentation anchor
    α::FScalar{T}       # effective augmentation penalty; owned by setaug!
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
    mulkkt!(rd, rp, Q, B, p, y)
    axpy!(-τ, c, rd)
    axpy!(-τ, g, rp)
    rτ = dot(c, p) + dot(g, y)
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
function solvepredictor!(s::HSDSolver{T}, gap::T; ptol::T, ytol::T, τtol::T) where {T}
    return solvepredictor!(
        s.wrk, s.kkt, s.settings, s.H, s.B, s.Q, s.c, s.g, s.ng[], s.d,
        s.τ[], s.κ[], gap;
        ptol, ytol, τtol,
    )
end

function solvepredictor!(
        w::HSDWorkspace{T},
        kkt::BorderedSolver{T},
        set::HSDSettings{T},
        H::BlockSparseMatrix{T},
        B::BlockSparseMatrix{T},
        Q::BlockSparseMatrix{T},
        c::AbstractVector{T},
        g::AbstractVector{T},
        ng::T,
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
    #   [  H          -Bᵀ             -c ] [ Δpa ]   [ rd - d ]
    #   [  B           0              -g ] [ Δya ] = [ rp     ]
    #   [ cᵀ - 2pᵀQ/τ  gᵀ  pᵀQp/τ² + κ/τ ] [ Δτa ]   [ rτ - κ ]
    #
    pbase, prefn, ppass, pstat, pfres, ppres, pdres, Δτa = solvekkt!(
        kkt, w.Δpa, w.Δya, H, B, c, g, ng, w.f, w.rp, gap;
        ptol, ytol, τtol, stall=set.refine_stall, itmax=set.refine_itmax,
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
    return pbase, prefn, ppass, pstat, Δτa, Δκa, pfres, ppres, pdres
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
function solvecorrector!(s::HSDSolver{T}, μ::T, gap::T, Δτa::T, Δκa::T; ptol::T, ytol::T, τtol::T) where {T}
    return solvecorrector!(
        s.wrk, s.kkt, s.settings, s.H, s.B, s.Q, s.c, s.g, s.ng[], s.K, s.p, s.d,
        s.caches, s.conewrk, s.ν, s.τ[], s.κ[],
        μ, gap, Δτa, Δκa;
        ptol, ytol, τtol,
    )
end

function solvecorrector!(
        w::HSDWorkspace{T},
        kkt::BorderedSolver{T},
        set::HSDSettings{T},
        H::BlockSparseMatrix{T},
        B::BlockSparseMatrix{T},
        Q::BlockSparseMatrix{T},
        c::AbstractVector{T},
        g::AbstractVector{T},
        ng::T,
        K::AbstractVector,
        p::AbstractVector{T},
        d::AbstractVector{T},
        caches::Caches{T},
        conewrk::ConeWorkspace{T},
        ν::Integer,
        τ::T,
        κ::T,
        μ::T,
        gap::T,
        Δτa::T,
        Δκa::T;
        ptol::T,
        ytol::T,
        τtol::T,
    ) where {T}
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
    #
    # warm-start the corrector from the predictor's column-free directions
    #
    #   Δp₀a = Δpa - Δτa Δp2,   Δy₀a = Δya - Δτa Δy2
    #
    # (the 2-row base's answer differs from these only by the corrector's second-order RHS) — seed both
    # into the output buffers and let the base refine from there.
    #
    copyto!(w.Δp, w.Δpa); axpy!(-Δτa, kkt.Δp2, w.Δp)
    copyto!(w.Δy, w.Δya); axpy!(-Δτa, kkt.Δy2, w.Δy)
    #
    # solve for the directions Δp, Δy, and Δτ
    #
    #   [ H           -Bᵀ               -c ] [ Δp ]   [ rd*                         ]
    #   [ B            0                -g ] [ Δy ] = [ rp                          ]
    #   [ cᵀ - 2pᵀQ/τ  gᵀ    pᵀQp/τ² + κ/τ ] [ Δτ ]   [ rτ - κ + (σμ - Δτa·Δκa) / τ ]
    #
    cbase, crefn, cpass, cstat, cfres, cpres, cdres, Δτ = solvekkt!(
        kkt, w.Δp, w.Δy, H, B, c, g, ng, w.f, w.rp, fτ;
        pwarm=true, ywarm=true, ptol, ytol, τtol, stall=set.refine_stall, itmax=set.refine_itmax,
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
    return cbase, crefn, cpass, cstat, Δτ, Δκ, cfres, cpres, cdres
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

    R, P, B, kkt = makekkt(B; elim=settings.elim)
    kkt = BorderedSolver(kkt, B)   # HSD solves the 3-row bordered system

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
    nB = FScalar{T}(undef)
    α = FScalar{T}(undef)

    ρ[] = resolve_rgmin(B)
    nc[] = norm(c)
    ng[] = norm(g)
    nB[] = norm(B)
    Δy0 = FVector{T}(undef, m); fill!(Δy0, false)

    τ[] = one(T)
    κ[] = one(T)

    solver = HSDSolver(Q, H, Hc, B, c, g, p, d, y, cones,
        scaling, P, hsdwrk, caches, conewrk, kkt,
        hist, ν, settings, ρ, τ, κ, Δy0, nc, ng, nB, α, TimerOutput()
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
    # Gondzio (IPM 2025) stopping test: relative complementarity μ/(1+|cᵀx|) ≤ ε_o. The embedded μ carries a
    # τ² factor to de-homogenize into original-problem units; pobj is already de-homogenized (÷τ).
    max(pres, dres) < s.settings.feas_tol && μ < s.settings.gap_tol * τ^2 * (1 + abs(pobj))
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
# initkkt! (HSD) / step!
############################################################################################

# HSD bordered do-once: factor F and solve/cache the Woodbury column + capacitance. The border w₂ solve
# needs the refinement tolerances, so step! computes them before calling this. Updates the running ρ-floor
# and the warm start s.Δy0. Returns (ok, ρ, wtuple) where wtuple carries the border solve's counts.
function initkkt!(s::HSDSolver{T}) where {T}
    flag, ρ, wtuple = initkkt!(
        s.kkt, s.H, s.B, s.c, s.g;
        α=s.α[], rgmin=s.ρ[], p=s.p, τ=s.τ[], κ=s.κ[], Qp=s.wrk.Qp, y0=s.Δy0,
    )
    s.ρ[] = max(s.ρ[], ρ)
    flag && copyto!(s.Δy0, s.kkt.Δy2)   # warm start for next iteration's w₂ solve
    return flag, ρ, wtuple
end

function step!(s::HSDSolver{T}) where {T}
    status = CONTINUE

    pbase = cbase = wbase = 0   # base-solve CRAIG per role (woodbury counted into craig1, archive-wins)
    prefn = crefn = wrefn = 0   # refinement CRAIG per role
    ppass = cpass = wpass = 0   # refinement passes per role
    pstat = cstat = wstat = KKT_SOLVED   # refinement exit status per role

    step = zero(T)
    pfres = ppres = pdres = T(NaN)
    cfres = cpres = cdres = T(NaN)
    wfres = wpres = wdres = T(NaN)
    αmin = αmax = T(NaN)
    ρ = zero(T)      # ρ-shift actually applied this step (0 = none / no factorization); recorded below

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
            #
            # choose augmentation parameter α
            #
            setaug!(s, T(CTRL_CAP_HSD))
            #
            # compute the KKT-solve tolerances (the border w₂ solve inside initkkt! needs them)
            #
            #   force: min(θ μ/μ₁, ceil)
            #
            μ1 = isempty(s.hist.μ) ? μ : first(s.hist.μ)
            # forcing (absolute μ-schedule) converted to the solver's absolute residual targets ptol
            # (primal) / ytol (dual). No accuracy floor here: the KKT solve computes its own per-row
            # roundoff floor from the terms it actually differences (refinekkt!/refinebdr!) — the floor is
            # a property of the linear solve, not the Newton step. The homogeneous τ is NOT folded in here —
            # it lives in the HSD residual scaling, not the KKT solve.
            #
            tol = min(s.settings.forcing_frac * μ / μ1, s.settings.forcing_ceil)
            ptol = tol * (1 + s.ng[])
            ytol = tol * (1 + s.nc[])
            τtol = tol * (1 + s.ng[] + s.nc[])   # border-row target (bordered predictor/corrector only)
            #
            # factor F and solve/cache the Woodbury column w₂ + capacitance (the bordered do-once)
            #
            #   [ H  -Bᵀ ] [ Δp2 ]   [ c ]
            #   [ B   0  ] [ Δy2 ] = [ g ] ,   S = Δp2ᵀ W Δp2 + (Δp2 - p/τ)ᵀ Q (Δp2 - p/τ) + κ/τ
            #
            initok, ρ, wtuple = @timeit s.timers "initkkt" initkkt!(s)
            if !initok
                if s.settings.verbose > 1
                    @warn "Failed to initialize KKT solver."
                end

                status = nearstatus(s, NUMERICAL_FAILURE)
            else
                wbase, wrefn, wpass, wstat, wfres, wpres, wdres = wtuple
                #
                # solve for the Mehrota predictor direction
                #
                #   [  H          -Bᵀ             -c ] [ Δpa ]   [ rd - d ]
                #   [  B           0              -g ] [ Δya ] = [ rp     ]
                #   [ cᵀ - 2pᵀQ/τ  gᵀ  pᵀQp/τ² + κ/τ ] [ Δτa ]   [ rτ - κ ]
                #
                pbase, prefn, ppass, pstat, Δτa, Δκa, pfres, ppres, pdres = @timeit s.timers "predictor" solvepredictor!(s, gap; ptol, ytol, τtol)

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
                cbase, crefn, cpass, cstat, Δτ, Δκ, cfres, cpres, cdres = @timeit s.timers "corrector" solvecorrector!(s, μ, gap, Δτa, Δκa; ptol, ytol, τtol)

                for v in vtxs(s.B)
                    if s.K[v] isa CofreeCone
                        fill!(view(w.Δd, colrange(s.B, v)), zero(T))
                    end
                end

                if s.settings.verbose > 1 && (pstat !== KKT_SOLVED || cstat !== KKT_SOLVED || wstat !== KKT_SOLVED)
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
                #   
                # compute optimal augmentation window
                #
                #   α* ∈ [αmin, αmax]
                #
                pok = pstat === KKT_SOLVED
                cok = cstat === KKT_SOLVED
                wok = wstat === KKT_SOLVED
                state = pok && cok && wok

                if s.settings.policy == 1
                    # Tier 1: aggregate the predictor, corrector, and Woodbury solves (worst case each end).
                    αmin, αmax = augwindow1(s.α[], _nanmax(pfres, cfres, wfres), _nanmax(pdres, cdres, wdres), ptol, ytol)
                else
                    αmin = augmin(s.α[], pfres, ptol, state, pbase, max(ppass, cpass, wpass), T(CTRL_BDG_HSD))
                    αmax = augmax(s.α[], pdres, ytol, state, pbase, T(CTRL_GAP_HSD))
                end

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

    push!(s.hist, (; μ, step, pres, dres, gap, α=s.α[], ρ, τ=s.τ[], κ=s.κ[],
        pbase, prefn, ppass, pstat, cbase, crefn, cpass, cstat, wbase, wrefn, wpass, wstat,
        pfres, ppres, pdres, cfres, cpres, cdres, wfres, wpres, wdres, αmin, αmax))

    return status
end

function reinit!(solver::HSDSolver{T}, prob::IPMProblem{T}; frac::Real=0.1) where {T}
    @assert ncols(prob.B) == ncols(solver.B)
    @assert nrows(prob.B) == nrows(solver.B)
    @assert nvtxs(prob.B) == nvtxs(solver.B)
    @assert nouts(prob.B) == nouts(solver.B)

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

    solver.ρ[] = resolve_rgmin(solver.B)

    return solver
end
