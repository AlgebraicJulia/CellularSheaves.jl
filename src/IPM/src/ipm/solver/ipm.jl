struct IPMSolver{T, I, V} <: AbstractSolver{T}
    Q::BlockSparseMatrix{T, I}
    H::BlockSparseMatrix{T, I}
    B::BlockSparseMatrix{T, I}
    f::FVector{T}
    g::FVector{T}
    p::FVector{T}
    d::FVector{T}
    y::FVector{T}
    K::FVector{V}
    scaling::IPMScaling{T}
    C::FPermutation{I}
    R::FPermutation{I}
    wrk::IPMWorkspace{T}
    caches::Caches{T, I}
    sched::ConeSchedule{T, I}
    kkt::UzawaSolver{:L, T, I}
    hist::IPMHistory{T}
    ν::Int
    settings::IPMSettings{T}
    ρ::FScalar{T}
    nf::FScalar{T}
    ng::FScalar{T}
    sg::FScalar{T}     # ‖g‖ in original (unscaled) units — stopping-test primal denominator
    sf::FScalar{T}     # ‖f‖ in original (unscaled) units — stopping-test dual denominator
    nB::FScalar{T}      # ‖B‖ — fixed for the solver's lifetime; the cold-start augmentation anchor
    δ::FScalar{T}       # reciprocal augmentation 1/α; owned by setaug!
    timers::TimerOutput
end

function result(s::IPMSolver{T}, status::IPMStatus) where {T}
    pu = copy(s.p)
    du = copy(s.d)
    yu = copy(s.y)

    unscale!(pu, du, yu, s.scaling)

    p = Vector{T}(undef, length(s.p))
    d = Vector{T}(undef, length(s.d))
    y = Vector{T}(undef, length(s.y))

    ldiv!(p, s.C, pu)
    ldiv!(d, s.C, du)
    ldiv!(y, s.R, yu)

    niter = 0
    nsolve = 0

    for row in s.hist
        niter += 1
        nsolve += row.piter + row.citer
    end

    w = s.wrk; scl = s.scaling
    if isempty(s.hist)
        mu = pres = dres = pobj = dobj = T(NaN)
    else
        residuals!(s)
        pQp  = dot(s.p, s.Q, s.p)
        pobj = pQp / 2 - dot(s.f, s.p)
        dobj = dot(s.g, s.y) - pQp / 2
        mu   = iszero(s.ν) ? T(NaN) : dot(s.p, s.d) / s.ν
        pres = scalenorm(w.rp, scl.rscl) / (1 + s.sg[])
        dres = scalenorm(w.rd, scl.cscl) / (1 + s.sf[])
    end

    return IPMResult{T}(p, d, y, status, niter, nsolve, s.hist, s.timers,
                        mu, pres, dres, pobj, dobj)
end

############################################################################################
# residuals!
############################################################################################

#
# compute negated residuals
#
#   [ rd ] = [ d + f ] - [  Q  -Bᵀ ] [ p ]
#   [ rp ]   [   g   ]   [  B   0  ] [ y ]
#
function residuals!(s::IPMSolver{T}) where {T}
    w = s.wrk
    mulkkt!(w.rd, w.rp, s.Q, s.B, s.p, s.y)

    @inbounds for i in eachindex(w.rd, s.d, s.f)
        w.rd[i] = s.d[i] + s.f[i] - w.rd[i]
    end

    @inbounds for i in eachindex(w.rp, s.g)
        w.rp[i] = s.g[i] - w.rp[i]
    end

    return w.rd, w.rp
end

############################################################################################
# solvepredictor! / solvecorrector!
############################################################################################

#
# solve for the Mehrotra predictor direction
#
#   [ H  -Bᵀ ] [ Δpa ]   [ rd - d ]
#   [ B   0  ] [ Δya ] = [ rp     ]
#
function solvepredictor!(s::IPMSolver{T}; ptol::T, ytol::T) where {T}
    return solvepredictor!(
        s.wrk, s.kkt, s.settings, s.H, s.B, s.Q, s.d;
        ptol, ytol,
    )
end

function solvepredictor!(
        w::IPMWorkspace{T},
        kkt::KKTSolver{T},
        set::IPMSettings{T},
        H::BlockSparseMatrix{T},
        B::BlockSparseMatrix{T},
        Q::BlockSparseMatrix{T},
        d::AbstractVector{T};
        ptol::T,
        ytol::T,
    ) where {T}
    axpby!(-1, d, 0, w.f)
    axpby!(1, w.rd, 1, w.f)
    #
    # solve for the directions Δpa and Δya to force_tol (base + internal refinement)
    #
    #   [ H  -Bᵀ ] [ Δpa ]   [ rd - d ]
    #   [ B   0  ] [ Δya ] = [ rp     ]
    #
    piter, ppass, pstat, dmin, dmax = solvekkt!(
        kkt, w.Δpa, w.Δya, H, B, w.f, w.rp;
        gtol=ptol, ftol=ytol, stall=set.refine_stall_tol, rfmax=set.refine_max_iter, cgmax=set.newton_max_iter,
    )
    #
    # recover Δda:
    #
    #   Δda ← Q Δpa - Bᵀ Δya - rd
    #
    copyto!(w.Δda, w.rd)
    mul!(w.Δda, B', w.Δya, -1, -1)
    mul!(w.Δda, Q, w.Δpa, 1, 1)

    return piter, ppass, pstat, dmin, dmax
end

#
# solve for the Mehrotra combined direction
#
#   [ H  -Bᵀ ] [ Δp ]   [ rd* ]
#   [ B   0  ] [ Δy ] = [ rp  ]
#
# where rd* is the corrected dual residual
#
function solvecorrector!(s::IPMSolver{T}, μ::T; ptol::T, ytol::T) where {T}
    return solvecorrector!(
        s.wrk, s.kkt, s.settings, s.H, s.B, s.Q, s.K, s.p, s.d,
        s.caches, s.sched, s.ν, μ;
        ptol, ytol,
    )
end

function solvecorrector!(
        w::IPMWorkspace{T},
        kkt::KKTSolver{T},
        set::IPMSettings{T},
        H::BlockSparseMatrix{T},
        B::BlockSparseMatrix{T},
        Q::BlockSparseMatrix{T},
        K::AbstractVector,
        p::AbstractVector{T},
        d::AbstractVector{T},
        caches::Caches{T},
        sched::ConeSchedule{T},
        ν::Integer,
        μ::T;
        ptol::T,
        ytol::T,
    ) where {T}
    #
    # compute the largest step length τa ∈ (0, 1]
    # such that the perturbed iterates
    #
    #   p + τa Δpa ∈ K
    #   d + τa Δda ∈ K*
    #
    # lie within their respective cones
    #
    τa = maxsteps(sched, K, p, d, w.Δpa, w.Δda, caches, B, w.step; step_frac=one(T))
    #
    # compute the centering parameter
    #
    #   σμ ← clamp(μa (μa / μ)², 0, μ)
    #
    # where
    #
    #   μa  = ⟨p + τa Δpa, d + τa Δda⟩ / ν
    #
    σμ = zero(T)

    for j in cols(B)
        σμ += (p[j] + τa * w.Δpa[j]) * (d[j] + τa * w.Δda[j])
    end

    if iszero(ν)
        σμ = zero(T)          # no cones ⇒ no centering; the affine step is exact
    else
        σμ /= ν
        σμ = clamp(σμ * (σμ / μ)^2, zero(T), μ)
    end
    #
    # set f to the Mehrota corrector term:
    #
    #   f ← -d + (σμ e - Δpa ∘ Δda) / p
    #
    # where e is the Jordan identity element e ∈ K.
    #
    initcorrector!(sched, K, w.f, caches, p, d, w.Δpa, w.Δda, σμ, B)

    axpy!(1, w.rd, w.f)
    #
    # solve for the directions Δp and Δy
    #
    #   [ H  -Bᵀ ] [ Δp ]   [ rd* ]
    #   [ B   0  ] [ Δy ] = [ rp  ]
    #
    copyto!(w.Δp, w.Δpa)
    copyto!(w.Δy, w.Δya)

    citer, cpass, cstat, _, _ = solvekkt!(
        kkt, w.Δp, w.Δy, H, B, w.f, w.rp;
        xwarm=true, ywarm=true, gtol=ptol, ftol=ytol, stall=set.refine_stall_tol, rfmax=set.refine_max_iter, cgmax=set.newton_max_iter,
    )
    #
    # recover Δd:
    #
    #   Δd ← Q Δp - Bᵀ Δy - rd
    #
    copyto!(w.Δd, w.rd)
    mul!(w.Δd, B', w.Δy, -1, -1)
    mul!(w.Δd, Q, w.Δp, 1, 1)

    return citer, cpass, cstat
end

############################################################################################
# startingpoint!
############################################################################################

function startingpoint!(
        p::AbstractVector{T},
        d::AbstractVector{T},
        B::BlockSparseMatrix{T},
        Q::BlockSparseMatrix{T},
        g::AbstractVector{T},
        f::AbstractVector{T},
        K::AbstractVector,
    ) where {T}
    identitypoint!(p, d, B, K)

    z = B * p
    w = Q * p

    np = norm(p)
    nz = norm(z)
    nw = norm(w)

    if nz > eps(T) * np
        sp = max(one(T), norm(g) / nz)
    else
        sp = one(T)
    end

    if np > eps(T)
        sd = max(one(T), (norm(f) + sp * nw) / np)
    else
        sd = one(T)
    end

    lmul!(sp, p)
    lmul!(sd, d)

    return p, d
end

############################################################################################
# constructor / reinit! / init
############################################################################################

function reinit!(s::IPMSolver; p0=nothing, d0=nothing, y0=nothing)
    return reinit!(s, p0, d0, y0)
end

function reinit!(s::IPMSolver{T}, p0, d0, y0) where {T}
    if isnothing(p0) && isnothing(d0)
        startingpoint!(s.p, s.d, s.B, s.Q, s.g, s.f, s.K)

        if isnothing(y0)
            fill!(s.y, zero(T))
        else
            mul!(s.y, s.R, y0)
            s.y ./= s.scaling.rscl
        end
    else
        isnothing(p0) || mul!(s.p, s.C, p0)
        isnothing(d0) || mul!(s.d, s.C, d0)

        if isnothing(d0)
            for v in vtxs(s.B)
                r = colrange(s.B, v)
                dualshadow!(view(s.d, r), view(s.p, r), cache(s.caches, v, s.K[v]), s.sched.large)
            end
        elseif isnothing(p0)
            for v in vtxs(s.B)
                r = colrange(s.B, v)
                primalshadow!(view(s.p, r), view(s.d, r), cache(s.caches, v, s.K[v]), s.sched.large)
            end
        end

        if isnothing(y0)
            fill!(s.y, zero(T))
        else
            mul!(s.y, s.R, y0)
        end

        scale!(s.p, s.d, s.y, s.scaling)
    end

    for v in vtxs(s.B)
        initcache!(cache(s.caches, v, s.K[v]))
    end

    empty!(s.hist)

    return s
end

function IPMSolver(prob::IPMProblem{T, I}, settings::IPMSettings{T}; p0=nothing, d0=nothing, y0=nothing) where {T, I}
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

    scaling = IPMScaling{T}(n, m)

    if settings.scale_max_iter > 0
        equilibrate!(scaling, B, Q, f, g; itmax=settings.scale_max_iter)
    end

    kkt = UzawaSolver(S, B)

    C = C * prob.C
    R = prob.R

    p = FVector{T}(undef, n)
    d = FVector{T}(undef, n)
    y = FVector{T}(undef, m)

    caches = Caches(cones, B)

    H = copy(Q)
    sched = ConeSchedule{T}(cones, B, nthreads())
    ipmwrk = IPMWorkspace{T}(m, n, nvtxs(B))
    hist = IPMHistory{T}()
    ρ = FScalar{T}(undef)
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

    solver = IPMSolver(Q, H, B, f, g, p, d, y, cones,
        scaling, C, R, ipmwrk, caches, sched, kkt,
        hist, ν, settings, ρ, nf, ng, sg, sf, nB, δ, TimerOutput()
    )

    return reinit!(solver; p0, d0, y0)
end

function IPMSolver(prob::IPMProblem{T}; p0=nothing, d0=nothing, y0=nothing, kw...) where {T}
    settings = IPMSettings{T}(; kw...)
    return IPMSolver(prob, settings; p0, d0, y0)
end

function CommonSolve.init(prob::IPMProblem{T}, settings::IPMSettings{T}; p0=nothing, d0=nothing, y0=nothing) where {T}
    return IPMSolver(prob, settings; p0, d0, y0)
end

function CommonSolve.init(prob::IPMProblem{T}; kw...) where {T}
    return IPMSolver(prob; kw...)
end

############################################################################################
# mu / isoptimal
############################################################################################
#
# compute the centrality parameter
#
#   μ = pᵀd / ν
#
function mu(s::IPMSolver{T}) where {T}
    if iszero(s.ν)
        μ = zero(T)
    else
        μ = dot(s.p, s.d) / s.ν
    end

    return μ
end

function nearstatus(s::IPMSolver, status::IPMStatus)
    if isnearoptimal(s)
        status = NEAR_OPTIMAL
    end

    return status
end

############################################################################################
# step!
############################################################################################

function step!(s::IPMSolver{T}) where {T}
    status = CONTINUE

    piter = citer = 0       # total CG per role (base + refinement)
    ppass = cpass = 0       # refinement passes per role
    pstat = cstat = KKT_SOLVED   # refinement exit status per role
    step = zero(T)
    dmin = dmax = T(NaN)   # per-solve min-cost δ-window (kktwindow!)
    ρ = zero(T)      # ρ-shift actually applied this step (0 = none / no factorization); recorded below

    w = s.wrk
    #
    # compute negated residuals
    #
    #   [ rd ]   [ d - f ]   [  Q  -Bᵀ ] [ p ]
    #   [ rp ] = [   g   ] - [  B   0  ] [ y ]
    #
    residuals!(s)
    #
    # compute the centrality parameter
    #
    #   μ = pᵀd / ν
    #
    μ = mu(s)
    pres = scalenorm(w.rp, s.scaling.rscl) / (1 + s.sg[])
    dres = scalenorm(w.rd, s.scaling.cscl) / (1 + s.sf[])

    pQp = dot(s.p, s.Q, s.p)
    pobj = pQp / 2 - dot(s.f, s.p)
    dobj = dot(s.g, s.y) - pQp / 2

    if isoptimal(s, pobj, dobj, pres, dres)
        status = OPTIMAL
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

            initok, ρ = @timeit s.timers "initkkt" initkkt!(s)

            if !initok
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

                tol = min(FORCING_FRAC * μ / μ1, FORCING_CEIL)
                ptol = tol * (1 + s.ng[])
                ytol = tol * (1 + s.nf[])
                #
                # solve for the Mehrotra predictor direction
                #
                #   [ H  -Bᵀ ] [ Δpa ]   [ rd - d ]
                #   [ B   0  ] [ Δya ] = [ rp     ]
                #
                piter, ppass, pstat, dmin, dmax = @timeit s.timers "predictor" solvepredictor!(s; ptol, ytol)

                for v in vtxs(s.B)
                    if s.K[v] isa CofreeCone
                        fill!(view(w.Δda, colrange(s.B, v)), zero(T))
                    end
                end
                #
                # solve for the Mehrotra combined direction
                #
                #   [ H  -Bᵀ ] [ Δp ]   [ rd* ]
                #   [ B   0  ] [ Δy ] = [ rp  ]
                #
                # where rd* is the corrected dual residual
                #
                citer, cpass, cstat = @timeit s.timers "corrector" solvecorrector!(s, μ; ptol, ytol)

                for v in vtxs(s.B)
                    if s.K[v] isa CofreeCone
                        fill!(view(w.Δd, colrange(s.B, v)), zero(T))
                    end
                end

                if s.settings.verbose > 1 && (pstat !== KKT_SOLVED || cstat !== KKT_SOLVED)
                    @info "KKT solve above target tolerance" pstat cstat
                end
                #
                # find the largest step sizes such that
                #
                #   p + αp Δp ∈ K
                #   d + αd Δd ∈ K*
                #
                step = @timeit s.timers "maxsteps" maxsteps(s, w.Δp, w.Δd; step_frac=s.settings.step_frac)
                #
                # compute the updated iterates
                #
                #   p ← p + α Δp ∈ K
                #   d ← d + α Δd ∈ K*
                #
                axpy!(step, w.Δp, s.p)
                axpy!(step, w.Δd, s.d)
                axpy!(step, w.Δy, s.y)

                if isstalled(s)
                    if s.settings.verbose > 1
                        @warn "Stalling detected."
                    end

                    status = nearstatus(s, STALLED)
                end
            end
        end
    end

    push!(s.hist, (; μ, step, pres, dres, pobj, dobj, δ=s.δ[], ρ, piter, ppass, pstat, citer, cpass, cstat,
        dmin, dmax))

    return status
end

############################################################################################
# solve / reinit!
############################################################################################

"""
    solve(problem::IPMProblem;
        verbose=false,
        step_frac=0.99,
        feas_tol=1e-8,
        gap_tol=1e-8,
        max_iter=100,
        near_factor=1000.0,
        stall_tol=1e-3,
        refine_stall_tol=0.5,
        scale_max_iter=10,
        refine_max_iter=10,
        newton_max_iter=100,
        aug_tol=1e-7,
    )

Solve an [`IPMProblem`](@ref).
"""
function CommonSolve.solve(prob::IPMProblem{T}, settings::IPMSettings{T}; p0=nothing, d0=nothing, y0=nothing) where {T}
    return solve!(init(prob, settings; p0, d0, y0))
end

function CommonSolve.solve(prob::IPMProblem{T}; kw...) where {T}
    return solve!(init(prob; kw...))
end

