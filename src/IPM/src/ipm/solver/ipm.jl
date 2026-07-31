struct IPMSolver{T, I, W, C} <: AbstractSolver{T}
    Q::BlockSparseMatrix{T, I}
    H::BlockSparseMatrix{T, I}
    B::BlockSparseMatrix{T, I}
    c::FVector{T}
    g::FVector{T}
    p::FVector{T}
    d::FVector{T}
    y::FVector{T}
    K::FVector{C}
    scaling::IPMScaling{T}
    P::FPermutation{I}
    wrk::IPMWorkspace{T}
    caches::Caches{T, I}
    conewrk::ConeWorkspace{T}
    kkt::W
    hist::IPMHistory{T}
    ν::Int
    settings::IPMSettings{T}
    ρ::FScalar{T}
    nc::FScalar{T}
    ng::FScalar{T}
    α::FScalar{T}       # effective augmentation penalty; fixed at construction (raug); no controller
    timers::TimerOutput
end

function result(s::IPMSolver{T}, status::IPMStatus) where {T}
    p = Vector{T}(undef, length(s.p))
    d = Vector{T}(undef, length(s.d))
    y = Vector{T}(undef, length(s.y))

    ldiv!(p, s.P, s.p)
    ldiv!(d, s.P, s.d)
    copyto!(y, s.y)

    unscale!(p, d, y, s.scaling)

    niter = 0
    nsolve = 0

    for row in s.hist
        niter += 1
        nsolve += row.pbase + row.prefn + row.cbase + row.crefn
    end

    return IPMResult{T}(p, d, y, status, niter, nsolve, s.hist, s.timers)
end

############################################################################################
# mulkkt!
############################################################################################

#
# compute the matrix-vector product
#
#   [ u ] = [ A  -Bᵀ ] [ x ]
#   [ v ]   [ B   0  ] [ y ]
#
function mulkkt!(
        u::AbstractVector{T},
        v::AbstractVector{T},
        A::AbstractMatrix{T},
        B::BlockSparseMatrix{T},
        x::AbstractVector{T},
        y::AbstractVector{T},
    ) where {T}
    mul!(u, Symmetric(A, :L), x)
    mul!(u, B', y, -one(T), one(T))
    mul!(v, B, x)
    return
end

############################################################################################
# residuals!
############################################################################################

#
# compute negated residuals
#
#   [ rd ] = [ d + c ] - [  Q  -Bᵀ ] [ p ]
#   [ rp ]   [   g   ]   [  B   0  ] [ y ]
#
function residuals!(s::IPMSolver{T}) where {T}
    w = s.wrk
    mulkkt!(w.rd, w.rp, s.Q, s.B, s.p, s.y)

    @inbounds for i in eachindex(w.rd, s.d, s.c)
        w.rd[i] = s.d[i] + s.c[i] - w.rd[i]
    end

    @inbounds for i in eachindex(w.rp, s.g)
        w.rp[i] = s.g[i] - w.rp[i]
    end

    return w.rd, w.rp
end

############################################################################################
# refinekkt! — solver-owned refinement loop for the 2-row system
############################################################################################

function refinekkt!(
    Δp::AbstractVector{T},
    Δy::AbstractVector{T},
    wrk::KKTWorkspace{T},
    H::BlockSparseMatrix{T},
    B::BlockSparseMatrix{T},
    f::AbstractVector{T},
    rp::AbstractVector{T},
    sp::AbstractVector{T},  # scratch for primal residual
    sd::AbstractVector{T},  # scratch for dual residual
    dp::AbstractVector{T},
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
    res0 = res1 = T(NaN)       # residual entering pass 1 (post base solve) and pass 2
    res0_d = res0_p = T(NaN)   # its dual /(1+nc) and primal /(1+ng) components (the ceiling lives in the dual)
    res_exit = T(NaN)          # last residual before returning (any status) — graded severity of the exit

    for i in 1:itmax
        #
        # compute the residuals
        #
        #   [ sd ] = [ f  ] - [ H -Bᵀ ] [ Δp ]
        #   [ sp ]   [ rp ]   [ B  0  ] [ Δy ]
        #
        mulkkt!(sd, sp, H, B, Δp, Δy)
        axpby!(one(T), f,  -one(T), sd)
        axpby!(one(T), rp, -one(T), sp)

        dres = norm(sd, Inf) / (one(T) + nc)
        pres = norm(sp, Inf) / (one(T) + ng)
        res = max(dres, pres)
        res_exit = res
        i == 1 && (res0 = res; res0_d = dres; res0_p = pres)
        i == 2 && (res1 = res)

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
        # update Δp and Δy:
        #
        #   Δp ← Δp + dp
        #   Δy ← Δy + dy
        #
        axpy!(one(T), dp, Δp)
        axpy!(one(T), dy, Δy)
    end

    return npass, niter, status, res0, res1, res0_d, res0_p, res_exit
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
function solvepredictor!(s::IPMSolver{T}; force_tol::T, floor_tol::T) where {T}
    return solvepredictor!(
        s.wrk, s.kkt, s.settings, s.H, s.B, s.Q, s.d, s.nc[], s.ng[];
        force_tol, floor_tol,
    )
end

function solvepredictor!(
        w::IPMWorkspace{T},
        kkt::KKTWorkspace{T},
        set::IPMSettings{T},
        H::BlockSparseMatrix{T},
        B::BlockSparseMatrix{T},
        Q::BlockSparseMatrix{T},
        d::AbstractVector{T},
        nc::T,
        ng::T;
        force_tol::T,
        floor_tol::T,
    ) where {T}
    atol = max(force_tol, floor_tol)

    axpby!(-1, d, 0, w.f)
    axpby!(1, w.rd, 1, w.f)
    #
    # solve for the directions Δpa and Δya
    #
    #   [ H  -Bᵀ ] [ Δpa ]   [ rd - d ]
    #   [ B   0  ] [ Δya ] = [ rp     ]
    #
    pbase = solve_kkt!(kkt, w.Δpa, w.Δya, H, B, w.f, w.rp; atol)
    r0_p = kkt.r0[]; r1_p = kkt.r1[]; craig_p = kkt.itrwrk.stats.status
    s2min_p = kkt.s2min[]; s2max_p = kkt.s2max[]   # σ̂² of B on the predictor base-solve rhs subspace
    #
    # refine Δpa, Δya to force_tol
    #
    ppass, prefn, pstat, pres0, pres1, pres0_d, pres0_p, pres_exit = refinekkt!(
        w.Δpa, w.Δya, kkt, H, B,
        w.f, w.rp, w.sy, w.sp, w.dp, w.dy, nc, ng;
        itmax=set.refine_itmax, force_tol, floor_tol, stall=set.refine_stall,
    )
    #
    # recover Δda:
    #
    #   Δda ← Q Δpa - Bᵀ Δya - rd
    #
    copyto!(w.Δda, w.rd)
    mul!(w.Δda, B', w.Δya, -1, -1)
    mul!(w.Δda, Symmetric(Q, :L), w.Δpa, 1, 1)

    return pbase, prefn, ppass, pstat, pres0, pres1, r0_p, craig_p, r1_p, pres0_d, pres0_p, pres_exit, s2min_p, s2max_p
end

#
# solve for the Mehrotra combined direction
#
#   [ H  -Bᵀ ] [ Δp ]   [ rd* ]
#   [ B   0  ] [ Δy ] = [ rp  ]
#
# where rd* is the corrected dual residual
#
function solvecorrector!(s::IPMSolver{T}, μ::T; force_tol::T, floor_tol::T) where {T}
    return solvecorrector!(
        s.wrk, s.kkt, s.settings, s.H, s.B, s.Q, s.K, s.p, s.d,
        s.caches, s.conewrk, s.ν, s.nc[], s.ng[], μ;
        force_tol, floor_tol,
    )
end

function solvecorrector!(
        w::IPMWorkspace{T},
        kkt::KKTWorkspace{T},
        set::IPMSettings{T},
        H::BlockSparseMatrix{T},
        B::BlockSparseMatrix{T},
        Q::BlockSparseMatrix{T},
        K::AbstractVector,
        p::AbstractVector{T},
        d::AbstractVector{T},
        caches::Caches{T},
        conewrk::ConeWorkspace{T},
        ν::Integer,
        nc::T,
        ng::T,
        μ::T;
        force_tol::T,
        floor_tol::T,
    ) where {T}
    atol = max(force_tol, floor_tol)
    #
    # compute the largest step length τa ∈ (0, 1]
    # such that the perturbed iterates
    #
    #   p + τa Δpa ∈ K
    #   d + τa Δda ∈ K*
    #
    # lie within their respective cones
    #
    τa = one(T)

    for v in vtxs(B)
        τpv, τdv = maxsteps(K[v], v, p, d, w.Δpa, w.Δda, caches, B, conewrk)
        τa = min(τa, τpv, τdv)
    end
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
    for v in vtxs(B)
        initcorrector!(K[v], v, w.f, caches, p, d, w.Δpa, w.Δda, σμ, B, conewrk)
    end

    axpy!(1, w.rd, w.f)
    #
    # solve for the directions Δp and Δy
    #
    #   [ H  -Bᵀ ] [ Δp ]   [ rd* ]
    #   [ B   0  ] [ Δy ] = [ rp  ]
    #
    cbase = solve_kkt!(kkt, w.Δp, w.Δy, H, B, w.f, w.rp, w.Δya; atol)
    r0_c = kkt.r0[]; r1_c = kkt.r1[]; craig_c = kkt.itrwrk.stats.status
    #
    # use iterative refinement to improve
    # the solutions Δp and Δy
    #
    cpass, crefn, cstat, cres0, cres1, cres0_d, cres0_p, cres_exit = refinekkt!(
        w.Δp, w.Δy, kkt, H, B,
        w.f, w.rp, w.sy, w.sp, w.dp, w.dy, nc, ng;
        itmax=set.refine_itmax, force_tol, floor_tol, stall=set.refine_stall,
    )
    #
    # recover Δd:
    #
    #   Δd ← Q Δp - Bᵀ Δy - rd
    #
    copyto!(w.Δd, w.rd)
    mul!(w.Δd, B', w.Δy, -1, -1)
    mul!(w.Δd, Symmetric(Q, :L), w.Δp, 1, 1)

    return cbase, crefn, cpass, cstat, cres0, cres1, r0_c, craig_c, r1_c, cres0_d, cres0_p, cres_exit
end

############################################################################################
# startingpoint
############################################################################################

function startingpoint(B::BlockSparseMatrix{T}, g::AbstractVector{T}, c::AbstractVector{T}, p::AbstractVector{T}) where {T}
    z = B * p

    np = norm(p)
    nz = norm(z)

    if nz > eps(T) * np
        sp = max(one(T), norm(g) / nz)
    else
        sp = one(T)
    end

    if np > eps(T)
        sd = max(one(T), norm(c) / np)
    else
        sd = one(T)
    end

    return sp, sd
end

############################################################################################
# init
############################################################################################

function CommonSolve.init(prob::IPMProblem{T}; kw...) where {T}
    return init(prob, IPMSettings{T}(; kw...))
end

function CommonSolve.init(prob::IPMProblem{T, I}, settings::IPMSettings{T}) where {T, I}
    n = size(prob.B, 2)
    m = size(prob.B, 1)
    ν = conedegree(prob.K, prob.B)

    scaling = IPMScaling{T}(n, m)

    if settings.scale_itmax > 0
        B = copy(prob.B)
        Q = copy(prob.Q)
        c = copy(prob.c)
        g = copy(prob.g)

        equilibrate!(scaling, B, Q, c, g; itmax=settings.scale_itmax)
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
    sp, sd = startingpoint(B, g, c, p)
    lmul!(sp, p)
    lmul!(sd, d)

    caches = Caches(cones, B)

    for v in vtxs(B)
        initcache!(cache(caches, v, cones[v]))
    end

    H = allocblockdiag(B)
    conewrk = ConeWorkspace{T}(cones, B)
    ipmwrk = IPMWorkspace{T}(m, n)
    hist = IPMHistory{T}()
    ρ = FScalar{T}(undef)
    nc = FScalar{T}(undef)
    ng = FScalar{T}(undef)
    α = FScalar{T}(undef)

    ρ[] = settings.rgmin
    nc[] = norm(c)
    ng[] = norm(g)

    nH = (sd / sp) * sqrt(n)
    nQ = norm(Symmetric(Q, :L))
    nB = norm(B)
    α[] = settings.aaug + settings.raug * (nH + nQ) / nB^2

    return IPMSolver(Q, H, B, c, g, p, d, y, cones,
        scaling, P, ipmwrk, caches, conewrk, kkt,
        hist, ν, settings, ρ, nc, ng, α, TimerOutput()
    )
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

function isoptimal(s::IPMSolver{T}, μ::T, pres::T, dres::T) where {T}
    pQp = dot(s.p, Symmetric(s.Q, :L), s.p)
    pobj = pQp / 2 - dot(s.c, s.p)
    dobj = dot(s.g, s.y) - pQp / 2
    return max(pres, dres) < s.settings.feas_tol && (iszero(s.ν) || μ < s.settings.gap_tol || pobj - dobj < s.settings.gap_tol * (1 + abs(pobj) + abs(dobj)))
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

    pbase = cbase = 0       # base-solve CRAIG per role
    prefn = crefn = 0       # refinement CRAIG per role
    ppass = cpass = 0       # refinement passes per role
    pstat = cstat = REACHED_FORCE   # refinement exit status per role
    pres0 = pres1 = cres0 = cres1 = T(NaN)   # pass-0/pass-1 residuals per role
    pres0_d = pres0_p = cres0_d = cres0_p = T(NaN)   # pass-1 residual split into dual/primal per role
    pres_exit = cres_exit = T(NaN)           # refinement exit residual per role (graded severity)
    r0_p = r0_c = T(NaN)                     # pre-CRAIG base residual per role
    r1_p = r1_c = T(NaN)                     # post-CRAIG base residual per role
    s2min_p = s2max_p = T(NaN)               # σ̂²min/max of B on the predictor base-solve rhs subspace
    craig_p = craig_c = ""                   # CRAIG termination status per role
    bar_hdiag_med = bar_hdiag_frac_mid = T(NaN)   # pre-Q barrier-Hessian diagonal stats (cone coords)
    step = zero(T)

    w = s.wrk
    #
    # compute negated residuals
    #
    #   [ rd ]   [ d - c ]   [  Q  -Bᵀ ] [ p ]
    #   [ rp ] = [   g   ] - [  B   0  ] [ y ]
    #
    residuals!(s)
    #
    # compute the centrality parameter
    #
    #   μ = pᵀd / ν
    #
    μ = mu(s)
    pres = norm(w.rp) / (1 + s.ng[])
    dres = norm(w.rd) / (1 + s.nc[])

    if isoptimal(s, μ, pres, dres)
        status = OPTIMAL
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
        # barrier-Hessian diagonal stats (pre-Q: the dⱼ/pⱼ degeneracy signature, before Q masks it)
        #
        bar_hdiag_med, bar_hdiag_frac_mid = barrier_hdiag_stats(s)
        #
        # add the quadratic term
        #
        #   H ← H + Q
        #
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
                # solve for the Mehrotra predictor direction
                #
                #   [ H  -Bᵀ ] [ Δpa ]   [ rd - d ]
                #   [ B   0  ] [ Δya ] = [ rp     ]
                #
                pbase, prefn, ppass, pstat, pres0, pres1, r0_p, craig_p, r1_p, pres0_d, pres0_p, pres_exit, s2min_p, s2max_p = @timeit s.timers "predictor" solvepredictor!(s; force_tol, floor_tol)

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
                cbase, crefn, cpass, cstat, cres0, cres1, r0_c, craig_c, r1_c, cres0_d, cres0_p, cres_exit = @timeit s.timers "corrector" solvecorrector!(s, μ; force_tol, floor_tol)

                for v in vtxs(s.B)
                    if s.K[v] isa CofreeCone
                        fill!(view(w.Δd, colrange(s.B, v)), zero(T))
                    end
                end

                if s.settings.verbose > 1 && (pstat != REACHED_FORCE || cstat != REACHED_FORCE)
                    @info "KKT solve above target tolerance" pstat cstat
                end
                #
                # find the largest step sizes such that
                #
                #   p + αp Δp ∈ K
                #   d + αd Δd ∈ K*
                #
                @timeit s.timers "maxsteps" pstep, dstep = maxsteps(s, w.Δp, w.Δd; step_frac=s.settings.step_frac)
                #
                # take a common step in the primal and dual spaces
                #
                #   α = min(αp, αd)
                #
                step = min(pstep, dstep)
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
                elseif isnumfail(s)
                    if s.settings.verbose > 1
                        @warn "Step collapse detected."
                    end

                    status = nearstatus(s, NUMERICAL_FAILURE)
                end
            end
        end
    end

    push!(s.hist, (; μ, step, pres, dres, α=s.α[], ρ=s.ρ[], pbase, prefn, ppass, pstat, cbase, crefn, cpass, cstat, pres0, pres1, cres0, cres1, r0_p, r0_c, craig_p, craig_c, r1_p, r1_c, pres0_d, pres0_p, cres0_d, cres0_p, pres_exit, cres_exit, bar_hdiag_med, bar_hdiag_frac_mid, s2min_p, s2max_p))
    if status == CONTINUE && atfloor(s.hist; patience=s.settings.floor_patience)
        if s.settings.verbose > 1
            @warn "Refinement floor reached $(s.settings.floor_patience) consecutive times"
        end

        status = nearstatus(s, NUMERICAL_FAILURE)
    end

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
        itmax=100,
        near_factor=1000.0,
        stall_tol=1e-6,
        forcing_frac=0.1,
        forcing_ceil=0.3,
        refine_itmax=10,
        refine_stall=0.5,
        floor_patience=3,
        scale_itmax=10,
        rgmin=1e-9,
        rgmax=1e-6,
        aaug=0.0,
        raug=1e7,
    )

Solve an [`IPMProblem`](@ref).
"""
function CommonSolve.solve(prob::IPMProblem{T}; kw...) where {T}
    settings = IPMSettings{T}(; kw...)
    return solve!(init(prob, settings))
end

"""
    reinit!(solver, prob; frac=0.1)

Reinitialize an [`IPMSolver`](@ref), updating the vectors `b` and `g`.
"""
function reinit!(solver::IPMSolver{T}, prob::IPMProblem{T}; frac::Real=0.1) where {T}
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

    mul!(solver.c, solver.P, c)
    copyto!(solver.g, g)

    solver.nc[] = norm(solver.c)
    solver.ng[] = norm(solver.g)

    p, d, y = identitypoint(solver.B, solver.K)
    sp, sd = startingpoint(solver.B, solver.g, solver.c, p)
    lmul!(sp, p)
    lmul!(sd, d)
    axpby!(frac, p, 1 - frac, solver.p)
    axpby!(frac, d, 1 - frac, solver.d)
    axpby!(frac, y, 1 - frac, solver.y)

    empty!(solver.hist)

    solver.ρ[] = solver.settings.rgmin

    nH = (sd / sp) * sqrt(T(n))
    nQ = norm(Symmetric(solver.Q, :L))
    nB = norm(solver.B)
    solver.α[] = solver.settings.aaug + solver.settings.raug * (nH + nQ) / nB^2

    return solver
end
