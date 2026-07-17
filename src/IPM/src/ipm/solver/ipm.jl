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
    scaling::Scaling{T}
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
    nB::FScalar{T}
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
    npred = 0
    ncorr = 0

    for row in s.hist
        niter += 1
        npred += row.npred
        ncorr += row.ncorr
    end

    return IPMResult{T}(p, d, y, status, niter, npred, ncorr, s.hist, s.timers)
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
# refinekkt! — A2/A3: solver-owned refinement loop for 2-row system
############################################################################################

function refinekkt!(
    Δp::AbstractVector{T},
    Δy::AbstractVector{T},
    Δd::AbstractVector{T},
    wrk::KKTWorkspace{T},
    set::KKTSettings{T},
    H::BlockSparseMatrix{T},
    B::BlockSparseMatrix{T},
    Q::BlockSparseMatrix{T},
    f::AbstractVector{T},
    rp::AbstractVector{T},
    rd::AbstractVector{T},
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
    nH::T,
    nB::T,
) where {T}
    niter = 0
    status = REACHED_FORCE
    prv = typemax(T)

    nf = norm(f)
    nrp = norm(rp)

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

        if res ≤ force_tol
            status = REACHED_FORCE
            break
        end
        #
        # compute backward-error threshold:
        #
        #   100ε max(σd/(1+‖c‖), σp/(1+‖g‖))
        #
        # where
        #
        #   σd = ‖H‖ ‖Δp‖ + ‖B‖ ‖Δy‖ + ‖f‖
        #   σp = ‖B‖ ‖Δp‖            + ‖rp‖
        #
        np, ny = norm(Δp), norm(Δy)

        σd = nH * np + nB * ny + nf
        σp = nB * np           + nrp

        dynam_tol = 100eps(T) * max(σd / (1 + nc), σp / (1 + ng))

        if res ≤ max(dynam_tol, floor_tol)
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
        #   [ H -Bᵀ ] [ dp ] = [ sp ]
        #   [ B  0  ] [ dy ]   [ sd ]
        #
        niter += solve_kkt!(wrk, set, dp, dy, H, B, sd, sp; atol=max(force_tol, dynam_tol, floor_tol))
        #
        # update Δp and Δy:
        #
        #   Δp ← Δp + dp
        #   Δy ← Δy + dy
        #
        axpy!(one(T), dp, Δp)
        axpy!(one(T), dy, Δy)
    end
    #
    # re-compute Δd:
    #
    #   Δd = Q Δp - Bᵀ Δy - rd 
    #
    copyto!(Δd, rd)
    mul!(Δd,           B',     Δy, one(T),  one(T))
    mul!(Δd, Symmetric(Q, :L), Δp, one(T), -one(T))

    return niter, status
end

############################################################################################
# newton!
############################################################################################

#
# solve for the directions Δp and Δy
#
#   [ H -Bᵀ ] [ Δp ] = [ f  ]
#   [ B  0  ] [ Δy ]   [ rp ]
#
function newton!(
        Δp::AbstractVector{T},
        Δy::AbstractVector{T},
        Δd::AbstractVector{T},
        wrk::KKTWorkspace{T},
        set::KKTSettings{T},
        H::BlockSparseMatrix{T},
        B::BlockSparseMatrix{T},
        f::AbstractVector{T},
        rp::AbstractVector{T},
        rd::AbstractVector{T},
        Q::BlockSparseMatrix{T},
        y0 = nothing;
        atol::T = T(0.1),
    ) where {T}
    #
    # solve for Δp and Δy:
    #
    #   [ H -Bᵀ ] [ Δp ] = [ f  ]
    #   [ B  0  ] [ Δy ]   [ rp ]
    #
    niter = solve_kkt!(wrk, set, Δp, Δy, H, B, f, rp, y0; atol)
    #
    # recover Δd:
    #
    #   Δd ← Q Δp - Bᵀ Δy - rd
    #
    copyto!(Δd, rd)
    mul!(Δd, B', Δy, -1, -1)
    mul!(Δd, Symmetric(Q, :L), Δp, 1, 1)

    return niter
end

############################################################################################
# steplengths
############################################################################################

function steplengths(
        pmax::T,
        dmax::T,
        rd::AbstractVector{T},
        Δp::AbstractVector{T},
        Δy::AbstractVector{T},
        Δd::AbstractVector{T},
        Q::BlockSparseMatrix{T},
        B::BlockSparseMatrix{T},
        g::AbstractVector{T},
        h::AbstractVector{T}
    ) where {T}
    mul!(g, Symmetric(Q, :L), Δp)
    mul!(h,           B',     Δy)

    cgg = dot( g,  g)
    chh = dot( h,  h)
    cuu = dot(Δd, Δd)

    cgh = dot( g,  h)
    cgu = dot( g, Δd)
    chu = dot( h, Δd)

    crg = -dot(rd,  g)
    crh = -dot(rd,  h)
    cru = -dot(rd, Δd)

    aA = cuu - chu^2 / chh
    aB = cgg - cgh^2 / chh

    bA = -2(cru + pmax * cgu - chu * (crh + pmax * cgh) / chh)
    bB =  2(crg - dmax * cgu - cgh * (crh - dmax * chu) / chh)

    tol = sqrt(eps(T))

    flagA = aA > tol * cuu

    if !flagA
        dA = dmax
    else
        dA = clamp(-bA / 2aA, zero(T), dmax)
    end

    flagB = aB > tol * cgg

    if !flagB
        pB = pmax
    else
        pB = clamp(-bB / 2aB, zero(T), pmax)
    end

    yA = (crh + pmax * cgh - dA   * chu) / chh
    yB = (crh + pB *   cgh - dmax * chu) / chh

    fA = cgg * pmax^2 + chh * yA^2 + cuu * dA^2 +
         2(crg * pmax - crh * yA - cru * dA - cgh * pmax * yA - cgu * pmax * dA + chu * yA * dA)
    fB = cgg * pB^2 + chh * yB^2 + cuu * dmax^2 +
         2(crg * pB - crh * yB - cru * dmax - cgh * pB * yB - cgu * pB * dmax + chu * yB * dmax)

    if !flagB || fA <= fB
        pstep, dstep = pmax, dA
    else
        pstep, dstep = pB, dmax
    end

    pstep = max(min(pmax, dmax), pstep)
    ystep = (crh + pstep * cgh - dstep * chu) / chh

    return pstep, ystep, dstep
end

function steplengths(s::IPMSolver{T}, pmax::T, dmax::T, Δp::AbstractVector{T}, Δy::AbstractVector{T}, Δd::AbstractVector{T}) where {T}
    w = s.wrk
    return steplengths(pmax, dmax, w.rd, Δp, Δy, Δd, s.Q, s.B, w.g, w.h)
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
    w = s.wrk
    atol = max(force_tol, floor_tol)

    axpby!(-1, s.d, 0, w.f)
    axpby!(1, w.rd, 1, w.f)
    #
    # solve for the directions Δpa and Δya
    #
    #   [ H  -Bᵀ ] [ Δpa ]   [ rd - d ]
    #   [ B   0  ] [ Δya ] = [ rp     ]
    #
    return newton!(w.Δpa, w.Δya, w.Δda, s.kkt, s.settings.kkt, s.H, s.B, w.f, w.rp, w.rd, s.Q; atol)
end

#
# solve for the Mehrotra combined direction
#
#   [ H  -Bᵀ ] [ Δp ]   [ rd* ]
#   [ B   0  ] [ Δy ] = [ rp  ]
#
# where rd* is the corrected dual residual
#
function solvecorrector!(s::IPMSolver{T}, μ::T; force_tol::T, floor_tol::T, nH::T, nB::T) where {T}
    w = s.wrk
    atol = max(force_tol, floor_tol)
    #
    # compute the largest step lengths τpa, τda ∈ (0, 1]
    # such that the perturbed iterates
    #
    #   p + τpa Δpa ∈ K
    #   d + τda Δda ∈ K*
    #
    # lie within their respective cones
    #
    τpa = one(T)
    τda = one(T)

    for v in vtxs(s.B)
        τpv, τdv = maxsteps(s, v, w.Δpa, w.Δda)
        τpa = min(τpa, τpv)
        τda = min(τda, τdv)
    end
    #
    # compute the centering parameter
    #
    #   σμ ← clamp(μa (μa / μ)², 0, μ)
    #
    # where
    #
    #   μa  = ⟨p + τpa Δpa, d + τda Δda⟩ / ν
    #
    σμ = zero(T)

    for j in cols(s.B)
        σμ += (s.p[j] + τpa * w.Δpa[j]) * (s.d[j] + τda * w.Δda[j])
    end

    σμ /= s.ν
    σμ = clamp(σμ * (σμ / μ)^2, zero(T), μ)
    #
    # set f to the Mehrota corrector term:
    #
    #   f ← -d + (σμ e - Δpa ∘ Δda) / p
    #
    # where e is the Jordan identity element e ∈ K.
    #
    for v in vtxs(s.B)
        initcorrector!(s.K[v], v, w.f, s.caches, s.p, s.d, w.Δpa, w.Δda, σμ, s.B, s.conewrk)
    end

    axpy!(1, w.rd, w.f)
    #
    # solve for the directions Δp and Δy
    #
    #   [ H  -Bᵀ ] [ Δp ]   [ rd* ]
    #   [ B   0  ] [ Δy ] = [ rp  ]
    #
    ncorr = newton!(w.Δp, w.Δy, w.Δd, s.kkt, s.settings.kkt, s.H, s.B, w.f, w.rp, w.rd, s.Q, w.Δya; atol)
    #
    # use iterative refinement to improve
    # the solutions Δp and Δy 
    #
    nrefn, refstat = refinekkt!(
        w.Δp, w.Δy, w.Δd, s.kkt, s.settings.kkt, s.H, s.B, s.Q,
        w.f, w.rp, w.rd, w.sy, w.sp, w.dp, w.dy, s.nc[], s.ng[];
        itmax=s.settings.refine_itmax, force_tol, floor_tol, stall=s.settings.refine_stall,
        nH, nB
    )

    ncorr += nrefn

    return ncorr, refstat
end

############################################################################################
# startingpoint
############################################################################################

function startingpoint(B::BlockSparseMatrix{T}, g::AbstractVector{T}, c::AbstractVector{T}, cones::AbstractVector) where {T}
    p, d, y = identitypoint(B, cones)

    z = B * p

    np = norm(p)
    nz = norm(z)

    if nz > eps(T) * np
        ξ = max(one(T), norm(g) / nz)
    else
        ξ = one(T)
    end

    if np > eps(T)
        η = max(one(T), norm(c) / np)
    else
        η = one(T)
    end

    lmul!(ξ, p)
    lmul!(η, d)

    return p, d, y
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

    scaling = Scaling{T}(n, m)

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

    R, P, B, kkt = make_kkt(settings.kkt, B)

    c = P * c
    Q = halfselectvtxs(halfselectvtxs(Q, R.perm), R.perm)
    cones = tounion(prob.K, R.perm)

    p, d, y = startingpoint(B, g, c, cones)

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
    nB = FScalar{T}(undef)

    ρ[] = settings.rgmin
    nc[] = norm(c)
    ng[] = norm(g)
    nB[] = norm(B)

    return IPMSolver(Q, H, B, c, g, p, d, y, cones,
        scaling, P, ipmwrk, caches, conewrk, kkt,
        hist, ν, settings, ρ, nc, ng, nB, TimerOutput()
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
    pobj = dot(s.c, s.p) + pQp / 2
    dobj = dot(s.g, s.y) - pQp / 2
    return max(pres, dres) < s.settings.feas_tol && (iszero(s.ν) || μ < s.settings.gap_tol || pobj - dobj < s.settings.gap_tol * (1 + abs(pobj) + abs(dobj)))
end

############################################################################################
# step!
############################################################################################

function step!(s::IPMSolver{T}) where {T}
    status = CONTINUE

    npred = 0
    ncorr = 0
    pstep = zero(T)
    dstep = zero(T)
    refstat = REACHED_FORCE

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

            if isnearoptimal(s)
                status = NEAR_OPTIMAL
            else
                status = NUMERICAL_FAILURE
            end
        else
            nH = norm(Symmetric(s.H, :L))

            if !@timeit s.timers "initkkt" initkkt!(s, nH)
                if s.settings.verbose > 1
                    @warn "Failed to initialize KKT solver."
                end

                if isnearoptimal(s)
                    status = NEAR_OPTIMAL
                else
                    status = NUMERICAL_FAILURE
                end
            else
                #
                # compute tolerances for predictor and corrector solves
                #
                #   force: min(θμ, ceil)
                #   floor: 100ϵ (1 + max(‖rp‖, ‖rd‖))
                #
                force_tol = min(s.settings.forcing_frac * μ, s.settings.forcing_ceil)
                floor_tol = 100eps(T) * (1 + max(norm(w.rp, Inf), norm(w.rd, Inf)))
                #
                # solve for the Mehrotra predictor direction
                #
                #   [ H  -Bᵀ ] [ Δpa ]   [ rd - d ]
                #   [ B   0  ] [ Δya ] = [ rp     ]
                #
                npred = @timeit s.timers "predictor" solvepredictor!(s; force_tol, floor_tol)
                #
                # solve for the Mehrotra combined direction
                #
                #   [ H  -Bᵀ ] [ Δp ]   [ rd* ]
                #   [ B   0  ] [ Δy ] = [ rp  ]
                #
                # where rd* is the corrected dual residual
                #
                ncorr, refstat = @timeit s.timers "corrector" solvecorrector!(s, μ; force_tol, floor_tol, nH, nB=s.nB[])

                if s.settings.verbose > 1 && refstat != REACHED_FORCE
                    @info "KKT solve above target tolerance" refstat
                end
                #
                # find the largest step sizes such that
                #
                #   p + αp Δp ∈ K
                #   d + αd Δd ∈ K*
                #
                @timeit s.timers "maxsteps" pstep, dstep = maxsteps(s, w.Δp, w.Δd; step_frac=s.settings.step_frac)
                #
                # adjust αp and αd using Meszaros' heuristic
                #
                pstep, ystep, dstep = steplengths(s, pstep, dstep, w.Δp, w.Δy, w.Δd)
                #
                # compute the updated iterates
                #
                #   p ← p + α Δp ∈ K
                #   d ← d + α Δd ∈ K*
                #
                axpy!(pstep, w.Δp, s.p)
                axpy!(dstep, w.Δd, s.d)
                axpy!(ystep, w.Δy, s.y)

                if isstalled(s)
                    if s.settings.verbose > 1
                        @warn "Stalling detected."
                    end

                    if isnearoptimal(s)
                        status = NEAR_OPTIMAL
                    else
                        status = STALLED
                    end
                elseif isnumfail(s)
                    if s.settings.verbose > 1
                        @warn "Step collapse detected."
                    end

                    if isnearoptimal(s)
                        status = NEAR_OPTIMAL
                    else
                        status = NUMERICAL_FAILURE
                    end
                end
            end
        end
    end

    push!(s.hist, (; μ, pstep, dstep, pres, dres, npred, ncorr, refstat))

    if status == CONTINUE && atfloor(s.hist; patience=s.settings.floor_patience)
        if s.settings.verbose > 1
            @warn "Refinement floor reached $(s.settings.floor_patience) consecutive times"
        end

        if isnearoptimal(s)
            status = NEAR_OPTIMAL
        else
            status = NUMERICAL_FAILURE
        end
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
        kkt=UzawaSettings(),
    )

Solve an [`IPMProblem`](@ref).
"""
function CommonSolve.solve(prob::IPMProblem{T}; kw...) where {T}
    settings = IPMSettings{T}(; kw...)
    return solve!(init(prob, settings))
end

"""
    reinit!(solver, prob; frac=0.1, rgfrac=0.0)

Reinitialize an [`IPMSolver`](@ref), updating the vectors `b` and `g`.
"""
function reinit!(solver::IPMSolver{T}, prob::IPMProblem{T}; frac::Real=0.1, rgfrac::Real=0.0) where {T}
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

    mul!(solver.c, solver.P, c)
    copyto!(solver.g, g)

    solver.nc[] = norm(solver.c)
    solver.ng[] = norm(solver.g)

    p, d, y = startingpoint(solver.B, solver.g, solver.c, solver.K)
    axpby!(frac, p, 1 - frac, solver.p)
    axpby!(frac, d, 1 - frac, solver.d)
    axpby!(frac, y, 1 - frac, solver.y)

    empty!(solver.hist)

    ρprev = solver.ρ[]
    ρcold = solver.settings.rgmin
    solver.ρ[] = rgfrac * ρprev + (1 - rgfrac) * ρcold

    return solver
end
