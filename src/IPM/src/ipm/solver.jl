function IPMResult(s::IPMSolver{T}, status::IPMStatus) where {T}
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

function conedegree(cones::AbstractVector, B::BlockSparseMatrix)
    ν = 0

    for v in vtxs(B)
        ν += degree(cones[v], ncols(B, v))
    end

    return ν
end

function residuals!(
        rp::AbstractVector,
        rd::AbstractVector,
        B::BlockSparseMatrix,
        p::AbstractVector,
        d::AbstractVector,
        y::AbstractVector,
        c::AbstractVector,
        g::AbstractVector,
        Q::BlockSparseMatrix,
    )
    #
    # compute the primal residual:
    #
    #   rp = g - B p
    #
    copyto!(rp, g)
    mul!(rp, B, p, -1, 1)
    #
    # compute the dual residual:
    #
    #   rd =  c - d + Q p - Bᵀ y
    #
    copyto!(rd, c)
    mul!(rd, Symmetric(Q, :L), p,  1, 1)
    mul!(rd,           B',     y, -1, 1)
    axpy!(-1, d, rd)

    return rp, rd
end

function scale!(
        cone::AbstractCone,
        v::Integer,
        H::BlockSparseMatrix,
        caches::Caches,
        p::AbstractVector,
        d::AbstractVector,
        B::BlockSparseMatrix,
        Q::BlockSparseMatrix,
        conewrk::ConeWorkspace,
    )
    r = colrange(B, v)
    Hv = block(H, v, v, v)
    cv = cache(caches, v, cone)
    ok = scale!(Hv, view(p, r), view(d, r), cv, conewrk)
    axpy!(true, block(Q, v, v, v), Hv)
    return ok
end

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
        sp::AbstractVector{T},
        sy::AbstractVector{T},
        dp::AbstractVector{T},
        dy::AbstractVector{T},
        y0 = nothing;
        itmax::Integer = 0,
        rtol::T = T(0.3),
        stall::T = T(0.5),
    ) where {T}
    kkt_iters = solve_kkt!(wrk, set, Δp, Δy, H, B, f, rp, y0; rtol)

    if itmax > 0
        kkt_iters += refine_kkt!(
            Δp, Δy, wrk, set, H, B, f, rp, sp, sy, dp, dy;
            itmax, rtol, stall,
        )
    end

    copyto!(Δd, rd)
    mul!(Δd, B', Δy, -1, 1)
    mul!(Δd, Symmetric(Q, :L), Δp, 1, 1)

    return kkt_iters
end

function initcorrector!(
        cone::AbstractCone,
        v::Integer,
        f::AbstractVector,
        caches::Caches,
        p::AbstractVector,
        d::AbstractVector,
        Δp::AbstractVector,
        Δd::AbstractVector,
        σμ::Real,
        B::BlockSparseMatrix,
        conewrk::ConeWorkspace,
    )
    r = colrange(B, v)
    cv = cache(caches, v, cone)
    corr!(view(f, r), view(p, r), view(d, r), view(Δp, r), view(Δd, r), σμ, cv, conewrk)
    return
end

function maxsteps(
        cone::AbstractCone,
        v::Integer,
        p::AbstractVector{T},
        d::AbstractVector{T},
        Δp::AbstractVector{T},
        Δd::AbstractVector{T},
        caches::Caches{T},
        B::BlockSparseMatrix{T},
        conewrk::ConeWorkspace{T};
        step_frac::T = one(T),
    ) where {T}
    r = colrange(B, v)
    τp, τd = maxsteps(view(p, r), view(Δp, r), view(d, r), view(Δd, r), cache(caches, v, cone), conewrk)

    if !(cone isa CofreeCone)
        τp *= step_frac
        τd *= step_frac
    end

    return τp, τd
end

# Steplengths in Interior-Point Algorithms of Quadratic Programming
# Meszaros
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
    #
    # compute coefficients:
    #
    #   g = Q  Δp
    #   h = Bᵀ Δy
    #
    mul!(g, Symmetric(Q, :L), Δp)
    mul!(h,           B',     Δy)

    cgg = dot( g,  g)
    chh = dot( h,  h)
    cuu = dot(Δd, Δd)

    cgh = dot( g,  h)
    cgu = dot( g, Δd)
    chu = dot( h, Δd)

    crg = dot(rd,  g)
    crh = dot(rd,  h)
    cru = dot(rd, Δd)
    #
    # candidate A: minimize
    #
    #   aA d² + bA d + cA
    #
    # over 0 ≤ d ≤ dmax, where
    #
    #   aA = ⟨Δd, Δd⟩ - ⟨h, Δd⟩² / ⟨h, h⟩
    #   bA = -2(⟨rd, Δd⟩ + pmax ⟨g, Δd⟩ - ⟨h, Δd⟩ (⟨rd, h⟩ + pmax ⟨g, h⟩) / ⟨h, h⟩)
    #
    # candidate B: minimize
    #
    #   aB p² + bB p + cB
    #
    # over 0 ≤ p ≤ pmax, where
    #
    #   aB = ⟨g, g⟩ - ⟨g, h⟩² / ⟨h, h⟩
    #   bB = 2(⟨rd, g⟩ - dmax ⟨g, Δd⟩ - ⟨g, h⟩ (⟨rd, h⟩ - dmax ⟨h, Δd⟩) / ⟨h, h⟩)
    #
    aA = cuu - chu^2 / chh
    aB = cgg - cgh^2 / chh

    bA = -2(cru + pmax * cgu - chu * (crh + pmax * cgh) / chh)
    bB =  2(crg - dmax * cgu - cgh * (crh - dmax * chu) / chh)


    tol = sqrt(eps(T))
    #
    # if aA is zero, then the reduced objective is flat in d, and the
    # the efficient-set solution dmax.
    #
    flagA = aA > tol * cuu
    
    if !flagA
        dA = dmax
    else
        dA = clamp(-bA / 2aA, zero(T), dmax)
    end
    #
    # if aB is zero, then the reduced objective is flat in p, and the
    # efficient-set solution is pmax.
    #
    flagB = aB > tol * cgg

    if !flagB
        pB = pmax
    else
        pB = clamp(-bB / 2aB, zero(T), pmax)
    end
    #
    # select the candidate with smaller dual residual
    #
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
    #
    # compute the primal floor (method 2).
    #
    pstep = max(min(pmax, dmax), pstep)
    ystep = (crh + pstep * cgh - dstep * chu) / chh

    return pstep, ystep, dstep
end

function solvepredictor!(s::IPMSolver{T}; rtol::T) where {T}
    w = s.wrk

    axpby!(-1, s.d, 0, w.f)
    axpby!(-1, w.rd, 1, w.f)

    return newton!(w.Δpa, w.Δya, w.Δda, s.kkt, s.settings.kkt, s.H, s.B, w.f, w.rp, w.rd, s.Q, w.sp, w.sy, w.dp, w.dy;
                   itmax=s.settings.refine_itmax, rtol, stall=s.settings.refine_stall)
end

function solvecorrector!(s::IPMSolver{T}, μ::T; rtol::T) where {T}
    w = s.wrk

    τpa = one(T)
    τda = one(T)

    for v in vtxs(s.B)
        τpv, τdv = maxsteps(s.K[v], v, s.p, s.d, w.Δpa, w.Δda, s.caches, s.B, s.conewrk)
        τpa = min(τpa, τpv)
        τda = min(τda, τdv)
    end

    σμ = zero(T)

    for j in cols(s.B)
        σμ += (s.p[j] + τpa * w.Δpa[j]) * (s.d[j] + τda * w.Δda[j])
    end

    σμ /= s.ν
    σμ = clamp(σμ * (σμ / μ)^2, zero(T), μ)

    for v in vtxs(s.B)
        initcorrector!(s.K[v], v, w.f, s.caches, s.p, s.d, w.Δpa, w.Δda, σμ, s.B, s.conewrk)
    end

    axpy!(-1, w.rd, w.f)

    kkt_iters = newton!(w.Δp, w.Δy, w.Δd, s.kkt, s.settings.kkt, s.H, s.B, w.f, w.rp, w.rd, s.Q, w.sp, w.sy, w.dp, w.dy, w.Δya;
                        itmax=s.settings.refine_itmax, rtol, stall=s.settings.refine_stall)

    return kkt_iters
end

function startingpoint(B::BlockSparseMatrix{T}, g::AbstractVector{T}, c::AbstractVector{T}, cones::AbstractVector) where {T}
    m, n = size(B)

    p = FVector{T}(undef, n)
    d = FVector{T}(undef, n)
    y = FVector{T}(undef, m)
    z = FVector{T}(undef, m)

    for v in vtxs(B)
        r = colrange(B, v)
        identity!(view(p, r), cones[v])
        identity!(view(d, r), cones[v])
    end

    mul!(z, B, p)

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

    fill!(y, zero(T))

    return p, d, y
end

function isstalled(hist::IPMHistory{T}; window=6, threshold=0.5) where {T}
    n = length(hist)

    if n < 2 * window
        return false
    end

    floor = √eps(T)

    prev_start = n - 2 * window + 1
    prev_end = n - window
    curr_start = n - window + 1
    curr_end = n

    μ_best_prev = minimum(view(hist.μ, prev_start:prev_end))
    μ_best_curr = minimum(view(hist.μ, curr_start:curr_end))

    if μ_best_curr < threshold * μ_best_prev
        return false
    end

    pres_best_prev = minimum(view(hist.pres, prev_start:prev_end))
    pres_best_curr = minimum(view(hist.pres, curr_start:curr_end))

    if pres_best_prev > floor && pres_best_curr < threshold * pres_best_prev
        return false
    end

    dres_best_prev = minimum(view(hist.dres, prev_start:prev_end))
    dres_best_curr = minimum(view(hist.dres, curr_start:curr_end))

    if dres_best_prev > floor && dres_best_curr < threshold * dres_best_prev
        return false
    end

    return true
end

function isstalled(s::IPMSolver)
    return isstalled(s.hist)
end

function isnearoptimal(hist::IPMHistory; feas_tol, gap_tol, near_factor)
    if isempty(hist)
        return false
    end

    μ  = hist.μ[end]
    rp = hist.pres[end]
    rd = hist.dres[end]

    return rp < near_factor * feas_tol && rd < near_factor * feas_tol && μ < near_factor * gap_tol
end

function isnearoptimal(s::IPMSolver)
    return isnearoptimal(s.hist; feas_tol=s.settings.feas_tol, gap_tol=s.settings.gap_tol, near_factor=s.settings.near_factor)
end

function isnumfail(hist::IPMHistory; window=3, threshold=1e-6)
    if length(hist.pstep) < window
        return false
    end

    τavg = sum(hist.pstep[end-window+1:end]) / window
    τavg = min(τavg, sum(hist.dstep[end-window+1:end]) / window)

    if τavg > threshold
        return false
    end

    if length(hist.pres) < window + 1
        return true
    end

    return hist.pres[end] > 0.9 * hist.pres[end - window] || hist.dres[end] > 0.9 * hist.dres[end - window]
end

function isnumfail(s::IPMSolver)
    return isnumfail(s.hist; threshold=s.settings.stall_tol)
end

function initkkt!(s::IPMSolver{T}) where {T}
    flag, s.ρ[] = initkkt!(s.kkt, s.settings.kkt, s.H, s.ρ[], s.settings.rgmax)
    return flag
end

function CommonSolve.init(prob::IPMProblem{T}; kw...) where {T}
    return init(prob, IPMSettings{T}(; kw...))
end

function CommonSolve.init(prob::IPMProblem{T, I}, settings::IPMSettings{T}) where {T, I}
    n = size(prob.B, 2)
    m = size(prob.B, 1)
    ν = conedegree(prob.K, prob.B)
    #
    # equilibrate problem data
    #
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
    #
    # initialize kkt solver
    #
    R, P, B, kkt = make_kkt(settings.kkt, B)
    #
    # permute problem data
    #
    c = P * c
    Q = halfselectvtxs(halfselectvtxs(Q, R.perm), R.perm)
    cones = tounion(prob.K, R.perm)
    #
    # compute starting point
    #
    p, d, y = startingpoint(B, g, c, cones)
    #
    # initialize per-cone caches
    #
    caches = Caches(cones, B)

    for v in vtxs(B)
        initcache!(cache(caches, v, cones[v]))
    end

    H = allocblockdiag(B)
    conewrk = ConeWorkspace{T}(cones, B)
    ipmwrk = IPMWorkspace{T}(m, n)
    hist = IPMHistory{T}()
    ρ = FScalar{T}(undef)

    ρ[] = settings.rgmin

    return IPMSolver(Q, H, B, c, g, p, d, y, cones,
        scaling, P, ipmwrk, caches, conewrk, kkt,
        hist, ν, settings, ρ, TimerOutput()
    )
end

function step!(s::IPMSolver{T}) where {T}
    status = CONTINUE

    npred = 0
    ncorr = 0
    pstep = zero(T)
    dstep = zero(T)

    w = s.wrk
    #
    # compute the primal and dual residuals:
    #
    #   rp = g - B p
    #   rd = c - d + Q p - Bᵀ y
    #
    residuals!(w.rp, w.rd, s.B, s.p, s.d, s.y, s.c, s.g, s.Q)

    if iszero(s.ν)
        μ = zero(T)
    else
        μ = dot(s.p, s.d) / s.ν
    end
    #
    # compute the forcing term
    #
    #   η = min(η₀, C μ / μ₀),
    #
    # which modifies the tolerance
    # of the KKT system solver.
    #
    if isempty(s.hist)
        μ0 = μ
    else
        μ0 = s.hist.μ[1]
    end

    if iszero(μ0)
        rtol = s.settings.forcing_ceil
    else
        rtol = min(s.settings.forcing_ceil, s.settings.forcing_frac * μ / μ0)
    end

    #
    # check the quantities
    #
    #   - ‖rp‖ / (1 + ‖g‖)
    #   - ‖rd‖ / (1 + ‖c‖)
    #   - μ
    #
    # for convergence
    #
    pres = norm(w.rp) / (1 + norm(s.g))
    dres = norm(w.rd) / (1 + norm(s.c))

    if pres < s.settings.feas_tol && dres < s.settings.feas_tol && (iszero(s.ν) || μ < s.settings.gap_tol)
        status = OPTIMAL
    else
        #
        # compute the sum
        #
        #   H = Hf(w) + Q
        #
        # where f is the barrier function and
        # w is the scaling point.
        #
        flag = true

        @timeit s.timers "scale" for v in vtxs(s.B)
            flag = flag && scale!(s.K[v], v, s.H, s.caches, s.p, s.d, s.B, s.Q, s.conewrk)
        end

        if !flag
            if s.settings.verbose
                @warn "Scaling failed."
            end

            if isnearoptimal(s)
                status = NEAR_OPTIMAL
            else
                status = NUMERICAL_FAILURE
            end
        elseif !@timeit s.timers "initkkt" initkkt!(s)
            if s.settings.verbose
                @warn "Failed to initialize KKT solver."
            end

            if isnearoptimal(s)
                status = NEAR_OPTIMAL
            else
                status = NUMERICAL_FAILURE
            end
        else
            npred = @timeit s.timers "predictor" solvepredictor!(s; rtol)
            ncorr = @timeit s.timers "corrector" solvecorrector!(s, μ; rtol)
            #
            # take a step in the direction
            #
            #   (Δp, Δy, Δd)
            #
            pstep = one(T)
            dstep = one(T)

            @timeit s.timers "maxsteps" for v in vtxs(s.B)
                τpv, τdv = maxsteps(s.K[v], v, s.p, s.d, w.Δp, w.Δd, s.caches, s.B, s.conewrk; step_frac=s.settings.step_frac)
                pstep = min(pstep, τpv)
                dstep = min(dstep, τdv)
            end

            pstep, ystep, dstep = steplengths(
                pstep, dstep, w.rd, w.Δp, w.Δy, w.Δd,
                s.Q, s.B, w.g, w.h
            )

            axpy!(pstep, w.Δp, s.p)
            axpy!(dstep, w.Δd, s.d)
            axpy!(ystep, w.Δy, s.y)

            if isstalled(s)
                if s.settings.verbose
                    @warn "Stalling detected."
                end

                if isnearoptimal(s)
                    status = NEAR_OPTIMAL
                else
                    status = STALLED
                end
            elseif isnumfail(s)
                if s.settings.verbose
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

    push!(s.hist, (; μ, pstep, dstep, pres, dres, npred, ncorr))
    return status
end

function CommonSolve.solve!(s::IPMSolver)
    status = CONTINUE; i = 0

    if s.settings.verbose
        showsettings(stdout, s.settings)
        println(stdout)
        showtop(stdout)
    end

    while status == CONTINUE
        if i + 1 ≥ s.settings.itmax
            status = ITERATION_LIMIT
        else
            status = step!(s); i += 1

            if s.settings.verbose
                showrow(stdout, i, s.hist[i])
            end
        end

        if status != CONTINUE
            break
        end
    end

    if s.settings.verbose
        showbot(stdout)
    end

    return IPMResult(s, status)
end

"""
    solve(problem::IPMProblem;
        verbose=true,
        step_frac=0.99,
        feas_tol=1e-8,
        gap_tol=1e-8,
        itmax=100,
        near_factor=1000.0,
        stall_tol=1e-6,
        forcing_ceil=0.3,
        forcing_frac=1.0,
        refine_itmax=10,
        refine_stall=0.5,
        scale_itmax=10,
        kkt=UzawaSettings(),
    )

Solve an [`IPMProblem`](@ref).
"""
function CommonSolve.solve(prob::IPMProblem{T}; kw...) where {T}
    settings = IPMSettings{T}(; kw...)
    return solve!(init(prob, settings))
end

function CommonSolve.solve(prob::IPMProblem{T}, settings::IPMSettings{T}) where {T}
    return solve!(init(prob, settings))
end

"""
    reinit!(solver, prob; frac=0.1, rgfrac=0.0)

Reinitialize an [`IPMSolver`](@ref) for a new problem with the same structure.

Only `c` and `g` may differ from the original problem; `Q`, `B`, and cones must match.

# Arguments
- `frac::Real=0.1`: blend fraction for starting point. `frac=0` uses warm start from
  previous solution; `frac=1` uses cold start.
- `rgfrac::Real=0.0`: carry fraction for regularization ρ. `rgfrac=0` resets
  to rgmin (cold); `rgfrac=1` keeps previous ρ; values in between blend:
  `ρ_new = rgfrac * ρ_old + (1 - rgfrac) * rgmin`. Carrying ρ across time-stepping
  chains can save factorization ratchet iterations when consecutive problems have
  similar conditioning.
"""
function reinit!(solver::IPMSolver{T}, prob::IPMProblem{T}; frac::Real=0.1, rgfrac::Real=0.0) where {T}
    @assert ncols(prob.B) == ncols(solver.B)
    @assert nrows(prob.B) == nrows(solver.B)
    @assert nvtxs(prob.B) == nvtxs(solver.B)
    @assert nouts(prob.B) == nouts(solver.B)
    #
    # scale c and g
    #
    c = copy(prob.c)
    g = copy(prob.g)

    for j in cols(solver.B)
        c[j] *= solver.scaling.cscl[j]
    end

    for i in rows(solver.B)
        g[i] *= solver.scaling.rscl[i]
    end
    #
    # update solver state
    #
    mul!(solver.c, solver.P, c)
    copyto!(solver.g, g)
    #
    # choose a starting point by interpolating
    # between a cold-start and the converged
    # iterate
    #
    p, d, y = startingpoint(solver.B, solver.g, solver.c, solver.K)
    axpby!(frac, p, 1 - frac, solver.p)
    axpby!(frac, d, 1 - frac, solver.d)
    axpby!(frac, y, 1 - frac, solver.y)
    #
    # clear history
    #
    empty!(solver.hist)
    #
    # update regularization: blend previous ρ with cold-start rgmin
    #
    ρ_old = solver.ρ[]
    ρ_cold = solver.settings.rgmin
    solver.ρ[] = rgfrac * ρ_old + (1 - rgfrac) * ρ_cold

    return solver
end
