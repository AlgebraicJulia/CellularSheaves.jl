######################################################################
# mle_spline.jl
#
# Maximum-likelihood spline estimation — the corpus's first genuinely
# NON-QUADRATIC objective, exercising ExponentialCone leaves.
#
# Two likelihoods, one assembly (mode ∈ {:density, :intensity}):
#
#   :density    minimize  −Σ_j w_j log f(x_j)  + ½λr Σ_s PᵀGP
#               s.t.      C^k continuity,  ∫f = 1 (pin edge),  P ≥ 0
#
#   :intensity  minimize  ∫λ − Σ_j w_j log λ(t_j)  + ½λr Σ_s PᵀGP
#               s.t.      C^k continuity,  P ≥ 0
#               (inhomogeneous-Poisson log-likelihood; ∫λ is LINEAR in
#                the control points, so it rides in c — no pin needed)
#
# Each observation x_j hangs one 3-dim ExponentialCone leaf (u,y,t) off
# the segment containing it, wired by a single 2-row edge:
#
#       φ(x_j)ᵀ P − u = 0        (u is the spline value at x_j)
#       y = 1                    (the homogenizer pin)
#
# with objective coefficient −w_j on t. The library's cone is
#     ExponentialCone = { (x,y,z) : x,y > 0,  y·log(x/y) ≥ z },
# so with y pinned to 1 the leaf enforces t ≤ log u, and minimizing
# −Σ w t drives t to log f(x_j) at the optimum (epigraph tight). This
# is duality.md §5 verbatim: Move 1 epigraph scalar t, Move 2
# homogenizer y = 1, pins folded into the co-boundary through g.
# NOTE the coordinate order (u, y, t) — value first, bound LAST — the
# reverse of MOSEK's exponential cone; see mle_spline.md.
#
# The segment stalks keep the PositiveCone: the exp leaves certify
# positivity only AT the observations, the orthant certifies f ≥ 0
# everywhere (and remains the shape prior). Q carries the roughness
# Gram, so one problem mixes Q + orthant + exponential cones.
#
# Weights w_j support binned data (see `binned`): leaves scale with the
# number of DISTINCT observation sites, not the raw sample size.
#
# Literature: Papp & Alizadeh, Shape-Constrained Estimation Using
# Nonnegative Splines (JCGS 2014); Alizadeh & Papp, arrival-rate
# estimation by ML nonnegative splines (Ann. Oper. Res. 2013).
#
# NOTE: written against the CellularSheaves.IPM PR-67 API; not executed
# here (no Julia toolchain). mle_spline_oracle.py is the numerical
# ground truth for every value quoted in mle_spline.md.
######################################################################

include("nonnegative_spline_lp.jl")   # Bernstein/jet numerics + IPM imports

# ---- samplers (no Distributions.jl dependency) -----------------------

"""Poisson(λ) via cumulative exponentials (log-space Knuth; safe for large λ)."""
function rand_poisson(rng, λ::Float64)
    k = -1; t = 0.0
    while t < λ
        t += -log(rand(rng)); k += 1
    end
    return k
end

# ---- instance --------------------------------------------------------

struct MLESplineInstance
    M::Int; n::Int; k::Int; ρ::Int
    mode::Symbol                  # :density | :intensity
    knots::Vector{Float64}        # length M+1, equal-width on [0,1]
    obs::Vector{Float64}          # observation locations
    w::Vector{Float64}            # per-observation weights (binned counts)
    λr::Float64; ε::Float64       # roughness weight, strict-convexity floor
end

"""
    generate_mle_instance(; mode, M, n, k, ρ, N, λr, ε, seed)

Sample observations from a bimodal ground truth on [0,1]:
:density draws N points by rejection from the normalized shape;
:intensity draws an inhomogeneous Poisson process with E[#events] = N
by thinning. Same shape function as the regression instances.
"""
function generate_mle_instance(;
        mode::Symbol = :density, M::Int = 6, n::Int = 4, k::Int = 2, ρ::Int = 2,
        N::Int = 150, λr::Float64 = 1e-3, ε::Float64 = 1e-8, seed::Int = 1)
    @assert n ≥ k + 1 && n ≥ ρ
    @assert mode in (:density, :intensity)
    rng = MersenneTwister(seed)
    knots = collect(range(0.0, 1.0; length = M + 1))
    shape(x) = exp(-((x - 0.35) / 0.15)^2) + 0.5 * exp(-((x - 0.75) / 0.10)^2)

    grid = 0.0:1e-4:1.0
    smax = maximum(shape, grid)
    sint = 1e-4 * (sum(shape, grid) - (shape(0.0) + shape(1.0)) / 2)   # trapezoid ∫shape

    obs = Float64[]
    if mode === :density
        while length(obs) < N                     # rejection sampling
            x = rand(rng)
            rand(rng) * smax ≤ shape(x) && push!(obs, x)
        end
    else                                          # thinning at rate λmax
        λmax = N * smax / sint
        for _ in 1:rand_poisson(rng, λmax)
            x = rand(rng)
            rand(rng) * λmax ≤ N * shape(x) / sint && push!(obs, x)
        end
    end
    sort!(obs)
    return MLESplineInstance(M, n, k, ρ, mode, knots, obs,
                             ones(length(obs)), λr, ε)
end

"""
    binned(inst; nbins) :: MLESplineInstance

Coalesce observations into nbins equal-width bins: one leaf per occupied
bin at its midpoint, weight = count. An approximation of the likelihood
(evaluation moves to midpoints) that caps the leaf count at nbins.
"""
function binned(inst::MLESplineInstance; nbins::Int = 64)
    edges = collect(range(0.0, 1.0; length = nbins + 1))
    mids = Float64[]; w = Float64[]
    for b in 1:nbins
        lo, hi = edges[b], edges[b + 1]
        cnt = count(x -> lo ≤ x < hi, inst.obs) +
              (b == nbins ? count(==(1.0), inst.obs) : 0)
        cnt > 0 && (push!(mids, (lo + hi) / 2); push!(w, Float64(cnt)))
    end
    return MLESplineInstance(inst.M, inst.n, inst.k, inst.ρ, inst.mode,
                             inst.knots, mids, w, inst.λr, inst.ε)
end

# ---- native model ----------------------------------------------------

"""
    build_mle_problem(inst) :: (IPMProblem, info)

Segment stalks: PositiveCone, roughness Gram in Q, (intensity) ∫λ in c.
One ExponentialCone leaf (u,y,t) per observation, wired by a 2-row edge
[φᵀP − u = 0; y = 1] with rhs (0,1); objective −w_j on t. C^k continuity
on the spine; (density) one ∫f = 1 pin edge.
"""
function build_mle_problem(inst::MLESplineInstance)
    n, k, M = inst.n, inst.k, inst.M
    ns = n + 1
    L = jet_operator(n, 1, k, 1.0, :left)
    R = jet_operator(n, 1, k, 1.0, :right)

    row_ids = Int[]; col_ids = Int[]; blocks = Matrix{Float64}[]
    col_cone = AbstractCone[]; col_dim = Int[]
    rhs = Dict{Int, Vector{Float64}}()
    new_vertex!(dim, cone) = (push!(col_dim, dim); push!(col_cone, cone); length(col_dim))
    ec = Ref(0); new_edge!() = (ec[] += 1; ec[])
    place!(e, v, mat) = (push!(row_ids, e); push!(col_ids, v); push!(blocks, Matrix{Float64}(mat)))

    # segment vertices (orthant: f ≥ 0 everywhere, not just at the data)
    seg = [new_vertex!(ns, PositiveCone()) for _ in 1:M]

    # observation leaves: stalk (u, y, t) ∈ ExponentialCone, y log(u/y) ≥ t
    leaf = Int[]; obsmap = Tuple{Int, Vector{Float64}}[]
    for x in inst.obs
        s, tloc = locate(inst.knots, x)
        φ = bernstein_eval(n, tloc)
        push!(obsmap, (s, φ))
        lv = new_vertex!(3, ExponentialCone()); push!(leaf, lv)
        e = new_edge!()
        place!(e, seg[s], vcat(φ', zeros(1, ns)))          # φᵀP − u = 0
        place!(e, lv, [-1.0 0.0 0.0; 0.0 1.0 0.0])         # …and  y = 1
        rhs[e] = [0.0, 1.0]
    end

    # C^k continuity between consecutive pieces
    for s in 1:(M - 1)
        e = new_edge!(); place!(e, seg[s], R); place!(e, seg[s + 1], -L)
    end

    # density: ∫ f = 1  (piece integral = width/(n+1) · Σ_i P_i)
    if inst.mode === :density
        e = new_edge!()
        for s in 1:M
            wint = (inst.knots[s + 1] - inst.knots[s]) / (n + 1)
            place!(e, seg[s], reshape(fill(wint, ns), 1, ns))
        end
        rhs[e] = [1.0]
    end

    # assemble
    B = blocksparse(row_ids, col_ids, blocks)
    nvtx = length(col_dim)
    @assert nvtxs(B) == nvtx "column blocks not contiguous 1..nvtx"
    g = zeros(size(B, 1)); for (e, v) in rhs; g[rowrange(B, e)] .= v; end

    # Intensity is emitted NORMALIZED: λ = s·f with s = Σw. Exactly
    # equivalent (F_raw = F_scaled − (Σw)·log s) and it brings the exp-cone
    # arguments to O(1); the raw form stalled Mosek in the field. See
    # mle_spline.md §3.1.
    sc = inst.mode === :intensity ? sum(inst.w) : 1.0
    Q = IPM.allocblockdiag(B); fill!(Q, 0)
    G = rough_gram(n, inst.ρ)
    Qs = sc^2 .* (inst.λr .* G .+ inst.ε .* Matrix{Float64}(I, ns, ns))
    for s in 1:M; block(Q, seg[s], seg[s], seg[s]) .= Qs; end

    c = zeros(size(B, 2))
    if inst.mode === :intensity                            # +s·∫f (linear)
        for s in 1:M
            wint = (inst.knots[s + 1] - inst.knots[s]) / (n + 1)
            c[colrange(B, seg[s])] .+= sc * wint
        end
    end
    for (j, lv) in enumerate(leaf)                         # −Σ w_j t_j
        c[colrange(B, lv)[3]] = -inst.w[j]
    end

    cones = AbstractCone[col_cone[v] for v in 1:nvtx]
    prob = IPMProblem(Q, B, c, g, cones)
    info = (segcols = [colrange(B, seg[s]) for s in 1:M],
            leafcols = [colrange(B, lv) for lv in leaf],
            obsmap = obsmap, nvtx = nvtx, sc = sc)
    return prob, info
end

mle_settings() = IPMSettings{Float64}(
    kkt = UzawaSettings{Float64}(raug = 1e5), feas_tol = 1e-8, gap_tol = 1e-8, itmax = 200)

"""Build+solve demo; reports likelihood, epigraph tightness, min coeff, ∫f."""
function mle_demo(; mode = :density, kwargs...)
    inst = generate_mle_instance(; mode, kwargs...)
    prob, info = build_mle_problem(inst)
    solve(prob, mle_settings())                            # JIT warmup
    t = @elapsed (res = solve(prob, mle_settings()))

    Pseg(s) = res.p[info.segcols[s]]
    tight = maximum(abs(res.p[info.leafcols[j]][3] -
                        log(dot(info.obsmap[j][2], Pseg(info.obsmap[j][1]))))
                    for j in eachindex(inst.obs))
    mincoef = minimum(minimum(Pseg(s)) for s in 1:inst.M)
    ∫f = sum((inst.knots[s + 1] - inst.knots[s]) / (inst.n + 1) * sum(Pseg(s))
             for s in 1:inst.M)
    nll = -sum(inst.w[j] * log(dot(info.obsmap[j][2], Pseg(info.obsmap[j][1])))
               for j in eachindex(inst.obs))
    @printf("MLE %-9s  M=%d n=%d N=%d  vtx=%d  %.1f ms  iters=%d  status=%s\n",
            string(mode), inst.M, inst.n, length(inst.obs), info.nvtx,
            1e3t, res.niter, res.status)
    @printf("    −loglik=%.6f  epigraph-tight=%.2e  min-coeff=%.2e  ∫f=%.6f\n",
            nll, tight, mincoef, ∫f)
    inst.mode === :intensity &&
        @printf("    normalized: λ = %.0f·f  (∫f ≈ 1 at the MLE; raw obj = scaled − %.0f·log %.0f)\n",
                info.sc, sum(inst.w), info.sc)
    return inst, res, info
end

# =====================================================================
# JuMP reference model (Mosek exponential-cone comparison)
# =====================================================================

"""Build the MLE spline as a JuMP model. NOTE MOI's exponential cone is
(x,y,z): y·e^{x/y} ≤ z, so t ≤ log u reads [t, 1, u] — bound FIRST,
value LAST: the reverse of the sheaf leaf's (u, y, t)."""
function build_jump_mle(inst::MLESplineInstance, optimizer)
    n, k, M = inst.n, inst.k, inst.M
    ns = n + 1
    L = jet_operator(n, 1, k, 1.0, :left)
    R = jet_operator(n, 1, k, 1.0, :right)
    G = rough_gram(n, inst.ρ)
    sc = inst.mode === :intensity ? sum(inst.w) : 1.0     # normalized intensity
    Qs = sc^2 .* (inst.λr .* G .+ inst.ε .* Matrix{Float64}(I, ns, ns))
    N = length(inst.obs)

    model = Model(optimizer); set_silent(model)
    P = [@variable(model, [1:ns], lower_bound = 0.0) for _ in 1:M]
    @variable(model, t[1:N])

    for (j, x) in enumerate(inst.obs)
        s, tloc = locate(inst.knots, x)
        φ = bernstein_eval(n, tloc)
        @constraint(model, [t[j], 1.0, φ' * P[s]] in MOI.ExponentialCone())
    end
    for s in 1:(M - 1)
        @constraint(model, R * P[s] .== L * P[s + 1])
    end

    intw = [(inst.knots[s + 1] - inst.knots[s]) / (n + 1) for s in 1:M]
    if inst.mode === :density
        @constraint(model, sum(intw[s] * sum(P[s]) for s in 1:M) == 1.0)
    end
    obj = sum(0.5 * P[s]' * Qs * P[s] for s in 1:M) - sum(inst.w[j] * t[j] for j in 1:N)
    if inst.mode === :intensity
        obj += sc * sum(intw[s] * sum(P[s]) for s in 1:M)
    end
    @objective(model, Min, obj)
    return model, P, t
end

# =====================================================================
# Benchmark: Sheaf IPM vs Mosek (exp-cone path)
# =====================================================================

function run_mle_benchmark(; optimizer = nothing, dual_optimizer = nothing, solver_name::String = "Mosek", nwarmup::Int = 2, nruns::Int = 5)
    optimizer === nothing && error("Pass optimizer, e.g. run_mle_benchmark(optimizer=Mosek.Optimizer)")
    # (mode, M, n, N, nbins, raug, tol) — tuned per case from sweep
    # intensity uses tol=1e-6 (dual stalls at ~1e-6 anyway)
    cases = [
        (:density,   8, 4, 200,  0,  1e5, 1e-8),
        (:density,   16, 4, 1000, 0,  1e5, 1e-8),
        (:density,   16, 4, 1000, 64, 1e5, 1e-8),
        (:intensity, 8, 4, 200,  0,  1e3, 1e-6),
        (:intensity, 16, 4, 1000, 0,  1e4, 1e-6),
        (:intensity, 16, 4, 1000, 64, 1e1, 1e-6),
    ]
    sname = rpad(solver_name, 9)
    sname_d = rpad(solver_name * "-D", 9)
    println("="^95)
    println("MLE-Spline Benchmark (exp-cone leaves): Sheaf IPM vs $solver_name")
    println("="^95)
    if dual_optimizer !== nothing
        @printf("%-11s %6s %4s %3s %6s %6s %6s %9s %9s %9s %7s %7s\n",
                "Mode", "raug", "M", "n", "N", "leaves", "|V|", "IPM(ms)", sname, sname_d, "P/IPM", "D/IPM")
    else
        @printf("%-11s %6s %4s %3s %6s %6s %6s %9s %9s %8s\n",
                "Mode", "raug", "M", "n", "N", "leaves", "|V|", "IPM(ms)", sname, "speedup")
    end
    println("-"^95)
    for (mode, M, n, N, nbins, raug, tol) in cases
        inst = generate_mle_instance(; mode, M, n, N)
        nbins > 0 && (inst = binned(inst; nbins))
        prob, info = build_mle_problem(inst)
        settings = IPMSettings{Float64}(
            kkt = UzawaSettings{Float64}(raug = raug),
            feas_tol = tol, gap_tol = tol, itmax = 200)
        for _ in 1:nwarmup; solve(prob, settings); end
        t_ipm = minimum([@elapsed solve(prob, settings) for _ in 1:nruns])
        for _ in 1:nwarmup
            m, _, _ = build_jump_mle(inst, optimizer); optimize!(m)
        end
        t_mosek = minimum([begin
            m, _, _ = build_jump_mle(inst, optimizer)
            @elapsed optimize!(m)
        end for _ in 1:nruns])

        if dual_optimizer !== nothing
            for _ in 1:nwarmup
                m, _, _ = build_jump_mle(inst, dual_optimizer); optimize!(m)
            end
            t_dual = minimum([begin
                m, _, _ = build_jump_mle(inst, dual_optimizer)
                @elapsed optimize!(m)
            end for _ in 1:nruns])
            @printf("%-11s %6.0e %4d %3d %6d %6d %6d %9.1f %9.1f %9.1f %6.2fx %6.2fx\n",
                    string(mode), raug, M, n, N, length(inst.obs), info.nvtx,
                    t_ipm * 1000, t_mosek * 1000, t_dual * 1000,
                    t_mosek / t_ipm, t_dual / t_ipm)
        else
            @printf("%-11s %6.0e %4d %3d %6d %6d %6d %9.1f %9.1f %7.2fx\n",
                    string(mode), raug, M, n, N, length(inst.obs), info.nvtx,
                    t_ipm * 1000, t_mosek * 1000, t_mosek / t_ipm)
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    include(joinpath(@__DIR__, "..", "benchmark_utils.jl"))
    opts = parse_benchmark_args(ARGS)

    if opts.mosek
        using MosekTools
    else
        using Clarabel
    end
    optimizer, dual_optimizer = get_optimizers(opts; lp_only = false)
    if !opts.mosek
        dual_optimizer = nothing  # skip dualization (ExponentialCone not supported)
    end
    solver_name = opts.mosek ? "Mosek" : "Clarabel"
    print_benchmark_config(opts; lp_only = false)

    run_mle_benchmark(; optimizer, dual_optimizer, solver_name,
                      nwarmup = opts.nwarmup, nruns = opts.nruns)
end
