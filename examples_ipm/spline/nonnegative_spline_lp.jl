######################################################################
# nonnegative_spline_lp.jl
#
# Shape-constrained Bézier-spline estimation — ORTHANT (LP) variant.
# The clean "cone-on-vertex" instance: the shape cone sits on the
# segment control-point stalk itself, with no reification for the
# nonnegativity case.
#
#   shape = :nonneg    ->  PositiveCone directly on each segment stalk
#                          (P_i >= 0  ⟹  f >= 0 by the convex-hull property)
#   shape = :monotone  ->  ΔP  >= 0   via a PositiveCone difference-leaf
#   shape = :convex    ->  Δ²P >= 0   via a PositiveCone difference-leaf
#
# This is a SUFFICIENT (strict-subcone) relaxation; it is a dense sieve,
# tightenable by degree elevation. The exact cone is the companion file
# nonnegative_spline_exact.jl. See nonnegative_spline.md for the theory.
#
# Objective per piece: data-fit (least squares to samples) + λ·roughness
# (trace-normalized min-ρ-derivative Bernstein Gram) + tiny ε for strict
# convexity on data-sparse pieces. Pieces are joined by C^k continuity.
#
# Standalone: the small Bernstein/jet closure this file needs is inlined
# below, so there is no dependency on the corridor benchmark files.
######################################################################

using CellularSheaves.IPM
using CellularSheaves.BlockSparseArrays: colrange, rowrange, blocksparse, block, nvtxs
using LinearAlgebra
using Random
using Printf

# ---- inlined Bernstein/jet closure (from the corridor numerics) -----

"""falling(n,r) = n·(n−1)···(n−r+1); the scalar from differentiating a degree-n Bézier r times."""
falling(n::Integer, r::Integer) = prod(n - j for j in 0:(r - 1); init = 1.0)

"""Degree-m Bernstein Gram 𝒢_{ik} = ∫₀¹ b_i b_k = C(m,i)C(m,k) / [(2m+1) C(2m,i+k)], (m+1)×(m+1)."""
function bernstein_gram(m::Integer)
    G = Matrix{Float64}(undef, m + 1, m + 1)
    for i in 0:m, k in 0:m
        G[i + 1, k + 1] = binomial(m, i) * binomial(m, k) / ((2m + 1) * binomial(2m, i + k))
    end
    return G
end

"""Endpoint r-th-derivative coefficients (length n+1): f⁽ʳ⁾(endpoint) = Σ_i w[i]·P_i. side=:left/:right."""
function endpoint_jet_coeffs(n::Integer, r::Integer, T::Real, side::Symbol)
    w = zeros(n + 1); scale = falling(n, r) / T^r
    if side === :left
        for l in 0:r; w[l + 1] = scale * (-1)^(r - l) * binomial(r, l); end
    elseif side === :right
        for l in 0:r; w[(n - r + l) + 1] = scale * (-1)^(r - l) * binomial(r, l); end
    else
        error("side must be :left or :right")
    end
    return w
end

"""Restriction map reading a segment's boundary jet (value…k-th derivative), (k+1)d × (n+1)d."""
function jet_operator(n::Integer, d::Integer, k::Integer, T::Real, side::Symbol)
    R = zeros((k + 1) * d, (n + 1) * d)
    for r in 0:k
        w = endpoint_jet_coeffs(n, r, T, side)
        for i in 0:n, a in 1:d
            R[r * d + a, i * d + a] = w[i + 1]
        end
    end
    return R
end

# ---- local Bernstein numerics (spline-specific) ---------------------

"""m-th forward-difference operator, (n+1-m)×(n+1): (Δ^m P)_i = Σ_j C(m,j)(-1)^{m-j} P_{i+j}."""
function diffop(n::Int, m::Int)
    D = zeros(n + 1 - m, n + 1)
    for i in 0:(n - m), j in 0:m
        D[i + 1, i + j + 1] = (-1)^(m - j) * binomial(m, j)
    end
    return D
end

"""Bernstein evaluation vector at t∈[0,1]: φ_i = C(n,i) t^i (1−t)^{n−i}."""
function bernstein_eval(n::Int, t::Float64)
    φ = zeros(n + 1)
    for i in 0:n
        φ[i + 1] = binomial(n, i) * t^i * (1 - t)^(n - i)
    end
    return φ
end

"""Trace-normalized minimum-ρ-derivative Bernstein Gram (T=1), (n+1)×(n+1)."""
function rough_gram(n::Int, ρ::Int)
    Δ = diffop(n, ρ)
    G = (falling(n, ρ)^2) * (Δ' * bernstein_gram(n - ρ) * Δ)
    G ./= tr(G)
    return G
end

# ---- instance -------------------------------------------------------

struct ShapeSplineInstance
    M::Int; n::Int; k::Int; ρ::Int
    shape::Symbol                      # :nonneg | :monotone | :convex
    density::Bool                      # add ∫f = 1 edge
    knots::Vector{Float64}             # length M+1, domain partition
    Qpiece::Vector{Matrix{Float64}}    # per-piece objective Gram
    cpiece::Vector{Vector{Float64}}    # per-piece linear term
    λ::Float64; ε::Float64
end

"""Map a global x to (piece index s∈1..M, local t∈[0,1])."""
function locate(knots, x)
    M = length(knots) - 1
    s = clamp(searchsortedlast(knots, x), 1, M)
    t = (x - knots[s]) / (knots[s + 1] - knots[s])
    return s, clamp(t, 0.0, 1.0)
end

"""
    generate_shape_instance(; M, n, k, ρ, shape, density, N, σ, λ, ε, seed)

Sample a shaped ground-truth function, generate noisy data, and accumulate
the per-piece data-fit + roughness objective. Domain is [0,1].
"""
function generate_shape_instance(;
        M::Int = 8, n::Int = 6, k::Int = 3, ρ::Int = 2,
        shape::Symbol = :nonneg, density::Bool = false,
        N::Int = 200, σ::Float64 = 0.05, λ::Float64 = 1e-2, ε::Float64 = 1e-8,
        seed::Int = 1)
    @assert n ≥ k + 1 && n ≥ ρ
    rng = MersenneTwister(seed)
    knots = collect(range(0.0, 1.0; length = M + 1))
    truth(x) = shape === :monotone ? 0.5 + 0.5 * tanh(6 * (x - 0.5)) :
               shape === :convex   ? (x - 0.5)^2 + 0.1 :
               exp(-((x - 0.35) / 0.15)^2) + 0.5 * exp(-((x - 0.75) / 0.10)^2)

    G = rough_gram(n, ρ)
    Qp = [λ .* G .+ ε .* Matrix{Float64}(I, n + 1, n + 1) for _ in 1:M]
    cp = [zeros(n + 1) for _ in 1:M]

    xs = sort(rand(rng, N))
    for x in xs
        y = truth(x) + σ * randn(rng)
        s, t = locate(knots, x)
        φ = bernstein_eval(n, t)
        Qp[s] .+= 2.0 .* (φ * φ')       # ½·2·φφ' ⇒ (φ'P)² in the objective
        cp[s] .-= 2.0 .* y .* φ          # −2 y φ'P
    end
    return ShapeSplineInstance(M, n, k, ρ, shape, density, knots, Qp, cp, λ, ε)
end

# ---- native model ---------------------------------------------------

"""
    build_shape_lp_problem(inst) :: (IPMProblem, info)

Segment stalks carry the objective; the shape cone is a PositiveCone on
the stalk (:nonneg) or on a difference-leaf (:monotone/:convex). Pieces
are joined by C^k continuity. Optional ∫f = 1 edge for density.
"""
function build_shape_lp_problem(inst::ShapeSplineInstance)
    n, k, M = inst.n, inst.k, inst.M
    d = 1; ns = n + 1
    L = jet_operator(n, d, k, 1.0, :left)
    R = jet_operator(n, d, k, 1.0, :right)
    mdiff = inst.shape === :monotone ? 1 : inst.shape === :convex ? 2 : 0

    row_ids = Int[]; col_ids = Int[]; blocks = Matrix{Float64}[]
    col_cone = AbstractCone[]; col_dim = Int[]
    rhs = Dict{Int, Vector{Float64}}()
    new_vertex!(dim, cone) = (push!(col_dim, dim); push!(col_cone, cone); length(col_dim))
    ec = Ref(0); new_edge!() = (ec[] += 1; ec[])
    place!(e, v, mat) = (push!(row_ids, e); push!(col_ids, v); push!(blocks, Matrix{Float64}(mat)))

    # segment vertices
    seg = Int[]
    for _ in 1:M
        cone = (inst.shape === :nonneg) ? PositiveCone() : CofreeCone()
        push!(seg, new_vertex!(ns, cone))
    end

    # shape leaves for monotone/convex: Δ^m P − leaf = 0, leaf ≥ 0
    if mdiff > 0
        D = diffop(n, mdiff)
        for s in 1:M
            lv = new_vertex!(ns - mdiff, PositiveCone()); e = new_edge!()
            place!(e, seg[s], D); place!(e, lv, -Matrix{Float64}(I, ns - mdiff, ns - mdiff))
        end
    end

    # C^k continuity between consecutive pieces
    for s in 1:(M - 1)
        e = new_edge!(); place!(e, seg[s], R); place!(e, seg[s + 1], -L)
    end

    # optional ∫ f = 1  (piece integral = width/(n+1) · Σ_i P_i)
    if inst.density
        e = new_edge!()
        for s in 1:M
            w = (inst.knots[s + 1] - inst.knots[s]) / (n + 1)
            place!(e, seg[s], reshape(fill(w, ns), 1, ns))
        end
        rhs[e] = [1.0]
    end

    # assemble
    B = blocksparse(row_ids, col_ids, blocks)
    nvtx = length(col_dim)
    @assert nvtxs(B) == nvtx "column blocks not contiguous 1..nvtx"
    g = zeros(size(B, 1)); for (e, v) in rhs; g[rowrange(B, e)] .= v; end
    Q = IPM.allocblockdiag(B); fill!(Q, 0)
    for s in 1:M; block(Q, seg[s], seg[s], seg[s]) .= inst.Qpiece[s]; end
    c = zeros(size(B, 2))
    for s in 1:M; c[colrange(B, seg[s])] .= inst.cpiece[s]; end
    cones = AbstractCone[col_cone[v] for v in 1:nvtx]

    prob = IPMProblem(Q, B, c, g, cones)
    info = (segcols = [colrange(B, seg[s]) for s in 1:M], nvtx = nvtx)
    return prob, info
end

lp_settings() = IPMSettings{Float64}(
    kkt = UzawaSettings{Float64}(raug = 1e5), feas_tol = 1e-8, gap_tol = 1e-8, itmax = 100)

"""Quick build+solve demo; reports objective and worst shape-violation of the fit."""
function shape_lp_demo(; shape = :nonneg, kwargs...)
    inst = generate_shape_instance(; shape = shape, kwargs...)
    prob, info = build_shape_lp_problem(inst)
    solve(prob, lp_settings())                         # JIT warmup
    t = @elapsed (res = solve(prob, lp_settings()))
    n = inst.n
    worst = Inf
    for s in 1:inst.M
        P = res.p[info.segcols[s]]
        v = shape === :monotone ? minimum(diffop(n,1) * P) :
            shape === :convex   ? minimum(diffop(n,2) * P) : minimum(P)
        worst = min(worst, v)
    end
    @printf("LP %-9s  M=%d n=%d k=%d  vtx=%d  %.1f ms  iters=%d  status=%s  min-coeff=%.2e\n",
            string(shape), inst.M, n, inst.k, info.nvtx, 1e3t, res.niter, res.status, worst)
    return inst, res, info
end

# =====================================================================
# JuMP reference model (for Mosek comparison)
# =====================================================================

using JuMP

"""Build the LP shape-spline as a JuMP model."""
function build_jump_lp(inst::ShapeSplineInstance, optimizer)
    n, k, M = inst.n, inst.k, inst.M
    ns = n + 1
    L = jet_operator(n, 1, k, 1.0, :left)
    R = jet_operator(n, 1, k, 1.0, :right)
    mdiff = inst.shape === :monotone ? 1 : inst.shape === :convex ? 2 : 0

    model = Model(optimizer); set_silent(model)

    # segment variables
    if inst.shape === :nonneg
        P = [@variable(model, [1:ns], lower_bound = 0.0) for _ in 1:M]
    else
        P = [@variable(model, [1:ns]) for _ in 1:M]
    end

    # objective: Σ_s ½ P_s' Q_s P_s + c_s' P_s
    obj = sum(0.5 * P[s]' * inst.Qpiece[s] * P[s] + inst.cpiece[s]' * P[s] for s in 1:M)
    @objective(model, Min, obj)

    # shape constraints for monotone/convex
    if mdiff > 0
        D = diffop(n, mdiff)
        for s in 1:M
            @constraint(model, D * P[s] .>= 0)
        end
    end

    # C^k continuity
    for s in 1:(M - 1)
        @constraint(model, R * P[s] .== L * P[s + 1])
    end

    # optional ∫ f = 1
    if inst.density
        int_coeffs = [(inst.knots[s + 1] - inst.knots[s]) / (n + 1) for s in 1:M]
        @constraint(model, sum(int_coeffs[s] * sum(P[s]) for s in 1:M) == 1.0)
    end

    return model, P
end

# =====================================================================
# Benchmark: Sheaf IPM vs Mosek
# =====================================================================

function run_lp_benchmark(; optimizer = nothing, dual_optimizer = nothing, solver_name::String = "Mosek", nwarmup::Int = 2, nruns::Int = 5)
    optimizer === nothing && error("Pass optimizer, e.g. run_lp_benchmark(optimizer=Mosek.Optimizer)")

    cases = [
        (:nonneg,   50, 6,  1e5),   # M=50 pieces, degree 6
        (:nonneg,   100, 6, 1e5),
        (:monotone, 50, 6,  1e5),
        (:monotone, 100, 6, 1e5),
        (:convex,   50, 6,  1e5),
        (:convex,   100, 6, 1e5),
    ]

    sname = rpad(solver_name, 9)
    sname_d = rpad(solver_name * "-D", 9)
    println("="^90)
    println("Shape-Spline LP Benchmark: Sheaf IPM vs $solver_name")
    println("="^90)
    println()
    if dual_optimizer !== nothing
        @printf("%-12s %6s %7s %6s %5s %5s %9s %9s %9s %7s %7s\n",
                "Shape", "raug", "splines", "degree", "|V|", "|E|", "IPM(ms)", sname, sname_d, "P/IPM", "D/IPM")
    else
        @printf("%-12s %6s %7s %6s %5s %5s %9s %9s %8s\n",
                "Shape", "raug", "splines", "degree", "|V|", "|E|", "IPM(ms)", sname, "speedup")
    end
    println("-"^90)

    for (shape, M, n, raug) in cases
        inst = generate_shape_instance(; shape, M, n, N = 20 * M)
        prob, info = build_shape_lp_problem(inst)
        nV = info.nvtx
        # |E|: continuity edges (M-1) + shape-leaf edges (M if monotone/convex, 0 if nonneg)
        nE = (M - 1) + (shape === :nonneg ? 0 : M)

        settings = IPMSettings{Float64}(
            kkt = UzawaSettings{Float64}(raug = raug),
            feas_tol = 1e-8, gap_tol = 1e-8, itmax = 100)
        # warmup
        for _ in 1:nwarmup; solve(prob, settings); end
        t_ipm = minimum([@elapsed solve(prob, settings) for _ in 1:nruns])

        for _ in 1:nwarmup
            m, _ = build_jump_lp(inst, optimizer)
            optimize!(m)
        end
        t_mosek = minimum([begin
            m, _ = build_jump_lp(inst, optimizer)
            @elapsed optimize!(m)
        end for _ in 1:nruns])

        if dual_optimizer !== nothing
            for _ in 1:nwarmup
                m, _ = build_jump_lp(inst, dual_optimizer)
                optimize!(m)
            end
            t_dual = minimum([begin
                m, _ = build_jump_lp(inst, dual_optimizer)
                @elapsed optimize!(m)
            end for _ in 1:nruns])

            @printf("%-12s %6.0e %7d %6d %5d %5d %9.1f %9.1f %9.1f %6.2fx %6.2fx\n",
                    string(shape), raug, M, n, nV, nE, t_ipm * 1000, t_mosek * 1000, t_dual * 1000,
                    t_mosek / t_ipm, t_dual / t_ipm)
        else
            @printf("%-12s %6.0e %7d %6d %5d %5d %9.1f %9.1f %7.2fx\n",
                    string(shape), raug, M, n, nV, nE, t_ipm * 1000, t_mosek * 1000,
                    t_mosek / t_ipm)
        end
    end
end

# =====================================================================
# Script entry point
# =====================================================================

if abspath(PROGRAM_FILE) == @__FILE__
    include(joinpath(@__DIR__, "..", "benchmark_utils.jl"))

    opts = parse_benchmark_args(ARGS)

    if opts.mosek
        using MosekTools
    else
        using Clarabel
    end
    optimizer, dual_optimizer = get_optimizers(opts)
    solver_name = opts.mosek ? "Mosek" : "Clarabel"
    print_benchmark_config(opts; lp_only = true)

    run_lp_benchmark(; optimizer, dual_optimizer, solver_name,
                     nwarmup = opts.nwarmup, nruns = opts.nruns)
end
