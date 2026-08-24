# =============================================================================
#
# # Isotonic Regression
#
# Let G = (V, E) be a directed acyclic graph with vertices V = {1, …, n}.
# A vector x ∈ ℝⁿ is *isotonic* with respect to G if
#
#   xᵢ ≤ xⱼ for all (i, j) ∈ E.
#
# Given vertex weights vertex weights w ∈ ℝⁿ and an m × n matrix
#
#   A = [ a₁ a₂ … aₙ ]
#
# of observations, the weighted ℓ₂ isotonic regression of A is the
# projection of A onto the set of isotonic vectors:
#
#           n
#   min 1/2 ∑ wᵢ ‖ xᵢ - aᵢ ‖²   (1)
#    x     i=1
#
#   s.t. xⱼ - xᵢ - sᵢⱼ = 0,  (i, j) ∈ E
#                  sᵢⱼ ≥ 0,  (i, j) ∈ E  
#
# This is a quadratic program.
#
# # References
#
#  - Kyng, Rao & Sachdeva,
#    Fast, Provable Algorithms for Isotonic Regression in all ℓₚ-norms,
#    NIPS 2015, 2719–2727. arXiv:1507.00710.
#
# =============================================================================

using AppleAccelerate
using LinearAlgebra, Random, Printf
using Graphs

using CellularSheaves.IPM
using CellularSheaves.IPM: CofreeCone, PositiveCone, IPMSettings, HSDSettings
using CellularSheaves.BlockSparseArrays: colrange, block, blocksparse

include("../utils.jl")

# -----------------------------------------------------------------------------
# Partial orders
# -----------------------------------------------------------------------------

# A random d-regular graph, oriented by a random permutation — the random-graph
# family in the reference experiments.
function randreg_edges(n::Int, d::Int; rng = Random.default_rng())
    G = random_regular_graph(n, d; rng)
    pos = invperm(randperm(rng, n))
    E = Tuple{Int, Int}[]

    for edge in edges(G)
        u = src(edge)
        v = dst(edge)
        pos[u] < pos[v] ? push!(E, (u, v)) : push!(E, (v, u))
    end

    return E
end

# -----------------------------------------------------------------------------
# Data
# -----------------------------------------------------------------------------

# A uniformly random linear extension of the order, by Kahn's algorithm with
# random tie-breaking.
function random_linear_extension(n::Int, E::Vector{Tuple{Int, Int}}; rng = Random.default_rng())
    indeg = zeros(Int, n)
    succ = [Int[] for _ in 1:n]

    for (u, v) in E
        push!(succ[u], v)
        indeg[v] += 1
    end

    ready = [v for v in 1:n if indeg[v] == 0]
    value = zeros(Float64, n)

    for k in 1:n
        isempty(ready) && error("not a DAG")
        i = rand(rng, eachindex(ready))
        ready[i], ready[end] = ready[end], ready[i]
        v = pop!(ready)
        value[v] = k

        for u in succ[v]
            indeg[u] -= 1
            indeg[u] == 0 && push!(ready, u)
        end
    end

    return value
end

# Observations: a random linear extension perturbed by i.i.d. Gaussian noise,
# following the reference experiments, which use noise levels 1 and 10.
function isotonic_data(n::Int, E::Vector{Tuple{Int, Int}}; noise::Float64 = 1.0, seed::Int = 0)
    rng = Random.MersenneTwister(seed)
    return random_linear_extension(n, E; rng) .+ noise .* randn(rng, n)
end

# -----------------------------------------------------------------------------
# Build
# -----------------------------------------------------------------------------

function build_isotonic(n::Int, E::Vector{Tuple{Int, Int}}, a::AbstractVector;
                        weights::AbstractVector = ones(n))
    m = length(E)
    # Solved in the deviation variables δ = x - a (recover x = a + δ). This keeps
    # the objective equal to the fit residual ½δᵀWδ instead of ½xᵀWx - (Wa)ᵀx,
    # whose magnitude ∝ Σaᵢ² would dwarf the residual and make the relative gap
    # test μ < gap_tol·(1+|pobj|) stop while the gap is still wide open. Only the
    # linear term (c → 0) and the constraint RHS (g → aᵥ - aᵤ) move; Q and B, and
    # hence all sparsity and the condensed matrix W + AᵀHA, are unchanged.
    # stalks: 1..n are the δ_i, then n+e is s_e for edge e

    one11 = ones(1, 1)

    row_ids = Int[]
    col_ids = Int[]
    blks = Matrix{Float64}[]
    g = zeros(m)

    for (e, (u, v)) in enumerate(E)
        push!(row_ids, e); push!(col_ids, u);     push!(blks,  copy(one11))
        push!(row_ids, e); push!(col_ids, v);     push!(blks, -copy(one11))
        push!(row_ids, e); push!(col_ids, n + e); push!(blks,  copy(one11))
        g[e] = a[v] - a[u]                          # δ_u - δ_v + s_e = a_v - a_u
    end

    # every stalk is one-dimensional, so the block dimensions are given
    # explicitly rather than inferred — a vertex incident to no edge would
    # otherwise be assigned a zero-dimensional stalk
    B = blocksparse(row_ids, col_ids, blks, ones(Int, m), ones(Int, n + m))

    # Q: w_i on each δ_i stalk (the fidelity term), zero on the slack stalks
    Q = IPM.allocblockdiag(B)
    fill!(Q, 0)

    for i in 1:n
        block(Q, i, i, i)[1, 1] = weights[i]
    end

    # ½ δᵀWδ has no linear term; the data lives entirely in g now.
    c = zeros(size(B, 2))

    cones = vcat(IPM.AbstractCone[CofreeCone() for _ in 1:n],
                 IPM.AbstractCone[PositiveCone() for _ in 1:m])

    return IPMProblem(Q, B, c, g, cones)
end

function build_iso_jump(n::Int, E::Vector{Tuple{Int, Int}}, a::AbstractVector;
                        weights::AbstractVector, optimizer)
    model = Model(optimizer)
    @variable(model, x[1:n])

    for (u, v) in E
        @constraint(model, x[u] <= x[v])
    end

    @objective(model, Min, 0.5 * sum(weights[i] * (x[i] - a[i])^2 for i in 1:n))
    return model
end

# -----------------------------------------------------------------------------
# Benchmark
# -----------------------------------------------------------------------------

const PROBLEMS = [
    ("RandReg-1",    () -> (1000, randreg_edges(1000, 6; rng = Random.MersenneTwister(0)), 1.0)),
    ("RandReg-10",   () -> (1000, randreg_edges(1000, 6; rng = Random.MersenneTwister(0)), 10.0)),
    ("RandReg-1-lg", () -> (2500, randreg_edges(2500, 6; rng = Random.MersenneTwister(0)), 1.0)),
]

function benchmark(; tol = 1e-6, clarabel = false, mosek = false)
    @printf("%-11s %6s %6s %9s %9s %9s %9s %8s %8s\n",
            "problem", "n", "m", "IPM", "HSD", "Clarabel", "Mosek", "Cla/IPM", "Mos/IPM")

    for (name, make) in PROBLEMS
        n, E, noise = make()
        a = isotonic_data(n, E; noise)
        w = ones(n)
        prob = build_isotonic(n, E, a; weights = w)

        mi = measure_ipm(prob, IPMSettings{Float64}(feas_tol = tol, gap_tol = tol))
        mh = measure_ipm(prob, HSDSettings{Float64}(feas_tol = tol, gap_tol = tol))

        if clarabel
            mc = measure_jump(() -> build_iso_jump(n, E, a; weights = w, optimizer = clarabel_opt(; tol)))
        else
            mc = (t = NaN, status = "—", obj = NaN)
        end

        if mosek
            mk = measure_jump(() -> build_iso_jump(n, E, a; weights = w, optimizer = mosek_opt(; tol)))
        else
            mk = (t = NaN, status = "—", obj = NaN)
        end

        cratio = (isfinite(mc.t) && isfinite(mi.t)) ? @sprintf("%.2f", mc.t / mi.t) : "—"
        mratio = (isfinite(mk.t) && isfinite(mi.t)) ? @sprintf("%.2f", mk.t / mi.t) : "—"

        @printf("%-11s %6d %6d %s %s %s %s %8s %8s\n",
                name, n, length(E), fmt_time(mi), fmt_time(mh), fmt_time(mc), fmt_time(mk), cratio, mratio)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    opts = parse_args()
    tol = any(a -> startswith(a, "--tol="), ARGS) ? opts.tol : 1e-6
    benchmark(; clarabel = opts.clarabel, mosek = opts.mosek, tol)
end

# =============================================================================
# Sample run: 2026-08-24 (--clarabel --mosek), tol = 1e-6.
#
#   problem          n      m       IPM       HSD  Clarabel     Mosek  Cla/IPM  Mos/IPM
#   RandReg-1     1000   3000   34.8ms   64.7ms  135.6ms   69.2ms*     3.90     1.99
#   RandReg-10    1000   3000   38.1ms   55.4ms  135.8ms   74.7ms*     3.56     1.96
#   RandReg-1-lg  2500   7500  216.4ms  427.0ms 1904.6ms  205.1ms*     8.80     0.95
# =============================================================================
