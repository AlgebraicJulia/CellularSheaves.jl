######################################################################
# num_sheaf.jl
#
# NETWORK UTILITY MAXIMIZATION (Kelly) as a cellular sheaf on the
# source–link bipartite graph — the example where dual decomposition of
# the sheaf IS a deployed protocol: the capacity duals λ_l are link
# congestion prices, and solving the dual by gradient ascent is
# literally distributed TCP congestion control (Kelly–Maulloo–Tan 1998;
# Low–Lapsley 1999; Chiang–Low–Calderbank–Doyle, "layering as
# optimization decomposition", Proc. IEEE 2007).
#
#   max Σ_s U_s(x_s)  s.t.  Σ_{s: l∈path(s)} x_s + σ_l = c_l,  x, σ ≥ 0
#
# THE FAIRNESS LADDER IS A CONE LADDER (Mo–Walrand α-fairness):
#   :throughput   α=0  U = w·x       LP          starves long paths
#   :proportional α=1  U = w·log x   EXP leaves  (u,1,t): t ≤ log u —
#                                     the mle_spline pattern verbatim
#   :delay        α=2  U = −w/x      SOC leaves  q·x ≥ w ⟺
#                                     (q+x, q−x, 2√w) ∈ Q³ (bound first)
#   :maxmin       α=∞  shared floor  LP (Cofree t cell + slack leaves)
#
# MEASURED (num_oracle.py, CLARABEL; Kelly's linear network, 1 long
# flow over L=3 unit-capacity links + 3 shorts — closed forms at every
# rung):
#   x_long: throughput 0, proportional 1/(L+1) = 0.25 (2.98e-6),
#           delay 1/(1+√L) = 0.366025 (9.3e-6), maxmin 0.5 (1.5e-10);
#   assembled sheaf ≡ direct model in VALUE at every rung (≤ 5.8e-9;
#   the LP rungs have non-unique rates — value is the invariant).
#   KELLY'S IDENTITY from the duals: w_s/x_s = Σ_{l∈path} λ_l to
#   9.3e-5, with λ = 4/3 = (L+1)/L analytic; α=2 version w/x² = path
#   price to 1.4e-3. Jain's index on a random mesh: throughput 0.597 <
#   {delay 0.683, proportional 0.681} < maxmin 0.787 — measured NOT
#   monotone between α=1 and 2 (Jain ≠ the α-fair order). H¹ = 0.
#
# Written against the CellularSheaves.IPM PR-67 API; not executed here.
# num_oracle.py is the numerical ground truth.
######################################################################

using CellularSheaves.IPM
using CellularSheaves.BlockSparseArrays: colrange, rowrange, blocksparse, block, nvtxs
using LinearAlgebra
using Printf
using Random

# ---- instance -------------------------------------------------------------

struct NUMInstance
    routes::Vector{Vector{Int}}      # links used by each source
    cap::Vector{Float64}
    w::Vector{Float64}
    mode::Symbol                     # :throughput | :proportional | :delay | :maxmin
    ε::Float64
end

num_instance(routes, cap; w = ones(length(routes)), mode = :proportional,
             ε = 1e-9) = NUMInstance(routes, cap, w, mode, ε)

"""Kelly's linear network: source 1 crosses all L links; sources
2..L+1 each use one link. Closed forms per mode in the header."""
kelly_linear(; L = 3, mode = :proportional) =
    num_instance(vcat([collect(1:L)], [[l] for l in 1:L]), ones(L); mode)

function random_mesh(; ns = 8, nl = 6, seed = 0, mode = :proportional)
    rng = MersenneTwister(seed)
    routes = [sort(shuffle(rng, collect(1:nl))[1:rand(rng, 1:3)]) for _ in 1:ns]
    num_instance(routes, 0.5 .+ rand(rng, nl); mode)
end

# ---- builder -----------------------------------------------------------------

"""Cells: source rates x_s (PositiveCone, dim 1), link slacks σ_l
(PositiveCone), plus per-mode utility leaves. Edges: one capacity
hyperedge per link (capacities in g) touching its sources' cells and
the slack; leaf wiring per mode. The capacity duals in res.y are the
LINK PRICES."""
function build_num_problem(inst::NUMInstance)
    ns, nl = length(inst.routes), length(inst.cap)
    mode = inst.mode
    vσ(l) = ns + l
    nv = ns + nl +
         (mode === :proportional ? ns :
          mode === :delay ? 2ns :
          mode === :maxmin ? 1 + ns : 0)
    row_ids = Int[]; col_ids = Int[]; blks = Matrix{Float64}[]
    place!(e, v, m) = (push!(row_ids, e); push!(col_ids, v); push!(blks, Matrix{Float64}(m)))
    ec = Ref(0); gval = Dict{Int, Vector{Float64}}()

    for l in 1:nl                                       # capacity hyperedges
        id = (ec[] += 1)
        for (s, path) in enumerate(inst.routes)
            l in path && place!(id, s, ones(1, 1))
        end
        place!(id, vσ(l), ones(1, 1))
        gval[id] = [inst.cap[l]]
    end
    if mode === :proportional                           # exp leaf (u, 1, t)
        for s in 1:ns
            v = ns + nl + s
            id = (ec[] += 1)
            place!(id, v, [1.0 0.0 0.0; 0.0 1.0 0.0])
            place!(id, s, [-1.0; 0.0][:, :])
            gval[id] = [0.0, 1.0]
        end
    elseif mode === :delay                              # q cell + SOC leaf
        for s in 1:ns
            q = ns + nl + s
            v = ns + nl + ns + s
            id = (ec[] += 1)
            place!(id, v, Matrix{Float64}(I, 3, 3))
            place!(id, q, [-1.0; -1.0; 0.0][:, :])
            place!(id, s, [-1.0; 1.0; 0.0][:, :])
            gval[id] = [0.0, 0.0, 2sqrt(inst.w[s])]
        end
    elseif mode === :maxmin                             # floor + slack leaves
        tcell = ns + nl + 1
        for s in 1:ns
            id = (ec[] += 1)
            place!(id, s, ones(1, 1))
            place!(id, tcell, -ones(1, 1))
            place!(id, ns + nl + 1 + s, -ones(1, 1))
            gval[id] = [0.0]
        end
    end

    B = blocksparse(row_ids, col_ids, blks)
    g = zeros(size(B, 1))
    for (e, v) in gval
        g[rowrange(B, e)] .= v
    end
    @assert nvtxs(B) == nv
    Q = IPM.allocblockdiag(B); fill!(Q, 0)
    c = zeros(size(B, 2))
    for v in 1:nv
        d = length(colrange(B, v))
        block(Q, v, v, v) .= inst.ε .* Matrix{Float64}(I, d, d)
    end
    if mode === :throughput
        for s in 1:ns; c[colrange(B, s)[1]] = -inst.w[s]; end
    elseif mode === :proportional
        for s in 1:ns; c[colrange(B, ns + nl + s)[3]] = -inst.w[s]; end
    elseif mode === :delay
        for s in 1:ns; c[colrange(B, ns + nl + s)[1]] = 1.0; end
    else
        c[colrange(B, ns + nl + 1)[1]] = -1.0
    end
    cones = AbstractCone[]
    for v in 1:nv
        if v ≤ ns + nl
            push!(cones, PositiveCone())
        elseif mode === :proportional
            push!(cones, ExponentialCone())
        elseif mode === :delay
            push!(cones, v ≤ ns + nl + ns ? CofreeCone() : IPM.SecondOrderCone())
        else
            push!(cones, v == ns + nl + 1 ? CofreeCone() : PositiveCone())
        end
    end
    prob = IPMProblem(Q, B, c, g, cones)
    info = (xcols = Dict(s => colrange(B, s) for s in 1:ns), ns = ns, nl = nl,
            nv = nv, h1 = 0)
    return prob, info
end

num_rates(res, info) = [res.p[info.xcols[s][1]] for s in 1:info.ns]

"""Link congestion prices from the capacity-hyperedge duals (rows l in
edge order). Kelly's identity check: w/x vs path price."""
function num_prices(inst::NUMInstance, prob, res)
    [res.y[rowrange(prob.B, l)[1]] for l in 1:length(inst.cap)]
end

# ---- demos ----------------------------------------------------------------------

num_settings() = IPMSettings{Float64}(kkt = UzawaSettings{Float64}(raug = 1e5))

"""The fairness–cone ladder on Kelly's linear network, all four closed
forms: 0, 1/(L+1), 1/(1+√L), 1/2."""
function kelly_ladder_demo(; L = 3)
    forms = Dict(:throughput => 0.0, :proportional => 1 / (L + 1),
                 :delay => 1 / (1 + sqrt(L)), :maxmin => 0.5)
    for mode in (:throughput, :proportional, :delay, :maxmin)
        inst = kelly_linear(; L, mode)
        prob, info = build_num_problem(inst)
        res = solve(prob, num_settings())
        x0 = num_rates(res, info)[1]
        @printf("  %-13s x_long = %.6f   closed form %.6f   it=%d %s\n",
                mode, x0, forms[mode], res.niter, res.status)
    end
end

"""Kelly's identity — the dual is TCP's congestion signal: at the
proportional-fair optimum, w_s/x_s equals the sum of link prices on
the source's path (oracle: 9.3e-5; λ = (L+1)/L analytic)."""
function price_demo(; L = 3)
    inst = kelly_linear(; L, mode = :proportional)
    prob, info = build_num_problem(inst)
    res = solve(prob, num_settings())
    x = num_rates(res, info)
    λ = abs.(num_prices(inst, prob, res))    # sign per library convention
    for (s, path) in enumerate(inst.routes)
        @printf("  source %d: w/x = %.5f   path price = %.5f\n",
                s, inst.w[s] / x[s], sum(λ[path]))
    end
end

function mesh_demo(; seed = 0)
    for mode in (:throughput, :delay, :proportional, :maxmin)
        inst = random_mesh(; seed, mode)
        prob, info = build_num_problem(inst)
        res = solve(prob, num_settings())
        x = max.(num_rates(res, info), 1e-12)
        J = sum(x)^2 / (length(x) * sum(abs2, x))
        @printf("  %-13s Jain = %.4f   it=%d %s\n", mode, J, res.niter, res.status)
    end
    println("  (oracle: 0.597 / 0.683 / 0.681 / 0.787 — not monotone α=1→2)")
end

# =====================================================================
# JuMP reference (:proportional)
# =====================================================================

using JuMP

function build_jump_num(inst::NUMInstance, optimizer)
    ns, nl = length(inst.routes), length(inst.cap)
    model = Model(optimizer); set_silent(model)
    @variable(model, x[1:ns] >= 0)
    @variable(model, t[1:ns])
    for s in 1:ns
        @constraint(model, [t[s], 1.0, x[s]] in MOI.ExponentialCone())
    end
    for l in 1:nl
        @constraint(model, sum(x[s] for s in 1:ns if l in inst.routes[s]) <= inst.cap[l])
    end
    @objective(model, Max, sum(inst.w[s] * t[s] for s in 1:ns))
    return model, x
end

# =====================================================================
# Benchmark: Sheaf IPM vs Mosek
# =====================================================================

"""Benchmark IPM vs Mosek on NUM problems across the fairness ladder."""
function run_num_benchmark(; optimizer = nothing, dual_optimizer = nothing,
                           solver_name::String = "Mosek", nwarmup::Int = 2, nruns::Int = 5)
    optimizer === nothing && error("Pass optimizer")

    # (inst_fn, label, raug) - raug around 1e1-1e4 works best
    cases = [
        # Kelly's linear network (L=3)
        (() -> kelly_linear(; L = 3, mode = :throughput), "kelly L=3 tput", 1e3),
        (() -> kelly_linear(; L = 3, mode = :proportional), "kelly L=3 prop", 1e1),
        (() -> kelly_linear(; L = 3, mode = :delay), "kelly L=3 delay", 1e2),
        (() -> kelly_linear(; L = 3, mode = :maxmin), "kelly L=3 maxmin", 1e3),
        # Larger Kelly
        (() -> kelly_linear(; L = 10, mode = :proportional), "kelly L=10 prop", 1e1),
        (() -> kelly_linear(; L = 20, mode = :proportional), "kelly L=20 prop", 1e1),
        # Random mesh
        (() -> random_mesh(; ns = 20, nl = 10, mode = :proportional), "mesh 20×10 prop", 1e2),
        (() -> random_mesh(; ns = 50, nl = 25, mode = :proportional), "mesh 50×25 prop", 1e2),
        (() -> random_mesh(; ns = 100, nl = 50, mode = :proportional), "mesh 100×50 prop", 1e1),
    ]

    sname = rpad(solver_name, 9)
    sname_d = rpad(solver_name * "-D", 9)
    println("="^100)
    println("Network Utility Maximization Benchmark: Sheaf IPM vs $solver_name")
    println("="^100)
    if dual_optimizer !== nothing
        @printf("%-18s %6s %5s %4s %5s %9s %9s %9s %7s %7s\n",
                "Config", "raug", "DOF", "|V|", "H1", "IPM(ms)", sname, sname_d, "P/IPM", "D/IPM")
    else
        @printf("%-18s %6s %5s %4s %5s %9s %9s %8s\n",
                "Config", "raug", "DOF", "|V|", "H1", "IPM(ms)", sname, "speedup")
    end
    println("-"^100)

    for (inst_fn, label, raug) in cases
        inst = inst_fn()
        # Skip non-proportional modes for JuMP (only proportional reference implemented)
        if inst.mode !== :proportional
            prob, info = build_num_problem(inst)
            settings = IPMSettings{Float64}(
                kkt = UzawaSettings{Float64}(raug = raug),
                feas_tol = 1e-8, gap_tol = 1e-8, itmax = 200)
            for _ in 1:nwarmup; solve(prob, settings); end
            t_ipm = minimum([@elapsed solve(prob, settings) for _ in 1:nruns])
            @printf("%-18s %6.0e %5d %4d %5d %9.1f %9s %9s %7s %7s\n",
                    label, raug, size(prob.B, 2), info.nv, info.h1,
                    t_ipm * 1000, "—", "—", "—", "—")
            continue
        end

        prob, info = build_num_problem(inst)
        settings = IPMSettings{Float64}(
            kkt = UzawaSettings{Float64}(raug = raug),
            feas_tol = 1e-8, gap_tol = 1e-8, itmax = 200)
        for _ in 1:nwarmup; solve(prob, settings); end
        t_ipm = minimum([@elapsed solve(prob, settings) for _ in 1:nruns])

        for _ in 1:nwarmup
            m, _ = build_jump_num(inst, optimizer); optimize!(m)
        end
        t_mosek = minimum([@elapsed begin
            m, _ = build_jump_num(inst, optimizer); optimize!(m)
        end for _ in 1:nruns])

        if dual_optimizer !== nothing
            for _ in 1:nwarmup
                m, _ = build_jump_num(inst, dual_optimizer); optimize!(m)
            end
            t_dual = minimum([@elapsed begin
                m, _ = build_jump_num(inst, dual_optimizer); optimize!(m)
            end for _ in 1:nruns])
            @printf("%-18s %6.0e %5d %4d %5d %9.1f %9.1f %9.1f %6.2fx %6.2fx\n",
                    label, raug, size(prob.B, 2), info.nv, info.h1,
                    t_ipm * 1000, t_mosek * 1000, t_dual * 1000,
                    t_mosek / t_ipm, t_dual / t_ipm)
        else
            @printf("%-18s %6.0e %5d %4d %5d %9.1f %9.1f %7.2fx\n",
                    label, raug, size(prob.B, 2), info.nv, info.h1,
                    t_ipm * 1000, t_mosek * 1000, t_mosek / t_ipm)
        end
    end
end

# =====================================================================
# Command-line interface
# =====================================================================

if abspath(PROGRAM_FILE) == @__FILE__
    include(joinpath(@__DIR__, "..", "benchmark_utils.jl"))
    using Dualization

    opts = parse_benchmark_args(ARGS)
    println("Network Utility Maximization benchmark")
    println("Solver: $(opts.mosek ? "Mosek" : "Clarabel (open-source)")")
    println()

    if opts.mosek
        using MosekTools
        optimizer = Mosek.Optimizer
        dual_optimizer = Dualization.dual_optimizer(Mosek.Optimizer)
        solver_name = "Mosek"
    else
        using Clarabel
        optimizer = Clarabel.Optimizer
        dual_optimizer = nothing  # skip dualization (ExponentialCone issues)
        solver_name = "Clarabel"
    end

    run_num_benchmark(; optimizer, dual_optimizer, solver_name,
                      nwarmup = opts.nwarmup, nruns = opts.nruns)
end
