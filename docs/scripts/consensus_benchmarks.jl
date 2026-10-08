# Benchmarks for optimal consensus of planar double integrators
# (ControlSheaves.PredictiveConsensus).
#
# Open loop: one horizon-T plan for the fleet, from
#   - the joint Riccati recursion (exact, unbounded only),
#   - one centralized conic QP (Mumblebee IPM),
#   - predicted-trajectory sweeps, multicolor (Gauss–Seidel) and damped
#     parallel (Jacobi), each agent solving its own trajectory QP; with the
#     bound also cold starts (μ′ = 0 and μ′ = 1e-6) and an exact polish,
# reporting time, sweeps, neighbour-exchange rounds and the cost gap to the
# optimum.
#
# Closed loop: the fleet driven for STEPS steps by
#   - the centralized infinite-horizon LQR (clipped to the bound if bounded),
#   - the second-order diffusion law u = -kp η - kv v, gains tuned by grid
#     search for the lowest cost (the best case for diffusion),
#   - distributed receding-horizon control with 1, 3 and 10 sweeps per step
#     (interior warm starts; bounded cases also run 3 sweeps cold),
# scored by the team cost of the run against the open-loop optimum over the
# same window (a lower bound for every causal controller).
#
# Run with:  julia -t 8 --project=docs docs/scripts/consensus_benchmarks.jl
# Writes docs/figures/consensus/benchmarks.csv.

using CellularSheaves
using CellularSheaves.ControlSheaves.CoordinationBenchmarks: coordination_scenario
using CellularSheaves.ControlSheaves.PredictiveConsensus
using CellularSheaves.ControlSheaves.PredictiveConsensus: init, solve, solve!
using CellularSheaves.NetworkSheaves.SchwarzMethods: MulticolorSweep, ParallelSweep
using LinearAlgebra
using Printf
using Random

const OUT = joinpath(@__DIR__, "..", "figures", "consensus")
mkpath(OUT)
const DT = 0.1
const HORIZON = 40
const STEPS = 120
const BOUND = 2.0

# A fleet: a named agent graph with its pinned targets.
struct FleetCase
    name::String
    family::Symbol
    size_parameter::Int
end

# One measured method on one case.
struct ConsensusRow
    case::String
    agents::Int
    bounded::Bool
    mode::String          # "open loop" or "closed loop"
    method::String
    seconds::Float64
    sweeps::Int
    rounds::Int
    ipm_iterations::Int   # interior-point iterations, all local solves
    cost::Float64
    gap::Float64          # cost / optimum - 1
end

# Wall time and result of f(); runs under two seconds are repeated three times.
function timed(f)
    GC.gc()
    t = @elapsed result = f()
    t < 2 || return t, result
    for _ in 1:2
        GC.gc()
        t = min(t, @elapsed f())
    end
    return t, result
end

function fleet(c::FleetCase, horizon; bounded)
    scenario = coordination_scenario(c.family; size_parameter = c.size_parameter)
    span = sqrt(scenario.nagents)
    targets = reshape([span * (k - 1) / max(1, scenario.ntargets - 1) for k in 1:scenario.ntargets for _ in 1:2],
        2, scenario.ntargets)
    lq = ConsensusLQ(scenario, PlanarDoubleIntegrator(DT), horizon; targets,
        control_bound = bounded ? BOUND : Inf)
    rng = Xoshiro(7)
    x0 = [span .* randn(rng, 2, scenario.nagents); randn(rng, 2, scenario.nagents)]
    return lq, x0
end

function measure(c::FleetCase; bounded)
    rows = ConsensusRow[]
    lq, x0 = fleet(c, HORIZON; bounded)
    N = lq.scenario.nagents
    prob = ConsensusProblem(lq, x0)
    add!(mode, method, t, sweeps, rounds, ipm, cost, opt) =
        push!(rows, ConsensusRow(c.name, N, bounded, mode, method, t, sweeps, rounds, ipm, cost, cost / opt - 1))

    t_qp, central = timed(() -> solve(prob, CentralizedQP()))
    optimum = central.cost
    if !bounded
        t, plan = timed(() -> solve(prob, JointRiccati()))
        add!("open loop", "joint Riccati", t, 0, 0, 0, plan.cost, optimum)
    end
    add!("open loop", "centralized QP (IPM)", t_qp, 0, 0, central.local_iterations, central.cost, optimum)
    colors = length(init(prob, PredictedTrajectorySweeps()).colors)
    variants = [
        ("sweeps, multicolor", PredictedTrajectorySweeps(tol = 1e-7, maxiter = 2000), colors),
        ("sweeps, damped parallel", PredictedTrajectorySweeps(sweep = ParallelSweep(), tol = 1e-7, maxiter = 5000), 1)]
    if bounded
        append!(variants, [
            ("sweeps, multicolor, cold, μ′ = 0",
                PredictedTrajectorySweeps(tol = 1e-7, maxiter = 2000, barrier = 0.0, warm_start = false), colors),
            ("sweeps, multicolor, cold, μ′ = 1e-6",
                PredictedTrajectorySweeps(tol = 1e-7, maxiter = 2000, warm_start = false), colors),
            ("sweeps, multicolor + 20 exact polish",
                PredictedTrajectorySweeps(tol = 1e-7, maxiter = 2000, polish = 20), colors)])
    end
    for (name, alg, per_sweep) in variants
        t, plan = timed(() -> solve(prob, alg))
        add!("open loop", name, t, plan.converged ? plan.iterations : -1, plan.iterations * per_sweep,
            plan.local_iterations, plan.cost, optimum)
    end

    # Closed loop, scored against the open-loop optimum over the whole window.
    window, _ = fleet(c, STEPS; bounded)
    best = solve(ConsensusProblem(window, x0), bounded ? CentralizedQP() : JointRiccati()).cost
    t, run = timed(() -> rollout(lq, RiccatiFeedback(lq), x0, STEPS))
    add!("closed loop", bounded ? "centralized LQR (clipped)" : "centralized LQR", t, 0, 0, 0, run.cost, best)
    tuned = argmin(((kp, kv) for kp in (0.25, 0.5, 1.0, 2.0, 4.0, 8.0), kv in (0.5, 1.0, 2.0, 4.0, 8.0))) do (kp, kv)
        cost = rollout(lq, SecondOrderDiffusion(lq, kp, kv), x0, STEPS).cost
        isfinite(cost) ? cost : Inf
    end
    t, run = timed(() -> rollout(lq, SecondOrderDiffusion(lq, tuned...), x0, STEPS))
    add!("closed loop", @sprintf("diffusion (kp = %g, kv = %g)", tuned...), t, 0, STEPS, 0, run.cost, best)
    horizon_variants = [(s, true) for s in (1, 3, 10)]
    bounded && push!(horizon_variants, (3, false))
    for (sweeps, warm) in horizon_variants
        alg = PredictedTrajectorySweeps(tol = 1e-7, warm_start = warm)
        t, (run, ctrl) = timed() do
            ctrl = RecedingHorizon(prob, alg; sweeps)
            rollout(lq, ctrl, x0, STEPS), ctrl
        end
        add!("closed loop", "receding horizon, $sweeps sweep" * (sweeps == 1 ? "" : "s") * "/step" * (warm ? "" : ", cold"),
            t, sweeps, STEPS * sweeps * colors, ctrl.local_iterations, run.cost, best)
    end
    return rows
end

cases = [
    FleetCase("grid 3×3", :grid, 3),
    FleetCase("grid 5×5", :grid, 5),
    FleetCase("ring 32", :ring, 32),
    FleetCase("grid 8×8", :grid, 8),
]

measure(FleetCase("warm-up", :grid, 2); bounded = false)
measure(FleetCase("warm-up", :grid, 2); bounded = true)
rows = reduce(vcat, [measure(c; bounded) for c in cases for bounded in (false, true)])

open(joinpath(OUT, "benchmarks.csv"), "w") do io
    println(io, "case,agents,bounded,mode,method,seconds,sweeps,rounds,ipm_iterations,cost,gap")
    for r in rows
        println(io, join((repr(r.case), r.agents, r.bounded, repr(r.mode), repr(r.method), r.seconds,
            r.sweeps, r.rounds, r.ipm_iterations, r.cost, r.gap), ","))
    end
end

for c in cases, bounded in (false, true)
    selected = filter(r -> r.case == c.name && r.bounded == bounded, rows)
    @printf("\n%s: %d agents, horizon %d, %s, threads = %d\n", c.name, first(selected).agents, HORIZON,
        bounded ? "‖u‖ ≤ $BOUND" : "unbounded", Threads.nthreads())
    @printf("%-12s %-40s %9s %7s %7s %9s %12s %10s\n", "", "method", "seconds", "sweeps", "rounds", "IPM its", "cost", "gap")
    for r in selected
        @printf("%-12s %-40s %9.4f %7s %7d %9d %12.4f %10.2e\n", r.mode, r.method, r.seconds,
            r.sweeps < 0 ? "no conv" : string(r.sweeps), r.rounds, r.ipm_iterations, r.cost, r.gap)
    end
end
println("\nwrote ", joinpath(OUT, "benchmarks.csv"))
