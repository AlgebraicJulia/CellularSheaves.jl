"""
    CoordinationProfiling

Stage-level accounting for the two coordination laws benchmarked in
[`CoordinationBenchmarks`](@ref).

Where `step_profile` answers "how long did each stage take", this module answers
the questions that decide whether a law can actually be deployed on a fleet:
how much memory each stage allocates and holds, how much arithmetic it performs,
how many matrix entries it touches, and what the latency would be if every agent
had its own processor instead of the whole benchmark running on one core.

The parallel-latency model is the important one. Wall-clock here is measured as
total arithmetic on a single core for *both* laws, which is symmetric but
describes neither deployment. [`critical_path_profile`](@ref) models the
deployed system: Diffusion's per-tick latency is the busiest agent, not the sum
over agents, and Direct's solve latency is the heaviest root-to-leaf chain of the
elimination tree, not the whole factor.
"""
module CoordinationProfiling

using ..CoordinationBenchmarks
using ..CoordinationBenchmarks: CoordinationScenario, DirectPlan,
    diffusion_residual!, diffusion_control!, direct_control!,
    _integrate!, _target_state, _initial_state, _permute!, _invpermute!,
    stable_gain_ceiling, tree_statistics
using CellularSheaves.NetworkSheaves.DistributedSolve
using CliqueTrees.Multifrontal
using Graphs
using LinearAlgebra
using SparseArrays

export StageMeasurement, stage_profile, stage_totals,
    CriticalPathProfile, critical_path_profile

# ===========================================================================
# measurement
# ===========================================================================

"""
    StageMeasurement

One stage of one control law, measured along every axis that matters for
deployment.

# Fields
- `stage`, `description`, `cadence`: what the stage is and how often it runs.
- `calls`: how many times the mission invokes it.
- `t_min`, `t_median`, `t_p90`: per-call time in seconds. The spread between
  `t_min` and `t_p90` is the honest measure of how much garbage collection and
  scheduler noise affect the stage; a minimum alone hides it.
- `bytes`, `allocs`: heap allocated per call. A stage that allocates cannot run
  in a hard-real-time control loop, independent of how fast it is.
- `flops`: arithmetic operations per call, derived from the sheaf and the factor
  structure rather than timed, so it is machine independent.
- `nonzeros`: matrix entries the stage reads per call.
- `resident`: bytes one agent keeps between calls. This is the deployable
  figure. Diffusion holds its own state and its neighbours', so O(degree); a
  Direct worker holds its own frontal panel, so O(treewidth^2). Neither holds
  the fleet-wide object.
- `resident_fleet`: the same quantity summed over the fleet, reported only so
  that the per-agent figure cannot be mistaken for it.
"""
struct StageMeasurement
    stage::Symbol
    description::String
    cadence::String
    calls::Int
    t_min::Float64
    t_median::Float64
    t_p90::Float64
    bytes::Int
    allocs::Int
    flops::Int
    nonzeros::Int
    resident::Int
    resident_fleet::Int
end

"""Per-call time distribution and allocation behaviour of `op`."""
function _measure(op, reps)
    op()
    ts = Vector{Float64}(undef, reps)
    for r in 1:reps
        t0 = time_ns()
        op()
        ts[r] = (time_ns() - t0) * 1e-9
    end
    sort!(ts)
    ## measured inside a function: at global scope the boxing of non-const
    ## globals is attributed to the call and reports allocations that the stage
    ## never performs
    bytes = @allocated op()
    allocs = @allocations op()
    idx90 = clamp(ceil(Int, 0.9 * reps), 1, reps)
    return (t_min = ts[1], t_median = ts[(reps + 1) ÷ 2], t_p90 = ts[idx90],
            bytes = Int(bytes), allocs = Int(allocs))
end

# ===========================================================================
# arithmetic models
# ===========================================================================

"""Residual and separator width of each supernode of the factor."""
function _supernode_dims(L)
    S = L.S
    return [(length(collect(neighbors(S.res, j))), length(collect(neighbors(S.sep, j))))
            for j in 1:nv(S.res)]
end

"""Dense Cholesky cost of one supernode: factor the residual block, solve the
off-diagonal panel against it, then update the separator with its outer product."""
_supernode_factor_flops(r, s) = (r^3) ÷ 3 + r^2 * s + r * s^2

"""One triangular substitution over a supernode: the residual block against
itself, plus the off-diagonal panel applied to the separator."""
_supernode_solve_flops(r, s) = r^2 + 2r * s

"""Directed agent-to-agent and agent-to-target incidences of the sheaf."""
function _incidences(s::CoordinationScenario)
    agent = sum(length(s.agent_nbrs[i]) for i in 1:s.nagents; init = 0)
    target = sum(length(s.target_nbrs[i]) for i in 1:s.nagents; init = 0)
    return (agent, target)
end

# ===========================================================================
# stage inventories, measured
# ===========================================================================

"""
    stage_profile(method, scenario; ticks, solves, reps) -> Vector{StageMeasurement}

Measure every stage of one control law, from the first setup step to the last
plant integration.

`ticks` should be the tick count the law actually needs to reach the formation,
and `solves` the number of harmonic solves the mission performs, so that the
per-stage costs sum to the end-to-end cost rather than to an arbitrary horizon.
With targets held fixed Direct solves once, so `solves = 1`.
"""
function stage_profile(method::Symbol, s::CoordinationScenario;
        ticks::Integer, solves::Integer = 1, reps::Integer = 200,
        dt = 0.005, safety = 0.4)
    method in (:diffusion, :direct) ||
        throw(ArgumentError("unknown method: $method"))
    n = size(s.H, 1)
    d = s.dim
    k = safety * stable_gain_ceiling(method, s, dt)
    p = _target_state(s)
    q = _initial_state(s)
    qtmp = copy(q)
    u, eta, qstar = zeros(n), zeros(n), zeros(n)
    n_agent, n_target = _incidences(s)
    out = StageMeasurement[]

    push_stage!(sym, desc, cadence, calls, op, flops, nz, resident, fleet) = begin
        m = _measure(op, reps)
        push!(out, StageMeasurement(sym, desc, cadence, calls, m.t_min, m.t_median,
            m.t_p90, m.bytes, m.allocs, flops, nz, resident, fleet))
    end

    if method === :diffusion
        ## an agent's whole working set: its own state, one buffer per
        ## neighbour, and the targets it senses
        per_agent = maximum(8d * (1 + length(s.agent_nbrs[i]) + length(s.target_nbrs[i]))
                            for i in 1:s.nagents; init = 0)
        ## the radio exchange is charged in slots by the communication model, not
        ## in seconds; it is listed so the inventory is complete
        push!(out, StageMeasurement(:exchange, "exchange neighbour states", "per tick",
            ticks, 0.0, 0.0, 0.0, 0, 0, 0, 0,
            maximum(8d * length(s.agent_nbrs[i]) for i in 1:s.nagents; init = 0),
            8d * n_agent))
        ## one subtraction and one accumulation per incidence per component
        push_stage!(:residual, "accumulate disagreement", "per tick", ticks,
            () -> diffusion_residual!(eta, s, q, p),
            2d * (n_agent + n_target), n_agent + n_target, per_agent, 8n)
        push_stage!(:control, "scale by the gain", "per tick", ticks,
            () -> diffusion_control!(u, eta, k), n, 0, 8d, 8n)
        push_stage!(:integrate, "integrate the plant", "per tick", ticks,
            () -> _integrate!(qtmp, u, dt), 2n, 0, 8d, 8n)
        return out
    end

    Hc = copy(s.H)
    plan0 = direct_plan(s)
    dims0 = _supernode_dims(plan0.factor.L)
    ## a worker holds its own frontal panel: the dense residual block plus the
    ## off-diagonal panel against its separator
    panel(r, sep) = 8 * (r * r + r * sep)
    per_worker = maximum(panel(r, sep) for (r, sep) in dims0; init = 0)
    fleet_factor = sum(panel(r, sep) for (r, sep) in dims0; init = 0)

    push_stage!(:assemble, "assemble the operator blocks", "once", 1,
        () -> copy(s.H), 0, nnz(s.H), per_worker,
        Int(Base.summarysize(s.H) + Base.summarysize(s.Bmat)))

    sym_only = _measure(() -> ChordalCholesky(copy(Hc)), max(reps ÷ 20, 3))
    push!(out, StageMeasurement(:symbolic, "order and build the clique tree", "once", 1,
        sym_only.t_min, sym_only.t_median, sym_only.t_p90, sym_only.bytes,
        sym_only.allocs, 0, nnz(s.H), 0, 0))

    plan = plan0
    dims = dims0
    factor_flops = sum(_supernode_factor_flops(r, sep) for (r, sep) in dims; init = 0)
    solve_flops = sum(_supernode_solve_flops(r, sep) for (r, sep) in dims; init = 0)
    stats = tree_statistics(plan)

    ## the numeric factorization is measured together with the symbolic analysis
    ## it depends on, then the symbolic part subtracted; two independent minima
    ## can cross, so clamp
    both = _measure(() -> cholesky!(ChordalCholesky(copy(Hc)), NoPivot()), max(reps ÷ 20, 3))
    ## one baseline is subtracted from all three quantiles. Subtracting the
    ## symbolic p90 from the combined p90 would difference two independent
    ## distributions and can invert the ordering, reporting a p90 below the min.
    push!(out, StageMeasurement(:factor, "numeric factorization", "once", 1,
        max(both.t_min - sym_only.t_min, 0.0), max(both.t_median - sym_only.t_min, 0.0),
        max(both.t_p90 - sym_only.t_min, 0.0), both.bytes, both.allocs,
        factor_flops, nnz(s.H), per_worker, fleet_factor))

    push_stage!(:workspace, "allocate the solve workspace", "once", 1,
        () -> TreeWorkspace(plan.factor.L, Float64), 0, 0,
        maximum(8 * sep * sep for (_, sep) in dims0; init = 0),
        Int(Base.summarysize(plan.workspace)))

    rhs, scratch, perm = plan.rhs, plan.scratch, plan.perm
    push_stage!(:rhs, "form the right-hand side", "per solve", solves,
        () -> mul!(rhs, s.Bmat, vec(p)), 2nnz(s.Bmat), nnz(s.Bmat), 8d, 8n)
    push_stage!(:permute, "apply the fill-reducing permutation", "per solve", solves,
        () -> _permute!(scratch, rhs, perm), 0, 0, 8d, 8n)
    push_stage!(:forward, "forward substitution up the tree", "per solve", solves,
        () -> tree_forward_ldiv!(plan.factor.L, scratch, plan.workspace),
        solve_flops, stats.fill, per_worker, fleet_factor)
    push_stage!(:backward, "backward substitution down the tree", "per solve", solves,
        () -> tree_backward_ldiv!(plan.factor.L, scratch, plan.workspace),
        solve_flops, stats.fill, per_worker, fleet_factor)
    push_stage!(:unpermute, "undo the permutation", "per solve", solves,
        () -> _invpermute!(qstar, scratch, perm), 0, 0, 8d, 8n)
    push_stage!(:control, "scale the tracking error", "per tick", ticks,
        () -> direct_control!(u, q, qstar, k), 2n, 0, 8d, 8n)
    push_stage!(:integrate, "integrate the plant", "per tick", ticks,
        () -> _integrate!(qtmp, u, dt), 2n, 0, 8d, 8n)
    return out
end

"""
    stage_totals(stages) -> NamedTuple

Sum a stage inventory over the whole mission, weighting each stage by its call
count. `resident` is a maximum rather than a sum, since the stages hold their
memory concurrently.
"""
function stage_totals(stages::Vector{StageMeasurement})
    seconds = sum(st.calls * st.t_min for st in stages; init = 0.0)
    bytes = sum(st.calls * st.bytes for st in stages; init = 0)
    allocs = sum(st.calls * st.allocs for st in stages; init = 0)
    flops = sum(st.calls * st.flops for st in stages; init = 0)
    resident = maximum((st.resident for st in stages); init = 0)
    resident_fleet = maximum((st.resident_fleet for st in stages); init = 0)
    return (; seconds, bytes, allocs, flops, resident, resident_fleet)
end

# ===========================================================================
# parallel latency
# ===========================================================================

"""
    CriticalPathProfile

What one law would cost on a fleet where every agent carries its own processor,
as opposed to the single-core totals reported elsewhere.

# Fields
- `serial_flops`: total arithmetic, the quantity a one-core benchmark measures.
- `parallel_flops`: arithmetic along the critical path, which is what a
  distributed execution actually waits for.
- `speedup`: the ratio, an upper bound on what parallel hardware could buy.
- `depth`: hops along the critical path. For Diffusion this is the tick count,
  since each tick is one synchronous round; for Direct it is the height of the
  elimination tree traversed twice.
- `busiest`: the heaviest single unit of work on the path, the agent or
  supernode that sets the pace.
"""
struct CriticalPathProfile
    method::Symbol
    serial_flops::Int
    parallel_flops::Int
    speedup::Float64
    depth::Int
    busiest::Int
end

"""Heaviest root-to-leaf chain of the elimination tree under a per-supernode weight."""
function _critical_chain(L, weight)
    S = L.S
    ns = nv(S.res)
    memo = Dict{Int, Int}()
    function down(j)
        haskey(memo, j) && return memo[j]
        ch = collect(neighbors(S.chd, j))
        best = isempty(ch) ? 0 : maximum(down(c) for c in ch)
        return (memo[j] = best + weight(j))
    end
    roots = [j for j in 1:ns if iszero(S.pnt[j])]
    return maximum((down(r) for r in roots); init = 0)
end

"""
    critical_path_profile(method, scenario; ticks, solves) -> CriticalPathProfile

Model the deployed cost of one law rather than its single-core cost.

Diffusion parallelizes perfectly across agents: every agent forms its own
disagreement from data it already holds, so a tick costs the busiest agent, not
the sum over the fleet. Direct parallelizes across the elimination tree, but only
down to its critical path: independent subtrees factor and substitute
concurrently, while a chain of supernodes must run in sequence.

Both numbers are derived from the sheaf and the factor structure, so neither
depends on this machine.
"""
function critical_path_profile(method::Symbol, s::CoordinationScenario;
        ticks::Integer, solves::Integer = 1)
    method in (:diffusion, :direct) ||
        throw(ArgumentError("unknown method: $method"))
    d = s.dim
    n = size(s.H, 1)

    if method === :diffusion
        n_agent, n_target = _incidences(s)
        per_tick_serial = 2d * (n_agent + n_target) + 3n
        busiest = maximum(2d * (length(s.agent_nbrs[i]) + length(s.target_nbrs[i])) + 3d
                          for i in 1:s.nagents; init = 0)
        return CriticalPathProfile(:diffusion, ticks * per_tick_serial,
            ticks * busiest, per_tick_serial / max(busiest, 1), Int(ticks), busiest)
    end

    plan = direct_plan(s)
    dims = _supernode_dims(plan.factor.L)
    serial_factor = sum(_supernode_factor_flops(r, sep) for (r, sep) in dims; init = 0)
    serial_solve = sum(_supernode_solve_flops(r, sep) for (r, sep) in dims; init = 0)
    serial = serial_factor + 2solves * serial_solve + solves * 2nnz(s.Bmat) + ticks * 4n

    par_factor = _critical_chain(plan.factor.L,
        j -> _supernode_factor_flops(dims[j][1], dims[j][2]))
    par_solve = _critical_chain(plan.factor.L,
        j -> _supernode_solve_flops(dims[j][1], dims[j][2]))
    ## the per-agent stages stay parallel; only the tree passes serialize
    parallel = par_factor + 2solves * par_solve + solves * 2d + ticks * 4d
    stats = tree_statistics(plan)
    busiest = maximum(_supernode_factor_flops(r, sep) for (r, sep) in dims; init = 0)
    return CriticalPathProfile(:direct, serial, parallel,
        serial / max(parallel, 1), 2stats.depth, busiest)
end

end # module
