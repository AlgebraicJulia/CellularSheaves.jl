# Shared closed-loop simulation driver for the `nested/` literate examples.
#
# Not a standalone documentation page: `docs/make.jl`'s literate walk skips files whose name
# starts with `_`, so this produces no generated markdown/notebook. Both examples in this
# directory `include` this file (located via `pkgdir(CellularSheaves)`, not `@__DIR__` -- see
# their own comments) to pull in
# `NestedEscortProblem`/`NestedEscortResult`/`run_nested_escort_simulation`.
#
# Mirrors `Layered.run_layered_escort_simulation` (`src/ControlSheaves/Layered.jl`) structurally,
# built on `NestedSystems`'s tree-shaped primitives (`SheafTower`, `solve_hierarchical`,
# `resolve_dynamics`) instead of `LayeredEscortSpec`'s flat ones. Kept here rather than promoted
# into `src/` for now -- no precompile cost while iterating -- but shaped as a Problem/Result
# pair so a future `ext/CellularSheavesPlots` recipe could dispatch on `NestedEscortResult` the
# same way `LayeredEscortPlots.jl` dispatches on `LayeredSimulationResult`.

using CellularSheaves.ControlSheaves.NestedSystems
using CellularSheaves.ControlSheaves.AgentControllers

"""
    NestedEscortProblem(tower, bindings, target_trajectories; target_velocities=nothing,
                        target_accelerations=nothing, dt=0.05, steps=200)

A closed-loop tracking problem over a [`SheafTower`](@ref): `tower` supplies the topology and
per-timestep [`solve_hierarchical`](@ref) reference solve, `bindings` resolves per-agent
dynamics/gain/initial state via [`resolve_dynamics`](@ref), and `target_trajectories[t]` is a
function `time -> Vector{Float64}` giving target `t`'s `D`-dimensional world position.
`target_velocities`/`target_accelerations`, if supplied, enable feedforward tracking (see
[`run_nested_escort_simulation`](@ref)).
"""
struct NestedEscortProblem
    tower::SheafTower
    bindings::SystemBinding
    target_trajectories::Vector{Any}
    target_velocities::Union{Nothing,Vector{Any}}
    target_accelerations::Union{Nothing,Vector{Any}}
    dt::Float64
    steps::Int
end

NestedEscortProblem(tower::SheafTower, bindings::SystemBinding, target_trajectories::Vector;
                     target_velocities=nothing, target_accelerations=nothing,
                     dt::Float64=0.05, steps::Int=200) =
    NestedEscortProblem(tower, bindings, Vector{Any}(target_trajectories),
                        target_velocities === nothing ? nothing : Vector{Any}(target_velocities),
                        target_accelerations === nothing ? nothing : Vector{Any}(target_accelerations),
                        dt, steps)

"""
    NestedEscortResult

Output of [`run_nested_escort_simulation`](@ref). `sim_data[step][agent]` is the agent's raw
dynamics state (length `state_dim(dyn)`); `qstar_history[step][agent]` is the `D`-dimensional
hierarchical reference the agent was tracking that step (`tower.agent_vertices` order);
`target_history[step][target]` is the target's `D`-dimensional world position that step.
"""
struct NestedEscortResult
    problem::NestedEscortProblem
    sim_data::Vector{Vector{Vector{Float64}}}
    qstar_history::Vector{Vector{Vector{Float64}}}
    target_history::Vector{Vector{Vector{Float64}}}
end

"""
    run_nested_escort_simulation(prob::NestedEscortProblem; use_feedforward::Bool=false) -> NestedEscortResult

Run a closed-loop simulation over `prob.tower`'s spec: at each step, solve the tower
hierarchically for the current target positions, lift the result to per-agent references, and
step each agent's [`AgentState`](@ref) toward that reference.

`use_feedforward` toggles two things, exactly as in `Layered.run_layered_escort_simulation`:
whether `AgentState` uses a joint (position+velocity) Tikhonov filter, and whether velocity/
acceleration references are additionally solved for (via the same hierarchical machinery, on
`prob.target_velocities`/`prob.target_accelerations`) and passed to [`step_agent!`](@ref) for
differential-flatness feedforward.
"""
function run_nested_escort_simulation(prob::NestedEscortProblem; use_feedforward::Bool=false)
    tower = prob.tower
    resolved = resolve_dynamics(tower.spec, prob.bindings)   # tower.agent_vertices order (Issue 012)
    agent_states = [AgentState(initial_state(r.dynamics, r.initial_position), r.dynamics,
                                prob.dt, r.K_lqr, 0.02; use_velocity=use_feedforward)
                    for r in resolved]

    time_grid = 0:prob.dt:(prob.steps * prob.dt)
    sim_data = Vector{Vector{Vector{Float64}}}()
    qstar_history = Vector{Vector{Vector{Float64}}}()
    target_history = Vector{Vector{Vector{Float64}}}()

    for step in 1:prob.steps
        t = time_grid[step]
        p_targets = [traj(t) for traj in prob.target_trajectories]
        q = solve_hierarchical(tower, p_targets)[end]
        qstar_agents = [q[v] for v in tower.agent_vertices]
        push!(qstar_history, qstar_agents)
        push!(target_history, [q[v] for v in tower.target_vertices])

        qdot_agents = qddot_agents = nothing
        if use_feedforward && prob.target_velocities !== nothing
            v = solve_hierarchical(tower, [vt(t) for vt in prob.target_velocities])[end]
            qdot_agents = [v[vtx] for vtx in tower.agent_vertices]
        end
        if use_feedforward && prob.target_accelerations !== nothing
            a = solve_hierarchical(tower, [at(t) for at in prob.target_accelerations])[end]
            qddot_agents = [a[vtx] for vtx in tower.agent_vertices]
        end

        step_states = Vector{Vector{Float64}}()
        for i in eachindex(agent_states)
            pd = length(position_indices(agent_states[i].dyn))
            q_ref = qstar_agents[i][1:pd]
            if use_feedforward && qdot_agents !== nothing && qddot_agents !== nothing
                step_agent!(agent_states[i], q_ref, qdot_agents[i][1:pd], qddot_agents[i][1:pd], prob.dt)
            elseif use_feedforward && qdot_agents !== nothing
                step_agent!(agent_states[i], q_ref, qdot_agents[i][1:pd], prob.dt)
            else
                step_agent!(agent_states[i], q_ref, prob.dt)
            end
            push!(step_states, copy(agent_states[i].x))
        end
        push!(sim_data, step_states)
    end

    return NestedEscortResult(prob, sim_data, qstar_history, target_history)
end

"""
    agent_index_ranges(team_sizes::Vector{Int}) -> Vector{UnitRange{Int}}

Contiguous index ranges into `tower.agent_vertices` for a sequence of leaf teams of the given
sizes, in the order those teams were added as children. Agent indices are assigned depth-first
over the spec's tree (see `NestedSystems._assign_agent_ranges!`), so each leaf team occupies a
contiguous block in visitation order -- this reconstructs those blocks from the team sizes alone,
without reaching into `NestedSystems`'s private bookkeeping.
"""
function agent_index_ranges(team_sizes::Vector{Int})
    ranges = Vector{UnitRange{Int}}()
    off = 0
    for m in team_sizes
        push!(ranges, (off + 1):(off + m))
        off += m
    end
    return ranges
end

"""
    translation_pin(n_members::Int, D::Int, k::Int) -> RawRestriction

A restriction map pinning direct member `k`'s **translation** components (the first `D-1`
coordinates) to a target, deliberately leaving the homogeneous row (`D`) unconstrained.

Used for every *redundant* pin beyond the first when several agents around a ring each observe
the same target (see `redundant_pin` below). `project(k)` materializes as a full `D×D` identity
block, homogeneous row included; a single such pin is fine (it unambiguously forces that row to
`1`), but *several* pins to the same target are mutually inconsistent for a rigid body by
construction -- that inconsistency is the whole point, it's what rebalances a ring's vote count
against competing edges elsewhere in the tower. Left unmodified, the least-squares compromise that
inconsistency produces doesn't stay confined to the intended translation components: it drags the
homogeneous row along with it, away from `1`, which -- because these affine restriction maps only
represent *pure* translation correctly when the homogeneous row is exactly `1` -- rescales the
entire rigid body's geometry. `translation_pin` avoids this by only ever contesting translation:
only the first, unmodified `project(1)` pin is left responsible for anchoring the homogeneous row.
"""
function translation_pin(n_members::Int, D::Int, k::Int)
    M = zeros(D, n_members * D)
    for i in 1:(D - 1)
        M[i, (k - 1) * D + i] = 1.0
    end
    return RawRestriction(M)
end

"""
    redundant_pin(n_members::Int, D::Int, k::Int) -> RestrictionSpec

`project(1)` for the first (`k == 1`) pin, [`translation_pin`](@ref) for every pin after that --
the system-map to use for member `k` of a `1:2:n_members`-style redundant-pin scheme.
"""
redundant_pin(n_members::Int, D::Int, k::Int) = k == 1 ? project(1) : translation_pin(n_members, D, k)

"""
    snapshot_steps(steps::Int, n_snapshots::Int) -> Vector{Int}

`n_snapshots` step indices evenly spaced across `1:steps` (endpoints included). The temporal
resolution a simulation needs for accurate dynamics (every `dt`) is far finer than a static plot
needs for a legible composite of a formation's shape over time -- this picks the sparse subset
worth actually drawing.
"""
snapshot_steps(steps::Int, n_snapshots::Int) = unique(round.(Int, range(1, steps, length=n_snapshots)))

"""
    fade_alpha(i::Int, n::Int; lo=0.15, hi=0.85) -> Float64

Linear alpha ramp from `lo` (first of `n` items, `i == 1`) to `hi` (last, `i == n`) -- used so a
sequence of overlaid snapshots reads as "earlier is fainter, later is more solid" rather than a
uniform, hard-to-parse overlay.
"""
fade_alpha(i::Int, n::Int; lo::Real=0.15, hi::Real=0.85) = n <= 1 ? hi : lo + (hi - lo) * (i - 1) / (n - 1)
