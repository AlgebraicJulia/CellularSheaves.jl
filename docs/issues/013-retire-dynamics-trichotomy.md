# Issue 013 — Retire the Dynamics Trichotomy in Favour of the Inheritance Cascade

**Parent design**: `007-nested-layered-systems-design.md` (§4.7)
**Prerequisites**: 012 (`SystemBinding` / `resolve_dynamics` must exist first)
**Blocks**: nothing

---

## Background

`Layered.jl` currently expresses dynamics assignment as a closed trichotomy of three
`AbstractDynamicsSpec` subtypes:

- `HomogeneousDynamics` — one model for every agent
- `TeamHomogeneousDynamics` — one model per team
- `IndividualizedDynamics` — one model per agent

Issue 012 replaces this with an arbitrary-depth inheritance cascade (`SystemBinding` /
`AgentBinding` / `resolve_dynamics`, most specific wins). The three fixed tiers become degenerate
cases of declaring at the root, at each leaf team, and at each agent respectively — so the
trichotomy is now redundant and should be removed rather than maintained alongside.

**This is a deletion plus a port, not a pure deletion.** `LayeredEscortProblem` accepts a flat
`LayeredEscortSpec`, whereas `resolve_dynamics` walks a `NestedSystemSpec` tree. A flat spec is
structurally a depth-1 tree — rings and support pods as children of a single root — so the port is
mechanical, but it must be written.

---

## Codebase Orientation

| File | Why |
|---|---|
| `src/ControlSheaves/Layered.jl:33–68` | The types being removed, plus `get_agent_dynamics_config`'s three dispatch methods. |
| `src/ControlSheaves/Layered.jl:400–447` | `LayeredEscortProblem`'s `dynamics_spec::AbstractDynamicsSpec` field (line 406) and constructor (line 429). |
| `src/ControlSheaves/Layered.jl:484–506` | `run_layered_escort_simulation`'s per-agent setup. Line 502 calls `get_agent_dynamics_config(prob.dynamics_spec, i, team_idx)`; lines 486–501 compute `team_idx` by scanning `ring_node_ranges` / `support_node_ranges`. That scan is exactly the flat tree walk the adapter replaces. |
| `src/ControlSheaves/Layered.jl:8–14` | Export list. |
| `src/ControlSheaves/NestedSystems.jl` | Issue 012's `SystemBinding`, `AgentBinding`, `ResolvedAgent`, `resolve_dynamics`. |
| `test/network_sheaves/LayeredEscort.jl:96, 126, 172` | Three `HomogeneousDynamics(dyn, K)` call sites to migrate. |
| `docs/literate/layered/multilayer_escort.jl:77` | Literate example using `HomogeneousDynamics`. **Missing this breaks the docs build**, which is a separate command from the test suite — see `CLAUDE.md`. |
| `Agents.md:193` | Cites `get_agent_dynamics_config(spec::HomogeneousDynamics, ...)` as the canonical example of the repo's "multiple dispatch over `isa` branching" convention. Removing the function orphans this reference; it must be repointed at a surviving example. |

---

## Requested Implementation

### 1. Flat-spec adapter

```julia
"""
    resolve_dynamics(spec::LayeredEscortSpec, ctx::SystemBinding) -> Vector{ResolvedAgent}

Resolve dynamics bindings for a flat layered escort specification, treating it as a depth-1
tree: each escort ring and each support pod is a child of a single implicit root.

Children are addressed by the conventional names `:ring1, :ring2, …` and `:support1, :support2, …`,
matching the ordering of `spec.rings` and `spec.supports`. Per-agent overrides use the agent's
**local** index within its ring or pod, consistent with `RingSpec.observers`.

Agents come back in global agent-index order, matching `spec.ring_node_ranges` followed by
`spec.support_node_ranges` — the same ordering `run_layered_escort_simulation` iterates.
"""
function resolve_dynamics(spec::LayeredEscortSpec, ctx::SystemBinding)
```

Implement by reusing Issue 012's cascade logic rather than duplicating the precedence rules —
extract the fold into a shared helper if needed. The precedence must be identical.

### 2. Port `LayeredEscortProblem`

- Change the field to `bindings::SystemBinding`.
- Change the constructor's positional `dynamics_spec::AbstractDynamicsSpec` argument to
  `bindings::SystemBinding`.
- In `run_layered_escort_simulation`, call `resolve_dynamics(spec, prob.bindings)` **once**
  before the time loop and index the resulting vector, replacing both the per-agent
  `get_agent_dynamics_config` call at line 502 and the `team_idx` scan at lines 486–501. This
  removes an O(n_agents × n_teams) scan from setup as a side benefit.
- The `initial_positions` field overlaps with `AgentBinding.initial_position`. **Keep it, at the
  lowest priority.** Initial positions resolve in three tiers, most specific first:

  ```
  AgentBinding.initial_position   (per-agent binding)
  → LayeredEscortProblem.initial_positions[i]   (problem-level vector)
  → _default_initial_position(dyn, i)           (airstrip fallback, Layered.jl:464)
  ```

  Document this ordering explicitly in the constructor docstring. Existing callers passing
  `initial_positions` and no per-agent bindings must see unchanged behaviour.

### 3. Delete

- `AbstractDynamicsSpec`, `HomogeneousDynamics`, `TeamHomogeneousDynamics`,
  `IndividualizedDynamics`, and all three `get_agent_dynamics_config` methods.
- Their entries in `Layered.jl`'s export list.

### 4. Migrate call sites

- `test/network_sheaves/LayeredEscort.jl` — `HomogeneousDynamics(dyn, K)` becomes
  `SystemBinding(dynamics=dyn, K_lqr=K)`.
- `docs/literate/layered/multilayer_escort.jl` — same substitution. Rebuild the docs to confirm.
- `Agents.md:193` — repoint the dispatch-convention example at a surviving one. Good candidates:
  `initial_state` in `AgentControllers.jl` (line 119, dispatching on dynamics type) or
  `materialize_restriction` from Issue 011. Keep the convention's wording; only change the example.

---

## Tests to Write

```julia
@testset "flat adapter: root binding applies to every agent" begin
    spec = LayeredEscortSpec([RingSpec(1, 6, 0.3), RingSpec(2, 4, 0.25)], [SupportSpec(1, 2, 3)])
    r = resolve_dynamics(spec, SystemBinding(dynamics=QuadrotorDynamics()))
    @test length(r) == spec.n_agents
    @test all(x -> x.dynamics isa QuadrotorDynamics, r)
end

@testset "flat adapter: per-ring override beats the root default" begin
    spec = LayeredEscortSpec([RingSpec(1, 6, 0.3), RingSpec(2, 4, 0.25)], [SupportSpec(1, 2, 3)])
    ctx = SystemBinding(dynamics=QuadrotorDynamics(),
                        children=Dict(:ring2 => SystemBinding(dynamics=PlanarQuadrotorDynamics())))
    r = resolve_dynamics(spec, ctx)
    @test all(x -> x.dynamics isa QuadrotorDynamics, r[spec.ring_node_ranges[1]])
    @test all(x -> x.dynamics isa PlanarQuadrotorDynamics, r[spec.ring_node_ranges[2]])
end

@testset "flat adapter: support pods are addressable" begin
    spec = LayeredEscortSpec([RingSpec(1, 6, 0.3), RingSpec(2, 4, 0.25)], [SupportSpec(1, 2, 3)])
    ctx = SystemBinding(dynamics=QuadrotorDynamics(),
                        children=Dict(:support1 => SystemBinding(dynamics=PlanarQuadrotorDynamics())))
    r = resolve_dynamics(spec, ctx)
    @test all(x -> x.dynamics isa PlanarQuadrotorDynamics, r[spec.support_node_ranges[1]])
end

@testset "simulation is unchanged by the port" begin
    # Same scenario as the pre-port test at LayeredEscort.jl:96, now via SystemBinding.
    prob = LayeredEscortProblem(spec, SystemBinding(dynamics=dyn_test, K_lqr=K_test), target_trajs;
                                dt=0.05, steps=10)
    res = run_layered_escort_simulation(prob)
    @test length(res.sim_data) == 10
    @test res.sim_data[1][1][1] ≈ 0.0 atol=0.05     # airstrip default preserved
end

@testset "AgentBinding.initial_position beats the initial_positions fallback" begin
    ctx = SystemBinding(dynamics=dyn_test,
                        children=Dict(:ring1 => SystemBinding(
                            agents=Dict(1 => AgentBinding(initial_position=[7.0, 7.0, 7.0])))))
    prob = LayeredEscortProblem(spec, ctx, target_trajs;
                                initial_positions=fill([0.0, 0.0, 1.5], spec.n_agents))
    res = run_layered_escort_simulation(prob)
    @test res.sim_data[1][1][1:3] ≈ [7.0, 7.0, 7.0] atol=0.05
end

@testset "removed names are gone" begin
    @test !isdefined(CellularSheaves.ControlSheaves.Layered, :HomogeneousDynamics)
    @test !isdefined(CellularSheaves.ControlSheaves.Layered, :get_agent_dynamics_config)
end
```

Also **remove** Issue 012's `"root-only binding matches HomogeneousDynamics semantics"` testset —
it references a type that no longer exists. Its intent is preserved by the flat-adapter tests above.

---

## Verification Checklist

- [ ] All three types, `AbstractDynamicsSpec`, and all `get_agent_dynamics_config` methods deleted.
- [ ] Export list in `Layered.jl` updated.
- [ ] `resolve_dynamics(::LayeredEscortSpec, ::SystemBinding)` reuses Issue 012's precedence logic
      rather than reimplementing it — precedence must not be able to drift between the two.
- [ ] `run_layered_escort_simulation` resolves bindings **once** before the time loop, and the
      `team_idx` scan at lines 486–501 is gone.
- [ ] `initial_positions` still works as a documented lower-precedence fallback.
- [ ] `test/network_sheaves/LayeredEscort.jl` migrated; numerical results unchanged.
- [ ] `docs/literate/layered/multilayer_escort.jl` migrated and **the docs build passes**
      (`julia --project=docs docs/make.jl`) — this is a separate check from the test suite.
- [ ] `Agents.md:193`'s dispatch-convention example repointed at a surviving function.
- [ ] Issue 012's `HomogeneousDynamics` comparison testset removed.
- [ ] Full test suite passes.

---

## Out of Scope

- Removing or deprecating `Layered.jl` itself. The flat path keeps its simulation loop, which the
  nested path does not yet have; this issue only changes how it receives dynamics.
- Changing `AgentControllers.jl`'s `AbstractAgentDynamics` hierarchy — only the *assignment*
  mechanism is being replaced, not the dynamics models.
- Porting `LayeredEscortProblem` to `NestedSystemSpec`. That is a larger migration for later.
