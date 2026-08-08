# Issue 012 — Nested Dynamics Context with Most-Specific-Wins Cascade

**Parent design**: `007-nested-layered-systems-design.md` (§4.7)
**Prerequisites**: 009 (needs the spec tree types)
**Independent of**: 010, 011 — may be implemented in parallel with them
**Blocks**: DSL resolver (Issue 6)

---

## Mathematical / Design Background

Topology and dynamics are orthogonal concerns. A `NestedSystemSpec` (Issue 009) describes *who
agrees with whom*; it says nothing about what an agent flies like. This issue supplies the
second half: binding an `AbstractAgentDynamics`, an LQR gain, and an initial position to every
leaf agent in the tree.

Doing this with one flat entry per leaf agent is intractable for deep trees — a three-level tree
of six-agent rings already has dozens of leaves, nearly all of which share one dynamics model.
Two design decisions (design doc §4.7) address this:

### Cascade with most-specific-wins

Dynamics declared at any node are inherited by every leaf beneath it. Where both an inherited and
a local value apply, the local one wins:

```
per-agent  >  leaf team  >  nearer ancestor  >  …  >  root default
```

This generalizes the existing flat trichotomy in `Layered.jl` — `HomogeneousDynamics` (one model
for all), `TeamHomogeneousDynamics` (per team), `IndividualizedDynamics` (per agent) — from three
fixed tiers to an arbitrary-depth chain. Those three become the degenerate cases of declaring at
the root, at each leaf team, and at each agent respectively.

### Nested context, not flat dotted keys

The binding context mirrors the system tree structurally rather than using flat string keys like
`"systemA.child1.agent[2]"`. It is built from custom types with keyword-argument constructors, so
malformed input fails at construction rather than silently binding nothing — a typo'd child name
is an error, not a missing binding that surfaces as a confusing solve-time failure.

**Target trajectories do not cascade.** They are inherently per-target and are supplied
separately at solve time; this issue does not touch them.

---

## Codebase Orientation

| File | Why |
|---|---|
| `src/ControlSheaves/AgentControllers.jl` | `AbstractAgentDynamics` (line 14), `initial_state(dyn, position)` (line 119), `AgentState(x0, dyn, dt, K_lqr, eps; use_velocity)` (line 202). The resolved bindings must be directly usable to construct these. |
| `src/ControlSheaves/Layered.jl:33–68` | `AbstractDynamicsSpec` and the `HomogeneousDynamics`/`TeamHomogeneousDynamics`/`IndividualizedDynamics` trichotomy, plus `get_agent_dynamics_config` dispatch. This issue generalizes that pattern; read it before designing. |
| `src/ControlSheaves/NestedSystems.jl` | Issue 009's `LeafTeam`, `RefinedSystem`, `NestedSystemSpec`. The context must mirror this tree and be validated against it. |
| `src/ControlSheaves/TrackingDSL/Resolver.jl` | The repo's existing late-binding precedent — read for the staging discipline, though it uses a flat `env::Dict` which §4.7 explicitly moves away from here. |

---

## Requested Implementation

Extend `src/ControlSheaves/NestedSystems.jl` (or a sibling `DynamicsBinding.jl` included from it).

```julia
"""
    AgentBinding(; dynamics=nothing, K_lqr=nothing, initial_position=nothing)

Per-agent dynamics override. Any field left `nothing` falls back to the enclosing team's
binding, then to successively more distant ancestors — see [`resolve_dynamics`](@ref).
Fields are resolved independently: an agent may override only its initial position while
inheriting its dynamics model.
"""
Base.@kwdef struct AgentBinding
    dynamics::Union{Nothing,AbstractAgentDynamics} = nothing
    K_lqr::Union{Nothing,Matrix{Float64}} = nothing
    initial_position::Union{Nothing,Vector{Float64}} = nothing
end

"""
    SystemBinding(; dynamics=nothing, K_lqr=nothing, children=Dict(), agents=Dict())

Dynamics bindings for one node of the system tree, mirroring the structure of the
`NestedSystemSpec` it is resolved against.

`dynamics`/`K_lqr` apply to every leaf agent in this subtree unless a descendant overrides them.
`children` maps a child system's name to its own `SystemBinding`. `agents` maps a local agent
index within a `LeafTeam` to an `AgentBinding`.

Child names are validated against the spec during [`resolve_dynamics`](@ref): a name present here
but absent from the spec is an error, so typos surface immediately.
"""
Base.@kwdef struct SystemBinding
    dynamics::Union{Nothing,AbstractAgentDynamics} = nothing
    K_lqr::Union{Nothing,Matrix{Float64}} = nothing
    children::Dict{Symbol,SystemBinding} = Dict{Symbol,SystemBinding}()
    agents::Dict{Int,AgentBinding} = Dict{Int,AgentBinding}()
end

"""
    ResolvedAgent

A fully-bound leaf agent: its global agent index, its dotted path in the tree (for error
messages and debugging), and the dynamics, gain, and initial position that won the cascade.
"""
struct ResolvedAgent
    agent_index::Int
    path::Vector{Symbol}
    dynamics::AbstractAgentDynamics
    K_lqr::Matrix{Float64}
    initial_position::Vector{Float64}
end

"""
    resolve_dynamics(spec::NestedSystemSpec, ctx::SystemBinding) -> Vector{ResolvedAgent}

Walk `spec`'s tree carrying the inherited binding, overriding it wherever `ctx` supplies a more
specific value, and return one `ResolvedAgent` per leaf agent in global agent-index order.

Precedence, most specific first: per-agent, leaf team, nearer ancestor, root default. Each field
resolves independently.

Throws if any agent ends with no dynamics bound anywhere up its chain, or if `ctx` names a child
or agent index that does not exist in `spec`. Both errors report the full dotted path.
"""
function resolve_dynamics(spec::NestedSystemSpec, ctx::SystemBinding)
```

### Algorithm sketch

1. Recurse from `spec.root` alongside `ctx`, carrying an `inherited::AgentBinding`.
2. At each node, fold the node's own `dynamics`/`K_lqr` over `inherited` — a non-`nothing` local
   field replaces the inherited one, a `nothing` field leaves it intact. Fold fields
   independently; do not replace the whole binding wholesale.
3. **Validate as you descend**: every key in `ctx.children` must name a child of the current
   `RefinedSystem`, and every key in `ctx.agents` must be within the current `LeafTeam`'s
   `1:n_agents`. Throw with the dotted path on a mismatch.
4. At a `LeafTeam`, fold each `agents[i]` binding over the accumulated inherited one and emit a
   `ResolvedAgent`.
5. If `dynamics` is still `nothing` at emission, throw naming the agent's full path. If `K_lqr`
   is `nothing`, default to `zeros(0,0)`, matching `HomogeneousDynamics`'s existing convention
   (`Layered.jl:45`). If `initial_position` is `nothing`, fall back to the "airstrip" default in
   `Layered.jl:464` — reuse that function rather than duplicating the layout.
6. Emit in global agent-index order so the result indexes directly against the tower's
   `agent_vertices` from Issue 009.

### Convenience constructors

Provide shorthands so the common cases stay terse:

```julia
homogeneous_binding(dyn; K_lqr=zeros(0,0)) = SystemBinding(dynamics=dyn, K_lqr=K_lqr)
```

A root-level `SystemBinding(dynamics=dyn)` must reproduce `HomogeneousDynamics(dyn)` semantics
exactly — that equivalence is a test below.

---

## Tests to Write

Add to `test/ControlSheaves/NestedSystems.jl`:

```julia
@testset "root default cascades to every leaf" begin
    spec = three_level_spec()
    ctx = SystemBinding(dynamics=QuadrotorDynamics())
    resolved = resolve_dynamics(spec, ctx)
    @test length(resolved) == n_agents(spec)
    @test all(r -> r.dynamics isa QuadrotorDynamics, resolved)
end

@testset "most specific wins: agent over team over ancestor" begin
    spec = two_level_spec()
    ctx = SystemBinding(
        dynamics = QuadrotorDynamics(),
        children = Dict(:teamA => SystemBinding(
            dynamics = PlanarQuadrotorDynamics(),
            agents = Dict(2 => AgentBinding(dynamics = QuadrotorDynamics())))))
    r = resolve_dynamics(spec, ctx)
    @test r[1].dynamics isa PlanarQuadrotorDynamics   # team default
    @test r[2].dynamics isa QuadrotorDynamics         # per-agent override
    @test last(r).dynamics isa QuadrotorDynamics      # root default, other subtree
end

@testset "fields resolve independently" begin
    # override only the initial position; dynamics must still be inherited
    spec = two_level_spec()
    ctx = SystemBinding(dynamics=QuadrotorDynamics(),
                        children=Dict(:teamA => SystemBinding(
                            agents=Dict(1 => AgentBinding(initial_position=[9.0, 9.0, 9.0])))))
    r = resolve_dynamics(spec, ctx)
    @test r[1].dynamics isa QuadrotorDynamics
    @test r[1].initial_position == [9.0, 9.0, 9.0]
end

@testset "unbound dynamics throws naming the path" begin
    spec = two_level_spec()
    err = try resolve_dynamics(spec, SystemBinding()); nothing catch e; e end
    @test err !== nothing
    @test occursin("teamA", sprint(showerror, err))
end

@testset "typo'd child name is rejected" begin
    spec = two_level_spec()
    ctx = SystemBinding(dynamics=QuadrotorDynamics(),
                        children=Dict(:teamTypo => SystemBinding()))
    @test_throws Exception resolve_dynamics(spec, ctx)
end

@testset "out-of-range agent index is rejected" begin
    spec = two_level_spec()   # teamA has 6 agents
    ctx = SystemBinding(dynamics=QuadrotorDynamics(),
                        children=Dict(:teamA => SystemBinding(
                            agents=Dict(99 => AgentBinding()))))
    @test_throws Exception resolve_dynamics(spec, ctx)
end

@testset "root-only binding matches HomogeneousDynamics semantics" begin
    spec = two_level_spec()
    dyn = QuadrotorDynamics()
    r = resolve_dynamics(spec, SystemBinding(dynamics=dyn))
    old = HomogeneousDynamics(dyn)
    @test all(i -> r[i].dynamics === get_agent_dynamics_config(old, i, 1)[1], eachindex(r))
end

@testset "resolved order matches the tower's agent ordering" begin
    spec = three_level_spec()
    tower = build_sheaf_tower(spec)
    r = resolve_dynamics(spec, SystemBinding(dynamics=QuadrotorDynamics()))
    @test [x.agent_index for x in r] == collect(eachindex(tower.agent_vertices))
end
```

---

## Verification Checklist

- [ ] `Base.@kwdef` used so every field has a documented default and construction validates input.
- [ ] Fields resolve **independently** — overriding `initial_position` must not clear `dynamics`.
- [ ] Unknown child names and out-of-range agent indices throw, with the dotted path in the message.
- [ ] Unbound dynamics throws rather than silently defaulting.
- [ ] `K_lqr` defaults to `zeros(0,0)`, matching `Layered.jl:45`.
- [ ] Initial-position fallback reuses `Layered.jl`'s `_default_initial_position`, not a copy.
- [ ] Resolved agents come back in global agent-index order, aligned with Issue 009's
      `agent_vertices`.
- [ ] Root-only binding is behaviourally identical to `HomogeneousDynamics`.
- [ ] `@argcheck` at public boundaries; docstrings explain the precedence rule explicitly.
- [ ] Full test suite passes.

---

## Out of Scope

- Target trajectories — they do not cascade and are supplied at solve time.
- Any DSL surface syntax for declaring bindings; this is the Julia API the DSL resolver will
  target.
- Running simulations or stepping agents — this issue produces bindings, not trajectories.
- Removing `HomogeneousDynamics`/`TeamHomogeneousDynamics`/`IndividualizedDynamics`. They stay in
  place for the duration of this issue so the flat `Layered` path keeps working; **Issue 013**
  retires them once this cascade exists to replace them. The "root-only binding matches
  `HomogeneousDynamics` semantics" testset below is therefore temporary — Issue 013 deletes it.
