"""
    NestedDSLLowering

Turns a validated [`SystemFragment`](@ref) into the concrete `NestedSystems` API: a
[`NestedSystemSpec`](@ref), its compiled [`SheafTower`](@ref), the [`SystemBinding`](@ref)
cascade, and the name → index tables that let downstream code keep talking in names.

Lowering is where names become indices. It walks the fragment tree once, carrying the path taken
so far, and in that single pass:

* orders each node's children (via [`fragment_children`](@ref)) to fix their `RefinedSystem`
  indices;
* rewrites `@link` endpoints into `SystemEdge` index pairs;
* rewrites `@observe` paths — declared *relative* to whatever node they appear in — into
  absolute `Observation.system_path` index vectors, which is what lets one fragment be spliced
  in at different depths without editing;
* hoists every `@target` to the top level and numbers it;
* folds `@bind` terms into the nested `SystemBinding` structure;
* records each node's contiguous block of agent indices, mirroring `NestedSystems`'s own
  depth-first assignment.

That last table is what [`CompiledNestedSystem`](@ref) exposes as [`agent_range`](@ref): the
DSL's answer to hand-maintaining an `agent_index_ranges([4, 5, 2, …])` call whose argument has
to be kept in sync with the tree by hand.
"""
module NestedDSLLowering

using ...NestedSystems: AbstractSystemNode, LeafTeam, RefinedSystem, SystemEdge, TargetSpec,
                        Observation, NestedSystemSpec, SheafTower, build_sheaf_tower,
                        SystemBinding, AgentBinding, NestedEscortProblem,
                        solve_hierarchical, solve_direct, run_nested_escort_simulation
using ..NestedDSLTerm
using ..NestedDSLValidator

export CompiledNestedSystem, compile_nested_system, nested_spec, nested_tower, nested_bindings,
       agent_range, agent_vertices, target_index, target_vector, escort_problem

"""
    CompiledNestedSystem

Everything a [`SystemFragment`](@ref) compiles to, with the name tables kept alongside so that
callers never have to fall back to raw indices.

| Field | Contents |
|---|---|
| `fragment` | the source fragment, retained for introspection |
| `spec` | the [`NestedSystemSpec`](@ref) |
| `tower` | the compiled [`SheafTower`](@ref) |
| `bindings` | the [`SystemBinding`](@ref) cascade built from `@bind` |
| `targets` | target names, in the order they index `spec.targets` |
| `ranges` | dotted path → that node's contiguous block of agent indices |

Use [`agent_range`](@ref), [`agent_vertices`](@ref), [`target_index`](@ref), and
[`target_vector`](@ref) rather than reaching into `ranges`/`targets` directly.
"""
struct CompiledNestedSystem
    fragment::SystemFragment
    spec::NestedSystemSpec
    tower::SheafTower
    bindings::SystemBinding
    targets::Vector{Symbol}
    ranges::Dict{Vector{Symbol},UnitRange{Int}}
end

function Base.show(io::IO, c::CompiledNestedSystem)
    print(io, "CompiledNestedSystem(", length(c.tower.agent_vertices), " agents, ",
          length(c.targets), " targets, tower depth ", c.tower.depth, ")")
end

# ─────────────────────────────────────────────────────────────────────────────
# Entry points
# ─────────────────────────────────────────────────────────────────────────────

"""
    compile_nested_system(f::SystemFragment) -> CompiledNestedSystem

Validate `f`, lower it to a [`NestedSystemSpec`](@ref), compile the [`SheafTower`](@ref), and
resolve its `@bind` declarations — the one call that takes a fragment all the way to something
solvable.

```julia
c = compile_nested_system(@nested_system begin
    @team ring = ring(5; radius=1.0)
    @target t1
    @observe ring => t1
    @bind dynamics=QuadrotorDynamics()
end)

q = solve_hierarchical(c.tower, target_vector(c, Dict(:t1 => [0.0, 0.0, 1.5, 1.0])))
```

Use [`nested_spec`](@ref) instead when the tower is not wanted (it is the expensive part).
"""
function compile_nested_system(f::SystemFragment)
    validate_fragment(f)
    spec, targets, ranges = _lower_spec(f)
    return CompiledNestedSystem(f, spec, build_sheaf_tower(spec), _lower_bindings(f), targets,
                                ranges)
end

"""
    nested_spec(f::SystemFragment) -> NestedSystemSpec

Validate `f` and lower it to a [`NestedSystemSpec`](@ref) without compiling a tower.
"""
function nested_spec(f::SystemFragment)
    validate_fragment(f)
    return first(_lower_spec(f))
end

"""
    nested_tower(f::SystemFragment) -> SheafTower

Shorthand for `build_sheaf_tower(nested_spec(f))`.
"""
nested_tower(f::SystemFragment) = build_sheaf_tower(nested_spec(f))

"""
    nested_bindings(f::SystemFragment) -> SystemBinding

Resolve `f`'s `@bind` declarations into the nested [`SystemBinding`](@ref) cascade, without
building a spec or a tower.
"""
function nested_bindings(f::SystemFragment)
    validate_fragment(f)
    return _lower_bindings(f)
end

# ─────────────────────────────────────────────────────────────────────────────
# Spec lowering
# ─────────────────────────────────────────────────────────────────────────────

function _lower_spec(f::SystemFragment)
    target_names = collect_target_names(f)
    tindex = Dict{Symbol,Int}(n => i for (i, n) in enumerate(target_names))

    root = _lower_node(f, :root)

    observations = Observation[]
    _collect_observations!(observations, f, Int[], tindex)

    spec = NestedSystemSpec(root, [TargetSpec(n) for n in target_names], observations,
                            fragment_dim(f), fragment_affine(f))

    ranges = Dict{Vector{Symbol},UnitRange{Int}}()
    _assign_ranges!(ranges, f, Symbol[], Ref(1))

    return spec, target_names, ranges
end

"""
    _lower_node(f, name) -> RefinedSystem

Build the [`RefinedSystem`](@ref) for one fragment: its children in declaration order, and its
`@link` terms as [`SystemEdge`](@ref)s between those children's indices.
"""
function _lower_node(f::SystemFragment, name::Symbol)
    kids = fragment_children(f)
    index_of = Dict{Symbol,Int}(p.first => i for (i, p) in enumerate(kids))

    children = AbstractSystemNode[_lower_child(node) for (_, node) in kids]

    edges = SystemEdge[]
    for t in terms(f)
        t isa LinkTerm || continue
        push!(edges, SystemEdge(index_of[t.src.path[1]], index_of[t.dst.path[1]],
                                t.src.map, t.dst.map))
    end

    return RefinedSystem(name, children, edges)
end

_lower_child(t::TeamTerm) = LeafTeam(t.name, t.kind, t.n_agents, t.radius, copy(t.observers))
_lower_child(t::SubsystemTerm) = _lower_node(t.body, t.name)

"""
    _collect_observations!(out, f, prefix, tindex)

Walk the tree accumulating [`Observation`](@ref)s. `prefix` is the index path of the node being
visited, so an `@observe` written relative to a deeply nested fragment becomes an absolute path
here — the mechanism that makes a fragment position-independent.
"""
function _collect_observations!(out::Vector{Observation}, f::SystemFragment, prefix::Vector{Int},
                                tindex::Dict{Symbol,Int})
    for t in terms(f)
        t isa ObserveTerm || continue
        _, indices = resolve_fragment_path(f, t.system.path; context="@observe")
        push!(out, Observation([prefix; indices], tindex[t.target]; system_map=t.system.map))
    end
    for (i, (_, node)) in enumerate(fragment_children(f))
        node isa SubsystemTerm && _collect_observations!(out, node.body, [prefix; i], tindex)
    end
end

"""
    _assign_ranges!(ranges, f, path, next)

Record every node's contiguous block of agent indices, walking children depth-first in
declaration order — the same order `NestedSystems` itself assigns agents in, so these ranges
index directly into `SheafTower.agent_vertices`.
"""
function _assign_ranges!(ranges::Dict{Vector{Symbol},UnitRange{Int}}, f::SystemFragment,
                         path::Vector{Symbol}, next::Ref{Int})
    start = next[]
    for (name, node) in fragment_children(f)
        child_path = [path; name]
        if node isa TeamTerm
            ranges[child_path] = next[]:(next[] + node.n_agents - 1)
            next[] += node.n_agents
        else
            _assign_ranges!(ranges, node.body, child_path, next)
        end
    end
    ranges[path] = start:(next[] - 1)
    return ranges
end

# ─────────────────────────────────────────────────────────────────────────────
# Binding lowering
# ─────────────────────────────────────────────────────────────────────────────

"""
    _lower_bindings(f) -> SystemBinding

Flatten every `@bind` into an absolute-path table, then rebuild the tree-shaped
[`SystemBinding`](@ref) from it.

Going through a flat intermediate is what lets a binding be declared anywhere and still land in
the right place: `@bind mid.ringA K_lqr=K` written at the root and `@bind ringA K_lqr=K` written
inside `mid` produce the same absolute path and therefore the same binding.

Where two `@bind` terms set the *same* field of the *same* node, the later declaration wins,
matching how a plain Julia assignment would read.
"""
function _lower_bindings(f::SystemFragment)
    flat = Dict{Tuple{Vector{Symbol},Union{Nothing,Int}},BindTerm}()
    _flatten_binds!(flat, f, Symbol[])
    return _binding_for(f, Symbol[], flat)
end

function _flatten_binds!(flat, f::SystemFragment, prefix::Vector{Symbol})
    for t in terms(f)
        t isa BindTerm || continue
        key = ([prefix; t.path], t.agent)
        flat[key] = haskey(flat, key) ? _fold_bind(flat[key], t) : t
    end
    for (name, node) in fragment_children(f)
        node isa SubsystemTerm && _flatten_binds!(flat, node.body, [prefix; name])
    end
end

_fold_bind(old::BindTerm, new::BindTerm) = BindTerm(
    old.path, old.agent,
    new.dynamics === nothing ? old.dynamics : new.dynamics,
    new.K_lqr === nothing ? old.K_lqr : new.K_lqr,
    new.initial_position === nothing ? old.initial_position : new.initial_position)

function _binding_for(f::SystemFragment, path::Vector{Symbol}, flat)
    here = get(flat, (path, nothing), nothing)
    children = Dict{Symbol,SystemBinding}()
    for (name, node) in fragment_children(f)
        child_path = [path; name]
        sb = node isa TeamTerm ? _team_binding(node, child_path, flat) :
                                 _binding_for(node.body, child_path, flat)
        _is_trivial(sb) || (children[name] = sb)
    end
    return SystemBinding(dynamics = here === nothing ? nothing : here.dynamics,
                         K_lqr = here === nothing ? nothing : here.K_lqr,
                         children = children)
end

function _team_binding(node::TeamTerm, path::Vector{Symbol}, flat)
    here = get(flat, (path, nothing), nothing)
    agents = Dict{Int,AgentBinding}()
    for i in 1:node.n_agents
        b = get(flat, (path, i), nothing)
        b === nothing && continue
        agents[i] = AgentBinding(dynamics = b.dynamics, K_lqr = b.K_lqr,
                                 initial_position = b.initial_position)
    end
    return SystemBinding(dynamics = here === nothing ? nothing : here.dynamics,
                         K_lqr = here === nothing ? nothing : here.K_lqr,
                         agents = agents)
end

"""
    _is_trivial(sb::SystemBinding) -> Bool

Whether a binding declares nothing anywhere in its subtree. Trivial child bindings are dropped
rather than stored, so that `resolve_dynamics` never walks a tree of empty placeholders — and so
that a `SystemBinding` built by the DSL is as small and as readable as a hand-written one.
"""
_is_trivial(sb::SystemBinding) =
    sb.dynamics === nothing && sb.K_lqr === nothing && isempty(sb.agents) && isempty(sb.children)

# ─────────────────────────────────────────────────────────────────────────────
# Name-based lookups
# ─────────────────────────────────────────────────────────────────────────────

"""
    agent_range(c::CompiledNestedSystem, path...) -> UnitRange{Int}

The contiguous block of agent indices belonging to the system at `path`, given as symbols
(`agent_range(c, :mid, :ringA)`), a dotted string (`agent_range(c, "mid.ringA")`), or nothing at
all for the whole tree.

These are *agent* indices — the convention `NestedEscortResult.sim_data[step][a]` uses. To index
`solve_hierarchical`'s output instead, which is keyed by sheaf vertex, go through
[`agent_vertices`](@ref).
"""
function agent_range(c::CompiledNestedSystem, path::Vararg{Union{Symbol,AbstractString}})
    key = _path_key(path)
    haskey(c.ranges, key) || throw(NestedDSLError(
        "no system at `$(isempty(key) ? "<root>" : join(key, "."))` (known systems: " *
        "$(join(sort([join(k, ".") for k in keys(c.ranges) if !isempty(k)]), ", ")))"))
    return c.ranges[key]
end

"""
    agent_vertices(c::CompiledNestedSystem, path...) -> Vector{Int}

The finest-level sheaf vertices of the agents belonging to the system at `path` — that is,
`c.tower.agent_vertices[agent_range(c, path...)]`. These index directly into a
[`solve_hierarchical`](@ref) or [`solve_direct`](@ref) result.
"""
agent_vertices(c::CompiledNestedSystem, path::Vararg{Union{Symbol,AbstractString}}) =
    c.tower.agent_vertices[agent_range(c, path...)]

function _path_key(path)
    out = Symbol[]
    for p in path
        if p isa Symbol
            push!(out, p)
        else
            append!(out, Symbol.(split(String(p), ".")))
        end
    end
    return out
end

"""
    target_index(c::CompiledNestedSystem, name::Symbol) -> Int

The index of the target named `name`, i.e. its position in `c.spec.targets` and its vertex in
`c.tower.target_vertices`.
"""
function target_index(c::CompiledNestedSystem, name::Symbol)
    i = findfirst(==(name), c.targets)
    i === nothing && throw(NestedDSLError(
        "no target named `$name` (declared: $(join(c.targets, ", ")))"))
    return i
end

"""
    target_vector(c::CompiledNestedSystem, values::AbstractDict) -> Vector

Order a `name => value` dictionary into the positional vector that
[`solve_hierarchical`](@ref)/[`solve_direct`](@ref) expect, so that callers can keep referring to
targets by the names they declared.

Every declared target must be present; an unknown name is an error rather than being ignored.
"""
function target_vector(c::CompiledNestedSystem, values::AbstractDict)
    normalized = Dict{Symbol,Any}()
    for (k, v) in values
        name = Symbol(k)
        name in c.targets || throw(NestedDSLError(
            "unknown target `$k` (declared: $(join(c.targets, ", ")))"))
        haskey(normalized, name) && throw(NestedDSLError(
            "target `$name` given more than once (keys collide once normalized to symbols)"))
        normalized[name] = v
    end
    out = Vector{Any}(undef, length(c.targets))
    for (i, name) in enumerate(c.targets)
        haskey(normalized, name) || throw(NestedDSLError(
            "no value given for target `$name` (declared targets: $(join(c.targets, ", ")))"))
        out[i] = normalized[name]
    end
    return identity.(out)
end

"""
    escort_problem(c, trajectories; velocities=nothing, accelerations=nothing, dt=0.05, steps=200)

Build a [`NestedEscortProblem`](@ref) from a compiled system, with target trajectories given as a
`name => (t -> position)` dictionary rather than a positionally-ordered vector.

`velocities` and `accelerations` follow the same keying and enable feedforward tracking; see
[`run_nested_escort_simulation`](@ref).
"""
function escort_problem(c::CompiledNestedSystem, trajectories::AbstractDict;
                        velocities::Union{Nothing,AbstractDict}=nothing,
                        accelerations::Union{Nothing,AbstractDict}=nothing,
                        dt::Real=0.05, steps::Integer=200)
    return NestedEscortProblem(c.tower, c.bindings, target_vector(c, trajectories);
                               target_velocities = velocities === nothing ? nothing :
                                                   target_vector(c, velocities),
                               target_accelerations = accelerations === nothing ? nothing :
                                                      target_vector(c, accelerations),
                               dt = Float64(dt), steps = Int(steps))
end

end # module NestedDSLLowering
