"""
    NestedDSLValidator

Purely **symbolic** checks over a [`SystemFragment`](@ref): duplicate names, unresolvable paths,
member designators out of range, conflicting `@dim`/`@affine`, malformed bindings.

Nothing here computes a rank, a nullspace, or a restriction matrix. That is the same discipline
[`RestrictionSpec`](@ref) follows and it exists for the same reason: a specification error should
be reported against the *declaration the user wrote*, with a path and a name, long before the
tower compiler turns any of it into a matrix. By the time `NestedDSLLowering` runs, every name
in the fragment is known to resolve.

This module also owns the structural view of a fragment — [`fragment_children`](@ref),
[`node_arity`](@ref), [`resolve_fragment_path`](@ref) — which lowering reuses, so the two stages
can never disagree about what a path means or in what order children are numbered.
"""
module NestedDSLValidator

using ...NestedSystems: RestrictionSpec, ProjectMember, Centroid, project, centroid,
                        RefinedSystem, SystemEdge
using ..NestedDSLTerm

export validate_fragment, fragment_children, node_arity, resolve_fragment_path,
       collect_target_names, fragment_dim, fragment_affine

"""
    fragment_children(f::SystemFragment) -> Vector{Pair{Symbol,NestedTerm}}

The fragment's direct children — its [`TeamTerm`](@ref)s and [`SubsystemTerm`](@ref)s — paired
with their names, **in declaration order**.

Declaration order is load-bearing, not incidental: it fixes each child's index in the
[`RefinedSystem`](@ref) that lowering builds, which in turn fixes the depth-first agent
numbering that `SheafTower.agent_vertices` uses. Validation and lowering both go through this
function so the two can never disagree.
"""
function fragment_children(f::SystemFragment)
    out = Pair{Symbol,NestedTerm}[]
    for t in terms(f)
        t isa TeamTerm && push!(out, t.name => t)
        t isa SubsystemTerm && push!(out, t.name => t)
    end
    return out
end

"""
    node_arity(node::NestedTerm) -> Int

Number of a node's **direct members**: raw agents for a [`TeamTerm`](@ref), children for a
[`SubsystemTerm`](@ref). This is the arity a [`RestrictionSpec`](@ref) is declared against — a
refined child counts as one member however many agents it expands to.
"""
node_arity(node::TeamTerm) = node.n_agents
node_arity(node::SubsystemTerm) = length(fragment_children(node.body))

"""
    resolve_fragment_path(f, path; context) -> (node, indices)

Walk `path` from fragment `f`, returning the node it names and the child indices taken to get
there. Throws a [`NestedDSLError`](@ref) naming `context` and listing the available children
when a component does not resolve, or when the path tries to descend into a leaf team.
"""
function resolve_fragment_path(f::SystemFragment, path::AbstractVector{Symbol};
                               context::AbstractString="path")
    isempty(path) && throw(NestedDSLError("$context: empty path"))
    node::Union{Nothing,NestedTerm} = nothing
    frag = f
    indices = Int[]
    for (depth, name) in enumerate(path)
        if node !== nothing
            node isa SubsystemTerm || throw(NestedDSLError(
                "$context: `$(join(path, "."))` descends into `$(join(path[1:depth-1], "."))`, " *
                "which is a team and has no children"))
            frag = node.body
        end
        kids = fragment_children(frag)
        idx = findfirst(p -> p.first == name, kids)
        idx === nothing && throw(NestedDSLError(
            "$context: no system named `$name`" *
            (depth == 1 ? "" : " under `$(join(path[1:depth-1], "."))`") *
            " (available: $(isempty(kids) ? "none" : join([string(p.first) for p in kids], ", ")))"))
        push!(indices, idx)
        node = kids[idx].second
    end
    return node, indices
end

"""
    collect_target_names(f::SystemFragment) -> Vector{Symbol}

Every target declared anywhere in `f`'s tree, in declaration order (depth-first, following
[`fragment_children`](@ref)). Targets are global regardless of nesting depth; see
[`TargetTerm`](@ref).
"""
function collect_target_names(f::SystemFragment)
    names = Symbol[]
    _collect_targets!(names, f)
    return names
end

function _collect_targets!(names::Vector{Symbol}, f::SystemFragment)
    for t in terms(f)
        t isa TargetTerm && push!(names, t.name)
        t isa SubsystemTerm && _collect_targets!(names, t.body)
    end
end

"""
    fragment_dim(f::SystemFragment; default=4) -> Int
    fragment_affine(f::SystemFragment; default=true) -> Bool

The stalk dimension / affine flag declared anywhere in `f`'s tree, or `default` if never
declared. Two declarations of the *same* value are fine (a fragment may reasonably restate the
dimension it was written for); two *different* values are an error.
"""
fragment_dim(f::SystemFragment; default::Int=4) =
    _unique_setting(f, DimTerm, t -> t.D, "@dim", default)

fragment_affine(f::SystemFragment; default::Bool=true) =
    _unique_setting(f, AffineTerm, t -> t.affine, "@affine", default)

function _unique_setting(f::SystemFragment, ::Type{T}, get, label, default) where {T}
    vals = Any[]
    _collect_setting!(vals, f, T, get)
    isempty(vals) && return default
    allequal(vals) || throw(NestedDSLError(
        "$label is declared with conflicting values ($(join(unique(vals), ", "))) — a " *
        "specification has exactly one stalk dimension and one affine flag"))
    return vals[1]
end

function _collect_setting!(vals, f::SystemFragment, ::Type{T}, get) where {T}
    for t in terms(f)
        t isa T && push!(vals, get(t))
        t isa SubsystemTerm && _collect_setting!(vals, t.body, T, get)
    end
end

"""
    validate_fragment(f::SystemFragment) -> SystemFragment

Check `f` as a complete, root-level specification and return it unchanged.

Checks performed, all symbolic:

* every node has at least one child, and no two children of one node share a name;
* `@link` endpoints name **direct children** of the linking node, and differ from each other;
* `@observe` paths resolve, and name a target that is actually declared;
* every [`project`](@ref) designator is in range for the node it is declared against (an
  integer within the node's arity, or the name of an actual child);
* `@bind` paths resolve, agent indices lie within their team, and `[i]` is used only on a team;
* target names are unique and at least one target exists;
* `@dim`/`@affine` do not conflict.

Called automatically by `compile_nested_system`; call it directly to check a fragment without
building anything.
"""
function validate_fragment(f::SystemFragment)
    targets = collect_target_names(f)
    isempty(targets) && throw(NestedDSLError(
        "specification declares no targets — a nested system needs at least one `@target`"))
    _check_duplicate_targets(targets)

    fragment_dim(f)
    fragment_affine(f)

    _validate_node!(f, Symbol[], Set(targets), f)
    return f
end

function _check_duplicate_targets(targets::Vector{Symbol})
    seen = Set{Symbol}()
    for t in targets
        t in seen && throw(NestedDSLError("target `$t` is declared more than once"))
        push!(seen, t)
    end
end

function _validate_node!(f::SystemFragment, path::Vector{Symbol}, targets::Set{Symbol},
                         root::SystemFragment)
    where_ = isempty(path) ? "root" : "system `$(join(path, "."))`"
    kids = fragment_children(f)
    isempty(kids) && throw(NestedDSLError(
        "$where_ has no members — every system needs at least one `@team` or `@system`"))

    seen = Set{Symbol}()
    for (name, _) in kids
        name in seen && throw(NestedDSLError("$where_ declares two children named `$name`"))
        push!(seen, name)
    end

    for t in terms(f)
        t isa LinkTerm && _validate_link(t, f, where_)
        t isa ObserveTerm && _validate_observe(t, f, targets, where_)
        t isa BindTerm && _validate_bind(t, f, where_)
    end

    for (name, node) in kids
        node isa SubsystemTerm && _validate_node!(node.body, [path; name], targets, root)
    end
end

function _validate_link(t::LinkTerm, f::SystemFragment, where_::AbstractString)
    for (label, ref) in (("source", t.src), ("destination", t.dst))
        length(ref.path) == 1 || throw(NestedDSLError(
            "$where_: @link $label `$(join(ref.path, "."))` is not a direct child — an edge " *
            "connects two children of the same system. Give the two systems a common parent " *
            "and declare the link there."))
        node, _ = resolve_fragment_path(f, ref.path; context="$where_ @link $label")
        _validate_spec(ref.map, node, "$where_ @link $label `$(ref.path[1])`")
    end
    t.src.path == t.dst.path && throw(NestedDSLError(
        "$where_: @link connects `$(t.src.path[1])` to itself"))
end

function _validate_observe(t::ObserveTerm, f::SystemFragment, targets::Set{Symbol},
                           where_::AbstractString)
    node, _ = resolve_fragment_path(f, t.system.path; context="$where_ @observe")
    _validate_spec(t.system.map, node, "$where_ @observe `$(join(t.system.path, "."))`")
    t.target in targets || throw(NestedDSLError(
        "$where_: @observe references undeclared target `$(t.target)` (declared targets: " *
        "$(join(sort(collect(targets)), ", ")))"))
end

function _validate_bind(t::BindTerm, f::SystemFragment, where_::AbstractString)
    node = if isempty(t.path)
        nothing
    else
        first(resolve_fragment_path(f, t.path; context="$where_ @bind"))
    end
    if t.agent !== nothing
        node isa TeamTerm || throw(NestedDSLError(
            "$where_: @bind `$(join(t.path, "."))[$(t.agent)]` indexes an agent, but " *
            "$(node === nothing ? "this node" : "`$(join(t.path, "."))`") is not a team — " *
            "per-agent bindings belong on the team that owns the agent"))
        1 <= t.agent <= node.n_agents || throw(NestedDSLError(
            "$where_: @bind `$(join(t.path, "."))[$(t.agent)]` is out of range 1:$(node.n_agents)"))
    elseif t.initial_position !== nothing
        # `SystemBinding` carries dynamics and a gain, but an initial position belongs to one
        # agent — there is nothing sensible for a whole subtree to inherit it as.
        throw(NestedDSLError(
            "$where_: @bind $(isempty(t.path) ? "" : "`$(join(t.path, "."))` ")sets " *
            "initial_position without naming an agent — write `@bind team[i] " *
            "initial_position=…`, since an initial position applies to a single agent"))
    end
end

"""
    _validate_spec(spec, node, where_)

Check a [`RestrictionSpec`](@ref)'s member designator against `node`'s declared arity. Only
[`project`](@ref) carries a designator; [`centroid`](@ref) and `RawRestriction` are checked
numerically later, when they meet a concrete `D`.
"""
function _validate_spec(spec::ProjectMember, node::NestedTerm, where_::AbstractString)
    m = spec.member
    if m isa Int
        n = node_arity(node)
        1 <= m <= n || throw(NestedDSLError(
            "$where_: project($m) is out of range 1:$n" *
            (node isa TeamTerm ? " (the team has $n agents)" :
                                 " (the system has $n children)")))
    else
        node isa SubsystemTerm || throw(NestedDSLError(
            "$where_: project(:$m) selects a child by name, but this is a team whose agents are " *
            "unnamed — use an integer index in 1:$(node_arity(node))"))
        kids = fragment_children(node.body)
        any(p -> p.first == m, kids) || throw(NestedDSLError(
            "$where_: project(:$m) names no child (available: " *
            "$(join([string(p.first) for p in kids], ", ")))"))
    end
end

_validate_spec(::RestrictionSpec, ::NestedTerm, ::AbstractString) = nothing

end # module NestedDSLValidator
