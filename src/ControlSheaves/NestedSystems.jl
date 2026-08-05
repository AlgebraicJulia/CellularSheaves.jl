"""
    NestedSystems

Nested system specification and sheaf **tower** compiler (see
`docs/issues/007-nested-layered-systems-design.md` and
`docs/issues/009-nested-spec-and-tower-compiler.md`).

A nested specification is a tree of teams (raw agents in a rigid formation, see
[`build_escort_topology`](@ref)) and refined subsystems, together with a set of top-level
targets and an arbitrary many-to-many observation incidence between systems and targets.

[`build_sheaf_tower`](@ref) compiles this specification into a **tower** of cellular sheaves

```
H₀  (coarsest)
 ↑ f₁
H₁
 ↑ f₂
 ⋮
H_N (finest — raw agents and targets)
```

connected by graph homomorphisms `f_k`, each level the [`pushforward_sheaf`](@ref) of the one
below it. Because every homomorphism maps each target vertex to its own vertex at the coarser
level (targets are *singleton fibres* at every level — see the module docstring's parent design
doc, §4.2), a target's stalk stays `D`-dimensional and its fibre-section basis stays the identity
all the way to `H₀`. This is what lets the eventual harmonic-extension solve (Issue 010) pin
targets directly, with no pseudo-inverse anywhere in this module.
"""
module NestedSystems

using LinearAlgebra
using Graphs
using ArgCheck
using ...NetworkSheaves: vertex_stalks, get_restriction_map, EuclideanSheaf, underlying_graph,
                         coboundary_map
using ...NetworkSheaves.EuclideanSheaves: add_sheaf_edge!, harmonic_extension
using ...NetworkSheaves.Formations: build_escort_topology
using ...NetworkSheaves.GraphHomomorphisms: GraphHomomorphism, fiber_vertices
using ...NetworkSheaves.Pushforwards: pushforward_sheaf, all_fiber_bases

export AbstractSystemNode, LeafTeam, RefinedSystem, TargetSpec, Observation, NestedSystemSpec,
       SheafTower, build_sheaf_tower,
       solve_hierarchical, solve_direct, sheaf_energy, approximation_gap

# ---------------------------------------------------------------------------
# Specification types
# ---------------------------------------------------------------------------

"""
    AbstractSystemNode

Common supertype of the nested-system tree nodes: [`LeafTeam`](@ref) (raw agents) and
[`RefinedSystem`](@ref) (a subsystem whose members are themselves systems).
"""
abstract type AbstractSystemNode end

"""
    LeafTeam(name, kind, n_agents, radius; observers=[1])

A team of raw agents in a `kind` formation (see [`build_escort_topology`](@ref)). Leaves are
where actual agents live; every other node in the tree is structural.

`kind` is one of `:ring`, `:path`, `:star`, `:clique`. `observers` names which of the team's
own agents (local indices `1:n_agents`) are used to fix the team's internal geometry; it is
independent of which target(s) the team observes — that incidence is declared separately via
[`Observation`](@ref).
"""
struct LeafTeam <: AbstractSystemNode
    name::Symbol
    kind::Symbol
    n_agents::Int
    radius::Float64
    observers::Vector{Int}

    function LeafTeam(name::Symbol, kind::Symbol, n_agents::Int, radius::Real, observers::Vector{Int})
        @argcheck kind in (:ring, :path, :star, :clique) "kind must be one of :ring, :path, :star, :clique (got :$kind)"
        @argcheck n_agents >= 1 "n_agents must be positive (got $n_agents)"
        @argcheck all(1 <= o <= n_agents for o in observers) "observers must be within 1:n_agents"
        new(name, kind, n_agents, Float64(radius), observers)
    end
end

LeafTeam(name::Symbol, kind::Symbol, n_agents::Int, radius::Real; observers=[1]) =
    LeafTeam(name, kind, n_agents, radius, collect(Int, observers))

"""
    RefinedSystem(name, children, internal_edges=Tuple{Int,Int}[])

A subsystem whose vertices are themselves systems. `internal_edges` are consensus edges among
`children`, given as `(i, j)` index pairs into `children`.
"""
struct RefinedSystem <: AbstractSystemNode
    name::Symbol
    children::Vector{AbstractSystemNode}
    internal_edges::Vector{Tuple{Int,Int}}

    function RefinedSystem(name::Symbol, children::Vector{AbstractSystemNode}, internal_edges::Vector{Tuple{Int,Int}})
        @argcheck !isempty(children) "RefinedSystem $(name) must have at least one child"
        n = length(children)
        for (i, j) in internal_edges
            @argcheck 1 <= i <= n "internal_edges index $i out of range 1:$n for RefinedSystem $(name)"
            @argcheck 1 <= j <= n "internal_edges index $j out of range 1:$n for RefinedSystem $(name)"
            @argcheck i != j "internal_edges cannot connect child $i to itself"
        end
        new(name, children, internal_edges)
    end
end

RefinedSystem(name::Symbol, children::Vector{<:AbstractSystemNode}, internal_edges::Vector{Tuple{Int,Int}}=Tuple{Int,Int}[]) =
    RefinedSystem(name, Vector{AbstractSystemNode}(children), internal_edges)

"""
    TargetSpec(name)

An uncontrolled agent supplying a boundary condition at solve time. Targets are top-level
vertices and are never owned by a team.
"""
struct TargetSpec
    name::Symbol
end

"""
    Observation(system_path, target_index)

Declares that the system at `system_path` (a path of child indices from the root) observes
target `target_index`. Arbitrary many-to-many incidence is allowed: one system may observe
several targets and one target may be observed by several systems.
"""
struct Observation
    system_path::Vector{Int}
    target_index::Int
end

"""
    NestedSystemSpec(root, targets, observations, D, affine)

Full specification for a nested layered-control system: a tree of teams and subsystems
(`root`), a set of top-level [`TargetSpec`](@ref)s, the [`Observation`](@ref) incidence between
systems and targets, the stalk dimension `D`, and whether restriction maps are affine (see
[`build_escort_topology`](@ref)).
"""
struct NestedSystemSpec
    root::RefinedSystem
    targets::Vector{TargetSpec}
    observations::Vector{Observation}
    D::Int
    affine::Bool

    function NestedSystemSpec(root::RefinedSystem, targets::Vector{TargetSpec},
                               observations::Vector{Observation}, D::Int, affine::Bool)
        @argcheck D >= 1 "D must be positive (got $D)"
        @argcheck !isempty(targets) "at least one target is required"
        for obs in observations
            @argcheck !isempty(obs.system_path) "Observation.system_path must be non-empty"
            @argcheck 1 <= obs.target_index <= length(targets) "Observation.target_index $(obs.target_index) out of range 1:$(length(targets))"
            _resolve_path(root, obs.system_path)  # throws if the path does not resolve
        end
        new(root, targets, observations, D, affine)
    end
end

# ---------------------------------------------------------------------------
# Tree bookkeeping helpers
# ---------------------------------------------------------------------------

"""
    _resolve_path(root::RefinedSystem, path::Vector{Int}) -> AbstractSystemNode

Walk `path` (a sequence of child indices) from `root`, returning the node it addresses.
Throws an `ArgumentError` if the path is malformed (out-of-range index, or it tries to
descend into a `LeafTeam`, which has no children).
"""
function _resolve_path(root::RefinedSystem, path::Vector{Int})
    node::AbstractSystemNode = root
    for idx in path
        @argcheck node isa RefinedSystem "system_path descends past a LeafTeam (no children to index into)"
        @argcheck 1 <= idx <= length(node.children) "system_path index $idx out of range 1:$(length(node.children))"
        node = node.children[idx]
    end
    return node
end

"""
    _collapse_depth(node) -> Int

Number of pushforward steps needed for `node`'s subtree to collapse into a single vertex,
counting from the finest level `H_N` (where `_collapse_depth(::LeafTeam) == 1`, since raw
agents always collapse into their team's centre in one step). A `RefinedSystem` must wait for
its slowest child before it, too, can collapse: `1 + maximum(_collapse_depth, children)`.
Memoized in `S` (keyed by node identity, since two structurally-equal `LeafTeam`s are still
distinct tree nodes).
"""
function _collapse_depth!(node::LeafTeam, S::IdDict{AbstractSystemNode,Int})
    S[node] = 1
    return 1
end

function _collapse_depth!(node::RefinedSystem, S::IdDict{AbstractSystemNode,Int})
    haskey(S, node) && return S[node]
    s = 1 + maximum(_collapse_depth!(c, S) for c in node.children)
    S[node] = s
    return s
end

"""
    _build_parent_map!(node::RefinedSystem, parent_of)

Populate `parent_of[c] = node` for every child `c` of `node`, recursively. Root's direct
children get an entry too (pointing at the root), but it is never dereferenced: every root
child is guaranteed collapsed (a singleton slot) by `H₀`, so climbing (see
`_owner_slot_index`) never needs to pass through it.
"""
function _build_parent_map!(node::RefinedSystem, parent_of::IdDict{AbstractSystemNode,AbstractSystemNode})
    for c in node.children
        parent_of[c] = node
        c isa RefinedSystem && _build_parent_map!(c, parent_of)
    end
end

"""
    _assign_agent_ranges!(node, next_idx, agent_range, agent_owner)

Depth-first walk assigning consecutive `H_N` vertex indices to every raw agent. `next_idx` is a
running `Ref{Int}` counter (already offset past the target block). Records each node's
(contiguous) vertex range in `agent_range`, and appends each raw agent's owning `LeafTeam` to
`agent_owner` in visitation order.
"""
function _assign_agent_ranges!(node::LeafTeam, next_idx::Ref{Int},
                                agent_range::IdDict{AbstractSystemNode,UnitRange{Int}},
                                agent_owner::Vector{AbstractSystemNode})
    start = next_idx[]
    for _ in 1:node.n_agents
        push!(agent_owner, node)
    end
    next_idx[] += node.n_agents
    r = start:(next_idx[] - 1)
    agent_range[node] = r
    return r
end

function _assign_agent_ranges!(node::RefinedSystem, next_idx::Ref{Int},
                                agent_range::IdDict{AbstractSystemNode,UnitRange{Int}},
                                agent_owner::Vector{AbstractSystemNode})
    start = next_idx[]
    for c in node.children
        _assign_agent_ranges!(c, next_idx, agent_range, agent_owner)
    end
    r = start:(next_idx[] - 1)
    agent_range[node] = r
    return r
end

"""
    _level_owner_slots(node::RefinedSystem, level, depth, S) -> Vector{AbstractSystemNode}

Ordered list of the nodes that are represented as their own single vertex ("atomic") at
tower level `level` (coarser levels are smaller `level`), obtained by recursing into any
child that has not yet collapsed (`S[child] > depth - level`). Used to build the vertex set
of every pushforward level below `H_N`.
"""
function _level_owner_slots(node::RefinedSystem, level::Int, depth::Int, S::IdDict{AbstractSystemNode,Int})
    out = AbstractSystemNode[]
    for c in node.children
        if S[c] <= depth - level
            push!(out, c)
        else
            append!(out, _level_owner_slots(c, level, depth, S))
        end
    end
    return out
end

"""
    _owner_slot_index(node, pos_dict, parent_of) -> Int

Climb from `node` up through `parent_of` until reaching a node present in `pos_dict` (i.e. one
that is atomic at the target level), and return its position. Every node is guaranteed to find
such an ancestor-or-self without ever climbing past a root child (see `_build_parent_map!`).
"""
function _owner_slot_index(node::AbstractSystemNode,
                            pos_dict::IdDict{AbstractSystemNode,Int},
                            parent_of::IdDict{AbstractSystemNode,AbstractSystemNode})
    cur = node
    while !haskey(pos_dict, cur)
        cur = parent_of[cur]
    end
    return pos_dict[cur]
end

# ---------------------------------------------------------------------------
# Building H_N (the finest level)
# ---------------------------------------------------------------------------

"""
    _add_leafteam_edges!(F, leaf, D, affine, vertex_offset)

Build `leaf`'s formation via [`build_escort_topology`](@ref) and transplant its *consensus*
edges (agent-to-agent) into `F` at `vertex_offset`, mirroring `Layered.jl`'s ring transplant.
`build_escort_topology` also generates edges pinning `leaf.observers` to a synthetic local
target vertex (`n_agents + 1`); those are dropped here rather than transplanted, because a
`LeafTeam` has no target of its own in the nested architecture — teams are pinned to real
targets only via [`Observation`](@ref) edges. The consensus edges alone already give the team
an exact `D`-dimensional global-section space for any connected `kind` (see
`build_escort_topology`'s docstring), so dropping the synthetic pin does not affect rigidity.
"""
function _add_leafteam_edges!(F::EuclideanSheaf, leaf::LeafTeam, D::Int, affine::Bool, vertex_offset::Int)
    synthetic_target = leaf.n_agents + 1
    local_sheaf = build_escort_topology(leaf.kind, leaf.n_agents, synthetic_target, leaf.radius;
                                         observers=leaf.observers, D=D, affine=affine)
    for e in edges(underlying_graph(local_sheaf))
        u_loc, v_loc = src(e), dst(e)
        (u_loc > leaf.n_agents || v_loc > leaf.n_agents) && continue
        u, v = vertex_offset + u_loc, vertex_offset + v_loc
        add_sheaf_edge!(F, u, v,
                        get_restriction_map(local_sheaf, u_loc, v_loc),
                        get_restriction_map(local_sheaf, v_loc, u_loc))
    end
end

function _add_all_leafteam_edges!(F::EuclideanSheaf, node::RefinedSystem, D::Int, affine::Bool,
                                   agent_range::IdDict{AbstractSystemNode,UnitRange{Int}})
    for c in node.children
        if c isa LeafTeam
            _add_leafteam_edges!(F, c, D, affine, first(agent_range[c]) - 1)
        else
            _add_all_leafteam_edges!(F, c, D, affine, agent_range)
        end
    end
end

"""
    _add_internal_edges!(F, node, agent_range, D)

Add a consensus edge for every `(i, j)` in `node.internal_edges`, recursively over every
`RefinedSystem` in the tree (not just the root).

TODO(011): the representative vertex used for each endpoint is the child's first agent —
equivalent to `project(1)` — because declared per-edge restriction maps (`project`/`centroid`)
are Issue 011's territory, not this one's. This is the one place this issue touches that
territory; see the "Open detail" section of `docs/issues/009-nested-spec-and-tower-compiler.md`.
"""
function _add_internal_edges!(F::EuclideanSheaf, node::RefinedSystem,
                               agent_range::IdDict{AbstractSystemNode,UnitRange{Int}}, D::Int)
    I_D = Matrix{Float64}(I, D, D)
    for (i, j) in node.internal_edges
        u = first(agent_range[node.children[i]])  # TODO(011): project(1) placeholder
        v = first(agent_range[node.children[j]])  # TODO(011): project(1) placeholder
        add_sheaf_edge!(F, u, v, I_D, I_D)
    end
    for c in node.children
        c isa RefinedSystem && _add_internal_edges!(F, c, agent_range, D)
    end
end

"""
    _add_observation_edges!(F, spec, agent_range)

Add an edge from each [`Observation`](@ref)'s representative vertex (its system's first agent,
same TODO(011) caveat as `_add_internal_edges!`) to its target vertex.
"""
function _add_observation_edges!(F::EuclideanSheaf, spec::NestedSystemSpec,
                                  agent_range::IdDict{AbstractSystemNode,UnitRange{Int}})
    I_D = Matrix{Float64}(I, spec.D, spec.D)
    for obs in spec.observations
        node = _resolve_path(spec.root, obs.system_path)
        rep = first(agent_range[node])  # TODO(011): project(1) placeholder
        add_sheaf_edge!(F, rep, obs.target_index, I_D, I_D)
    end
end

# ---------------------------------------------------------------------------
# Tower type and compiler
# ---------------------------------------------------------------------------

"""
    SheafTower

Compiled tower of sheaves. `levels[1]` is the coarsest sheaf `H₀`; `levels[end]` is the finest
`H_N` containing one vertex per raw agent plus one per target. `homs[k] : levels[k+1] → levels[k]`,
and `bases[k][v]` is the fibre-section basis over vertex `v` of `levels[k]`, used to lift a value
at `v` down into `levels[k+1]`.

Targets occupy vertices `1:length(spec.targets)` at **every** level (including `H_N`), so
`target_vertices[t]` is target `t`'s vertex index and is valid at every level — this is what
makes the singleton-fibre invariant checkable and usable uniformly across the tower.
`agent_vertices[a]` is agent `a`'s vertex index in `levels[end]`.

Fibre bases (`bases`) are computed once during compilation and cached here; they are never
recomputed by downstream consumers (e.g. the Issue 010 solver).
"""
struct SheafTower
    spec::NestedSystemSpec
    levels::Vector{EuclideanSheaf{Float64}}
    homs::Vector{GraphHomomorphism}
    bases::Vector{Vector{Matrix{Float64}}}
    target_vertices::Vector{Int}
    agent_vertices::Vector{Int}
    depth::Int
end

"""
    build_sheaf_tower(spec::NestedSystemSpec) -> SheafTower

Compile a nested specification into a tower of sheaves connected by pushforwards.

Algorithm: flatten the tree to the finest level `H_N` (targets at vertices `1:n_targets`, raw
agents consecutively after, per team and per subsystem — see `_assign_agent_ranges!`), then
repeatedly [`pushforward_sheaf`](@ref) up the tower. At each step, every node whose entire
subtree has already collapsed becomes a single coarse vertex; targets always map to themselves
(singleton fibres, identity basis), which is asserted after every pushforward and fails loudly
if violated — this is the invariant the whole hierarchical architecture rests on (see the
module docstring). A fibre whose section space collapses to dimension 0 means the corresponding
team or subsystem is over-constrained and the tower cannot represent it; this is also rejected.
"""
function build_sheaf_tower(spec::NestedSystemSpec)
    root = spec.root
    D = spec.D
    n_targets = length(spec.targets)

    S = IdDict{AbstractSystemNode,Int}()
    _collapse_depth!(root, S)
    parent_of = IdDict{AbstractSystemNode,AbstractSystemNode}()
    _build_parent_map!(root, parent_of)

    depth = 1 + maximum(S[c] for c in root.children)

    agent_range = IdDict{AbstractSystemNode,UnitRange{Int}}()
    agent_owner = AbstractSystemNode[]
    next_idx = Ref(n_targets + 1)
    for c in root.children
        _assign_agent_ranges!(c, next_idx, agent_range, agent_owner)
    end
    n_agents = next_idx[] - 1 - n_targets
    total_HN = n_targets + n_agents

    # --- Build H_N ---
    F = EuclideanSheaf{Float64}(fill(D, total_HN))
    _add_all_leafteam_edges!(F, root, D, spec.affine, agent_range)
    _add_internal_edges!(F, root, agent_range, D)
    _add_observation_edges!(F, spec, agent_range)

    levels = Vector{EuclideanSheaf{Float64}}(undef, depth)
    homs = Vector{GraphHomomorphism}(undef, depth - 1)
    bases = Vector{Vector{Matrix{Float64}}}(undef, depth - 1)
    levels[depth] = F

    # Ordered non-target slot lists (and position lookups) for every coarser level.
    level_slots = Vector{Vector{AbstractSystemNode}}(undef, depth - 1)
    level_pos = Vector{IdDict{AbstractSystemNode,Int}}(undef, depth - 1)
    for k in 1:(depth - 1)
        lst = _level_owner_slots(root, k, depth, S)
        level_slots[k] = lst
        level_pos[k] = IdDict{AbstractSystemNode,Int}(n => i for (i, n) in enumerate(lst))
    end

    for k in (depth - 1):-1:1
        fine_level = levels[k + 1]
        n_fine = length(vertex_stalks(fine_level))
        vmap = Vector{Int}(undef, n_fine)

        for t in 1:n_targets
            vmap[t] = t
        end

        if k + 1 == depth
            for a in 1:n_agents
                owner = agent_owner[a]
                pos = _owner_slot_index(owner, level_pos[k], parent_of)
                vmap[n_targets + a] = n_targets + pos
            end
        else
            src_list = level_slots[k + 1]
            for (i, node) in enumerate(src_list)
                pos = _owner_slot_index(node, level_pos[k], parent_of)
                vmap[n_targets + i] = n_targets + pos
            end
        end

        n_target_coarse = n_targets + length(level_slots[k])
        hom = GraphHomomorphism(vmap, n_target_coarse)
        homs[k] = hom
        levels[k] = pushforward_sheaf(hom, fine_level)
        bases[k] = all_fiber_bases(hom, fine_level)

        # --- Assert the target invariant: singleton fibre, identity basis, stalk D. ---
        for t in 1:n_targets
            fverts = fiber_vertices(hom, t)
            length(fverts) == 1 || error(
                "SheafTower invariant violated: target $t is not a singleton fibre at level $k " *
                "(fibre vertices $fverts). Targets must map to their own vertex at every level.")
            stalk = vertex_stalks(levels[k])[t]
            stalk == D || error(
                "SheafTower invariant violated: target $t has stalk $stalk at level $k, expected $D.")
            B = bases[k][t]
            (size(B) == (D, D) && isapprox(B, Matrix{Float64}(I, D, D))) || error(
                "SheafTower invariant violated: target $t's fibre basis at level $k is not the " *
                "identity — a singleton fibre must have an identity basis.")
        end

        # --- Reject over-constrained (0-dimensional) fibres. ---
        for (v, B) in enumerate(bases[k])
            size(B, 2) == 0 && error(
                "SheafTower: fibre at coarse vertex $v (level $k) has a 0-dimensional global " *
                "section space — the corresponding team/subsystem is over-constrained and " *
                "cannot be represented by the tower.")
        end
    end

    target_vertices = collect(1:n_targets)
    agent_vertices = collect((n_targets + 1):(n_targets + n_agents))

    return SheafTower(spec, levels, homs, bases, target_vertices, agent_vertices, depth)
end

# ---------------------------------------------------------------------------
# Hierarchical & direct solves, and the approximation gap between them
# ---------------------------------------------------------------------------

"""
    _per_vertex(F::EuclideanSheaf, x::AbstractVector) -> Vector{Vector{Float64}}

Split a flat cochain `x` on `F` into one value per vertex, using each vertex's own stalk
dimension (stalks are **not** uniform above `H_N`).
"""
function _per_vertex(F::EuclideanSheaf, x::AbstractVector)
    dims = vertex_stalks(F)
    out = Vector{Vector{Float64}}(undef, length(dims))
    off = 0
    for v in eachindex(dims)
        out[v] = Vector{Float64}(x[(off + 1):(off + dims[v])])
        off += dims[v]
    end
    return out
end

"""
    _boundary_dict(tower, target_values) -> Dict{Int,Vector{Float64}}

Map each target's vertex index to its pinned value, checking arity and stalk dimension.
Targets share the same vertex index at every level of the tower, so the resulting dictionary
is a valid boundary condition for `tower.levels[1]` and `tower.levels[end]` alike.
"""
function _boundary_dict(tower::SheafTower, target_values::AbstractVector)
    @argcheck length(target_values) == length(tower.target_vertices) "expected $(length(tower.target_vertices)) target values, got $(length(target_values))"
    D = tower.spec.D
    for (t, val) in enumerate(target_values)
        @argcheck length(val) == D "target $t value has length $(length(val)), expected stalk dimension $D"
    end
    return Dict{Int,Vector{Float64}}(
        tower.target_vertices[t] => Vector{Float64}(target_values[t])
        for t in eachindex(target_values))
end

"""
    solve_hierarchical(tower::SheafTower, target_values) -> Vector{Vector{Vector{Float64}}}

Solve the tower top-down: one harmonic extension on the coarsest sheaf `tower.levels[1]` with
`target_values` pinned at the target vertices, then successive fibre-basis lifts down to each
finer level.

Because every fibre basis `B_v` spans the fibre's *exact* global-section space, each lift is
energy-preserving: the internal edges of a fibre contribute exactly zero, and the cross-edge
energy of the coarse sheaf equals that of the lifted cochain on the finer one. The coarse solve
therefore genuinely minimises finest-level energy — but only over configurations in which every
team stays rigidly in formation. See [`approximation_gap`](@ref) for what that restriction costs.

Returns one cochain per level, coarsest first, each as a vector of per-vertex values (per-vertex,
because stalk dimensions differ between vertices above the finest level). The last entry holds
the per-agent and per-target reference states on `tower.levels[end]`.
"""
function solve_hierarchical(tower::SheafTower, target_values::AbstractVector)
    boundary = _boundary_dict(tower, target_values)

    x0, _ = harmonic_extension(tower.levels[1], boundary)
    q = _per_vertex(tower.levels[1], Vector(x0))

    out = Vector{Vector{Vector{Float64}}}()
    push!(out, q)

    for k in eachindex(tower.homs)
        fine = tower.levels[k + 1]
        fine_dims = vertex_stalks(fine)
        q_fine = Vector{Vector{Float64}}(undef, length(fine_dims))

        for v in eachindex(q)
            lifted = tower.bases[k][v] * q[v]
            off = 0
            # `fiber_vertices` returns the fibre in ascending vertex order, which is exactly
            # the row ordering `fiber_section_basis` documents for the basis it returns.
            for u in fiber_vertices(tower.homs[k], v)
                d = fine_dims[u]
                q_fine[u] = lifted[(off + 1):(off + d)]
                off += d
            end
        end

        push!(out, q_fine)
        q = q_fine
    end

    return out
end

"""
    solve_direct(tower::SheafTower, target_values) -> Vector{Vector{Float64}}

Baseline: harmonic extension on the fully expanded finest sheaf `tower.levels[end]`, pinning the
same targets that [`solve_hierarchical`](@ref) pins. Every agent moves independently, so teams
may deform — the internal formation edges are penalised, not enforced.

Returns per-vertex values on the finest level. This is the unconstrained optimum against which
the hierarchical solution is measured; see [`approximation_gap`](@ref).
"""
function solve_direct(tower::SheafTower, target_values::AbstractVector)
    boundary = _boundary_dict(tower, target_values)
    x, _ = harmonic_extension(tower.levels[end], boundary)
    return _per_vertex(tower.levels[end], Vector(x))
end

"""
    sheaf_energy(F::EuclideanSheaf, x::AbstractVector) -> Float64

Dirichlet energy `‖δx‖²` of a flat cochain `x` on `F`, where `δ` is the
[`coboundary_map`](@ref). This is the quantity both solves minimise, and the common yardstick
that makes them comparable.
"""
function sheaf_energy(F::EuclideanSheaf, x::AbstractVector)
    d = coboundary_map(F)
    return sum(abs2, d * x)
end

sheaf_energy(F::EuclideanSheaf, q::AbstractVector{<:AbstractVector}) =
    sheaf_energy(F, reduce(vcat, q))

"""
    approximation_gap(tower::SheafTower, target_values)
        -> (; hierarchical, direct, gap, relative_gap)

Energy of both solutions measured on the *finest* sheaf, and their difference.

`gap` is guaranteed nonnegative up to floating-point tolerance, and this is a theorem rather
than an empirical observation. The hierarchical solution satisfies every constraint the direct
problem imposes plus the extra requirement that each team lie exactly on its space of global
sections; it is therefore a feasible point of the direct problem, and a feasible point can never
beat the optimum. Equality holds exactly when the direct optimum was already fibrewise-exact —
rigid teams pinned to a single target are the equality case. A positive gap quantifies what
insisting on rigid formations costs: it is the energy the direct solve recovers by letting teams
deform.

`relative_gap` is `gap / direct`, or `0.0` when both energies vanish.
"""
function approximation_gap(tower::SheafTower, target_values::AbstractVector)
    F = tower.levels[end]
    q_h = solve_hierarchical(tower, target_values)[end]
    q_d = solve_direct(tower, target_values)

    E_h = sheaf_energy(F, q_h)
    E_d = sheaf_energy(F, q_d)
    gap = E_h - E_d
    scale = max(E_d, eps(Float64))
    relative_gap = E_d <= eps(Float64) && abs(gap) <= eps(Float64) ? 0.0 : gap / scale

    return (hierarchical = E_h, direct = E_d, gap = gap, relative_gap = relative_gap)
end

end # module NestedSystems
