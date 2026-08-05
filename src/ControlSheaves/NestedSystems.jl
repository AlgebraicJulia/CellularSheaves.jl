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
using ..AgentControllers: AbstractAgentDynamics
using ..Layered: _default_initial_position

export AbstractSystemNode, LeafTeam, RefinedSystem, TargetSpec, Observation, NestedSystemSpec,
       SheafTower, build_sheaf_tower, n_agents,
       solve_hierarchical, solve_direct, sheaf_energy, approximation_gap,
       RestrictionSpec, ProjectMember, Centroid, RawRestriction, project, centroid,
       materialize_restriction, SystemEdge,
       AgentBinding, SystemBinding, ResolvedAgent, resolve_dynamics, homogeneous_binding

# ---------------------------------------------------------------------------
# Restriction-map declarations (Issue 011)
# ---------------------------------------------------------------------------

"""
    RestrictionSpec

Declarative description of the restriction map on one end of an edge: what a subsystem
*presents* to whatever it is wired to.

**Declaration is symbolic; lowering is numeric.** A `RestrictionSpec` is declared against the
subsystem's **raw joint state** — the plain concatenation of its direct members' values, of
dimension `n_members × D`. That dimension follows from declared arities alone, with no rank or
nullspace computation anywhere, which is what lets validation stay numeric-free. Only at tower
assembly is the declared map composed with the fibre-section basis discovered by the pushforward:

```
R_coarse = R · B_v        (D × total) · (total × k)  =  D × k
```

`EuclideanSheaf` permits non-square restriction maps, so `R_coarse` needs no special handling.

A node's **direct members** are its children if it is a [`RefinedSystem`](@ref), or its raw
agents if it is a [`LeafTeam`](@ref). Refinement is opaque: a refined child counts as exactly one
member no matter how many agents it eventually expands to.

Concrete specs: [`project`](@ref), [`centroid`](@ref), and [`RawRestriction`](@ref).
"""
abstract type RestrictionSpec end

"""
    ProjectMember(member)

Backing type for [`project`](@ref). `member` is a direct-member index, or a `Symbol` naming a
child of a [`RefinedSystem`](@ref).
"""
struct ProjectMember <: RestrictionSpec
    member::Union{Int,Symbol}
end

"""
    Centroid()

Backing type for [`centroid`](@ref).
"""
struct Centroid <: RestrictionSpec end

"""
    RawRestriction(M::Matrix{Float64})

Escape hatch: present an arbitrary declared matrix. `M` must be `D × (n_members · D)`, checked
when it is materialized against a specific node.
"""
struct RawRestriction <: RestrictionSpec
    M::Matrix{Float64}
end

"""
    project(i::Int) -> RestrictionSpec
    project(name::Symbol) -> RestrictionSpec

Present direct member `i`'s own block unchanged (or the child named `name`). The materialized
matrix is the selection matrix `[0 … I_D … 0]`.

This is the default on every edge, and it is the one spec that lowers all the way to a single
finest-level vertex — selecting a member, then that member's first member, and so on — so a
`project`-wired tower places its edges on raw agents in `H_N` exactly as it did before per-edge
maps existed.
"""
project(i::Union{Int,Symbol}) = ProjectMember(i)

"""
    centroid() -> RestrictionSpec

Present the unweighted average of the subsystem's **direct** members, `(1/N) Σ`.

Unlike [`project`](@ref), a centroid is a functional of several members at once, so it has no
representative vertex in `H_N`; it exists only at the level of abstraction where the subsystem
is a single vertex. See [`build_sheaf_tower`](@ref) for what that implies.
"""
centroid() = Centroid()

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
    SystemEdge(src, dst; src_map=project(1), dst_map=project(1))

A consensus edge between two children of a [`RefinedSystem`](@ref), given as index pairs into
`children`. `src_map` and `dst_map` declare what each endpoint *presents* on this edge (see
[`RestrictionSpec`](@ref)); a subsystem may present something different on each of its edges.

The defaults reproduce the plain `(i, j)` behaviour exactly, which is why an `internal_edges`
list of bare tuples is still accepted.
"""
struct SystemEdge
    src::Int
    dst::Int
    src_map::RestrictionSpec
    dst_map::RestrictionSpec
end

SystemEdge(src::Int, dst::Int; src_map::RestrictionSpec=project(1), dst_map::RestrictionSpec=project(1)) =
    SystemEdge(src, dst, src_map, dst_map)

SystemEdge(e::Tuple{Int,Int}) = SystemEdge(e[1], e[2])
SystemEdge(e::SystemEdge) = e

"""
    RefinedSystem(name, children, internal_edges=SystemEdge[])

A subsystem whose vertices are themselves systems. `internal_edges` are consensus edges among
`children`: either [`SystemEdge`](@ref)s, or bare `(i, j)` index pairs, which are promoted to
`SystemEdge`s with the default `project(1)` maps on both ends.
"""
struct RefinedSystem <: AbstractSystemNode
    name::Symbol
    children::Vector{AbstractSystemNode}
    internal_edges::Vector{SystemEdge}

    function RefinedSystem(name::Symbol, children::Vector{AbstractSystemNode}, internal_edges::Vector{SystemEdge})
        @argcheck !isempty(children) "RefinedSystem $(name) must have at least one child"
        n = length(children)
        for e in internal_edges
            @argcheck 1 <= e.src <= n "internal_edges index $(e.src) out of range 1:$n for RefinedSystem $(name)"
            @argcheck 1 <= e.dst <= n "internal_edges index $(e.dst) out of range 1:$n for RefinedSystem $(name)"
            @argcheck e.src != e.dst "internal_edges cannot connect child $(e.src) to itself"
            # Materializing here would need `D`, which the spec supplies later; the member
            # designators are checked against arity at tower-assembly time instead.
        end
        new(name, children, internal_edges)
    end
end

RefinedSystem(name::Symbol, children::Vector{<:AbstractSystemNode},
              internal_edges::AbstractVector=SystemEdge[]) =
    RefinedSystem(name, Vector{AbstractSystemNode}(children), SystemEdge[SystemEdge(e) for e in internal_edges])

"""
    TargetSpec(name)

An uncontrolled agent supplying a boundary condition at solve time. Targets are top-level
vertices and are never owned by a team.
"""
struct TargetSpec
    name::Symbol
end

"""
    Observation(system_path, target_index; system_map=project(1))

Declares that the system at `system_path` (a path of child indices from the root) observes
target `target_index`. Arbitrary many-to-many incidence is allowed: one system may observe
several targets and one target may be observed by several systems.

`system_map` declares what the observing system presents on this edge (see
[`RestrictionSpec`](@ref)); the target end is always the identity, since a target is a single
`D`-dimensional vertex at every level. The default `project(1)` reproduces the pre-Issue-011
behaviour of wiring the system's first agent to the target.
"""
struct Observation
    system_path::Vector{Int}
    target_index::Int
    system_map::RestrictionSpec
end

Observation(system_path::AbstractVector{Int}, target_index::Int; system_map::RestrictionSpec=project(1)) =
    Observation(collect(Int, system_path), target_index, system_map)

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
# Materializing restriction declarations against a node
# ---------------------------------------------------------------------------

"""
    _direct_member_names(node) -> Union{Vector{Symbol},Nothing}

Names of `node`'s direct members, or `nothing` when they are unnamed (a [`LeafTeam`](@ref)'s raw
agents are addressed by index only).
"""
_direct_member_names(node::RefinedSystem) = [c.name for c in node.children]
_direct_member_names(::LeafTeam) = nothing

"""
    _n_direct_members(node) -> Int

Number of direct members: children for a [`RefinedSystem`](@ref), raw agents for a
[`LeafTeam`](@ref).
"""
_n_direct_members(node::RefinedSystem) = length(node.children)
_n_direct_members(node::LeafTeam) = node.n_agents

"""
    _member_index(node, member) -> Int

Resolve a direct-member designator to an index, checking it against `node`'s arity. A `Symbol`
resolves against child names and throws a message naming the available children if unmatched.
"""
function _member_index(node::AbstractSystemNode, member::Int)
    n = _n_direct_members(node)
    @argcheck 1 <= member <= n "member index $member out of range 1:$n for system $(node.name)"
    return member
end

function _member_index(node::AbstractSystemNode, member::Symbol)
    names = _direct_member_names(node)
    names === nothing && throw(ArgumentError(
        "cannot select member :$member in $(node.name): a LeafTeam's agents are unnamed, " *
        "use an integer index in 1:$(_n_direct_members(node))"))
    idx = findfirst(==(member), names)
    idx === nothing && throw(ArgumentError(
        "system $(node.name) has no child named :$member (children: $(join(names, ", ")))"))
    return idx
end

"""
    materialize_restriction(r::RestrictionSpec, node::AbstractSystemNode, D::Int) -> Matrix{Float64}

Build the `D × (n_members · D)` matrix for `r` against `node`'s raw joint state.

Purely structural: it depends only on `node`'s declared arity, never on a rank, nullspace, or
pseudo-inverse computation. Composition with the fibre basis — the only numeric step — happens
later, at tower assembly.
"""
function materialize_restriction(r::ProjectMember, node::AbstractSystemNode, D::Int)
    @argcheck D >= 1 "D must be positive (got $D)"
    n = _n_direct_members(node)
    i = _member_index(node, r.member)
    R = zeros(Float64, D, n * D)
    R[:, ((i - 1) * D + 1):(i * D)] = Matrix{Float64}(I, D, D)
    return R
end

function materialize_restriction(::Centroid, node::AbstractSystemNode, D::Int)
    @argcheck D >= 1 "D must be positive (got $D)"
    n = _n_direct_members(node)
    R = zeros(Float64, D, n * D)
    for i in 1:n
        R[:, ((i - 1) * D + 1):(i * D)] = Matrix{Float64}(I, D, D) ./ n
    end
    return R
end

function materialize_restriction(r::RawRestriction, node::AbstractSystemNode, D::Int)
    @argcheck D >= 1 "D must be positive (got $D)"
    total = _n_direct_members(node) * D
    @argcheck size(r.M) == (D, total) "RawRestriction matrix is $(size(r.M)), expected ($D, $total) for system $(node.name)"
    return copy(r.M)
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
    n_agents(spec::NestedSystemSpec) -> Int
    n_agents(node::AbstractSystemNode) -> Int

Total number of raw agents in `spec`'s tree (or in `node`'s subtree).
"""
n_agents(spec::NestedSystemSpec) = n_agents(spec.root)
n_agents(node::LeafTeam) = node.n_agents
n_agents(node::RefinedSystem) = sum(n_agents, node.children)

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
    _fine_representative(node, r, agent_range) -> Union{Int,Nothing}

The single `H_N` vertex that `r` designates on `node`, or `nothing` when `r` designates no
single vertex.

Only [`project`](@ref) has a representative, and finding it is recursive: selecting a member of a
refined system leaves another system, whose own representative is *its* first member, and so on
down to a raw agent. Selecting a member of a [`LeafTeam`](@ref) lands on that agent directly.

An aggregate map such as [`centroid`](@ref) returns `nothing`: it is a functional of several
vertices at once and cannot be an edge endpoint in `H_N`, where those vertices are still
separate. Such an edge is placed at the coarsest level where its endpoints are single vertices —
see [`build_sheaf_tower`](@ref).
"""
function _fine_representative(node::LeafTeam, r::ProjectMember,
                               agent_range::IdDict{AbstractSystemNode,UnitRange{Int}})
    return first(agent_range[node]) + _member_index(node, r.member) - 1
end

function _fine_representative(node::RefinedSystem, r::ProjectMember,
                               agent_range::IdDict{AbstractSystemNode,UnitRange{Int}})
    child = node.children[_member_index(node, r.member)]
    return _fine_representative(child, ProjectMember(1), agent_range)
end

_fine_representative(::AbstractSystemNode, ::RestrictionSpec,
                     ::IdDict{AbstractSystemNode,UnitRange{Int}}) = nothing

"""
    _DeferredEdge

An edge that could not be expressed in `H_N` because at least one endpoint presents an aggregate
of several vertices. It is added at `level`, the finest level at which both endpoints are single
vertices. `u_node`/`v_node` are `nothing` for a target endpoint (targets are single vertices
everywhere and always present the identity).
"""
struct _DeferredEdge
    level::Int
    u_node::Union{AbstractSystemNode,Nothing}
    u_map::Union{RestrictionSpec,Nothing}
    u_target::Int
    v_node::Union{AbstractSystemNode,Nothing}
    v_map::Union{RestrictionSpec,Nothing}
    v_target::Int
end

"""
    _atomic_level(node, S, depth) -> Int

The finest tower level at which `node` is a single vertex. Below it, `node` has been split into
its direct members — which is precisely the level whose fibre basis an aggregate restriction map
must compose with.
"""
_atomic_level(node::AbstractSystemNode, S::IdDict{AbstractSystemNode,Int}, depth::Int) = depth - S[node]

"""
    _wire_declared_edges!(F, spec, agent_range, S, depth) -> Vector{_DeferredEdge}

Add every declared consensus and observation edge that can live in `H_N`, and return the ones
that cannot.

An edge lands in `H_N` exactly when **both** endpoints resolve to a single vertex there — which
is the case for `project`-wired edges, and so for every spec written before per-edge maps
existed. Keeping those edges at the finest level is what preserves the earlier behaviour and, more
importantly, what keeps [`solve_direct`](@ref)'s baseline meaningful: a target pinned through an
`H_N` edge still constrains individual agents when the tower's rigidity is relaxed.
"""
function _wire_declared_edges!(F::EuclideanSheaf, spec::NestedSystemSpec,
                                agent_range::IdDict{AbstractSystemNode,UnitRange{Int}},
                                S::IdDict{AbstractSystemNode,Int}, depth::Int)
    D = spec.D
    I_D = Matrix{Float64}(I, D, D)
    deferred = _DeferredEdge[]

    _wire_internal_edges!(F, spec.root, agent_range, S, depth, D, I_D, deferred)

    for obs in spec.observations
        node = _resolve_path(spec.root, obs.system_path)
        rep = _fine_representative(node, obs.system_map, agent_range)
        if rep === nothing
            push!(deferred, _DeferredEdge(_atomic_level(node, S, depth),
                                          node, obs.system_map, 0,
                                          nothing, nothing, obs.target_index))
        else
            add_sheaf_edge!(F, rep, obs.target_index, I_D, I_D)
        end
    end

    return deferred
end

function _wire_internal_edges!(F::EuclideanSheaf, node::RefinedSystem,
                                agent_range::IdDict{AbstractSystemNode,UnitRange{Int}},
                                S::IdDict{AbstractSystemNode,Int}, depth::Int, D::Int,
                                I_D::Matrix{Float64}, deferred::Vector{_DeferredEdge})
    for e in node.internal_edges
        cu, cv = node.children[e.src], node.children[e.dst]
        ru = _fine_representative(cu, e.src_map, agent_range)
        rv = _fine_representative(cv, e.dst_map, agent_range)
        if ru === nothing || rv === nothing
            # One aggregate endpoint forces *both* ends to be expressed at the coarse level,
            # since an edge cannot straddle two levels.
            lvl = min(_atomic_level(cu, S, depth), _atomic_level(cv, S, depth))
            push!(deferred, _DeferredEdge(lvl, cu, e.src_map, 0, cv, e.dst_map, 0))
        else
            add_sheaf_edge!(F, ru, rv, I_D, I_D)
        end
    end
    for c in node.children
        c isa RefinedSystem && _wire_internal_edges!(F, c, agent_range, S, depth, D, I_D, deferred)
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
    deferred = _wire_declared_edges!(F, spec, agent_range, S, depth)

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

    # Restriction matrices for deferred (aggregate-mapped) edge endpoints, filled in as the
    # loop below reaches each endpoint's own atomic level — see `_DeferredEdge`.
    R_cache = IdDict{AbstractSystemNode,Matrix{Float64}}()
    I_D = Matrix{Float64}(I, D, D)

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

        # --- Materialize any deferred (aggregate) endpoint whose node becomes atomic here. ---
        # Row ordering: `node`'s direct members occupy a contiguous block of `level_slots[k]`
        # in their own declared order (`_level_owner_slots` visits `node.children` in order and
        # `vmap` numbers vertices by list position), so `fiber_vertices(hom, cv)` — the row
        # order `bases[k][cv]` was built against — already matches `materialize_restriction`'s
        # column order. The dimension check below is a cheap guard against that invariant ever
        # breaking; a silent mismatch would misassign restriction blocks rather than error.
        for d in deferred
            for (node, spec_map) in ((d.u_node, d.u_map), (d.v_node, d.v_map))
                node === nothing && continue
                haskey(R_cache, node) && continue
                _atomic_level(node, S, depth) == k || continue
                cv = n_targets + level_pos[k][node]   # bases[k] is indexed by coarse vertex, targets first
                B = bases[k][cv]
                @argcheck size(B, 1) == _n_direct_members(node) * D (
                    "fibre basis at coarse vertex $cv (level $k, system $(node.name)) has " *
                    "$(size(B,1)) rows, expected $(_n_direct_members(node) * D) — row ordering " *
                    "does not match materialize_restriction's column ordering")
                R_cache[node] = materialize_restriction(spec_map, node, D) * B
            end
        end
        for d in deferred
            d.level == k || continue
            R_u = d.u_node === nothing ? I_D : R_cache[d.u_node]
            R_v = d.v_node === nothing ? I_D : R_cache[d.v_node]
            u_idx = d.u_node === nothing ? d.u_target : n_targets + level_pos[k][d.u_node]
            v_idx = d.v_node === nothing ? d.v_target : n_targets + level_pos[k][d.v_node]
            add_sheaf_edge!(levels[k], u_idx, v_idx, R_u, R_v)
        end

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

# ---------------------------------------------------------------------------
# Dynamics binding cascade (Issue 012)
# ---------------------------------------------------------------------------

"""
    AgentBinding(; dynamics=nothing, K_lqr=nothing, initial_position=nothing)

Per-agent dynamics override. Any field left `nothing` falls back to the enclosing team's
binding, then to successively more distant ancestors — see [`resolve_dynamics`](@ref). Fields
resolve independently: an agent may override only its initial position while inheriting its
dynamics model from an ancestor.
"""
Base.@kwdef struct AgentBinding
    dynamics::Union{Nothing,AbstractAgentDynamics} = nothing
    K_lqr::Union{Nothing,Matrix{Float64}} = nothing
    initial_position::Union{Nothing,Vector{Float64}} = nothing
end

"""
    SystemBinding(; dynamics=nothing, K_lqr=nothing, children=Dict(), agents=Dict())

Dynamics bindings for one node of the system tree, mirroring the structure of the
[`NestedSystemSpec`](@ref) it is resolved against.

`dynamics`/`K_lqr` apply to every leaf agent in this subtree unless a descendant overrides them.
`children` maps a child system's name to its own `SystemBinding`. `agents` maps a local agent
index within a [`LeafTeam`](@ref) to an [`AgentBinding`](@ref).

Child names and agent indices are validated against the spec during [`resolve_dynamics`](@ref):
a name or index present here but absent from the spec is an error, so a typo surfaces
immediately rather than as a silently-unbound agent.
"""
Base.@kwdef struct SystemBinding
    dynamics::Union{Nothing,AbstractAgentDynamics} = nothing
    K_lqr::Union{Nothing,Matrix{Float64}} = nothing
    children::Dict{Symbol,SystemBinding} = Dict{Symbol,SystemBinding}()
    agents::Dict{Int,AgentBinding} = Dict{Int,AgentBinding}()
end

"""
    homogeneous_binding(dyn; K_lqr=zeros(0,0)) -> SystemBinding

Shorthand for a root-level binding applying `dyn`/`K_lqr` to every agent in the tree — the
common case, and the degenerate case that reproduces the old flat `HomogeneousDynamics`.
"""
homogeneous_binding(dyn::AbstractAgentDynamics; K_lqr::Matrix{Float64}=zeros(0, 0)) =
    SystemBinding(dynamics=dyn, K_lqr=K_lqr)

"""
    ResolvedAgent

A fully-bound leaf agent: its global agent index (matching [`SheafTower`](@ref)'s
`agent_vertices`), its dotted path in the tree (for error messages and debugging), and the
dynamics, gain, and initial position that won the inheritance cascade.
"""
struct ResolvedAgent
    agent_index::Int
    path::Vector{Symbol}
    dynamics::AbstractAgentDynamics
    K_lqr::Matrix{Float64}
    initial_position::Vector{Float64}
end

"""
    _fold(inherited::AgentBinding, local_::AgentBinding) -> AgentBinding

Overlay `local_` on `inherited`, field by field: a non-`nothing` local field replaces the
inherited one, a `nothing` field leaves the inherited value intact. This independence — not
replacing the whole binding wholesale — is what lets an agent override just its initial
position while still inheriting its dynamics model.
"""
_fold(inherited::AgentBinding, local_::AgentBinding) = AgentBinding(
    dynamics = something(local_.dynamics, inherited.dynamics, Some(nothing)),
    K_lqr = something(local_.K_lqr, inherited.K_lqr, Some(nothing)),
    initial_position = something(local_.initial_position, inherited.initial_position, Some(nothing)),
)

_fold_system(inherited::AgentBinding, sb::SystemBinding) = _fold(inherited,
    AgentBinding(dynamics=sb.dynamics, K_lqr=sb.K_lqr))

"""
    resolve_dynamics(spec::NestedSystemSpec, ctx::SystemBinding) -> Vector{ResolvedAgent}

Walk `spec`'s tree carrying the inherited binding, overriding it wherever `ctx` supplies a more
specific value, and return one [`ResolvedAgent`](@ref) per leaf agent in global agent-index
order (matching [`SheafTower`](@ref)'s `agent_vertices`).

Precedence, most specific first: per-agent, leaf team, nearer ancestor, root default. Each field
resolves independently.

Throws if any agent ends with no `dynamics` bound anywhere up its chain, or if `ctx` names a
child or agent index that does not exist in `spec` — both errors report the full dotted path.
"""
function resolve_dynamics(spec::NestedSystemSpec, ctx::SystemBinding)
    resolved = ResolvedAgent[]
    next_idx = Ref(1)
    _resolve_dynamics!(resolved, next_idx, spec.root, ctx, AgentBinding(), Symbol[])
    return resolved
end

function _resolve_dynamics!(resolved::Vector{ResolvedAgent}, next_idx::Ref{Int},
                             node::RefinedSystem, ctx::SystemBinding,
                             inherited::AgentBinding, path::Vector{Symbol})
    !isempty(ctx.agents) && throw(ArgumentError(
        "SystemBinding declares `agents` at $(join([path; node.name], ".")), but $(node.name) " *
        "is a RefinedSystem — per-agent overrides belong on the LeafTeam that owns them"))
    for name in keys(ctx.children)
        any(c -> c.name == name, node.children) || throw(ArgumentError(
            "SystemBinding names child :$name at $(join([path; node.name], ".")), but " *
            "$(node.name) has no such child (children: $(join([c.name for c in node.children], ", ")))"))
    end

    here = _fold_system(inherited, ctx)
    for c in node.children
        child_ctx = get(ctx.children, c.name, SystemBinding())
        _resolve_dynamics!(resolved, next_idx, c, child_ctx, here, push!(copy(path), c.name))
    end
end

function _resolve_dynamics!(resolved::Vector{ResolvedAgent}, next_idx::Ref{Int},
                             node::LeafTeam, ctx::SystemBinding,
                             inherited::AgentBinding, path::Vector{Symbol})
    !isempty(ctx.children) && throw(ArgumentError(
        "SystemBinding declares children at $(join(path, ".")), but $(node.name) is a LeafTeam " *
        "with no children — use `agents` to bind individual agents by local index"))
    for idx in keys(ctx.agents)
        @argcheck 1 <= idx <= node.n_agents "SystemBinding names agent $idx at $(join(path, ".")), out of range 1:$(node.n_agents)"
    end

    here = _fold_system(inherited, ctx)
    for i in 1:node.n_agents
        local_binding = get(ctx.agents, i, AgentBinding())
        b = _fold(here, local_binding)
        b.dynamics === nothing && throw(ArgumentError(
            "no dynamics bound for agent $i at $(join([path; Symbol(i)], ".")) — declare a " *
            "`dynamics` at this agent, this team, or an ancestor"))
        agent_index = next_idx[]
        next_idx[] += 1
        K_lqr = something(b.K_lqr, Some(zeros(0, 0)))
        initial_position = something(b.initial_position, Some(_default_initial_position(b.dynamics, agent_index)))
        push!(resolved, ResolvedAgent(agent_index, [path; Symbol(i)], b.dynamics, K_lqr, initial_position))
    end
end

end # module NestedSystems
