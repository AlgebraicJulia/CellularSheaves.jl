"""
    NestedDSLTerm

AST node types for the nested-system DSL, plus the [`SystemFragment`](@ref) value that
collects them and the error type the later stages throw.

A fragment is an ordered list of **terms**, each describing one declaration made *at one node*
of the system tree. Fragments are immutable values: [`merge`](@ref) concatenates them, and
a [`SubsystemTerm`](@ref) nests one inside another. Nothing here is a `NestedSystemSpec` yet —
names are still names and no index has been assigned — which is what lets a fragment be built
in pieces, by ordinary Julia code, and assembled later.

Every term stores **already-evaluated Julia values** (a `Float64` radius, a
[`RestrictionSpec`](@ref), an `AbstractAgentDynamics`), not symbols to be resolved against a
context dict later. The DSL's macros expand to calls to the constructors in this module, so by
the time a term exists its numeric content has been computed by ordinary Julia evaluation in the
caller's own scope. Only *structural* names — child names, target names, member designators —
remain symbolic, and those are resolved by `NestedDSLLowering`.
"""
module NestedDSLTerm

using ArgCheck
using ...NestedSystems: RestrictionSpec, ProjectMember, project,
                        LeafTeam, RefinedSystem, SystemEdge, TargetSpec, Observation,
                        SystemBinding, AgentBinding
using ...AgentControllers: AbstractAgentDynamics

export NestedTerm, SystemFragment, MemberRef,
       TeamTerm, SubsystemTerm, LinkTerm, TargetTerm, ObserveTerm, BindTerm,
       DimTerm, AffineTerm,
       team_term, subsystem_term, link_term, target_term, observe_term, bind_term,
       dim_term, affine_term, member_ref, terms, NestedDSLError

"""
    NestedDSLError(msg)

Raised by every stage of the nested DSL — parsing, validation, and lowering — so that a caller
can distinguish a malformed specification from an unrelated `ArgumentError` thrown deep inside
the tower compiler.
"""
struct NestedDSLError <: Exception
    msg::String
end

Base.showerror(io::IO, e::NestedDSLError) = print(io, "NestedDSLError: ", e.msg)

"""
    NestedTerm

Supertype of every declaration a [`SystemFragment`](@ref) can hold.
"""
abstract type NestedTerm end

"""
    MemberRef(path, map=project(1))

A reference to a system in the tree, together with the [`RestrictionSpec`](@ref) it *presents*
on the edge being declared — the DSL's spelling of "this subsystem, seen through this map".

`path` is a dotted chain of child names relative to the node where the referring term was
declared: `[:mid, :ringA]` for `mid.ringA`. An empty `path` is not meaningful and is rejected
by the validator.

Endpoint syntax and the `MemberRef` it builds:

| Surface syntax | `path` | `map` |
|---|---|---|
| `ringA` | `[:ringA]` | `project(1)` |
| `mid.ringA` | `[:mid, :ringA]` | `project(1)` |
| `centroid(ringA)` | `[:ringA]` | `centroid()` |
| `project(hub, 3)` | `[:hub]` | `project(3)` |
| `project(mid, :ringA)` | `[:mid]` | `project(:ringA)` |
| `raw(ringA, M)` | `[:ringA]` | `RawRestriction(M)` |
| `via(ringA, s)` | `[:ringA]` | `s` |
"""
struct MemberRef
    path::Vector{Symbol}
    map::RestrictionSpec
end

MemberRef(path::AbstractVector{Symbol}) = MemberRef(collect(Symbol, path), project(1))

"""
    member_ref(path, map=project(1)) -> MemberRef

Functional constructor for [`MemberRef`](@ref); the DSL's endpoint syntax expands to a call to
this. `path` may be a single `Symbol` or a vector of them.
"""
member_ref(path::AbstractVector, map::RestrictionSpec=project(1)) =
    MemberRef(Symbol[Symbol(p) for p in path], map)
member_ref(name::Symbol, map::RestrictionSpec=project(1)) = MemberRef(Symbol[name], map)

"""
    TeamTerm(name, kind, n_agents, radius, observers)

Declares a leaf team of raw agents — the DSL's `@team` — lowering to a
[`LeafTeam`](@ref). See [`team_term`](@ref).
"""
struct TeamTerm <: NestedTerm
    name::Symbol
    kind::Symbol
    n_agents::Int
    radius::Float64
    observers::Vector{Int}
end

"""
    team_term(name, kind, n_agents, radius=1.0; observers=[1]) -> TeamTerm

Build a [`TeamTerm`](@ref), validating `kind` and `n_agents` up front so that a typo surfaces at
the declaration rather than at tower assembly.

`@team ringA = ring(4; radius=1.0)` expands to a call to this function.
"""
function team_term(name::Symbol, kind::Symbol, n_agents::Integer, radius::Real=1.0;
                   observers=[1])
    kind in (:ring, :path, :star, :clique) ||
        throw(NestedDSLError("team $name: kind must be one of ring, path, star, clique (got $kind)"))
    n_agents >= 1 || throw(NestedDSLError("team $name: n_agents must be positive (got $n_agents)"))
    obs = collect(Int, observers)
    all(1 <= o <= n_agents for o in obs) ||
        throw(NestedDSLError("team $name: observers $obs must lie within 1:$n_agents"))
    return TeamTerm(name, kind, Int(n_agents), Float64(radius), obs)
end

"""
    SubsystemTerm(name, body)

Declares a refined child system whose contents are the fragment `body` — the DSL's `@system` —
lowering to a [`RefinedSystem`](@ref).

Nesting is what makes the DSL's depth unbounded: `body` may itself contain `SubsystemTerm`s, to
any depth. Targets and bindings declared inside `body` are *not* rewritten here; lowering walks
the tree carrying the current path and prefixes them then, which is what lets the same fragment
be spliced in at different depths without change.
"""
struct SubsystemTerm <: NestedTerm
    name::Symbol
    body::Any   # SystemFragment; declared Any because SystemFragment is defined below
end

"""
    LinkTerm(src, dst)

A consensus edge between two **direct children** of the declaring node — the DSL's `@link` —
lowering to a [`SystemEdge`](@ref). Both [`MemberRef`](@ref)s must have single-component paths,
which the validator enforces: `SystemEdge` connects siblings only.
"""
struct LinkTerm <: NestedTerm
    src::MemberRef
    dst::MemberRef
end

"""
    TargetTerm(name)

Declares a target — the DSL's `@target` — lowering to a [`TargetSpec`](@ref). Targets are global
to the whole specification no matter how deeply the declaring fragment is nested, so that a
fragment can carry the target it observes along with it.
"""
struct TargetTerm <: NestedTerm
    name::Symbol
end

"""
    ObserveTerm(system, target)

Declares that `system` observes the target named `target` — the DSL's `@observe` — lowering to
an [`Observation`](@ref). `system.path` is relative to the declaring node and may be of any
depth; `system.map` is what the observing system presents on the edge.
"""
struct ObserveTerm <: NestedTerm
    system::MemberRef
    target::Symbol
end

"""
    BindTerm(path, agent, dynamics, K_lqr, initial_position)

Declares dynamics/gain/initial-position bindings — the DSL's `@bind` — lowering into the
[`SystemBinding`](@ref) cascade.

An empty `path` binds the declaring node itself (and hence every agent beneath it, unless
overridden); a non-empty `path` binds a descendant. `agent` is a local agent index within a
leaf team, or `nothing` to bind the node as a whole. Any field left `nothing` is simply not
declared here and is inherited, exactly as [`AgentBinding`](@ref) specifies.
"""
struct BindTerm <: NestedTerm
    path::Vector{Symbol}
    agent::Union{Nothing,Int}
    dynamics::Union{Nothing,AbstractAgentDynamics}
    K_lqr::Union{Nothing,Matrix{Float64}}
    initial_position::Union{Nothing,Vector{Float64}}
end

"""
    bind_term(path, agent; dynamics=nothing, K_lqr=nothing, initial_position=nothing) -> BindTerm

Functional constructor for [`BindTerm`](@ref). Rejects a binding that declares nothing at all,
since that is always a mistake rather than a no-op the user meant to write.
"""
function bind_term(path::AbstractVector, agent::Union{Nothing,Integer}=nothing;
                   dynamics=nothing, K_lqr=nothing, initial_position=nothing)
    dynamics === nothing && K_lqr === nothing && initial_position === nothing &&
        throw(NestedDSLError("@bind at $(_show_path(path)) declares nothing — give at least one " *
                             "of dynamics, K_lqr, initial_position"))
    return BindTerm(Symbol[Symbol(p) for p in path],
                    agent === nothing ? nothing : Int(agent),
                    dynamics,
                    K_lqr === nothing ? nothing : Matrix{Float64}(K_lqr),
                    initial_position === nothing ? nothing : Vector{Float64}(initial_position))
end

"""
    DimTerm(D)

Declares the stalk dimension — the DSL's `@dim`. Bubbles to the root during lowering no matter
where it is declared; declaring two *different* values anywhere in one specification is an error.
"""
struct DimTerm <: NestedTerm
    D::Int
end

"""
    AffineTerm(affine)

Declares whether restriction maps are affine — the DSL's `@affine`. Bubbles and conflicts exactly
like [`DimTerm`](@ref).
"""
struct AffineTerm <: NestedTerm
    affine::Bool
end

"""
    SystemFragment(terms=NestedTerm[])

An ordered, immutable bag of declarations describing the contents of **one node** of a system
tree — the value `@nested_system` returns and the unit of composition in this DSL.

A fragment is not tied to a position in any tree: the paths inside its terms are relative to
whatever node it ends up at. That is precisely what makes fragments composable —

* `merge(f, g)` concatenates two fragments describing the same node;
* `@include g` splices `g`'s declarations into the fragment being built;
* `@system name = g` makes `g` the body of a new child.

so a Julia function may build a fragment from ordinary arguments, a `for` loop may build one per
iteration, and the results merged — with no looping or conditional constructs in the DSL itself.
"""
struct SystemFragment
    terms::Vector{NestedTerm}
end

SystemFragment() = SystemFragment(NestedTerm[])

"""
    terms(f::SystemFragment) -> Vector{NestedTerm}

The fragment's declarations, in declaration order. Returns the internal vector; treat it as
read-only.
"""
terms(f::SystemFragment) = f.terms

"""
    merge(fragments::SystemFragment...) -> SystemFragment

Concatenate fragments describing the **same** node, preserving declaration order. This is the
DSL's composition operator: build pieces independently — in a loop, in a helper function, behind
an `if` — and merge them into the fragment for one node.

Merging does not check for duplicate names; that is the validator's job, and deferring it means
a partial fragment mid-composition is never spuriously invalid.
"""
Base.merge(f::SystemFragment, g::SystemFragment) = SystemFragment(vcat(f.terms, g.terms))
Base.merge(f::SystemFragment, g::SystemFragment, rest::SystemFragment...) =
    merge(merge(f, g), rest...)
Base.merge(f::SystemFragment) = f

Base.isempty(f::SystemFragment) = isempty(f.terms)
Base.length(f::SystemFragment) = length(f.terms)

function Base.show(io::IO, f::SystemFragment)
    n_team = count(t -> t isa TeamTerm, f.terms)
    n_sub = count(t -> t isa SubsystemTerm, f.terms)
    n_link = count(t -> t isa LinkTerm, f.terms)
    n_tgt = count(t -> t isa TargetTerm, f.terms)
    n_obs = count(t -> t isa ObserveTerm, f.terms)
    print(io, "SystemFragment($(length(f.terms)) terms: $n_team team, $n_sub system, ",
          "$n_link link, $n_tgt target, $n_obs observe)")
end

# Functional constructors the parser's macros expand into. Keeping the macro expansions as plain
# calls (rather than inlined `Expr`s building structs directly) means the same surface is usable
# from ordinary Julia without the macros at all.

"""
    subsystem_term(name, body::SystemFragment) -> SubsystemTerm

Functional constructor for [`SubsystemTerm`](@ref).
"""
function subsystem_term(name::Symbol, body::SystemFragment)
    return SubsystemTerm(name, body)
end

"""
    link_term(src::MemberRef, dst::MemberRef) -> LinkTerm

Functional constructor for [`LinkTerm`](@ref).
"""
link_term(src::MemberRef, dst::MemberRef) = LinkTerm(src, dst)

"""
    target_term(name::Symbol) -> TargetTerm

Functional constructor for [`TargetTerm`](@ref).
"""
target_term(name::Symbol) = TargetTerm(name)

"""
    observe_term(system::MemberRef, target::Symbol) -> ObserveTerm

Functional constructor for [`ObserveTerm`](@ref).
"""
observe_term(system::MemberRef, target::Symbol) = ObserveTerm(system, target)

"""
    dim_term(D::Integer) -> DimTerm

Functional constructor for [`DimTerm`](@ref).
"""
function dim_term(D::Integer)
    D >= 1 || throw(NestedDSLError("@dim must be positive (got $D)"))
    return DimTerm(Int(D))
end

"""
    affine_term(affine::Bool) -> AffineTerm

Functional constructor for [`AffineTerm`](@ref).
"""
affine_term(affine::Bool) = AffineTerm(affine)

_show_path(path) = isempty(path) ? "<this node>" : join(path, ".")

end # module NestedDSLTerm
