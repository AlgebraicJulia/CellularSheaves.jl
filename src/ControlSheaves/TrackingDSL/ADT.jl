"""
    TrackingDSLTerm

Abstract data type (AST) definitions for the name-driven Tracking DSL.

The DSL allows symbolic specification of multi-agent tracking problems with
late binding of numeric values via an external context dict. A `TrackingProgram`
is the root node holding a flat list of `TrackingStmt` nodes. Statements cover:

- Space and map declarations (`space X = R^6`, `map A : X -> X`)
- Agent/target declarations with optional dynamics (`agent a1 dynamics (A, B) period dt`)
- Time declarations (`horizon K`, `times Tall = 0:K`)
- Constraints (`consensus`, `track`, `consensus_sheaf`)
- Boundary conditions (`boundary agent a1 at t[begin] = x0`)
"""
module TrackingDSLTerm

using MLStyle: @data, @match

export TrackingProgram, TrackingStmt
export SpaceDecl, MapDecl, VectorDecl, AgentDecl, TargetDecl
export HorizonDecl, TimeAliasDecl, TimeSetDecl
export ConsensusConstraint, TrackConstraint, ConsensusSheafConstraint
export BoundaryConstraint
export TimeSpec, SingletonTime, TimeList, TimeRange, NamedTimeSet
export TimeRef, LiteralTime, NamedTime, BeginTime, EndTime
export SpaceExpr, BaseSpace, ProductSpace, DirectSumSpace

# ── Error hierarchy ─────────────────────────────────────────────────────────

"""
    TrackingDSLError

Supertype for all DSL-specific errors. Every subtype carries a human-readable
`msg::String` and, where applicable, a `context` field identifying the
offending statement or symbol.
"""
abstract type TrackingDSLError <: Exception end

"""
    TrackingSyntaxError

Raised when a statement cannot be parsed.
"""
struct TrackingSyntaxError <: TrackingDSLError
    msg::String
end
Base.showerror(io::IO, e::TrackingSyntaxError) =
    print(io, "TrackingDSL Syntax Error: ", e.msg)

"""
    TrackingDeclarationError

Raised for duplicate or missing name declarations.
"""
struct TrackingDeclarationError <: TrackingDSLError
    msg::String
end
Base.showerror(io::IO, e::TrackingDeclarationError) =
    print(io, "TrackingDSL Declaration Error: ", e.msg)

"""
    TrackingTypeError

Raised when a name is used with an incompatible type.
"""
struct TrackingTypeError <: TrackingDSLError
    msg::String
end
Base.showerror(io::IO, e::TrackingTypeError) =
    print(io, "TrackingDSL Type Error: ", e.msg)

"""
    TrackingUnboundSymbolError

Raised when lowering encounters a symbol without a bound value.
"""
struct TrackingUnboundSymbolError <: TrackingDSLError
    msg::String
end
Base.showerror(io::IO, e::TrackingUnboundSymbolError) =
    print(io, "TrackingDSL Unbound Symbol Error: ", e.msg)

"""
    TrackingTimeResolutionError

Raised when a time reference cannot be resolved.
"""
struct TrackingTimeResolutionError <: TrackingDSLError
    msg::String
end
Base.showerror(io::IO, e::TrackingTimeResolutionError) =
    print(io, "TrackingDSL Time Resolution Error: ", e.msg)

"""
    TrackingDimensionMismatchError

Raised when matrix/vector dimensions are inconsistent.
"""
struct TrackingDimensionMismatchError <: TrackingDSLError
    msg::String
end
Base.showerror(io::IO, e::TrackingDimensionMismatchError) =
    print(io, "TrackingDSL Dimension Mismatch Error: ", e.msg)

"""
    TrackingTemplateError

Raised when a consensus sheaf template cannot be attached.
"""
struct TrackingTemplateError <: TrackingDSLError
    msg::String
end
Base.showerror(io::IO, e::TrackingTemplateError) =
    print(io, "TrackingDSL Template Error: ", e.msg)

export TrackingDSLError, TrackingSyntaxError, TrackingDeclarationError,
    TrackingTypeError, TrackingUnboundSymbolError, TrackingTimeResolutionError,
    TrackingDimensionMismatchError, TrackingTemplateError

# ── Time references ──────────────────────────────────────────────────────────

"""
    TimeRef

A reference to a single point in time. Can be:
- `LiteralTime(n)` — an integer literal (e.g. `0`, `7`)
- `NamedTime(name)` — a named time alias (e.g. `tcapture`)
- `BeginTime()` — the special alias resolving to `0` (written `t[begin]` in the DSL)
- `EndTime()` — the special alias resolving to `k` from `horizon` (written `t[end]` in the DSL)
"""
abstract type TimeRef end
struct LiteralTime <: TimeRef; value::Int; end
struct NamedTime   <: TimeRef; name::Symbol; end
struct BeginTime   <: TimeRef; end
struct EndTime     <: TimeRef; end

"""
    TimeSpec

Describes the set of timesteps at which a constraint is active. Variants:

- `SingletonTime(t)` — activate at exactly one `TimeRef`
- `TimeList(ts)` — activate at an explicit list of `TimeRef`s
- `TimeRange(lo, hi)` — a range `lo:hi` of `TimeRef`s
- `NamedTimeSet(name)` — reference to a declared `times` alias
"""
abstract type TimeSpec end
struct SingletonTime <: TimeSpec; t::TimeRef; end
struct TimeList      <: TimeSpec; ts::Vector{TimeRef}; end
struct TimeRange     <: TimeSpec; lo::TimeRef; hi::TimeRef; end
struct NamedTimeSet  <: TimeSpec; name::Symbol; end

# ── Space expressions ─────────────────────────────────────────────────────────

"""
    SpaceExpr

A space expression describing the type of a stalk. Variants:

- `BaseSpace(name)` — a named scalar space (e.g. `X` declared as `space X = R^6`)
- `ProductSpace(a, b)` — Cartesian product `a × b`
- `DirectSumSpace(a, b)` — direct sum `a ⊕ b`
"""
abstract type SpaceExpr end
struct BaseSpace     <: SpaceExpr; name::Symbol; end
struct ProductSpace  <: SpaceExpr; a::SpaceExpr; b::SpaceExpr; end
struct DirectSumSpace <: SpaceExpr; a::SpaceExpr; b::SpaceExpr; end

# ── Boundary value reference ──────────────────────────────────────────────────

"""
    BoundaryRef

A reference to a boundary value in the context dict. Either a plain name `foo`
(unindexed) or an indexed form `x[a, t]` allowing late-bound agent/time indices.
"""
abstract type BoundaryRef end
struct PlainRef   <: BoundaryRef; name::Symbol; end
struct IndexedRef <: BoundaryRef
    name::Symbol
    agent_index::Symbol   # name of the agent index variable
    time_index::Symbol    # name of the time index variable
end

export BoundaryRef, PlainRef, IndexedRef

# ── Statement nodes ──────────────────────────────────────────────────────────

"""
    TrackingStmt

Supertype for all top-level statements in a `TrackingProgram`. Use `@data`
below for the concrete variants.
"""
abstract type TrackingStmt end

# --- Declarations ---

"""
    SpaceDecl(name, dim)

Declares a named Euclidean space, e.g. `space X = R^6`.
`dim` is the dimension or the name of a symbol bound to an integer.
"""
struct SpaceDecl <: TrackingStmt
    name::Symbol
    dim::Union{Int, Symbol}   # either literal or name to be resolved later
end

"""
    MapDecl(name, domain, codomain)

Declares a named linear map, e.g. `map A : X -> X`.
The domain and codomain are space names.
"""
struct MapDecl <: TrackingStmt
    name::Symbol
    domain::Symbol
    codomain::Symbol
end

"""
    VectorDecl(name, space)

Declares a named vector, e.g. `vector x0 : X`.
"""
struct VectorDecl <: TrackingStmt
    name::Symbol
    space::Symbol
end

"""
    AgentDecl(name, dynamics_A, dynamics_B, period)

Declares an agent with optional ZOH dynamics.
`dynamics_A` and `dynamics_B` are map names (or `nothing`).
`period` is the discretisation period name or literal.
"""
struct AgentDecl <: TrackingStmt
    name::Symbol
    dynamics_A::Union{Symbol, Nothing}
    dynamics_B::Union{Symbol, Nothing}
    period::Union{Symbol, Float64, Nothing}
end

"""
    TargetDecl(name)

Declares a tracking target with no dynamics (pinned by boundary).
"""
struct TargetDecl <: TrackingStmt
    name::Symbol
end

# --- Time declarations ---

"""
    HorizonDecl(name)

Declares the time horizon `horizon K`. `name` is the identifier (e.g. `:K`).
The value is supplied at resolve time via the context dict (`:K => 40`).
"""
struct HorizonDecl <: TrackingStmt
    name::Symbol
end

"""
    TimeAliasDecl(alias, ref)

Declares a named single-time alias, e.g. `time tcapture = 7`.
"""
struct TimeAliasDecl <: TrackingStmt
    alias::Symbol
    ref::TimeRef
end

"""
    TimeSetDecl(name, spec)

Declares a named time *set*, e.g. `times Tall = 0:K`.
"""
struct TimeSetDecl <: TrackingStmt
    name::Symbol
    spec::TimeSpec
end

# --- Constraints ---

"""
    ConsensusConstraint(name, agents, maps, time_spec)

Consensus edge between a pair of agents using named restriction maps,
active at the given `TimeSpec`.

`agents` is a length-2 tuple `(a1, a2)` of agent names.
`maps` is a length-2 tuple `(R1, R2)` of map names (may be equal).
"""
struct ConsensusConstraint <: TrackingStmt
    name::Symbol
    agents::NTuple{2, Symbol}
    maps::NTuple{2, Symbol}
    time_spec::TimeSpec
end

"""
    TrackConstraint(name, agent, target, maps, time_spec)

Tracking edge between an agent and a target using named restriction maps.
"""
struct TrackConstraint <: TrackingStmt
    name::Symbol
    agent::Symbol
    target::Symbol
    maps::NTuple{2, Symbol}
    time_spec::TimeSpec
end

"""
    ConsensusSheafConstraint(name, template, agents, time_spec)

Instantiate a named consensus sheaf template over a list of agents at
the given times.
"""
struct ConsensusSheafConstraint <: TrackingStmt
    name::Symbol
    template::Symbol
    agents::Vector{Symbol}
    time_spec::TimeSpec
end

# --- Boundary conditions ---

"""
    BoundaryConstraint(entity_kind, entity_name, time_ref, value_ref)

Pin the stalk of `entity` at `time` to the vector named `value_ref`.

`entity_kind` is `:agent` or `:target`.
`value_ref` is either a `PlainRef` (named vector) or an `IndexedRef` (late-bound index).
"""
struct BoundaryConstraint <: TrackingStmt
    entity_kind::Symbol          # :agent or :target
    entity_name::Symbol
    time_ref::TimeRef
    value_ref::BoundaryRef
end

# ── Program root ──────────────────────────────────────────────────────────────

"""
    TrackingProgram(statements)

Root node of a tracking DSL program. Holds a flat list of `TrackingStmt` nodes
in declaration order.
"""
struct TrackingProgram
    statements::Vector{TrackingStmt}
end

end # module TrackingDSLTerm
