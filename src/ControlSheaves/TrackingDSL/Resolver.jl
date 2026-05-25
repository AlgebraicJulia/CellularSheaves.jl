"""
    TrackingDSLResolver

Name resolution for `TrackingProgram` ASTs.

Resolution produces a `ResolvedProgram` that contains all numeric values needed
by the lowering pass.  The `t[begin]` and `t[end]` time references are
first-class: `t[begin]` resolves to `0`; `t[end]` resolves to the value of the
horizon `k`.

All numeric values (matrices, scalars, vectors) are supplied through an external
context dict passed to `resolve_tracking_program`.  Use [`set_indexed!`](@ref) to
register indexed boundary values.
"""
module TrackingDSLResolver

using LinearAlgebra
using ..TrackingDSLTerm

export resolve_tracking_program, ResolvedProgram, ResolvedAgent, ResolvedTarget,
    ResolvedConsensus, ResolvedTrack, ResolvedBoundary, set_indexed!

# ─────────────────────────────────────────────────────────────────────────────
# Resolved data structures
# ─────────────────────────────────────────────────────────────────────────────

"""
    ResolvedAgent

A fully resolved agent declaration with concrete numeric matrices.
"""
struct ResolvedAgent
    name::Symbol
    Ad::Matrix{Float64}
    Bd::Matrix{Float64}
    nx::Int
    nu::Int
end

"""
    ResolvedTarget

A fully resolved target declaration (no dynamics).
"""
struct ResolvedTarget
    name::Symbol
    stalk_dim::Int
end

"""
    ResolvedConsensus

A fully resolved consensus constraint.
"""
struct ResolvedConsensus
    name::Symbol
    agent_indices::NTuple{2, Int}
    restriction_maps::NTuple{2, Matrix{Float64}}
    timesteps::Vector{Int}
end

"""
    ResolvedTrack

A fully resolved tracking constraint.
"""
struct ResolvedTrack
    name::Symbol
    agent_index::Int
    target_index::Int
    restriction_maps::NTuple{2, Matrix{Float64}}
    timesteps::Vector{Int}
end

"""
    ResolvedBoundary

A single boundary condition: (vertex_kind, entity_index, time, stalk_vector).
"""
struct ResolvedBoundary
    entity_kind::Symbol   # :agent or :target
    entity_index::Int
    time::Int
    value::Vector{Float64}
end

"""
    ResolvedProgram

The output of the resolution phase.  All names have been replaced by concrete
numeric values.
"""
struct ResolvedProgram
    k::Int
    agents::Vector{ResolvedAgent}
    targets::Vector{ResolvedTarget}
    consensus::Vector{ResolvedConsensus}
    tracking::Vector{ResolvedTrack}
    boundaries::Vector{ResolvedBoundary}
    stalk_dim::Int      # common stalk dimension (nx + nu)
end

# ─────────────────────────────────────────────────────────────────────────────
# Resolution entry point
# ─────────────────────────────────────────────────────────────────────────────

"""
    resolve_tracking_program(prog::TrackingProgram, ctx::AbstractDict) -> ResolvedProgram

Resolve all symbolic names in `prog` to concrete numeric values using the
supplied context dict `ctx`.

`ctx` maps `Symbol` keys to their values (scalars, matrices, vectors).  Use
[`set_indexed!`](@ref) to register indexed boundary values:

```julia
ctx = Dict{Symbol,Any}(
    :K  => 5,
    :A  => [1.0 0.0; 0.0 1.0],
    :B  => [0.0; 1.0;;],
    :dt => 0.05,
)
set_indexed!(ctx, :x_ref, 1, 3, [1.0, 2.0, 0.0])  # x_ref[a=1, t=3]
resolved = resolve_tracking_program(prog, ctx)
```

Resolution steps:
1. Build value environment from `ctx`.
2. Resolve the horizon `k`.
3. Resolve space dimensions.
4. Resolve maps to concrete matrices.
5. Resolve agent dynamics.
6. Resolve time specs to integer vectors.
7. Resolve boundary conditions (including indexed references `x[a,t]`).

Raises `TrackingUnboundSymbolError` if any required name is missing from `ctx`.
Raises `TrackingDimensionMismatchError` for inconsistent matrix dimensions.

`t[begin]` resolves to `0`; `t[end]` resolves to the value of the horizon `K`.
"""
function resolve_tracking_program(prog::TrackingProgram, ctx::AbstractDict)
    env = Dict{Any,Any}(ctx)
    return _resolve(prog, env)
end

# ─────────────────────────────────────────────────────────────────────────────
# set_indexed! helper
# ─────────────────────────────────────────────────────────────────────────────

"""
    set_indexed!(ctx, name, agent_val, time_val, vec)

Store an indexed boundary value in `ctx` under the canonical key
`(name, agent_val, time_val)`.  Use this to supply values for boundary
conditions declared with an indexed reference like `boundary(:agent, a; at=t, value=x[a,t])`.

Because the key is a `Tuple`, `ctx` must support non-`Symbol` keys.
Use `Dict{Any,Any}` (not `Dict{Symbol,Any}`) when indexed values are needed:

```julia
ctx = Dict{Any,Any}(:K => 5, :A => ..., :B => ..., :dt => 0.05, :a => 1, :t_pin => 3)
set_indexed!(ctx, :x_ref, 1, 3, [0.0, 1.0, 0.0])
resolved = resolve_tracking_program(prog, ctx)
```
"""
function set_indexed!(ctx::AbstractDict, name::Symbol, agent_val::Int, time_val::Int, vec)
    ctx[(name, agent_val, time_val)] = vec
    return ctx
end

# ─────────────────────────────────────────────────────────────────────────────
# Internal: main resolution
# ─────────────────────────────────────────────────────────────────────────────

function _resolve(prog::TrackingProgram, env::Dict)
    # 1. Resolve horizon
    k = _resolve_horizon(prog, env)

    # 2. Collect space dimensions
    space_dims = Dict{Symbol, Union{Int, SpaceExpr}}()
    for stmt in prog.statements
        stmt isa SpaceDecl || continue
        d = _resolve_space_dim(stmt.dim, env, space_dims)
        space_dims[stmt.name] = d
    end

    # 3. Resolve map declarations
    map_matrices = Dict{Symbol, Matrix{Float64}}()
    for stmt in prog.statements
        stmt isa MapDecl || continue
        # map must be bound via bind statement
        mat = _require_matrix(stmt.name, env, "map '$(stmt.name)'")
        map_matrices[stmt.name] = mat
    end
    # Also allow maps that appear only in bind statements (without explicit MapDecl)
    for (k2, v) in env
        k2 isa Symbol || continue
        v isa Matrix && !haskey(map_matrices, k2) && (map_matrices[k2] = Matrix{Float64}(v))
    end

    # 4. Resolve agents
    agent_order = Symbol[]
    for stmt in prog.statements
        stmt isa AgentDecl && push!(agent_order, stmt.name)
    end
    resolved_agents = ResolvedAgent[]
    for stmt in prog.statements
        stmt isa AgentDecl || continue
        ra = _resolve_agent(stmt, map_matrices, env)
        push!(resolved_agents, ra)
    end

    # 5. Resolve targets
    resolved_targets = ResolvedTarget[]
    for stmt in prog.statements
        stmt isa TargetDecl || continue
        # Stalk dim = same as agent stalk (nx + nu)
        # Will be determined from agent dims at lowering
        push!(resolved_targets, ResolvedTarget(stmt.name, 0))
    end

    # 6. Common stalk dim
    stalk_dim = 0
    if !isempty(resolved_agents)
        stalk_dim = resolved_agents[1].nx + resolved_agents[1].nu
        for ra in resolved_agents
            ra.nx + ra.nu == stalk_dim ||
                throw(TrackingDimensionMismatchError(
                    "All agents must have the same stalk dimension (nx+nu), " *
                    "but agent '$(ra.name)' has $(ra.nx + ra.nu) ≠ $stalk_dim"))
        end
    end

    # Update target stalk dims
    resolved_targets = [ResolvedTarget(rt.name, stalk_dim) for rt in resolved_targets]

    # Index maps
    agent_index  = Dict(ra.name => i for (i, ra) in enumerate(resolved_agents))
    target_index = Dict(rt.name => i for (i, rt) in enumerate(resolved_targets))

    # 7. Resolve time specs
    # Pre-seed with horizon name (resolves both K => k and k => k lookups)
    time_aliases = Dict{Symbol,Int}()
    # Add the horizon variable name as an alias for k
    for stmt in prog.statements
        stmt isa HorizonDecl || continue
        time_aliases[stmt.name] = k
    end
    time_sets    = Dict{Symbol,Vector{Int}}()
    for stmt in prog.statements
        if stmt isa TimeAliasDecl
            time_aliases[stmt.alias] = _resolve_time_ref(stmt.ref, time_aliases, k)
        elseif stmt isa TimeSetDecl
            time_sets[stmt.name] = _resolve_time_spec(stmt.spec, time_aliases, time_sets, k)
        end
    end

    # 8. Resolve constraints
    resolved_consensus = ResolvedConsensus[]
    for stmt in prog.statements
        stmt isa ConsensusConstraint || continue
        ai = _require_index(stmt.agents[1], agent_index, "consensus '$(stmt.name)' agent")
        aj = _require_index(stmt.agents[2], agent_index, "consensus '$(stmt.name)' agent")
        Ri = _require_matrix(stmt.maps[1], map_matrices, "consensus '$(stmt.name)' map")
        Rj = _require_matrix(stmt.maps[2], map_matrices, "consensus '$(stmt.name)' map")
        ts = _resolve_time_spec(stmt.time_spec, time_aliases, time_sets, k)
        push!(resolved_consensus, ResolvedConsensus(stmt.name, (ai, aj), (Ri, Rj), ts))
    end

    resolved_tracking = ResolvedTrack[]
    for stmt in prog.statements
        stmt isa TrackConstraint || continue
        ai = _require_index(stmt.agent,  agent_index,  "track '$(stmt.name)' agent")
        ti = _require_index(stmt.target, target_index, "track '$(stmt.name)' target")
        Ra = _require_matrix(stmt.maps[1], map_matrices, "track '$(stmt.name)' agent map")
        Rt = _require_matrix(stmt.maps[2], map_matrices, "track '$(stmt.name)' target map")
        ts = _resolve_time_spec(stmt.time_spec, time_aliases, time_sets, k)
        push!(resolved_tracking, ResolvedTrack(stmt.name, ai, ti, (Ra, Rt), ts))
    end

    # 9. Resolve boundary conditions
    resolved_boundaries = ResolvedBoundary[]
    for stmt in prog.statements
        stmt isa BoundaryConstraint || continue
        rb = _resolve_boundary(stmt, agent_index, target_index, time_aliases, k, env)
        push!(resolved_boundaries, rb)
    end

    return ResolvedProgram(k, resolved_agents, resolved_targets,
        resolved_consensus, resolved_tracking, resolved_boundaries, stalk_dim)
end

# ─────────────────────────────────────────────────────────────────────────────
# Horizon resolution
# ─────────────────────────────────────────────────────────────────────────────

function _resolve_horizon(prog::TrackingProgram, env::Dict)
    for stmt in prog.statements
        stmt isa HorizonDecl || continue
        hname = stmt.name
        v = get(env, hname, nothing)
        v !== nothing || throw(TrackingUnboundSymbolError(
            "Horizon '$hname' is not bound. Add `:$hname => <integer>` to the context dict."))
        return Int(v)
    end
    throw(TrackingUnboundSymbolError(
        "No horizon declaration found. Add `horizon K` to the program and `:K => <integer>` to the context dict."))
end

# ─────────────────────────────────────────────────────────────────────────────
# Space dimension resolution
# ─────────────────────────────────────────────────────────────────────────────

function _resolve_space_dim(dim, env, space_dims)
    if dim isa Int
        return dim
    elseif dim isa Symbol
        v = get(env, dim, nothing)
        v !== nothing || throw(TrackingUnboundSymbolError(
            "Space dimension symbol '$dim' is not bound."))
        return Int(v)
    elseif dim isa DirectSumSpace
        a = _resolve_space_expr(dim.a, space_dims, env)
        b = _resolve_space_expr(dim.b, space_dims, env)
        return a + b
    elseif dim isa ProductSpace
        a = _resolve_space_expr(dim.a, space_dims, env)
        b = _resolve_space_expr(dim.b, space_dims, env)
        return a + b  # same as direct sum for euclidean spaces
    elseif dim isa BaseSpace
        return _resolve_space_expr(dim, space_dims, env)
    else
        throw(TrackingTypeError("Cannot resolve space dimension: $dim"))
    end
end

function _resolve_space_expr(expr, space_dims, env)
    if expr isa BaseSpace
        n = expr.name
        haskey(space_dims, n) ||
            throw(TrackingDeclarationError("Unknown space '$n'"))
        d = space_dims[n]
        d isa Int && return d
        return _resolve_space_dim(d, env, space_dims)
    elseif expr isa Int
        return expr
    else
        throw(TrackingTypeError("Cannot resolve space expression: $expr"))
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Agent resolution
# ─────────────────────────────────────────────────────────────────────────────

function _resolve_agent(stmt::AgentDecl, map_matrices, env)
    name = stmt.name
    A_name = stmt.dynamics_A
    B_name = stmt.dynamics_B
    (A_name === nothing || B_name === nothing) &&
        throw(TrackingUnboundSymbolError(
            "Agent '$name' has no dynamics. Add `dynamics (A, B) period dt`."))
    Ad = _require_matrix(A_name, map_matrices, "agent '$name' dynamics A")
    Bd = _require_matrix(B_name, map_matrices, "agent '$name' dynamics B")
    nx = size(Ad, 1)
    nu = size(Bd, 2)
    size(Ad, 2) == nx ||
        throw(TrackingDimensionMismatchError("Agent '$name': Ad must be square, got $(size(Ad))"))
    size(Bd, 1) == nx ||
        throw(TrackingDimensionMismatchError(
            "Agent '$name': Bd must have $nx rows (same as nx), got $(size(Bd, 1))"))
    return ResolvedAgent(name, Matrix{Float64}(Ad), Matrix{Float64}(Bd), nx, nu)
end

# ─────────────────────────────────────────────────────────────────────────────
# Time resolution
# ─────────────────────────────────────────────────────────────────────────────

function _resolve_time_ref(t::TimeRef, aliases::Dict{Symbol,Int}, k::Int)::Int
    if t isa BeginTime
        return 0
    elseif t isa EndTime
        return k
    elseif t isa LiteralTime
        return t.value
    elseif t isa NamedTime
        haskey(aliases, t.name) ||
            throw(TrackingTimeResolutionError("Undefined time alias '$(t.name)'"))
        return aliases[t.name]
    else
        throw(TrackingTimeResolutionError("Unknown TimeRef type: $(typeof(t))"))
    end
end

function _resolve_time_spec(spec::TimeSpec, aliases::Dict{Symbol,Int},
                             time_sets::Dict{Symbol,Vector{Int}}, k::Int)::Vector{Int}
    if spec isa SingletonTime
        return [_resolve_time_ref(spec.t, aliases, k)]
    elseif spec isa TimeList
        return [_resolve_time_ref(t, aliases, k) for t in spec.ts]
    elseif spec isa TimeRange
        a = _resolve_time_ref(spec.lo, aliases, k)
        b = _resolve_time_ref(spec.hi, aliases, k)
        return collect(a:b)
    elseif spec isa NamedTimeSet
        # NamedTimeSet can refer to either a declared time set OR a time alias (singleton)
        if haskey(time_sets, spec.name)
            return time_sets[spec.name]
        elseif haskey(aliases, spec.name)
            return [aliases[spec.name]]
        else
            throw(TrackingTimeResolutionError(
                "Undefined time set or alias '$(spec.name)'. " *
                "Declare it with `time($(spec.name) = ...)` or `times($(spec.name) = ...)`."))
        end
    else
        throw(TrackingTimeResolutionError("Unknown TimeSpec type: $(typeof(spec))"))
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Boundary condition resolution
# ─────────────────────────────────────────────────────────────────────────────

function _resolve_boundary(stmt::BoundaryConstraint, agent_index, target_index,
                            time_aliases, k, env)
    entity_kind = stmt.entity_kind
    entity_name = stmt.entity_name

    if entity_kind == :agent
        entity_idx = _require_index(entity_name, agent_index, "boundary agent")
    elseif entity_kind == :target
        entity_idx = _require_index(entity_name, target_index, "boundary target")
    else
        throw(TrackingDeclarationError("Unknown entity kind: $entity_kind"))
    end

    t = _resolve_time_ref_boundary(stmt.time_ref, time_aliases, k, env)

    val = _resolve_boundary_value(stmt.value_ref, env, entity_name, t)

    return ResolvedBoundary(entity_kind, entity_idx, t, val)
end

function _resolve_time_ref_boundary(t::TimeRef, aliases::Dict{Symbol,Int}, k::Int, env::Dict)::Int
    if t isa BeginTime
        return 0
    elseif t isa EndTime
        return k
    elseif t isa LiteralTime
        return t.value
    elseif t isa NamedTime
        # Could be a time alias or a late-bound scalar
        if haskey(aliases, t.name)
            return aliases[t.name]
        end
        v = get(env, t.name, nothing)
        v !== nothing && return Int(v)
        throw(TrackingTimeResolutionError(
            "Time reference '$(t.name)' is not bound. Add it to the context dict."))
    else
        throw(TrackingTimeResolutionError("Unknown TimeRef type: $(typeof(t))"))
    end
end

function _resolve_boundary_value(ref::PlainRef, env::Dict, entity_name::Symbol, t::Int)
    v = get(env, ref.name, nothing)
    v !== nothing || throw(TrackingUnboundSymbolError(
        "Boundary value '$(ref.name)' for entity '$entity_name' at time $t is not bound."))
    return Vector{Float64}(v)
end

function _resolve_boundary_value(ref::IndexedRef, env::Dict, entity_name::Symbol, t::Int)
    # Resolve agent and time index values from env, then look up canonical key
    a_val = get(env, ref.agent_index, nothing)
    t_val = get(env, ref.time_index,  nothing)
    if a_val !== nothing && t_val !== nothing
        key = (ref.name, Int(a_val), Int(t_val))
        v = get(env, key, nothing)
        v !== nothing && return Vector{Float64}(v)
    end
    throw(TrackingUnboundSymbolError(
        "Boundary indexed value '$(ref.name)[$(ref.agent_index),$(ref.time_index)]' " *
        "for entity '$entity_name' at time $t is not bound. " *
        "Use set_indexed!(ctx, :$(ref.name), agent_val, time_val, vec)."))
end

# ─────────────────────────────────────────────────────────────────────────────
# Utility helpers
# ─────────────────────────────────────────────────────────────────────────────

function _require_matrix(name::Symbol, source::Dict, ctx::String)::Matrix{Float64}
    v = get(source, name, nothing)
    v !== nothing || throw(TrackingUnboundSymbolError(
        "$ctx: matrix '$name' is not bound. Add `:$name => <matrix>` to the context dict."))
    v isa Expr && throw(TrackingUnboundSymbolError(
        "$ctx: matrix '$name' was bound to an unevaluated expression."))
    return Matrix{Float64}(v)
end

function _require_index(name::Symbol, index_map::Dict, ctx::String)::Int
    v = get(index_map, name, nothing)
    v !== nothing || throw(TrackingDeclarationError(
        "$ctx: '$name' is not declared."))
    return v
end

end # module TrackingDSLResolver
