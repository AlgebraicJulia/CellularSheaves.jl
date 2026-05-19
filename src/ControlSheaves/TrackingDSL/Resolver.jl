"""
    TrackingDSLResolver

Name resolution and late-binding semantics for `TrackingProgram` ASTs.

Resolution produces a `ResolvedProgram` that contains all numeric values needed
by the lowering pass.  The `initial` and `final` time aliases are first-class:
`initial` resolves to `0`; `final` resolves to the value of the horizon `k`.

Late binding is performed by `BindStmt` declarations processed in declaration
order. Indexed symbols like `x[a,t]` are resolved only after both `a` and `t`
are bound.
"""
module TrackingDSLResolver

using LinearAlgebra
using ..TrackingDSLTerm

export resolve_tracking_program, ResolvedProgram, ResolvedAgent, ResolvedTarget,
    ResolvedConsensus, ResolvedTrack, ResolvedBoundary, bind!

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
    resolve_tracking_program(prog::TrackingProgram) -> ResolvedProgram

Resolve all symbolic names in `prog` to concrete numeric values.

Resolution steps:
1. Process `BindStmt`s in declaration order to build a value environment.
2. Resolve the horizon `k`.
3. Resolve space dimensions.
4. Resolve maps to concrete matrices.
5. Resolve agent dynamics.
6. Resolve time specs to integer vectors.
7. Resolve boundary conditions (including indexed references `x[a,t]`).

Raises `TrackingUnboundSymbolError` if any required name is still unbound.
Raises `TrackingDimensionMismatchError` for inconsistent matrix dimensions.

# Example — symbolic + late binding

```julia
prog = parse_tracking_program(quote
    space X = R^2
    map A : X -> X
    agent a1 dynamics (A, B) period dt
    horizon K
    time initial = 0
    time final = K
    boundary agent a at t = x[a,t]
    bind a => 1
    bind t => 3
    bind K => 5
    bind A => [1.0 0.0; 0.0 1.0]
    bind B => [0.0; 1.0;;]
    bind dt => 0.05
    bind x[a,t] => [1.0, 0.0]
end)
resolved = resolve_tracking_program(prog)
```

`initial` resolves to `0`; `final` resolves to `K`.
"""
function resolve_tracking_program(prog::TrackingProgram)
    env = Dict{Any,Any}()
    _apply_binds!(env, prog.statements)
    return _resolve(prog, env)
end

# ─────────────────────────────────────────────────────────────────────────────
# bind! helper
# ─────────────────────────────────────────────────────────────────────────────

"""
    bind!(prog::TrackingProgram, name, value) -> TrackingProgram

Append a `BindStmt` for `name => value` to `prog.statements` in-place and
return `prog`.

This is a convenience wrapper so you can build programs incrementally:

```julia
prog = parse_tracking_program(quote
    agent a1 dynamics (A, B) period dt
    horizon K
    time initial = 0
    time final = K
    boundary agent a at t = x[a,t]
end)
bind!(prog, :K, 40)
bind!(prog, :a, 1)
bind!(prog, :t, 3)
bind!(prog, :A, [1.0 0.0; 0.0 1.0])
```

`initial` resolves to `0`; `final` resolves to the value of `K`.
"""
function bind!(prog::TrackingProgram, name::Symbol, value)
    push!(prog.statements, BindStmt(PlainRef(name), value))
    return prog
end

function bind!(prog::TrackingProgram, name::Symbol, agent_sym::Symbol, time_sym::Symbol, value)
    push!(prog.statements, BindStmt(IndexedRef(name, agent_sym, time_sym), value))
    return prog
end

# ─────────────────────────────────────────────────────────────────────────────
# Internal: apply binds
# ─────────────────────────────────────────────────────────────────────────────

function _apply_binds!(env::Dict, stmts)
    for stmt in stmts
        stmt isa BindStmt || continue
        key = _bind_key(stmt.lhs, env)
        # For indexed refs resolve the indices first
        if stmt.lhs isa IndexedRef
            # defer resolution until indices are available
            _store_indexed_bind!(env, stmt.lhs, stmt.rhs)
        else
            name = (stmt.lhs::PlainRef).name
            env[name] = stmt.rhs
        end
    end
    # Second pass: resolve any deferred indexed binds
    _resolve_deferred_indexed!(env)
end

function _bind_key(ref::PlainRef, env)
    return ref.name
end

function _bind_key(ref::IndexedRef, env)
    return ref  # used as placeholder
end

function _store_indexed_bind!(env::Dict, ref::IndexedRef, value)
    # Store raw indexed bind; key = (name, agent_val, time_val) or deferred
    push!(get!(env, :__indexed_binds__, []), (ref, value))
end

function _resolve_deferred_indexed!(env::Dict)
    binds = get(env, :__indexed_binds__, [])
    for (ref, value) in binds
        a_val = get(env, ref.agent_index, nothing)
        t_val = get(env, ref.time_index, nothing)
        if a_val !== nothing && t_val !== nothing
            key = (ref.name, Int(a_val), Int(t_val))
            env[key] = value
        else
            # Store with unresolved key for later use during boundary resolution
            key = (ref.name, ref.agent_index, ref.time_index)
            env[key] = value
        end
    end
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
    # Pre-seed with built-ins and horizon name
    time_aliases = Dict{Symbol,Int}(:initial => 0, :final => k)
    # Also add the horizon variable name as an alias for k
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
            "Horizon '$hname' is not bound. Add `bind $hname => <integer>`."))
        return Int(v)
    end
    # No HorizonDecl found — look for a direct :k bind
    v = get(env, :k, nothing)
    v !== nothing && return Int(v)
    throw(TrackingUnboundSymbolError(
        "No horizon declaration found. Add `horizon K` and `bind K => <integer>`."))
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
    if t isa InitialTime
        return 0
    elseif t isa FinalTime
        return k
    elseif t isa LiteralTime
        return t.value
    elseif t isa NamedTime
        # initial and final are built-in even if specified as NamedTime
        t.name == :initial && return 0
        t.name == :final   && return k
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
        # Also handle built-in initial/final
        if spec.name == :initial
            return [0]
        elseif spec.name == :final
            return [k]
        elseif haskey(time_sets, spec.name)
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
    if t isa InitialTime
        return 0
    elseif t isa FinalTime
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
            "Time reference '$(t.name)' is not bound. Add `bind $(t.name) => <int>`."))
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
    # Try to find (name, concrete_a, concrete_t) key first
    a_val = get(env, ref.agent_index, nothing)
    t_val = get(env, ref.time_index,  nothing)
    if a_val !== nothing && t_val !== nothing
        key = (ref.name, Int(a_val), Int(t_val))
        v = get(env, key, nothing)
        v !== nothing && return Vector{Float64}(v)
    end
    # Fall back to symbolic key lookup
    sym_key = (ref.name, ref.agent_index, ref.time_index)
    v = get(env, sym_key, nothing)
    v !== nothing && return Vector{Float64}(v)
    throw(TrackingUnboundSymbolError(
        "Boundary indexed value '$(ref.name)[$(ref.agent_index),$(ref.time_index)]' " *
        "for entity '$entity_name' at time $t is not bound."))
end

# ─────────────────────────────────────────────────────────────────────────────
# Utility helpers
# ─────────────────────────────────────────────────────────────────────────────

function _require_matrix(name::Symbol, source::Dict, ctx::String)::Matrix{Float64}
    v = get(source, name, nothing)
    v !== nothing || throw(TrackingUnboundSymbolError(
        "$ctx: matrix '$name' is not bound. " *
        "Use `bind!(prog, :$name, <matrix>)` or `@tracking_problem` with `bind($name => <matrix>)`."))
    v isa Expr && throw(TrackingUnboundSymbolError(
        "$ctx: matrix '$name' was bound to an unevaluated expression. " *
        "Use `@tracking_problem` (which evaluates bind values) or `bind!(prog, :$name, <matrix>)`."))
    return Matrix{Float64}(v)
end

function _require_index(name::Symbol, index_map::Dict, ctx::String)::Int
    v = get(index_map, name, nothing)
    v !== nothing || throw(TrackingDeclarationError(
        "$ctx: '$name' is not declared."))
    return v
end

end # module TrackingDSLResolver
