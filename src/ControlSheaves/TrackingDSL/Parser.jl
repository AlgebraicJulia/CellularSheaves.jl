"""
    TrackingDSLParser

MLStyle-driven parser for the name-driven Tracking DSL.

Provides the `@tracking_problem` macro and the functional `parse_tracking_program`
entry point.  All statement forms use **valid Julia syntax** so that Julia's own
parser tokenises the block before any DSL parsing occurs.

## Grammar

Statements inside a `@tracking_problem begin ... end` block must be valid Julia
expressions. The DSL interprets them by their leading function name:

```julia
@tracking_problem begin
    # Space declaration (method-definition syntax)
    space(X) = R^6
    space(U) = R^2
    space(Stalk) = X ⊕ U      # direct sum; also X * U or X × U

    # Map declaration  (use map_decl to avoid shadowing Base.map)
    map_decl(A, X, X)          # A : X → X
    map_decl(B, U, X)          # B : U → X
    map_decl(R_y, Stalk, R^1)

    # Agent / target declarations
    agent(a1; dynamics=(A,B), period=dt)
    agent(a2; dynamics=(A,B), period=dt)
    target(t1)
    target(t2)

    # Horizon
    horizon(K)

    # Time aliases (keyword-argument syntax)
    time(tcapture = 7)

    # Time sets
    times(Tall = 0:K)

    # Constraints
    consensus(c1; agents=(a1,a2), maps=(R_y,R_y), at=Tall)
    track(tr1; agent=a1, target=t1, maps=(R_z,R_z), at=Tall[end])
    consensus_sheaf(cS; template=CTemplate, over=(a1,a2,a3), at=tcapture)

    # Boundary conditions
    boundary(:agent, a1; at=Tall[begin], value=x0_a1)
    boundary(:target, t1; at=t, value=x_ref[t])
    boundary(:agent, a; at=t, value=x[a,t])   # indexed reference
end
```

### Time specifications

At any `at=` keyword, time can be given as:
- a named time set:      `at=Tall` (NamedTimeSet)
- begin/end of horizon:  `at=t[begin]` (BeginTime = 0), `at=t[end]` (EndTime = k)
- a literal integer:     `at=0`
- a range:               `at=(0:K)`
- an explicit list:      `at=[0,3]`

`t[begin]` and `t[end]` use Julia array-indexing syntax; the qualifying name `t`
is ignored and the references resolve to `0` and the horizon value `k`
respectively.
"""
module TrackingDSLParser

using MLStyle: @match

using ..TrackingDSLTerm

export @tracking_problem, parse_tracking_program

# ─────────────────────────────────────────────────────────────────────────────
# Public API
# ─────────────────────────────────────────────────────────────────────────────

"""
    @tracking_problem begin ... end

Parse a TrackingDSL program and return a `TrackingProgram` AST.

Statements are valid Julia expressions interpreted by their leading name.
Numeric values are **not** embedded in the program; pass them via the `ctx`
argument of `resolve_tracking_program` or `lower_tracking_program`.

The special time references `t[begin]` and `t[end]` are built-in: `t[begin]`
resolves to `0`; `t[end]` resolves to the horizon value `k` (from `horizon(K)` in
the program and `:K` in the context dict). The qualifying name `t` is syntactic
decoration and is ignored.

# Example

```julia
prog = @tracking_problem begin
    space(X) = R^2
    map_decl(A, X, X)
    map_decl(R_y, X, X)
    agent(a1; dynamics=(A,R_y), period=dt)
    agent(a2; dynamics=(A,R_y), period=dt)
    target(t1)
    horizon(K)
    times(Tall = 0:K)
    consensus(c1; agents=(a1,a2), maps=(R_y,R_y), at=Tall)
    track(tr1; agent=a1, target=t1, maps=(R_y,R_y), at=Tall[end])
    boundary(:agent, a1; at=Tall[begin], value=x0_a1)
end

ctx = Dict{Symbol,Any}(
    :K    => 40,
    :A    => [1.0 0.0; 0.0 1.0],
    :R_y  => [1.0 0.0; 0.0 0.0],
    :dt   => 0.05,
    :x0_a1 => zeros(2),
)
result = lower_tracking_program(prog, ctx)
```

`t[begin]` resolves to `0`; `t[end]` resolves to the value of `K`.
"""
macro tracking_problem(block::Expr)
    stmt_exprs = Any[]
    for line in block.args
        line isa LineNumberNode && continue
        stmt = _parse_stmt(line)
        stmt === nothing && continue
        push!(stmt_exprs, QuoteNode(stmt))
    end
    return :(TrackingProgram(TrackingStmt[$(stmt_exprs...)]))
end

"""
    parse_tracking_program(block::Expr) -> TrackingProgram

Functional entry point for parsing a TrackingDSL program from a Julia `Expr`
(as returned by `Meta.quot(quote ... end)`).

Returns a `TrackingProgram` AST without validating or resolving names.  Supply
numeric values via the `ctx` argument of `resolve_tracking_program` or
`lower_tracking_program`:

```julia
prog = parse_tracking_program(quote
    agent(a1; dynamics=(A,B), period=dt)
    horizon(K)
    times(Tall = 0:K)
end)
ctx = Dict{Symbol,Any}(:K => 40, :A => [1.0 0.0; 0.0 1.0], :B => [0.0; 1.0;;], :dt => 0.05)
resolved = resolve_tracking_program(prog, ctx)
```

# Example — indexed boundary reference

```julia
prog = parse_tracking_program(quote
    agent(a1; dynamics=(A,B), period=dt)
    horizon(K)
    boundary(:agent, a1; at=t_pin, value=x_ref[a,t_pin])
end)
ctx = Dict{Symbol,Any}(:K => 5, :A => ..., :B => ..., :dt => 0.05, :a => 1, :t_pin => 3)
set_indexed!(ctx, :x_ref, 1, 3, [0.0, 1.0, 0.0, 0.0])
resolved = resolve_tracking_program(prog, ctx)
```

`t[begin]` resolves to `0`; `t[end]` resolves to the horizon `K`.
"""
function parse_tracking_program(block::Expr)
    # Unwrap :quote and :block wrappings (from quote ... end passed as argument)
    inner = block
    if inner.head == :quote && length(inner.args) == 1
        inner = inner.args[1]
    end
    stmts = TrackingStmt[]
    for line in inner.args
        line isa LineNumberNode && continue
        stmt = _parse_stmt(line)
        stmt === nothing || push!(stmts, stmt)
    end
    return TrackingProgram(stmts)
end

# ─────────────────────────────────────────────────────────────────────────────
# Top-level statement dispatch
# ─────────────────────────────────────────────────────────────────────────────

function _parse_stmt(line)
    @match line begin
        # ── space(X) = R^n  (method-def form) ──────────────────────────────
        Expr(:(=), Expr(:call, :space, name::Symbol), rhs) =>
            SpaceDecl(name, _parse_space_dim(rhs))

        # ── map_decl(A, Dom, Cod) ───────────────────────────────────────────
        Expr(:call, :map_decl, name::Symbol, dom::Symbol, cod) =>
            MapDecl(name, dom, _extract_space_sym(cod))

        # ── agent(name; dynamics=(A,B), period=dt) ──────────────────────────
        Expr(:call, :agent, Expr(:parameters, kws...), name::Symbol) =>
            _parse_agent(name, kws)

        # ── target(name) ────────────────────────────────────────────────────
        Expr(:call, :target, name::Symbol) =>
            TargetDecl(name)

        # ── horizon(K) ──────────────────────────────────────────────────────
        Expr(:call, :horizon, name::Symbol) =>
            HorizonDecl(name)

        # ── time(alias = value)  [keyword-arg form] ─────────────────────────
        Expr(:call, :time, Expr(:parameters, kws...)) =>
            _parse_time_alias_kwarg(kws)
        Expr(:call, :time, kws...) =>
            _parse_time_alias_args(kws)

        # ── times(Name = lo:hi) ─────────────────────────────────────────────
        Expr(:call, :times, Expr(:parameters, kws...)) =>
            _parse_time_set_kwarg(kws)
        Expr(:call, :times, kws...) =>
            _parse_time_set_args(kws)

        # ── consensus(name; agents=(...), maps=(...), at=...) ───────────────
        Expr(:call, :consensus, Expr(:parameters, kws...), cname::Symbol) =>
            _parse_consensus(cname, kws)

        # ── track(name; agent=a, target=t, maps=(...), at=...) ──────────────
        Expr(:call, :track, Expr(:parameters, kws...), tname::Symbol) =>
            _parse_track(tname, kws)

        # ── consensus_sheaf(name; template=T, over=(...), at=...) ───────────
        Expr(:call, :consensus_sheaf, Expr(:parameters, kws...), sname::Symbol) =>
            _parse_consensus_sheaf(sname, kws)

        # ── boundary(:agent, name; at=time, value=expr) ─────────────────────
        Expr(:call, :boundary, Expr(:parameters, kws...), kind_lit, ename::Symbol) =>
            _parse_boundary(kind_lit, ename, kws)

        # ── Skip nothing ────────────────────────────────────────────────────
        ::Nothing => nothing

        _ => throw(TrackingSyntaxError("Unrecognised DSL statement: $(repr(line))"))
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Block unwrapper: `f(x) = y` in a quote block wraps `y` in Expr(:block, ...)
# ─────────────────────────────────────────────────────────────────────────────

function _unwrap_block(expr)
    if expr isa Expr && expr.head == :block
        # Filter out LineNumberNodes and grab the single non-trivial argument
        real_args = filter(a -> !(a isa LineNumberNode), expr.args)
        length(real_args) == 1 && return real_args[1]
    end
    return expr
end

# ─────────────────────────────────────────────────────────────────────────────
# Space dimension parsing
# ─────────────────────────────────────────────────────────────────────────────

function _parse_space_dim(rhs)
    # In a quote block, `space(X) = R^n` is a method def and the body is wrapped in
    # Expr(:block, LineNumberNode, actual_rhs). Unwrap that first.
    rhs2 = _unwrap_block(rhs)
    @match rhs2 begin
        Expr(:call, :^, :R, n::Int)    => n
        Expr(:call, :^, :R, n::Symbol) => n           # R^K (symbolic)
        Expr(:call, :⊕, a, b) => DirectSumSpace(_parse_space_ref(a), _parse_space_ref(b))
        Expr(:call, :×, a, b) => ProductSpace(_parse_space_ref(a), _parse_space_ref(b))
        Expr(:call, :*, a::Symbol, b::Symbol) => ProductSpace(BaseSpace(a), BaseSpace(b))
        s::Symbol => BaseSpace(s)
        n::Int    => n
        _ => throw(TrackingSyntaxError("Cannot parse space dimension: $(repr(rhs2))"))
    end
end

_parse_space_ref(s::Symbol) = BaseSpace(s)
_parse_space_ref(e) = throw(TrackingSyntaxError("Expected space name, got: $(repr(e))"))

function _extract_space_sym(cod)
    @match cod begin
        s::Symbol => s
        Expr(:call, :^, :R, n::Int)    => Symbol("R$n")  # normalise R^1 to R1
        Expr(:call, :^, :R, n::Symbol) => n
        _ => throw(TrackingSyntaxError("Cannot extract space symbol from codomain: $(repr(cod))"))
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Agent parsing
# ─────────────────────────────────────────────────────────────────────────────

function _parse_agent(name::Symbol, kws)
    dyn_A   = nothing
    dyn_B   = nothing
    period  = nothing
    for kw in kws
        @match kw begin
            Expr(:kw, :dynamics, Expr(:tuple, a::Symbol, b::Symbol)) => begin
                dyn_A = a; dyn_B = b
            end
            Expr(:kw, :period, p::Symbol) => (period = p)
            Expr(:kw, :period, p::Number) => (period = Float64(p))
            _ => throw(TrackingSyntaxError("agent '$name': unrecognised keyword: $(repr(kw))"))
        end
    end
    return AgentDecl(name, dyn_A, dyn_B, period)
end

# ─────────────────────────────────────────────────────────────────────────────
# Time alias parsing: time(alias = value)
# ─────────────────────────────────────────────────────────────────────────────

function _parse_time_alias_kwarg(kws)
    length(kws) == 1 || throw(TrackingSyntaxError(
        "time() takes exactly one keyword argument, got $(length(kws))"))
    kw = kws[1]
    @match kw begin
        Expr(:kw, alias::Symbol, rhs) => TimeAliasDecl(alias, _parse_time_ref(rhs))
        _ => throw(TrackingSyntaxError("time(): expected alias = value, got: $(repr(kw))"))
    end
end

function _parse_time_alias_args(args)
    length(args) == 1 || throw(TrackingSyntaxError(
        "time() takes exactly one keyword argument, got $(length(args))"))
    arg = args[1]
    @match arg begin
        Expr(:kw, alias::Symbol, rhs) => TimeAliasDecl(alias, _parse_time_ref(rhs))
        Expr(:(=), alias::Symbol, rhs) => TimeAliasDecl(alias, _parse_time_ref(rhs))
        _ => throw(TrackingSyntaxError("time(): expected alias = value, got: $(repr(arg))"))
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Time set parsing: times(Name = lo:hi)
# ─────────────────────────────────────────────────────────────────────────────

function _parse_time_set_kwarg(kws)
    length(kws) == 1 || throw(TrackingSyntaxError(
        "times() takes exactly one keyword argument, got $(length(kws))"))
    kw = kws[1]
    @match kw begin
        Expr(:kw, name::Symbol, rhs) => TimeSetDecl(name, _parse_time_spec_expr(rhs))
        _ => throw(TrackingSyntaxError("times(): expected Name = spec, got: $(repr(kw))"))
    end
end

function _parse_time_set_args(args)
    length(args) == 1 || throw(TrackingSyntaxError(
        "times() takes exactly one argument, got $(length(args))"))
    arg = args[1]
    @match arg begin
        Expr(:kw, name::Symbol, rhs) => TimeSetDecl(name, _parse_time_spec_expr(rhs))
        Expr(:(=), name::Symbol, rhs) => TimeSetDecl(name, _parse_time_spec_expr(rhs))
        _ => throw(TrackingSyntaxError("times(): expected Name = spec, got: $(repr(arg))"))
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Consensus parsing: consensus(name; agents=(a1,a2), maps=(R1,R2), at=Tall)
# ─────────────────────────────────────────────────────────────────────────────

function _parse_consensus(name::Symbol, kws)
    agents_val = nothing
    maps_val   = nothing
    at_val     = nothing
    for kw in kws
        @match kw begin
            Expr(:kw, :agents, v) => (agents_val = v)
            Expr(:kw, :maps,   v) => (maps_val   = v)
            Expr(:kw, :at,     v) => (at_val     = v)
            _ => throw(TrackingSyntaxError("consensus '$name': unrecognised keyword: $(repr(kw))"))
        end
    end
    agents_val !== nothing || throw(TrackingSyntaxError("consensus '$name': missing 'agents'"))
    maps_val   !== nothing || throw(TrackingSyntaxError("consensus '$name': missing 'maps'"))
    at_val     !== nothing || throw(TrackingSyntaxError("consensus '$name': missing 'at'"))
    agents = _parse_sym_tuple2(agents_val, "consensus '$name' agents")
    maps   = _parse_sym_tuple2(maps_val,   "consensus '$name' maps")
    ts     = _parse_at_spec(at_val)
    return ConsensusConstraint(name, agents, maps, ts)
end

# ─────────────────────────────────────────────────────────────────────────────
# Track parsing: track(name; agent=a, target=t, maps=(R1,R2), at=spec)
# ─────────────────────────────────────────────────────────────────────────────

function _parse_track(name::Symbol, kws)
    agent_val  = nothing
    target_val = nothing
    maps_val   = nothing
    at_val     = nothing
    for kw in kws
        @match kw begin
            Expr(:kw, :agent,  v) => (agent_val  = v)
            Expr(:kw, :target, v) => (target_val = v)
            Expr(:kw, :maps,   v) => (maps_val   = v)
            Expr(:kw, :at,     v) => (at_val     = v)
            _ => throw(TrackingSyntaxError("track '$name': unrecognised keyword: $(repr(kw))"))
        end
    end
    agent_val  !== nothing || throw(TrackingSyntaxError("track '$name': missing 'agent'"))
    target_val !== nothing || throw(TrackingSyntaxError("track '$name': missing 'target'"))
    maps_val   !== nothing || throw(TrackingSyntaxError("track '$name': missing 'maps'"))
    at_val     !== nothing || throw(TrackingSyntaxError("track '$name': missing 'at'"))
    agent_sym  = agent_val  isa Symbol ? agent_val  : throw(TrackingSyntaxError("track '$name': agent must be a symbol"))
    target_sym = target_val isa Symbol ? target_val : throw(TrackingSyntaxError("track '$name': target must be a symbol"))
    maps = _parse_sym_tuple2(maps_val, "track '$name' maps")
    ts   = _parse_at_spec(at_val)
    return TrackConstraint(name, agent_sym, target_sym, maps, ts)
end

# ─────────────────────────────────────────────────────────────────────────────
# Consensus sheaf: consensus_sheaf(name; template=T, over=(...), at=spec)
# ─────────────────────────────────────────────────────────────────────────────

function _parse_consensus_sheaf(name::Symbol, kws)
    tmpl_val  = nothing
    over_val  = nothing
    at_val    = nothing
    for kw in kws
        @match kw begin
            Expr(:kw, :template, v) => (tmpl_val = v)
            Expr(:kw, :over,     v) => (over_val = v)
            Expr(:kw, :at,       v) => (at_val   = v)
            _ => throw(TrackingSyntaxError("consensus_sheaf '$name': unrecognised keyword: $(repr(kw))"))
        end
    end
    tmpl_val !== nothing || throw(TrackingSyntaxError("consensus_sheaf '$name': missing 'template'"))
    over_val !== nothing || throw(TrackingSyntaxError("consensus_sheaf '$name': missing 'over'"))
    at_val   !== nothing || throw(TrackingSyntaxError("consensus_sheaf '$name': missing 'at'"))
    tmpl = tmpl_val isa Symbol ? tmpl_val : throw(TrackingSyntaxError("consensus_sheaf '$name': template must be a symbol"))
    agents = _parse_sym_tuple_n(over_val, "consensus_sheaf '$name' over")
    ts     = _parse_at_spec(at_val)
    return ConsensusSheafConstraint(name, tmpl, agents, ts)
end

# ─────────────────────────────────────────────────────────────────────────────
# Boundary: boundary(:agent, a1; at=t[begin], value=x0_a1)
#            boundary(:agent, a; at=t, value=x[a,t])
# ─────────────────────────────────────────────────────────────────────────────

function _parse_boundary(kind_lit, ename::Symbol, kws)
    # kind_lit should be :(:agent) or :(:target) (QuoteNode) or just :agent as Symbol
    kind = _extract_entity_kind(kind_lit)
    at_val    = nothing
    value_val = nothing
    for kw in kws
        @match kw begin
            Expr(:kw, :at,    v) => (at_val    = v)
            Expr(:kw, :value, v) => (value_val = v)
            _ => throw(TrackingSyntaxError("boundary: unrecognised keyword: $(repr(kw))"))
        end
    end
    at_val    !== nothing || throw(TrackingSyntaxError("boundary: missing 'at'"))
    value_val !== nothing || throw(TrackingSyntaxError("boundary: missing 'value'"))
    t   = _parse_time_ref(at_val)
    ref = _parse_boundary_ref(value_val)
    return BoundaryConstraint(kind, ename, t, ref)
end

function _extract_entity_kind(lit)
    @match lit begin
        QuoteNode(:agent)  => :agent
        QuoteNode(:target) => :target
        :agent             => :agent
        :target            => :target
        _ => throw(TrackingSyntaxError("boundary: entity kind must be :agent or :target, got: $(repr(lit))"))
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Time reference / spec helpers
# ─────────────────────────────────────────────────────────────────────────────

function _parse_time_ref(r)::TimeRef
    @match r begin
        Expr(:ref, _::Symbol, :begin) => BeginTime()
        Expr(:ref, _::Symbol, :end)   => EndTime()
        n::Int       => LiteralTime(n)
        s::Symbol    => NamedTime(s)
        _ => throw(TrackingSyntaxError("Cannot parse time reference: $(repr(r))"))
    end
end

function _parse_time_spec_expr(spec)::TimeSpec
    @match spec begin
        Expr(:ref, _::Symbol, :begin) => SingletonTime(BeginTime())
        Expr(:ref, _::Symbol, :end)   => SingletonTime(EndTime())
        n::Int        => SingletonTime(LiteralTime(n))
        s::Symbol     => NamedTimeSet(s)   # named time set reference
        # lo:hi range
        Expr(:call, :(:), lo, hi) => TimeRange(_parse_time_ref(lo), _parse_time_ref(hi))
        # [0, 3, t[end]] explicit list
        Expr(:vect, ts...)   => TimeList([_parse_time_ref(t) for t in ts])
        # (lo:hi) — parenthesised range
        Expr(:block, Expr(:call, :(:), lo, hi)) =>
            TimeRange(_parse_time_ref(lo), _parse_time_ref(hi))
        _ => throw(TrackingSyntaxError("Cannot parse time spec: $(repr(spec))"))
    end
end

function _parse_at_spec(spec)::TimeSpec
    # Same as _parse_time_spec_expr but a singleton symbol resolves to
    # NamedTimeSet for bare names expected to be time sets.
    @match spec begin
        Expr(:ref, _::Symbol, :begin) => SingletonTime(BeginTime())
        Expr(:ref, _::Symbol, :end)   => SingletonTime(EndTime())
        n::Int        => SingletonTime(LiteralTime(n))
        s::Symbol     => NamedTimeSet(s)
        Expr(:call, :(:), lo, hi) => TimeRange(_parse_time_ref(lo), _parse_time_ref(hi))
        Expr(:vect, ts...)        => TimeList([_parse_time_ref(t) for t in ts])
        _ => throw(TrackingSyntaxError("Cannot parse 'at' time spec: $(repr(spec))"))
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Boundary value reference helpers
# ─────────────────────────────────────────────────────────────────────────────

function _parse_boundary_ref(expr)::BoundaryRef
    @match expr begin
        s::Symbol => PlainRef(s)
        Expr(:ref, name::Symbol, a::Symbol, t::Symbol) => IndexedRef(name, a, t)
        Expr(:ref, name::Symbol, t::Symbol)            => IndexedRef(name, :_, t)
        _ => throw(TrackingSyntaxError("Cannot parse boundary value ref: $(repr(expr))"))
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Tuple helpers
# ─────────────────────────────────────────────────────────────────────────────

function _parse_sym_tuple2(expr, ctx::String)::NTuple{2,Symbol}
    @match expr begin
        Expr(:tuple, a::Symbol, b::Symbol) => (a, b)
        _ => throw(TrackingSyntaxError("$ctx: expected 2-tuple of symbols, got: $(repr(expr))"))
    end
end

function _parse_sym_tuple_n(expr, ctx::String)::Vector{Symbol}
    @match expr begin
        Expr(:tuple, args...) => begin
            all(a -> a isa Symbol, args) ||
                throw(TrackingSyntaxError("$ctx: all elements must be symbols"))
            collect(Symbol, args)
        end
        s::Symbol => [s]
        _ => throw(TrackingSyntaxError("$ctx: expected tuple of symbols, got: $(repr(expr))"))
    end
end

end # module TrackingDSLParser
