"""
    NestedDSLParser

MLStyle-driven surface syntax for the nested-system DSL: the [`@nested_system`](@ref) block and
the declaration macros (`@team`, `@system`, `@link`, `@target`, `@observe`, `@bind`, `@dim`,
`@affine`, `@include`) that go inside it.

## The design in one sentence

`@nested_system begin … end` is a **builder block**, not a quoted program: it binds a hidden
term list and then runs its body as ordinary Julia code, with each declaration macro appending
to that list. Everything that is not a declaration macro is left completely untouched.

That single choice is what buys composability. Because the block *executes*, a `for` loop is a
`for` loop, an `if` is an `if`, and a value computed three lines up is just a variable — the DSL
needs no looping or conditional constructs of its own, and none were added. And because the
block *returns a value* — a [`SystemFragment`](@ref) — a helper function can build one from
plain Julia arguments and hand it back to be `@include`d or made into a `@system`.

```julia
escort(name, m, r, tgt) = @nested_system begin
    @team \$(name) = ring(m; radius=r)
    @observe \$(name) => \$(tgt)
end

wing = @nested_system begin
    for i in 1:3                                  # a real Julia loop
        @include escort(Symbol(:ring, i), 5, 1.0, Symbol(:t, i))
    end
    @link centroid(ring1) => centroid(ring2)
end
```

## Names, and when a name is evaluated

Wherever the grammar expects a **name** — a team or system name, a path component, a target — a
bare identifier is taken *literally*, and `\$(expr)` means "evaluate `expr` now and use the
`Symbol` it produces". Nothing else is accepted in a name position, so a computed name is always
visibly marked rather than silently guessed at.

Every other position is an ordinary Julia expression, evaluated in the caller's scope at the
point the declaration runs: radii, agent counts, matrices, dynamics objects, `redundant_pin(…)`
calls, and so on.

## Statement reference

| Form | Meaning |
|---|---|
| `@dim 4` | stalk dimension (default `4`) |
| `@affine true` | affine restriction maps (default `true`) |
| `@team a = ring(5; radius=1.0, observers=[1])` | leaf team; `ring`/`path`/`star`/`clique`, or `team(kind_expr, n; …)` |
| `@system s begin … end` | refined child system, declared inline |
| `@system s = frag` | refined child system, from an existing fragment |
| `@include frag` | splice a fragment's declarations into *this* node |
| `@link centroid(a) => project(b, 2)` | consensus edge between two direct children |
| `@target t1 t2` | declare targets (global, however deeply nested) |
| `@observe centroid(a) => t1` | `a` observes target `t1` |
| `@bind dynamics=dyn K_lqr=K` | bind this node and everything below it |
| `@bind mid.ringA K_lqr=K_soft` | bind a descendant |
| `@bind ringA[2] initial_position=p` | bind one agent of a leaf team |

Endpoint syntax for `@link`/`@observe` is documented on [`MemberRef`](@ref).
"""
module NestedDSLParser

using MLStyle: @match

using ...NestedSystems: RestrictionSpec, Centroid, RawRestriction, project, centroid,
                        SystemEdge, SystemBinding, AgentBinding, LeafTeam, RefinedSystem,
                        TargetSpec, Observation, NestedSystemSpec
using ....NetworkSheaves.Formations: build_escort_topology
using ..NestedDSLTerm

export @nested_system, @team, @system, @link, @target, @observe, @bind, @dim, @affine, @include

# The block-local variable holding the term list under construction. Every declaration macro
# expands to a `push!` against this name in the caller's scope, which is why declarations work
# unchanged inside `for`/`if`/`let` bodies: they are ordinary code touching an ordinary local.
const BUILDER = :__nested_builder__

_builder() = esc(BUILDER)

_dsl_error(msg) = throw(NestedDSLError(msg))

# ─────────────────────────────────────────────────────────────────────────────
# The block macro
# ─────────────────────────────────────────────────────────────────────────────

"""
    @nested_system begin … end -> SystemFragment

Run a block of declarations, returning the [`SystemFragment`](@ref) they build.

The block is executed as ordinary Julia; only the declaration macros listed in the module
docstring are interpreted, and everything else — loops, conditionals, local variables, function
calls — behaves exactly as it would outside the block, with one exception: the whole block is
wrapped in a `let`, so it opens its own local scope. Variables from the enclosing scope are
still visible and mutable, but a bare assignment inside the block creates a new local rather than
leaking out (the same as writing a `let` block yourself at that spot). A nested `@nested_system`
(which `@system … begin … end` uses internally) shadows the outer builder, so declarations always
attach to the innermost enclosing node.

The returned fragment describes *one node*. It carries no notion of where in a tree it sits, so
it may be merged with another fragment for the same node ([`merge`](@ref)), spliced into one
(`@include`), or made the body of a child (`@system name = frag`).

To turn a fragment into a solvable problem, see `compile_nested_system`.
"""
macro nested_system(block)
    return _fragment_expr(block)
end

function _fragment_expr(block)
    B = _builder()
    return quote
        let $B = $(NestedTerm)[]
            $(esc(block))
            $(SystemFragment)($B)
        end
    end
end

_push(term_expr) = :(push!($(_builder()), $term_expr))

# ─────────────────────────────────────────────────────────────────────────────
# Name / path / endpoint parsing
# ─────────────────────────────────────────────────────────────────────────────

"""
    _name_expr(ex, what) -> Expr

Lower a name-position expression: a bare identifier becomes a literal `Symbol`, `\$(e)` becomes
a runtime `Symbol(e)`, and anything else is a syntax error.

Refusing the general case is deliberate. If `ring1` in a name position could mean "the variable
`ring1`", then every mistyped name would silently become a runtime `UndefVarError` far from the
declaration, and every DSL name would shadow-collide with the caller's locals.
"""
function _name_expr(ex, what::AbstractString)
    return @match ex begin
        s::Symbol => QuoteNode(s)
        Expr(:$, e) => :(Symbol($(esc(e))))
        _ => _dsl_error("$what expects a bare name (e.g. `ringA`) or `\$(expr)` evaluating to a " *
                        "Symbol; got `$ex`")
    end
end

_field_name_expr(ex) = @match ex begin
    QuoteNode(s::Symbol) => QuoteNode(s)
    Expr(:$, e) => :(Symbol($(esc(e))))
    _ => _dsl_error("path component must be a bare name or `\$(expr)`; got `$ex`")
end

"""
    _path_parts(ex) -> Vector

Flatten a dotted path expression (`mid.ringA`, `\$(nm).inner`) into its component name
expressions, outermost first.
"""
function _path_parts(ex)
    return @match ex begin
        Expr(:., inner, field) => Any[_path_parts(inner)..., _field_name_expr(field)]
        _ => Any[_name_expr(ex, "path")]
    end
end

_path_expr(ex) = Expr(:ref, Symbol, _path_parts(ex)...)

"""
    _restriction_expr(ex) -> Expr

Lower a member designator inside `project(path, i)`: `:name` selects a child by name (kept as a
literal), anything else is an ordinary expression evaluated for an integer index.
"""
_restriction_expr(ex) = ex isa QuoteNode ? :($(project)($ex)) : :($(project)($(esc(ex))))

"""
    _endpoint_expr(ex) -> Expr

Lower a `@link`/`@observe` endpoint into a [`member_ref`](@ref) call. See [`MemberRef`](@ref)
for the surface forms and the restriction map each one produces.
"""
function _endpoint_expr(ex)
    return @match ex begin
        Expr(:call, :centroid, p) => :($(member_ref)($(_path_expr(p)), $(centroid)()))
        Expr(:call, :project, p, i) => :($(member_ref)($(_path_expr(p)), $(_restriction_expr(i))))
        Expr(:call, :raw, p, M) => :($(member_ref)($(_path_expr(p)), $(RawRestriction)($(esc(M)))))
        Expr(:call, :via, p, s) => :($(member_ref)($(_path_expr(p)), $(esc(s))))
        Expr(:call, f::Symbol, _...) && if f in (:centroid, :project, :raw, :via) end =>
            _dsl_error("`$f` endpoint has the wrong number of arguments in `$ex` — expected " *
                       "centroid(path), project(path, i), raw(path, M), or via(path, spec)")
        _ => :($(member_ref)($(_path_expr(ex))))
    end
end

_pair_parts(ex, what) = @match ex begin
    Expr(:call, :(=>), lhs, rhs) => (lhs, rhs)
    _ => _dsl_error("$what expects `lhs => rhs`; got `$ex`")
end

# ─────────────────────────────────────────────────────────────────────────────
# Declaration macros
# ─────────────────────────────────────────────────────────────────────────────

"""
    @dim D

Declare the stalk dimension of the whole specification. Defaults to `4` (SE(3) homogeneous
coordinates) if never declared; declaring two different values anywhere is an error, caught by
the validator.
"""
macro dim(D)
    return _push(:($(dim_term)($(esc(D)))))
end

"""
    @affine flag

Declare whether restriction maps are affine (see [`build_escort_topology`](@ref)). Defaults to
`true`; conflicts are caught exactly as for [`@dim`](@ref).
"""
macro affine(flag)
    return _push(:($(affine_term)($(esc(flag)))))
end

"""
    @team name = ring(n; radius=r, observers=[…])
    @team name = path(n, r)
    @team name = team(kind_expr, n; radius=r)

Declare a leaf team of `n` raw agents in a `ring`, `path`, `star`, or `clique` formation.

`n` and `radius` are ordinary Julia expressions. Use the `team(kind_expr, …)` form when the
formation kind is itself computed. `radius` may be given positionally or by keyword and defaults
to `1.0`; `observers` defaults to `[1]`.
"""
macro team(ex)
    nm, call = @match ex begin
        Expr(:(=), nm, call) => (nm, call)
        _ => _dsl_error("@team expects `name = kind(n; …)`; got `$ex`")
    end
    return _push(_team_expr(nm, call))
end

function _team_expr(nm, call)
    head, positional, kwargs = @match call begin
        Expr(:call, f, args...) => (f, _split_call_args(args)...)
        _ => _dsl_error("@team expects a formation call on the right-hand side (e.g. " *
                        "`ring(5; radius=1.0)`); got `$call`")
    end

    kind_expr, rest = if head isa Symbol && head in (:ring, :path, :star, :clique)
        (QuoteNode(head), positional)
    elseif head == :team
        isempty(positional) && _dsl_error("@team … = team(kind, n; …) needs a formation kind")
        (:(Symbol($(esc(positional[1])))), positional[2:end])
    else
        _dsl_error("unknown formation `$head` — use ring, path, star, clique, or " *
                   "team(kind_expr, n; …)")
    end

    isempty(rest) && _dsl_error("@team `$call` needs an agent count")
    length(rest) > 2 && _dsl_error("@team `$call` takes at most an agent count and a radius " *
                                   "positionally; pass anything else by keyword")

    for (k, _) in kwargs
        k in (:radius, :observers) ||
            _dsl_error("@team `$call` got unknown keyword `$k` (expected radius or observers)")
    end

    args = Any[_name_expr(nm, "@team name"), kind_expr, esc(rest[1])]
    length(rest) == 2 && push!(args, esc(rest[2]))

    params = Any[]
    for (k, v) in kwargs
        if k == :radius
            length(rest) == 2 && _dsl_error("@team `$call` gives radius both positionally and " *
                                            "by keyword")
            push!(args, esc(v))
        else
            push!(params, Expr(:kw, k, esc(v)))
        end
    end

    return isempty(params) ? Expr(:call, team_term, args...) :
                             Expr(:call, team_term, Expr(:parameters, params...), args...)
end

"""
    _split_call_args(args) -> (positional, kwargs)

Split a call's argument expressions into positional expressions and `name => value` keyword
pairs, accepting both `f(x; k=v)` (a `:parameters` block) and `f(x, k=v)` (inline `:kw`).
"""
function _split_call_args(args)
    positional = Any[]
    kwargs = Pair{Symbol,Any}[]
    for a in args
        @match a begin
            Expr(:parameters, ps...) => for p in ps
                @match p begin
                    Expr(:kw, k::Symbol, v) => push!(kwargs, k => v)
                    _ => _dsl_error("unsupported keyword argument `$p`")
                end
            end
            Expr(:kw, k::Symbol, v) => push!(kwargs, k => v)
            _ => push!(positional, a)
        end
    end
    return positional, kwargs
end

"""
    @system name begin … end
    @system name = fragment

Declare a refined child system: either inline, with its own nested block of declarations, or
from a [`SystemFragment`](@ref) computed elsewhere.

Both forms nest to any depth — the inline block is itself a [`@nested_system`](@ref) block, and
a fragment may contain further `@system` declarations of its own.
"""
macro system(args...)
    if length(args) == 2
        nm, block = args
        return _push(:($(subsystem_term)($(_name_expr(nm, "@system name")),
                                         $(_fragment_expr(block)))))
    elseif length(args) == 1
        nm, body = @match args[1] begin
            Expr(:(=), nm, body) => (nm, body)
            _ => _dsl_error("@system expects `name begin … end` or `name = fragment`; got " *
                            "`$(args[1])`")
        end
        return _push(:($(subsystem_term)($(_name_expr(nm, "@system name")), $(esc(body)))))
    else
        _dsl_error("@system expects `name begin … end` or `name = fragment`")
    end
end

"""
    @include fragment…

Splice each fragment's declarations into the node currently being built, as though they had been
written out at this point.

This is the DSL's composition primitive. Combined with a Julia `for` loop it replaces any need
for iteration in the language itself, and combined with a helper function that returns a
fragment it replaces any need for DSL-level abstraction:

```julia
@nested_system begin
    for i in 1:n
        @include escort_ring(Symbol(:ring, i), m[i], radius(i))
    end
end
```

Paths inside a spliced fragment are interpreted relative to *this* node, since that is where the
declarations now live.
"""
macro include(frags...)
    isempty(frags) && _dsl_error("@include expects at least one fragment")
    B = _builder()
    pushes = [:(append!($B, $(terms)($(esc(f))))) for f in frags]
    return Expr(:block, pushes...)
end

"""
    @target name…

Declare one or more targets — uncontrolled vertices whose values are supplied as boundary
conditions at solve time.

Targets are **global**: a target declared inside a deeply nested fragment is still a top-level
target of the finished specification, which is what lets a reusable fragment carry the target it
tracks along with it. Declaring the same target name twice is an error.
"""
macro target(names...)
    isempty(names) && _dsl_error("@target expects at least one name")
    return Expr(:block, [_push(:($(target_term)($(_name_expr(n, "@target"))))) for n in names]...)
end

"""
    @link src => dst

Declare a consensus edge between two **direct children** of the node being built, with each
endpoint presenting whatever its designator says (see [`MemberRef`](@ref) for the endpoint
forms).

```julia
@link ringA => ringB                        # each presents its own first member
@link centroid(ringA) => centroid(ringB)    # each presents its members' average
@link pod => project(hub, 3)                # hub presents its third member
```

Edges connect siblings only — a deeper path is a validation error — because that is exactly what
[`SystemEdge`](@ref) expresses. Coupling two systems in different branches is done by giving
them a common parent and linking there.
"""
macro link(ex)
    lhs, rhs = _pair_parts(ex, "@link")
    return _push(:($(link_term)($(_endpoint_expr(lhs)), $(_endpoint_expr(rhs)))))
end

"""
    @observe system => target

Declare that `system` observes the target named `target`, presenting whatever its designator
says (see [`MemberRef`](@ref)).

Unlike [`@link`](@ref), the system may be at any depth below the declaring node, and the
incidence is unrestricted: one system may observe several targets and one target may be observed
by several systems.

```julia
@observe ringA => t1                                   # ringA's first agent tracks t1
@observe centroid(mid.ringA) => t1                     # its centroid tracks t1
@observe via(ringA, redundant_pin(5, 4, k)) => t1      # any RestrictionSpec you like
```
"""
macro observe(ex)
    lhs, rhs = _pair_parts(ex, "@observe")
    return _push(:($(observe_term)($(_endpoint_expr(lhs)), $(_name_expr(rhs, "@observe target")))))
end

"""
    @bind dynamics=dyn K_lqr=K
    @bind path.to.system K_lqr=K_soft
    @bind team[2] initial_position=p

Bind dynamics, LQR gain, and/or initial position, feeding the most-specific-wins cascade of
[`SystemBinding`](@ref).

With no path, the binding applies to the node being built and everything beneath it — the usual
way to give a whole tree its default dynamics. With a path, it applies to that descendant and
its subtree. With a trailing `[i]`, it applies to agent `i` of a leaf team alone.

Fields are independent: binding only `initial_position` on one agent leaves it inheriting its
dynamics and gain from above.
"""
macro bind(args...)
    isempty(args) && _dsl_error("@bind expects at least one of dynamics=, K_lqr=, " *
                                "initial_position=")
    first_is_kw = args[1] isa Expr && args[1].head in (:(=), :kw)
    path_ex = first_is_kw ? nothing : args[1]
    kw_args = first_is_kw ? args : args[2:end]

    isempty(kw_args) && _dsl_error("@bind at `$(path_ex)` declares nothing — give at least one " *
                                   "of dynamics=, K_lqr=, initial_position=")

    path_expr, agent_expr = if path_ex === nothing
        (Expr(:ref, Symbol), :nothing)
    else
        @match path_ex begin
            Expr(:ref, p, i) => (_path_expr(p), esc(i))
            _ => (_path_expr(path_ex), :nothing)
        end
    end

    params = Any[]
    for a in kw_args
        k, v = @match a begin
            Expr(:(=), k::Symbol, v) => (k, v)
            Expr(:kw, k::Symbol, v) => (k, v)
            _ => _dsl_error("@bind expects `key=value` arguments; got `$a`")
        end
        k in (:dynamics, :K_lqr, :initial_position) ||
            _dsl_error("@bind got unknown key `$k` (expected dynamics, K_lqr, or " *
                       "initial_position)")
        push!(params, Expr(:kw, k, esc(v)))
    end

    return _push(Expr(:call, bind_term, Expr(:parameters, params...), path_expr, agent_expr))
end

end # module NestedDSLParser
