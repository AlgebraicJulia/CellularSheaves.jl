"""
    TrackingDSL

Name-driven DSL for specifying multi-agent tracking problems as time-expanded
cellular sheaves.

This module provides a symbolic specification language (`@tracking_problem`,
`parse_tracking_program`) that decouples the structural description of a
tracking scenario from the numeric values of matrices and vectors.  Names are
declared first; concrete values are supplied via a context dict at resolve time.
The full pipeline is:

```
parse → validate → resolve(ctx) → lower → TrackingProblem
```

## Quickstart

```julia
using CellularSheaves
using CellularSheaves.ControlSheaves.TrackingDSL

prog = @tracking_problem begin
    space(X) = R^2
    space(U) = R^1
    map_decl(A, X, X)
    map_decl(B, U, X)
    map_decl(R_y, X, X)
    agent(a1; dynamics=(A,B), period=dt)
    agent(a2; dynamics=(A,B), period=dt)
    target(t1)
    horizon(K)
    time(initial = 0)
    time(final = K)
    times(Tall = initial:final)
    consensus(c1; agents=(a1,a2), maps=(R_y,R_y), at=Tall)
    track(tr1; agent=a1, target=t1, maps=(R_y,R_y), at=final)
    boundary(:agent, a1; at=initial, value=x0_a1)
end

ctx = Dict{Symbol,Any}(
    :K    => 10,
    :A    => [1.0 0.0; 0.0 1.0],
    :B    => [0.0; 1.0;;],
    :R_y  => [1.0 0.0; 0.0 0.0],
    :dt   => 0.05,
    :x0_a1 => [0.0, 0.0],
)
result = lower_tracking_program(prog, ctx)
sheaf  = build_time_expanded_tracking_sheaf(result.problem)
```

## Special aliases

`initial` always resolves to `0`.  `final` resolves to the value of the
horizon `k` declared with `horizon K` (after `K` is bound in `ctx`).

## Sub-modules

| Sub-module | Role |
|---|---|
| `TrackingDSLTerm` | AST node definitions and exception types |
| `TrackingDSLParser` | `@tracking_problem` macro and `parse_tracking_program` |
| `TrackingDSLValidator` | `validate_tracking_program` semantic checks |
| `TrackingDSLResolver` | `resolve_tracking_program`, `set_indexed!`, `ResolvedProgram` |
| `TrackingDSLLowering` | `lower_tracking_program`, `LoweredTrackingProblem` |
"""
module TrackingDSL

using Reexport

include("ADT.jl")
include("Parser.jl")
include("Validator.jl")
include("Resolver.jl")
include("Lowering.jl")

@reexport using .TrackingDSLTerm
@reexport using .TrackingDSLParser
@reexport using .TrackingDSLValidator
@reexport using .TrackingDSLResolver
@reexport using .TrackingDSLLowering

end # module TrackingDSL
