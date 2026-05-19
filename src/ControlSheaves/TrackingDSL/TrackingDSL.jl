"""
    TrackingDSL

Name-driven DSL for specifying multi-agent tracking problems as time-expanded
cellular sheaves.

This module provides a symbolic specification language (`@tracking_problem`,
`parse_tracking_program`) that decouples the structural description of a
tracking scenario from the numeric values of matrices and vectors.  Names are
declared first; concrete values are supplied later via `bind` statements or
the `bind!` helper.  The full pipeline is:

```
parse → validate → resolve → lower → TrackingProblem
```

## Quickstart

```julia
using CellularSheaves
using CellularSheaves.ControlSheaves.TrackingDSL

prog = @tracking_problem begin
    space X = R^2
    space U = R^1
    map A : X -> X
    map B : U -> X
    map R_y : X -> X
    agent a1 dynamics (A, B) period dt
    agent a2 dynamics (A, B) period dt
    target t1
    horizon K
    time initial = 0
    time final   = K
    times Tall   = initial:final
    consensus c1 between (a1, a2) using (R_y, R_y) at Tall
    track     tr1 agent a1 target t1 using (R_y, R_y) at final
    boundary  agent a1 at initial = x0_a1
    bind K    => 10
    bind A    => [1.0 0.0; 0.0 1.0]
    bind B    => [0.0; 1.0;;]
    bind R_y  => [1.0 0.0; 0.0 0.0]
    bind dt   => 0.05
    bind x0_a1 => [0.0, 0.0, 0.0]
end

result = prog |> validate_tracking_program |> resolve_tracking_program |> lower_tracking_program
sheaf  = build_time_expanded_tracking_sheaf(result.problem)
```

## Special aliases

`initial` always resolves to `0`.  `final` resolves to the value of the
horizon `k` declared with `horizon K` (after `K` is bound).

## Sub-modules

| Sub-module | Role |
|---|---|
| `TrackingDSLTerm` | AST node definitions and exception types |
| `TrackingDSLParser` | `@tracking_problem` macro and `parse_tracking_program` |
| `TrackingDSLValidator` | `validate_tracking_program` semantic checks |
| `TrackingDSLResolver` | `resolve_tracking_program`, `bind!`, `ResolvedProgram` |
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
