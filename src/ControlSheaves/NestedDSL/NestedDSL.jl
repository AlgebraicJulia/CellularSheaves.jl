"""
    NestedDSL

A composable specification language for the nested layered-control systems of
[`NestedSystems`](@ref) — teams, refined subsystems, targets, per-edge restriction maps, and the
dynamics cascade — written as ordinary Julia that *runs*.

## Why a builder block rather than a quoted program

`TrackingDSL`, the other DSL in this package, parses a block into a symbolic AST and resolves
numeric values later from a context dict. This one deliberately does the opposite:
[`@nested_system`](@ref) **executes** its block, and each declaration appends a term to the
fragment under construction.

The reason is the shape of the problems being specified. A nested system's interesting structure
is combinatorial — `n` escort rings, `n` support pods, a pin on every other agent around each
ring — and all of it is a function of parameters known at build time. A declarative language
would have to grow loops, conditionals, and functions to express that. Executing the block
instead means Julia's own `for`, `if`, and functions do that work, and the DSL is left to
express only what it is actually good at: the structure of one node.

```julia
spec = @nested_system begin
    @dim 4
    for i in 1:n                                     # Julia's loop, not the DSL's
        @team \$(Symbol(:ring, i)) = ring(m[i]; radius=r(i))
        @target \$(Symbol(:t, i))
        for k in 1:2:m[i]
            @observe via(\$(Symbol(:ring, i)), redundant_pin(m[i], 4, k)) => \$(Symbol(:t, i))
        end
    end
    @bind dynamics=QuadrotorDynamics() K_lqr=K
end
c = compile_nested_system(spec)
```

## Composition

[`@nested_system`](@ref) returns a [`SystemFragment`](@ref) — a first-class value describing the
contents of one node, with no commitment to where in a tree that node sits. Fragments compose
three ways:

| | |
|---|---|
| `merge(f, g)` | two fragments for the *same* node |
| `@include g` | splice `g`'s declarations into the node being built |
| `@system name = g` | make `g` the body of a new child |

Paths inside a fragment are relative, and lowering rewrites them against wherever the fragment
lands — so a helper function returning "an escort ring that tracks target `t`" can be dropped in
at the root, or three levels down, unchanged:

```julia
escort(name, m, r, tgt) = @nested_system begin
    @team \$(name) = ring(m; radius=r)
    @target \$(tgt)
    @observe \$(name) => \$(tgt)
end

fleet = @nested_system begin
    @system wingA begin
        @system flight1 begin
            @include escort(:alpha, 5, 1.0, :t1)
            @include escort(:beta, 5, 1.0, :t2)
            @link centroid(alpha) => centroid(beta)
        end
        @include escort(:gamma, 4, 0.8, :t3)
        @link centroid(flight1) => centroid(gamma)
    end
end
```

## Pipeline

```
@nested_system  →  SystemFragment  →  validate  →  lower  →  NestedSystemSpec + SheafTower
```

| Sub-module | Role |
|---|---|
| `NestedDSLTerm` | term types, [`SystemFragment`](@ref), [`NestedDSLError`](@ref) |
| `NestedDSLParser` | [`@nested_system`](@ref) and the declaration macros |
| `NestedDSLValidator` | symbolic checks: names, paths, arities, conflicts |
| `NestedDSLLowering` | [`compile_nested_system`](@ref) and the name→index tables |

## Coverage of the underlying API

Every feature of the `NestedSystems` type API has a surface form here: `@team` for
[`LeafTeam`](@ref) (all four formation kinds, `radius`, `observers`), `@system` for
[`RefinedSystem`](@ref) at unbounded depth, `@link` for [`SystemEdge`](@ref) with independent
per-endpoint restriction maps, `@target`/`@observe` for [`TargetSpec`](@ref)/[`Observation`](@ref)
with arbitrary many-to-many incidence, `@dim`/`@affine` for the remaining
[`NestedSystemSpec`](@ref) fields, and `@bind` for the whole
[`SystemBinding`](@ref)/[`AgentBinding`](@ref) cascade. The endpoint forms `project`, `centroid`,
`raw`, and `via` cover every [`RestrictionSpec`](@ref), with `via` accepting an arbitrary one —
including helpers like [`redundant_pin`](@ref) — so the DSL is never a bottleneck on what can be
expressed.
"""
module NestedDSL

using Reexport

import ..NestedSystems
using ..NestedSystems: LeafTeam, RefinedSystem, SystemEdge, TargetSpec, Observation,
                       NestedSystemSpec, SystemBinding, AgentBinding, RestrictionSpec,
                       redundant_pin

include("ADT.jl")
include("Parser.jl")
include("Validator.jl")
include("Lowering.jl")

@reexport using .NestedDSLTerm
@reexport using .NestedDSLParser
@reexport using .NestedDSLValidator
@reexport using .NestedDSLLowering

end # module NestedDSL
