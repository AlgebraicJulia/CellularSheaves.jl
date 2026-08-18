# NestedDSL

A composable specification language for the nested layered-control systems of
`NestedSystems`. `@nested_system` **executes** its block rather than quoting it, so
Julia's own loops, conditionals, and functions supply the abstraction, and the language itself
only has to describe the structure of one node. The block returns a
`SystemFragment` — a first-class value that can be merged, spliced, or nested — and
`compile_nested_system` takes a fragment all the way to a `NestedSystemSpec` and a `SheafTower`.

See `TrackingDSL` (the `tracking_dsl.md` page) for this package's other, deliberately
different DSL: that one parses to a symbolic AST and binds numbers later from a context dict.
See also `NetworkSheaves.CellularSheafParser` (the `dsl.md` page) for the unrelated
`@cellular_sheaf` macro.

```@autodocs
Modules = [CellularSheaves.ControlSheaves.NestedDSL,
           CellularSheaves.ControlSheaves.NestedDSL.NestedDSLTerm,
           CellularSheaves.ControlSheaves.NestedDSL.NestedDSLParser,
           CellularSheaves.ControlSheaves.NestedDSL.NestedDSLValidator,
           CellularSheaves.ControlSheaves.NestedDSL.NestedDSLLowering]
```
