# TrackingDSL

Name-driven DSL for specifying multi-agent tracking problems as time-expanded cellular
sheaves. `@tracking_problem` declares structure first; concrete numeric values are
supplied later via a context dict at resolve time. See `NetworkSheaves.CellularSheafParser`
(the `dsl.md` page) for the unrelated `@cellular_sheaf` macro.

```@autodocs
Modules = [CellularSheaves.ControlSheaves.TrackingDSL,
           CellularSheaves.ControlSheaves.TrackingDSL.TrackingDSLTerm,
           CellularSheaves.ControlSheaves.TrackingDSL.TrackingDSLParser,
           CellularSheaves.ControlSheaves.TrackingDSL.TrackingDSLValidator,
           CellularSheaves.ControlSheaves.TrackingDSL.TrackingDSLResolver,
           CellularSheaves.ControlSheaves.TrackingDSL.TrackingDSLLowering]
```
