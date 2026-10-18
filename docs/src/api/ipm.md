# IPM

The interior-point method for cone programs (positive, second-order, exponential,
and semidefinite cones) is provided by
[Mumblebee.jl](https://github.com/samuelsonric/mumblebee.jl), a dependency of
CellularSheaves. It solves

```math
\min_p \tfrac12 p^\top Q p - f^\top p \quad \text{s.t.} \quad Bp = g,\; p \in K,
```

supports warm starts (`init` with `p0`, `d0`, `y0`, and `Mumblebee.IPM.reinit!`
for new data `f`, `g` on a fixed structure), and differentiates the solution
through the KKT system (`frule!`, `frule2!`, `rrule!`).

`CellularSheaves.IPM` is an alias for `Mumblebee.IPM`, so

```julia
using CellularSheaves.IPM
```

brings `IPMProblem`, `IPMSettings`, `solve`, the cone types, and the derivative
rules into scope. Block-sparse problem data must use Mumblebee's own
`Mumblebee.BlockSparseArrays`, which is distinct from
`CellularSheaves.BlockSparseArrays`; the sparse-matrix constructor
`IPMProblem(Q, B, f, g, μ, K, s)` avoids the distinction. See the Mumblebee
documentation and examples for the API reference.
