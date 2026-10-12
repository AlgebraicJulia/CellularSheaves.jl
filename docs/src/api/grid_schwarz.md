# Grid Schwarz and Multigrid

Matrix-free solvers on structured grids: implicit axis-stencil operators,
boxes of a distributed grid with their halo exchanges (the overlap sheaf of the
box cover), red–black Gauss–Seidel, BiCGStab and an aggregation coarse space,
all as KernelAbstractions kernels (CPU threads and GPUs).

```@autodocs
Modules = [CellularSheaves.NetworkSheaves.GridSchwarz]
```

## Grid hierarchies

Cell-centred coarsening, with the maps between levels as the pullback and
pushforward of the constant sheaf along the aggregation homomorphism, and the
Galerkin coarse operator as their composite with the fine operator.

```@autodocs
Modules = [CellularSheaves.NetworkSheaves.GridMultigrid]
```

## Schwarz and multigrid commute

Restriction to a box cover commutes with the pushforward hierarchy when every box
is a union of whole blocks (base change), so the Galerkin coarse operator of a
local problem is the local problem of the Galerkin coarse operator, and additive
"multigrid in each box, then Schwarz" equals "Schwarz on every level of the
pushforward hierarchy". A dense reference implementation of both orders:

```@autodocs
Modules = [CellularSheaves.NetworkSheaves.MultilevelSchwarz]
```

## Coarse grids in the min-plus semiring

Deterministic optimal control on a state graph is linear in min-plus. The
pushforward along an aggregation is the quotient graph; its solution is a lower
bound on the cost-to-go and an exact A* heuristic, instead of a correction.

```@autodocs
Modules = [CellularSheaves.NetworkSheaves.TropicalAbstraction]
```

## Min-plus sheaves on a simplicial complex

Several agents that share a space: the clique complex of their interaction graph,
with the pushforward of the joint state graph of every simplex as its stalk
(symmetric in the agents), restriction maps that forget agents, and the
simplicial cover that glues the stalk values into a bound on the joint cost.

```@autodocs
Modules = [CellularSheaves.NetworkSheaves.TropicalCovers]
```
