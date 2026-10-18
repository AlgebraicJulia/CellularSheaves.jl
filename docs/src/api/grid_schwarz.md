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
