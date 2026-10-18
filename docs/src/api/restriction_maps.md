# Restriction Maps

A restriction map ``\mathcal F_{v \trianglelefteq e}`` need not be stored as a
dense matrix. [`AbstractRestrictionMap`](@ref) is a LinearMaps-style interface:
a map provides `size`, `mul!(y, R, x)` and `mul!(y, R', x)`, and gets
products with vectors and matrices, `Matrix(R)`, `sparse(R)` and a
`LinearOperator` from those.

| Map | Storage | Use |
|---|---|---|
| [`DenseRestriction`](@ref) | dense matrix | small stalks |
| [`SparseRestriction`](@ref) | sparse matrix | large, sparse maps |
| [`SelectionRestriction`](@ref) | index list | covers, overlaps, ghost layers |
| [`FunctionRestriction`](@ref) | two closures | maps known only by their action |

A sheaf chooses its storage through its second type parameter.
`EuclideanSheaf{T}(stalks)` keeps dense matrices (a
[`DenseEuclideanSheaf`](@ref)). `EuclideanSheaf{T,M}(stalks)` stores maps of
type `M`, for example `SelectionRestriction{T}`, or `AbstractRestrictionMap{T}`
to mix kinds. For such sheaves [`coboundary_map`](@ref) assembles a sparse
matrix, and [`coboundary_operator`](@ref) applies ``\delta`` and
``\delta^\mathsf{T}`` without assembling anything.

```@autodocs
Modules = [CellularSheaves.NetworkSheaves.RestrictionMaps]
```
