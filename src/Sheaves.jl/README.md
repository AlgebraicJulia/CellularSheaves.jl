# Sheaves.jl

Sheaves.jl is a submodule containing high-performance data structures and algorithms for
working with sheaves and their Laplacian matrices.

## Basic Use

Construct a sheaf using the function `sheaf`.

```julia-repl
using CellularSheaves.Sheaves

julia> I = [2, 1, 3, 2, 4, 3, 1, 4];

julia> J = [1, 2, 2, 3, 3, 4, 4, 1];

julia> V = [
           [-2.0;;              ],
           [-1.0  2.0           ],
           [ 1.0  1.0           ],
           [ 1.0;;              ],
           [-1.0;;              ],
           [ 1.0 -1.0           ],
           [ 0.0 -1.0; 1.0  0.0 ],
           [ 0.0; 1.0;;         ],
       ];

julia> S = sheaf(I, J, V)
{4, 8} Sheaf{Float64, Int64}:
 1  1  0  2
 1  2  1  0
 0  1  1  1
 2  0  1  2
```

Compute its coboundary matrix using the function `coboundary`.

```julia-repl
julia> coboundary(S)
5×6 SparseBlockMatrix{Float64, Int64} with 15 stored entries:
  2.0  -1.0   2.0  0.0  0.0   0.0
 -0.0   0.0   0.0  0.0  0.0  -1.0
 -1.0   0.0   0.0  0.0  1.0   0.0
  0.0  -1.0  -1.0  1.0  0.0   0.0
  0.0   0.0   0.0  1.0  1.0  -1.0
```

Compute its Laplacian matrix using the function `laplacian`.

```julia-repl
julia> laplacian(S)
6×6 SparseBlockMatrix{Float64, Int64} with 26 stored entries:
  5.0  -2.0   4.0   0.0  -1.0   0.0
 -2.0   2.0  -1.0  -1.0   0.0   0.0
  4.0  -1.0   5.0  -1.0   0.0   0.0
  0.0  -1.0  -1.0   2.0   1.0  -1.0
 -1.0   0.0   0.0   1.0   2.0  -1.0
  0.0   0.0   0.0  -1.0  -1.0   2.0
```
