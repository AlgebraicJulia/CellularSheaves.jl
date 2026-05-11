# Issue 001: Direct-Block Sheaf Laplacian Without Coboundary Matrix Materialization

**Origin:** Review comment on PR #15  
**Status:** Open

---

## Mathematical Background

Let `F` be a Euclidean sheaf on a graph `G = (V, E)`.  Each vertex `v ∈ V` has a stalk
`F(v) = ℝ^{d_v}` and each oriented edge `e = (u, v) ∈ E` has a stalk `F(e) = ℝ^{d_e}`.
The restriction maps are matrices

```
ρ_{u→e} : ℝ^{d_u} → ℝ^{d_e}
ρ_{v→e} : ℝ^{d_v} → ℝ^{d_e}
```

### Coboundary map

The coboundary map `d : C⁰(F) → C¹(F)` is the block-sparse matrix whose `(e, u)` block
is `ρ_{u→e}` and whose `(e, v)` block is `-ρ_{v→e}`.  Assembling the full `d` (an
`m × n` sparse block matrix where `m = Σ d_e` and `n = Σ d_v`) is an `O(m n)` operation
in the worst case and requires materializing a large intermediate matrix.

### Sheaf Laplacian — direct block formula

The sheaf Laplacian `L = dᵀ d` decomposes into vertex-indexed blocks.  For vertices
`u, v ∈ V` let `E(u)` be the set of edges incident to `u` and `E(u, v)` the set of
edges connecting `u` and `v` (at most one for simple graphs):

```
L_{u,u} =  Σ_{e ∈ E(u)}    ρ_{u→e}ᵀ ρ_{u→e}       (diagonal block)

L_{u,v} = -Σ_{e ∈ E(u,v)}  ρ_{u→e}ᵀ ρ_{v→e}       (off-diagonal block, u ≠ v)
```

These formulas allow the Laplacian (or any submatrix of it) to be assembled directly
from the restriction maps **without ever forming `d`**.  Each diagonal block costs
`O(d_e · d_u²)` and each off-diagonal block costs `O(d_e · d_u · d_v)`.

### Restricted blocks for boundary-value problems

`harmonic_extension` solves the BVP

```
L_II x_I = -L_IB x_B
```

where `I` and `B` are the index sets of interior and boundary vertices, respectively.
Today the implementation:

1. Builds the full coboundary map `d` (`O(Σ d_e · Σ d_v)` memory).
2. Computes `L = dᵀ d` (`O((Σ d_v)²)` memory for a dense result).
3. Extracts `L_II` and `L_IB` as submatrix copies.

Using the direct-block formula, `L_II` and `L_IB` can be assembled in a single pass over
the edges with at least one interior endpoint:

- For every edge `e = (u, v)`:
  - If `u ∈ I, v ∈ I`:  accumulate into `L_II[u,u]`, `L_II[v,v]`, `L_II[u,v]`.
  - If `u ∈ I, v ∈ B`:  accumulate into `L_II[u,u]`, `L_IB[u,v]`.
  - If `u ∈ B, v ∈ I`:  accumulate into `L_II[v,v]`, `L_IB[v,u]`.
  - If `u ∈ B, v ∈ B`:  skip (neither block is needed).

This is `O(|E_I| · d_max²)` in time and `O(n_I² + n_I · n_B)` in space, where `E_I`
is the set of edges with at least one interior endpoint, which is strictly cheaper than
materializing the full `(Σ d_v) × (Σ d_v)` Laplacian when the interior is a small
fraction of the graph.

---

## Codebase Orientation

| Symbol | Location | Notes |
|--------|----------|-------|
| `coboundary_map` | `src/network_sheaves/EuclideanSheaves.jl` | Builds the full block-sparse `d` |
| `sheaf_laplacian_matrix` | `src/network_sheaves/EuclideanSheaves.jl` | Returns `dᵀ d` |
| `sheaf_laplacian` | `src/network_sheaves/EuclideanSheaves.jl` | Returns a closure `x ↦ dᵀ(d x)` |
| `nullspace_ldlt(s)` | `src/network_sheaves/EuclideanSheaves.jl` | Calls `coboundary_map` → `dᵀ d` |
| `nearest_global_section_*` | `src/network_sheaves/EuclideanSheaves.jl` | Each builds `d` and forms `d dᵀ` or `dᵀ d` |
| `harmonic_extension` | `src/network_sheaves/EuclideanSheaves.jl` | Calls `sheaf_laplacian_matrix`, then slices `L[I,I]`, `L[I,B]` |
| `BlockSparseArrays` | `src/network_sheaves/BlockSparseArrays.jl` | Provides `blocksparse` for assembling `d` |

All callers that do `d = coboundary_map(s); L = d' * d` (or equivalently `B' * B`) are
candidates for the replacement.  Grep for `coboundary_map` and `sheaf_laplacian_matrix`
to find them.

---

## Requested Implementation

### 1. New function: `sheaf_laplacian_blocks`

```julia
"""
    sheaf_laplacian_blocks(s::EuclideanSheaf)
        -> (diag::Vector{Matrix}, offdiag::Dict{Tuple{Int,Int}, Matrix})

Compute the block decomposition of the sheaf Laplacian `L = dᵀ d` directly from
the restriction maps of `s`, **without forming the coboundary matrix `d`**.

Returns:
- `diag[v]`          — the `d_v × d_v` diagonal block `L_{v,v}`.
- `offdiag[(u,v)]`   — the `d_u × d_v` off-diagonal block `L_{u,v}` for each edge
                       `(u,v)` with `u < v`.  The block `L_{v,u} = L_{u,v}ᵀ`.
"""
```

### 2. New function: `restricted_laplacian_blocks`

```julia
"""
    restricted_laplacian_blocks(s::EuclideanSheaf,
                                interior::AbstractVector{Int},
                                boundary::AbstractVector{Int})
        -> (L_II::SparseMatrixCSC, L_IB::SparseMatrixCSC)

Compute the interior–interior and interior–boundary blocks of the sheaf Laplacian
**without assembling the full Laplacian**.

`L_II` is the `(Σ_{v ∈ I} d_v) × (Σ_{v ∈ I} d_v)` symmetric positive-semidefinite block.
`L_IB` is the `(Σ_{v ∈ I} d_v) × (Σ_{v ∈ B} d_v)` coupling block.
"""
```

### 3. Refactor `harmonic_extension`

Replace the lines

```julia
L     = sheaf_laplacian_matrix(s)
L_II  = L[I_idx, I_idx]
# …
L[I_idx, B_idx] * x_B
```

with a call to `restricted_laplacian_blocks` so the full Laplacian is never materialized.

### 4. Refactor remaining `d' * d` sites

Update every function that calls `coboundary_map` + `d' * d` or `sheaf_laplacian_matrix`:

- `nullspace_ldlt(s::EuclideanSheaf)`
- `nearest_global_section_ldl`
- `nearest_global_section_pinv`
- `energy_function(s::EuclideanSheaf)`

For functions that need the full Laplacian (e.g. `nullspace_ldlt`, `energy_function`),
provide a `sheaf_laplacian_matrix_direct` that assembles the sparse matrix from blocks
without the intermediate `d`.

### 5. Benchmarks

Add a file `benchmarks/sheaf_laplacian_benchmarks.jl` (using
[`BenchmarkTools.jl`](https://github.com/JuliaCI/BenchmarkTools.jl)) that compares:

- `sheaf_laplacian_matrix(s)` (current, via `dᵀ d`) vs. direct-block assembly.
- `harmonic_extension` (current) vs. refactored version with `restricted_laplacian_blocks`.
- Scale: graphs with 10, 100, 1 000 vertices, stalk dimension 3, and 10 % boundary
  vertices.

The benchmarks should demonstrate sub-linear scaling improvement of the restricted
approach as the fraction of interior vertices shrinks.

---

## Tests to Write

1. **Block correctness** — for a small hand-constructed sheaf, verify that every block
   returned by `sheaf_laplacian_blocks` matches the corresponding submatrix of the full
   `sheaf_laplacian_matrix`.

2. **Restricted blocks correctness** — for the same sheaf, verify that `L_II` and `L_IB`
   returned by `restricted_laplacian_blocks` match the exact submatrices of the full
   Laplacian.

3. **Regression on `harmonic_extension`** — verify that the refactored implementation
   returns the same `(x_p, null_basis)` as the current one on all existing test cases.

4. **Empty interior / all-interior edge cases** — ensure `restricted_laplacian_blocks`
   handles degenerate cases gracefully.

---

## Verification Checklist

- [ ] `sheaf_laplacian_blocks` returns blocks consistent with `sheaf_laplacian_matrix` for several test sheaves.
- [ ] `restricted_laplacian_blocks` returns blocks consistent with slicing the full Laplacian.
- [ ] `harmonic_extension` produces the same result before and after refactoring.
- [ ] No existing tests regress.
- [ ] Benchmark script runs without error and shows improvement for large interior sets.
- [ ] All new functions are exported from `EuclideanSheaves` and documented with docstrings.
- [ ] `coboundary_map` + `d' * d` pattern is eliminated from all non-test call sites.

---

## Out of Scope

- Changing the `coboundary_map` API or removing it (it remains useful for testing and for
  other algorithms that genuinely need `d`).
- Support for non-Euclidean (e.g. sheaves valued in arbitrary categories) stalks.
- Parallelising the block assembly loop.
- GPU/distributed extensions.
