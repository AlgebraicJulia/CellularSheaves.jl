# 001 — Discrete-time LTI control model on 0-cochains of a cellular sheaf

## Mathematical Background

Let `F` be a cellular sheaf on a graph `G = (V, E)`, with vertex stalks `F(v)`.
Its 0-cochain space is

```math
C^0(F) = \bigoplus_{v \in V} F(v).
```

This issue requests a **discrete-time linear time-invariant control model** whose
state lives on `C^0(F)`, not just on the global sections of `F`. The sheaf is
therefore used to organize the vertex-wise state decomposition and to retain the
existing network-sheaf geometry, but **sheaf consistency is not imposed as a hard
constraint** on the dynamics.

The requested system is

```math
x[k+1] = A x[k] + B u[k], \qquad y[k] = C x[k] + D u[k],
```

with

- `x[k] \in C^0(F)`,
- `u[k] \in \bigoplus_{v \in V} U_v`,
- `y[k] \in \bigoplus_{v \in V} Y_v`.

The matrices `B`, `C`, and `D` are assembled from **per-vertex** maps:

```math
B = \operatorname{diag}(B_v), \qquad
C = \operatorname{diag}(C_v), \qquad
D = \operatorname{diag}(D_v),
```

where

```math
B_v : U_v \to F(v), \qquad
C_v : F(v) \to Y_v, \qquad
D_v : U_v \to Y_v.
```

The state matrix `A : C^0(F) \to C^0(F)` is **user-supplied**. In particular,
this issue does **not** derive `A` from the coboundary map `d_F` or the sheaf
Laplacian `L_F = d_F^\top d_F`.

The intended nomenclature should match the codebase's existing use of
**cellular sheaf**, **network sheaf**, **0-cochain**, and **coboundary map**.
This keeps the feature aligned with the terminology used by Curry, Ghrist, and
the existing sheaf-Laplacian work in this repository.

## Codebase Orientation

| File | Why it matters |
| --- | --- |
| `src/network_sheaves/SheafInterface.jl:1-69` | Public abstract interface; add new exported operations here before implementing concrete methods elsewhere. |
| `src/network_sheaves/EuclideanSheaves.jl:47-179` | Defines `EuclideanSheaf`, vertex/edge stalk accessors, and `coboundary_map`; this is the concrete sheaf type the first implementation should target. |
| `src/network_sheaves/EuclideanSheaves.jl:252-303` | Shows current public API style for numerical routines (`nearest_global_section`, `nullspace_ldlt`). |
| `src/network_sheaves/NetworkSheaves.jl:1-25` | Module aggregator; new include and re-export must be wired here. |
| `test/runtests.jl:1-12` | One test file per source module; new tests must be included here. |
| `src/network_sheaves/Pushforwards.jl:39-227` | Good reference for docstring style, dimension checks, and returning sparse linear-algebra objects tied to sheaf structure. |

## Requested Implementation

Add a new source file `src/network_sheaves/LTIControl.jl` that exports a small,
self-contained API for a **discrete-time sheaf LTI system on 0-cochains**.

### Exact API

```julia
"""
    DiscreteSheafLTI(sheaf::EuclideanSheaf{T},
                     A::AbstractMatrix{T},
                     Bv::AbstractVector{<:AbstractMatrix{T}},
                     Cv::AbstractVector{<:AbstractMatrix{T}};
                     Dv::Union{Nothing,AbstractVector{<:AbstractMatrix{T}}}=nothing) where T

Construct a discrete-time LTI control system whose state space is the 0-cochain
space `C^0(sheaf) = ⨁_v sheaf(v)`.

The update and output equations are

`x[k+1] = A * x[k] + B * u[k]`

`y[k] = C * x[k] + D * u[k]`

where `B`, `C`, and `D` are assembled block-diagonally from the per-vertex maps
`Bv`, `Cv`, and `Dv`.

For each vertex `v`,
- `Bv[v]` must have size `(dim(F(v)), input_dim_v)`;
- `Cv[v]` must have size `(output_dim_v, dim(F(v)))`;
- `Dv[v]` must have size `(output_dim_v, input_dim_v)`.

If `Dv === nothing`, use the zero direct-feedthrough maps of the appropriate
sizes at each vertex.

The state matrix `A` is user-supplied and must have size
`(sum(vertex_stalks(sheaf)), sum(vertex_stalks(sheaf)))`.
"""
struct DiscreteSheafLTI{T}
    sheaf::EuclideanSheaf{T}
    A::SparseMatrixCSC{T,Int}
    B::SparseMatrixCSC{T,Int}
    C::SparseMatrixCSC{T,Int}
    D::SparseMatrixCSC{T,Int}
    input_dims::Vector{Int}
    output_dims::Vector{Int}
end

"""
    step(sys::DiscreteSheafLTI, x, u)

Return `(x_next, y)` for one time step of the discrete-time sheaf control
system.
"""
function step(sys::DiscreteSheafLTI, x, u)
```

### Implementation sketch

1. Add abstract declarations to `SheafInterface.jl`:
   - `state_dim`
   - `input_dim`
   - `output_dim`
   - `step`

   These should error by default, following the existing interface pattern.

2. In `LTIControl.jl`, implement a helper for block-diagonal assembly from
   vertex-local maps. Keep it sparse. Do **not** add a new dependency for this.
   A private helper such as

   ```julia
   _blockdiag(mats::AbstractVector{<:AbstractMatrix{T}}) where T
   ```

   is sufficient.

3. Validate inputs at the public constructor boundary using `@argcheck`.
   Required checks:
   - `length(Bv) == length(Cv) == nv(underlying_graph(sheaf))`
   - if `Dv !== nothing`, then `length(Dv) == nv(underlying_graph(sheaf))`
   - `size(A) == (n, n)` where `n = sum(vertex_stalks(sheaf))`
   - `size(Bv[v], 1) == vertex_stalks(sheaf)[v]`
   - `size(Cv[v], 2) == vertex_stalks(sheaf)[v]`
   - `size(Dv[v]) == (size(Cv[v], 1), size(Bv[v], 2))`

4. Compute

   ```math
   n = \sum_v \dim F(v), \qquad
   m = \sum_v \dim U_v, \qquad
   p = \sum_v \dim Y_v.
   ```

   Store `input_dims[v] = size(Bv[v], 2)` and `output_dims[v] = size(Cv[v], 1)`.

5. Assemble

   ```math
   B = \operatorname{diag}(B_v), \qquad
   C = \operatorname{diag}(C_v), \qquad
   D = \operatorname{diag}(D_v)
   ```

   as sparse matrices, and convert `A` to `sparse(A)` once in the constructor.

6. Implement:

   ```julia
   state_dim(sys) = size(sys.A, 1)
   input_dim(sys) = size(sys.B, 2)
   output_dim(sys) = size(sys.C, 1)
   ```

7. Implement

   ```julia
   step(sys, x, u) = (sys.A * x + sys.B * u, sys.C * x + sys.D * u)
   ```

   with dimension checks:
   - `length(x) == state_dim(sys)`
   - `length(u) == input_dim(sys)`

8. Wire the new module into `src/network_sheaves/NetworkSheaves.jl` and export
   the new public names there.

### Notes on representation

- This issue should target `EuclideanSheaf` first. Do not generalize to all
  possible `AbstractNetworkSheaf` implementations unless doing so is free.
- Keep the stored system matrices sparse even if a caller passes dense local
  blocks.
- Do not attempt to enforce `x[k]` being a global section. That is a different
  modeling choice.

## Tests to Write

Create `test/network_sheaves/LTIControl.jl` and include it from `test/runtests.jl`.

Use a concrete sheaf on a path graph `1 - 2 - 3` with vertex stalk dimensions
`[2, 1, 2]`. The edge restriction maps can be any small valid matrices; the
exact edge data does not affect the current state-space construction, but using
an actual sheaf instance ensures the decomposition comes from the package's
sheaf type.

Recommended test setup:

- `A = sparse(I, 5, 5)` or a small nontrivial sparse matrix on 5 state
  coordinates.
- `Bv = [ [1.0; 0.0], reshape([2.0], 1, 1), [0.0; 1.0] ]`
- `Cv = [ reshape([1.0, -1.0], 1, 2), reshape([3.0], 1, 1), reshape([2.0, 4.0], 1, 2) ]`
- `Dv = [ reshape([0.5], 1, 1), reshape([0.0], 1, 1), reshape([-1.0], 1, 1) ]`

The individual `@test` statements should include:

1. Constructor stores the expected state, input, and output dimensions.
2. `input_dims == [1, 1, 1]`.
3. `output_dims == [1, 1, 1]`.
4. `size(sys.A) == (5, 5)`.
5. `size(sys.B) == (5, 3)`.
6. `size(sys.C) == (3, 5)`.
7. `size(sys.D) == (3, 3)`.
8. `sys.B`, `sys.C`, and `sys.D` are block-diagonal with the expected nonzero entries.
9. `step(sys, x, u)` matches the hand-computed `(A*x + B*u, C*x + D*u)` for a
   concrete `x` and `u`.
10. Omitting `Dv` produces a zero `D` of the correct size.
11. A mismatched `A` size throws `ArgumentError`.
12. A mismatched `Bv[v]` row count throws `ArgumentError`.
13. A mismatched `Cv[v]` column count throws `ArgumentError`.
14. A mismatched `Dv[v]` size throws `ArgumentError`.
15. Calling `step` with a wrong-length `x` or `u` throws `ArgumentError`.

## Verification Checklist

- [ ] `docs/issues/001-discrete-time-lti-control.md` is committed with all required sections.
- [ ] `src/network_sheaves/LTIControl.jl` defines `DiscreteSheafLTI` and `step`.
- [ ] `SheafInterface.jl` declares the new public operations.
- [ ] `NetworkSheaves.jl` includes and re-exports the new module.
- [ ] `test/network_sheaves/LTIControl.jl` covers both valid and invalid inputs.
- [ ] The implementation keeps `A`, `B`, `C`, and `D` sparse.
- [ ] The implementation does not add new package dependencies.
- [ ] `julia --project=test test/runtests.jl` passes.

## Out of Scope

- Continuous-time systems `x'(t) = A x(t) + B u(t)`.
- Constraining the state to the global sections `H^0(F)`.
- Deriving `A` from `d_F`, `L_F`, diffusion, or consensus dynamics.
- Controllability, observability, stabilizability, or Kalman decomposition.
- Edge-space inputs or outputs.
- Model reduction or quotienting by fiber/global-section bases.
- Nonlinear, time-varying, or stochastic dynamics.
- A Catlab-based categorical control interface.

## References

- Justin Curry, *Sheaves, Cosheaves and Applications* — for the cochain-space
  viewpoint and cellular sheaf terminology.
- Robert Ghrist and collaborators on sheaf-theoretic signal processing and
  sheaf Laplacians — for the ambient network-sheaf language used here.
- The existing CellularSheaves.jl pushforward and Laplacian code paths, which
  already fix the package's conventions for `C^0(F)`, `d_F`, and sparse linear
  algebra.
