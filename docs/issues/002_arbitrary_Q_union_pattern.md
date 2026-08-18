# 002 — Arbitrary (non-block-diagonal) Q via a unioned factor pattern

## Summary

Today the quadratic `Q` in

```
min ½ pᵀQp − cᵀp   s.t.  Bp = g,  p ∈ K
```

must be **block diagonal** (one dense block per vertex stalk, `allocblockdiag(B)`). That
is a modeling limitation: `Q` can legitimately couple any two stalks. This issue lifts
`Q` to an **arbitrary symmetric** matrix.

The augmented Hessian factored each iteration is

```
F Fᵀ = (1/α)·A + BᵀB,   A = Q + H   (H = the block-diagonal cone/barrier Hessian)
```

so its sparsity is `pattern(Q) ∪ pattern(BᵀB)`. The elimination ordering and symbolic
factor are computed on that **union** graph; nothing constrains `Q ⊆ BᵀB`.

Block-diagonal `Q` is the special case `pattern(Q) ⊆ pattern(BᵀB)` ⟹ union = `BᵀB`, so
the change must be **bit-for-bit inert** on today's problems.

## The organizing invariant: one triangular matrix

**There is exactly one triangular object in the system — the Cholesky factor `F`
(`ChordalTriangular{:L}`).** Everything else is stored **full-symmetric**:

- `Q`, `H`, and `L = B'B` are full block-sparse symmetric matrices (`L` already is).
- Every matvec is a plain `mul!(u, M, x)` on the full matrix — **no `Symmetric(·,:L)`
  wrappers**. Correct because the stored matrix is symmetric; contiguous and
  cache-friendly.
- The factorization assembly into `F` (`axpy!`/`copyto!`, `utils.jl:153/324`) already
  reads only the **lower** triangle of its full source and skips the upper (line 182,
  `ulo < rlo && continue`). So each symmetric pair `(i,j)/(j,i)` contributes once.

This makes the whole thing clean: full-symmetric `H` factors correctly with **zero change
to the factorization glue**, and `Q`'s upper half is not wasted — the factorization
ignores it but the matvecs use it (that's what makes the full product fast).

Full-symmetric storage also makes the **elimination permutation trivial**: `P Q Pᵀ`
(`halfselectvtxs ∘ halfselectvtxs`) keeps a full-symmetric matrix full-symmetric. (A
lower-only store would cross the diagonal under permutation and silently corrupt — hence
full storage, deliberately.)

## Mathematical Background

`H` is the Nesterov–Todd / Tunçel cone Hessian, block-diagonal by construction (one block
per cone). `Q` is fixed problem data. `A = Q + H` ⟹ `pattern(A) = pattern(Q) ∪ diag`, and
`F Fᵀ = (1/α)A + BᵀB` ⟹ factor pattern `= pattern(Q) ∪ pattern(BᵀB)`.

Two things already generalize for free:

1. **Numeric hot paths** — all residual/objective/step matvecs go through `Symmetric(Q)`
   / `mulkkt!` today; switching to plain full `mul!` is a simplification, not new logic.
2. **The permutation** — `Q ← P Q Pᵀ` already permutes both block axes (init(), the
   `halfselectvtxs∘halfselectvtxs` line); off-diagonal blocks reorder correctly.

What is *not* general: the ordering/symbolic graph (built from `B` only), and equilibration
(hard-coded block-diagonal). And the per-iteration `H` assembly must be generalized.

## Codebase Orientation

| File | Role | Change |
|---|---|---|
| `src/BlockSparseArrays/src/…` | `weightedgraph`, block-sparse alloc | `weightedgraph(B, Q)` = union `Q`'s off-diagonal adjacency into `linegraph(B)` |
| `src/IPM/src/kkt/kkt.jl` | `makekkt(::Type{Solver}, B; elim)`; `mulkkt!` uses `Symmetric(A,:L)` | take `Q`, order on the union graph; drop the `Symmetric` wrapper |
| `src/IPM/src/ipm/solver/hsd.jl` | `H = allocblockdiag(B)` (449); step! `H += Q` loop (669); `Symmetric(s.Q/Q,:L)` (154/219/373) | `H = copy(Q)`; `fill!`/`scale!`/raw `axpy!`; drop wrappers |
| `src/IPM/src/ipm/solver/ipm.jl` | `H = allocblockdiag(B)` (317); step! loop (425); `Symmetric(·,:L)` (52/136/241/357) | same |
| `src/IPM/src/ipm/solver/solver.jl` | `scale!` writes `block(H,v,v,v)` | live O(deg) scan for the diagonal block (arc with `tgt == v`) |
| `src/IPM/src/scaling/scaling.jl` | `infnorms!` (57), `applyscaling!` (101); block-diagonal comment (39–40) | walk **all** of `Q`'s blocks (full store — no lower-only accumulation) |
| `src/IPM/src/ipm/problem.jl` | `IPMProblem` construction | `@argcheck` `Q` full-symmetric with all diagonal blocks present |

## Requested Implementation

### 1. Union ordering graph

```julia
# weightedgraph(B, Q): linegraph(B) (column-intersection = BᵀB pattern) UNIONed with
# Q's off-diagonal block adjacency {(v,w) : Q has block (v,w), v ≠ w}.
function weightedgraph(B::BlockSparseMatrix, Q::BlockSparseMatrix)
    weight, g = weightedgraph(B)
    return weight, union_offdiag_adjacency(g, Q)
end

function makekkt(::Type{Solver}, B, Q; elim = DEFAULT_ELIMINATION_ALGORITHM) where {Solver}
    weights, graph = weightedgraph(B, Q)
    R, P, S = symbolic(weights, graph; alg = elim)   # S ⊇ pattern(Q) ∪ pattern(BᵀB)
    B = selectvtxs(B, R.perm)
    F = FChordalTriangular{:N, :L, T, I}(S)
    return R, P, B, Solver(F, B)
end
```

Callers (ipm.jl:300, hsd.jl:435) pass the pre-permutation `Q`; `Q` is then permuted by
`R.perm` as it already is.

### 2. `H = copy(Q)` + per-iteration assembly

- init(): `H = allocblockdiag(B)` → `H = copy(Q)` (after `Q`'s permutation, so `H` shares
  the permuted full-symmetric pattern).
- step! `H += Q` loop → the three-line raw form:

```julia
fill!(s.H, 0)                        # zero all blocks
scale!(s)                            # barrier → diagonal blocks (overwrites the zeros)
axpy!(one(T), s.Q.val, s.H.val)      # H += Q, flat val-vector (H and Q share layout)
```

`scale!` stays an overwrite (the block is zero when it writes). Gives
`H_vv = barrier + Q_vv`, `H_vw = Q_vw`.

### 3. Drop the `Symmetric(·,:L)` wrappers

`mulkkt!` (kkt.jl:21) and every solver matvec: `mul!(u, Symmetric(A,:L), x)` →
`mul!(u, A, x)`; `Symmetric(s.Q,:L)` → `s.Q`. Sites: kkt.jl:21; hsd.jl 154/219/373;
ipm.jl 52/136/241/357 (and the `dot(p, Symmetric(Q,:L), p)` objective →
`dot(p, s.Q, p)`).

### 4. `scale!` diagonal-block find (solver.jl)

`Hv = block(H, v, v, v)` no longer holds — for a full-symmetric `H`, scan `v`'s arc range
for the arc with `tgt == v`:

```julia
diagblock(X, v) = (e = findfirst(a -> X.tgt[a] == v, srcrange(X, v)); block(X, v, v, e))
```

Live O(deg) scan; acceptable (`scale!` is already O(work) per vertex).

### 5. Equilibration (scaling.jl) — simpler with full storage

`Q̂ = D Q D`. With a full store, `infnorms!` and `applyscaling!` iterate **all** of `Q`'s
blocks directly:
- `infnorms!`: for block `(v,w)`, fold `‖Q_vw‖_∞` into column `v`'s ∞-norm (no
  transpose-accumulation — the `(w,v)` block is also stored).
- `applyscaling!`: scale block `(v,w)` by `s_v·s_w` (diagonal by `s_v²`), via `diagblock`
  for the diagonal.

### 6. Problem guard (problem.jl)

`@argcheck` `Q` is full-symmetric block-sparse with a diagonal block for every vertex.
No `Q ⊆ BᵀB` check.

## Tests to Write

```julia
# Regression: block-diagonal Q is bit-identical (union collapses to BᵀB).
@test solve(prob_blockdiag, set).pobj == baseline_pobj
@test solve(prob_blockdiag, set).history.μ == baseline_μ

# Full-Q objective agrees with a dense reference (Clarabel on the same program).
@test abs(solve(prob_denseQ, set).pobj - clarabel_obj) / (1 + abs(clarabel_obj)) < 1e-6

# Equilibration invariance for dense Q.
@test solve(prob_denseQ; scale=true).pobj ≈ solve(prob_denseQ; scale=false).pobj

# Union pattern: factor fill covers Q-only positions.
@test pattern(F) ⊇ pattern(Q) ∪ pattern(B'B)

# The full matvec (no wrapper) equals the old Symmetric(·,:L) matvec.
@test mul!(u1, Q_full, x) ≈ mul!(u2, Symmetric(Q_full, :L), x)
```

## Verification Checklist

- [ ] Block-diagonal `Q`: iterations and `pobj` **bit-identical** to pre-change.
- [ ] `weightedgraph(B, Q)`: unchanged for `Q ⊆ BᵀB`; gains exactly `Q`'s extra edges otherwise.
- [ ] `H = copy(Q)` after permutation; `fill!`/`scale!`/`axpy!` gives `barrier + Q` on the diagonal, `Q` off-diagonal.
- [ ] No `Symmetric(·,:L)` wrappers remain in the solver matvecs (grep); no `block(X,v,v,v)` "diagonal" access (grep) — all via the scan/`diagblock`.
- [ ] Dense-`Q` objective matches Clarabel; equilibration scale-invariant.
- [ ] `P Q Pᵀ` keeps `Q` full-symmetric (spot-check a permuted off-diagonal block equals its transpose partner).

## Out of Scope

- Non-symmetric `Q`.
- `Q` that changes across iterations (fixed problem data).
- Caching the diagonal-block arc (`diagarc`) — deferred; the live scan is cheap enough.
- Dynamic re-ordering as barrier fill grows — the union pattern is computed once.
