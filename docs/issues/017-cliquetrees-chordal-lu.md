# 017 — `ChordalLU`: supernodal multifrontal LU without pivoting (in CliqueTrees.jl)

**Target repository:** AlgebraicJulia/CliqueTrees.jl (`src/Multifrontal.jl`).
**Consumer:** CellularSheaves `SchwarzMethods` local solves and the
`DoubleIntegratorHJB` direct policy evaluation, both of which currently fall
back to UMFPACK for nonsymmetric matrices.

All file references below are to CliqueTrees v1.19.7
(`src/Multifrontal.jl/src/`, abbreviated `MF/`).

## Mathematical Background

`CliqueTrees.Multifrontal` factors symmetric matrices only (`ChordalCholesky`,
`ChordalLDLt`). The policy-evaluation matrices of a monotone (upwind)
discretization of a Hamilton–Jacobi–Bellman equation,
``A_u = \rho I + G_u``, are nonsymmetric: the upwind stencil follows the drift
and the control. They are, however, nonsingular M-matrices whose rows are
strictly diagonally dominant (by the discount ``\rho``).

For such matrices Gaussian elimination **without pivoting** is stable:

- every leading principal submatrix of a nonsingular M-matrix is a nonsingular
  M-matrix, so all pivots are positive and the LU factorization ``A = LU``
  exists; its Schur complements are again M-matrices and ``L``, ``U`` are
  M-matrices (Fiedler and Pták 1962; Funderlic and Plemmons, *LU decomposition
  of M-matrices by elimination without pivoting*, 1981);
- for a matrix diagonally dominant by rows (or columns) the growth factor of
  elimination without pivoting is at most 2 (Wilkinson 1961; Higham,
  *Accuracy and Stability of Numerical Algorithms*, §9.5).

Without pivoting, the elimination order may be chosen purely for sparsity.
For a matrix ``A`` with an unsymmetric pattern, the fill of ``L + U`` under a
symmetric permutation ``PAP^T`` is contained in the Cholesky fill of the
**symmetrized pattern** ``|A| + |A|^T`` (George and Ng 1987). So the existing
symbolic machinery applies unchanged: a fill-reducing ordering and supernodal
elimination tree of ``A + A^T`` give a static data structure that holds both
``L`` and ``U``, with ``\mathrm{struct}(U) = \mathrm{struct}(L)^T``. This is
the classical *symmetric-pattern multifrontal LU* (Duff and Reid 1984; Amestoy
and Duff, MA41).

UMFPACK's default threshold pivoting (tolerance 0.1) produced element growth of
``6\times10^{19}`` and a relative residual of 8.5, with no error raised, on a
``121^2`` HJB policy-evaluation matrix. CellularSheaves currently forces
partial pivoting to avoid this. A no-pivot factorization on a fixed symbolic
structure is both safer for this class and cheaper to refactor: a policy
iteration (semismooth Newton) refactors a matrix with the same pattern at
every step, so the ordering, the elimination tree and the storage can be
reused, and only the numeric phase reruns.

## Codebase Orientation

| File | Why it matters |
|---|---|
| `MF/chordal_symbolic.jl` | `ChordalSymbolic` (fields `res`, `sep`, `rel`, `chd`, `Dptr`, `Lptr`, `nMptr`, `nMval`, `nFval`, lines 15–28). `symbolic(A::SparseMatrixCSC)` (72–82) **rejects** a nonsymmetric matrix; `symmetric(graph, uplo)` in `MF/utils.jl:55–97` does not fully symmetrize. The LU must build the pattern of ``A + A^T`` itself and call `symbolic(pattern; check=false)`. |
| `MF/chordal_factorization.jl` | `ChordalFactorization{DIAG,UPLO}` (1–20) and the `ChordalCholesky`/`ChordalLDLt` aliases (45, 95). Storage: per-front full ``n_n\times n_n`` column-major diagonal block in `Dval`, off-diagonal block in `Lval` (``n_a\times n_n`` for `:L`, ``n_n\times n_a`` for `:U`). |
| `MF/abstract_factorization.jl` | `getproperty` (`F.L`, `F.U`, `F.P`, 23–41) and `copyto!` (56–71, via `sympermute`) assume symmetry: **do not** subtype `AbstractFactorization` for the LU. |
| `MF/cholesky.jl` | Numeric driver: `factorize!` (272), `chol_impl!` (411–441, postorder over fronts with a LIFO update stack), `chol_loop_snd!` (484–603: zero front, extend-add children, add the original entries, factor, push the Schur complement), `chol_loop_nod!` (610), `chol_send!` (1012). The LU driver is a copy of this loop with full squares instead of triangles. |
| `MF/blas/potrf.jl` | Recursive dense Cholesky (135) with a `potrf2!` base case: the model for a recursive no-pivot `getrf`. |
| `MF/blas/trsx.jl`, `MF/blas/gemx.jl`, `MF/blas/ger.jl` | `trsm!`/`trsv!`, `gemm!`/`gemv!`, `ger!`, reused as is. |
| `MF/chordal_triangular.jl` | `ChordalTriangular{DIAG,UPLO}` views over `(S, Dval, Lval)`; `copy_scatter!` (`:L` 858, `:U` 898) already scatters a full permuted matrix into the diagonal block and the `:L`/`:U` off-diagonal blocks. |
| `MF/divide.jl` | Solves: `div_impl!` (363) and the front loops (437, 590, 726, 880) are generic in `DIAG`/`UPLO` and read only `Dptr`/`Dval`/`Lptr`/`Lval`, so unit-lower and non-unit-upper views over a **shared** `Dval` solve correctly without changes. `DivisionWorkspace` (11) is large enough. |
| `MF/chordal_symbolic.jl:615–772` | `flatindices`: the value-to-slot map used for refactorization; the LU merges the `:L` map with the `:U` map of the separator columns. |
| `test/multifrontal.jl:3–115` | Model for the tests (loop over `UPLO`, refactorization through `flatindices`, BigFloat, `@inferred`/`@test_opt`). |

## Requested Implementation

New file `MF/lu.jl`, included after `cholesky.jl`, exporting `ChordalLU` and
`FChordalLU`.

Storage reuses the existing layout. ``L_{11}`` (unit lower) and ``U_{11}``
(upper) share one `Dval` buffer per front, as LAPACK `getrf` packs them;
`Lval` holds ``L_{21}`` (``n_a\times n_n``) and a new `Uval` holds ``U_{12}``
(``n_n\times n_a``). Then `ChordalTriangular{:U,:L}(S, Dval, Lval)` is ``L``
and `ChordalTriangular{:N,:U}(S, Dval, Uval)` is ``U``.

```julia
"""
    ChordalLU{T, I, ...} <: Factorization{T}

A supernodal multifrontal LU factorization without pivoting,

    P A Pᵀ = L U,

of a square sparse matrix `A` with a nonsymmetric pattern. `P` is a
fill-reducing symmetric permutation and the supernodal elimination tree is
computed from the pattern of A + Aᵀ, so struct(U) = struct(L)ᵀ and
`ChordalSymbolic` is reused unchanged. `L` is unit lower triangular and `U`
upper triangular.

No pivoting is done: the factorization exists and is stable when every
leading principal submatrix of P A Pᵀ is nonsingular with modest growth.
Nonsingular M-matrices and matrices diagonally dominant by rows or columns
qualify for every symmetric permutation P (growth ≤ 2 for diagonal
dominance). For other matrices use a pivoting LU.

    F = lu!(ChordalLU(A))           # analyse and factor
    x = F \\ b                       # also ldiv!(F, B), transpose(F) \\ b
    copyto!(F, A2); lu!(F)          # same pattern: numeric refactorization only

`lu!(F; check=true, tol=0)` stops (and sets `F.info` to the failing column)
at the first pivot with |u_jj| ≤ tol; `check=true` throws
`LinearAlgebra.ZeroPivotException`.
"""
struct ChordalLU{T, I, Dvl <: AbstractVector{T}, Lvl <: AbstractVector{T},
                 Uvl <: AbstractVector{T}, Prm, Ivp, Ifo} <: Factorization{T}
    S::ChordalSymbolic{I}
    Dval::Dvl
    Lval::Lvl
    Uval::Uvl
    perm::Prm
    invp::Ivp
    info::Ifo
end

const FChordalLU{T, I} = ChordalLU{T, I, FVector{T}, FVector{T}, FVector{T},
                                   FVector{I}, FVector{I}, FScalar{I}}

ChordalLU{T}(P::Permutation{I}, S::ChordalSymbolic{I}) where {T, I}   # allocate ndz(S), nlz(S), nlz(S)
ChordalLU(A::SparseMatrixCSC; alg=DEFAULT_ELIMINATION_ALGORITHM, snd=DEFAULT_SUPERNODE_TYPE)
ChordalLU(A::AbstractMatrix, P::Permutation, S::ChordalSymbolic)
Base.copyto!(F::ChordalLU, A::SparseMatrixCSC)        # C = A[perm, perm]; zero; scatter into D, L₂₁, U₁₂
Base.getproperty(F::ChordalLU, :L / :U / :P)
LinearAlgebra.lu!(F::ChordalLU; check::Bool=true, tol=zero(real(eltype(F))))
LinearAlgebra.lu!(W::FactorizationWorkspace, F::ChordalLU; check, tol)
LinearAlgebra.ldiv!(F::ChordalLU, B::AbstractVecOrMat)
LinearAlgebra.ldiv!(F::Transpose{<:Any, <:ChordalLU}, B::AbstractVecOrMat)
Base.:\(F::ChordalLU, b), LinearAlgebra.rdiv!(B, F)
LinearAlgebra.issuccess(F), LinearAlgebra.det(F), LinearAlgebra.logabsdet(F)
Base.size(F), SparseArrays.nnz(F), Base.show(io, mime, F)
flatindices(F::ChordalLU, A::SparseMatrixCSC)         # slots of A's entries in [Dval; Lval; Uval]
```

**Algorithm sketch.**

*Symbolic phase (reuse).*
1. `B = pattern(A)`, with ones on A's stored entries so that values cannot
   cancel.
2. `(P, S) = symbolic(B + transpose(B); check=false, alg, snd)`.
3. Allocate `Dval` (`ndz(S)`), `Lval` and `Uval` (`nlz(S)` each).

*Scatter (`copyto!`).*
1. Form `C = A[perm, perm]` with `colpermute`/`rowpermute` (`MF/utils.jl:225–250`).
   `sympermute` cannot be used, because it reads one triangle.
2. Zero the buffers, then `copy_scatter!` twice: the `:L` variant puts each
   column's `res` rows in the full diagonal block and its `sep` rows in
   ``L_{21}``; the `:U` variant puts the `res` rows of the separator columns in
   ``U_{12}``.

*Numeric phase (`lu_impl!`, a copy of `chol_impl!`).* Visit the fronts ``j``
in postorder with the same LIFO update stack (`Mptr`, `Mval`). Every
`FactorizationWorkspace` buffer is already sized for full squares: `Mval` holds
``\sum n_a^2`` and `Fval` holds ``n_{F}^2``. For a supernodal front, with
``n_n = |res(j)|``, ``n_a = |sep(j)|`` and the front ``F`` of size
``(n_n + n_a)^2``:

1. `fill!(F, 0)`, the full square rather than `zerotri!`.
2. Pop each child ``i`` and extend-add its full ``n_{a,i}\times n_{a,i}``
   update matrix into ``F`` through `rel(i)`. This needs a new `addscatter!`,
   the full-square counterpart of `addscattertri!`.
3. Add the original entries: ``F_{11} \mathrel{+}= D_j``, ``F_{21}
   \mathrel{+}= L_{21}``, ``F_{12} \mathrel{+}= U_{12}``.
4. `getrf_nopiv!(F₁₁; tol)`: recursive like `potrf!` (`MF/blas/potrf.jl:135`):
   factor the leading half, two `trsm!`, a `gemm!` Schur update, recurse,
   with an unblocked base case for ``n \le 64``. Record `info` at the first
   ``|u_{kk}| \le`` `tol`.
5. ``L_{21} \leftarrow F_{21} U_{11}^{-1}``: `trsm!` with side R, uplo U,
   non-unit.
6. ``U_{12} \leftarrow L_{11}^{-1} F_{12}``: `trsm!` with side L, uplo L,
   unit.
7. ``F_{22} \mathrel{-}= L_{21} U_{12}``: `gemm!` in place of `syrk!`, on the
   full square.
8. Copy ``F_{11}`` into `Dval`, then push ``F_{22}`` as the full
   ``n_a\times n_a`` update matrix.

For single-column fronts (`nod`), the pivot is ``f_{11}``,
``l = f_{21}/f_{11}``, ``u = f_{12}``, and the update is a `ger!`.

*Solve.* `ldiv!(F, B)`:
1. `mul!(C, F.P, B)`;
2. `ldiv!(W, F.L, C)`, the unit-lower forward pass;
3. `ldiv!(W, F.U, C)`, the upper backward pass;
4. `ldiv!(B, F.P, C)`.

The transpose solve is ``P^T (L^{-T} (U^{-T} (P b)))``, using the existing
`Transpose` wrappers of `ChordalTriangular`.

*Determinant.* ``\det A = \prod_j \det U_{11}^{(j)}`` (``\det P^2 = 1``), and
`logabsdet` sums ``\log|u_{kk}|``.

*Not reused:* `AbstractFactorization`'s `adjoint(F)=F`, `det = det(L)^2`,
`copyto!` via `sympermute`, and the regularization hooks (`initialize` and
`regularize` are sign-based and specific to the symmetric case).

## Tests to Write

In `test/multifrontal.jl`, add a `@testset "lu"`:

```julia
using SparseArrays, LinearAlgebra, Random
Random.seed!(1)
# Upwind M-matrix: rho*I + G with G a negated-off-diagonal row-stochastic generator.
function upwind(n; rho=0.5, density=4/n)
    G = sprand(n, n, density) .* 1.0
    G = G - Diagonal(G)
    return rho * I + Diagonal(vec(sum(G; dims=2))) - G
end
# SPD matrices are also safe without pivoting (growth ≤ 1).
for M in (upwind(500), upwind(2000), readmatrix("HB/bcsstk14"))
    b = randn(size(M, 1))
    F = lu!(ChordalLU(M))
    @test issuccess(F)
    @test M * (F \ b) ≈ b
    @test transpose(M) * (transpose(F) \ b) ≈ b
    @test M * (F \ [b b]) ≈ [b b]
    @test logabsdet(F)[1] ≈ logabsdet(lu(M))[1]
    @test det(F) ≈ det(lu(M)) rtol = 1e-8
    @test nnz(F) == nnz(F.L) + nnz(F.U) - size(M, 1)
    # numeric refactorization on the same pattern
    M2 = copy(M); nonzeros(M2) .*= 1 .+ rand(nnz(M2)) / 10
    copyto!(F, M2); lu!(F)
    @test M2 * (F \ b) ≈ b
end
# The factor of a symmetric positive definite matrix matches its Cholesky structure.
A = readmatrix("HB/bcsstk14")
F = lu!(ChordalLU(A))
@test F.L * F.U ≈ (F.P * A * F.P')        # dense check on a small matrix only
# Zero pivot without pivoting is reported, not hidden.
Z = sparse([0.0 1.0; 1.0 0.0])
@test_throws ZeroPivotException lu!(ChordalLU(Z))
@test !issuccess(lu!(ChordalLU(Z); check=false))
# Generic element types.
Mb = BigFloat.(upwind(60))
@test Mb * (lu!(ChordalLU(Mb)) \ ones(BigFloat, 60)) ≈ ones(BigFloat, 60)
@inferred lu!(ChordalLU(upwind(100)))
```

And in CellularSheaves, once a CliqueTrees release has `ChordalLU`, add to
`test/network_sheaves/SchwarzMethods.jl`, in the "nonsymmetric M-matrices"
testset, a check that local `ChordalLU` factors give the same solutions as the
UMFPACK ones (`@test r.u ≈ u rtol = 1e-8` for every sweep).

## Verification Checklist

- [ ] `symbolic` is called on the symmetrized pattern, and a nonsymmetric `A`
      is never passed to `symbolic` directly.
- [ ] `lu!` matches UMFPACK's `lu(A) \ b` to `1e-10` on the upwind test matrices.
- [ ] Refactoring through `copyto!` + `lu!` allocates nothing beyond the
      workspace (`@allocated` check after warm-up).
- [ ] The solve path reuses `div_impl!` unchanged; only factorization-level
      wiring is new.
- [ ] A zero or tiny pivot raises `ZeroPivotException` (or sets `info` with
      `check=false`); the code never silently returns a garbage factor.
- [ ] The docstring states the class of matrices for which no pivoting is safe
      and cites Funderlic–Plemmons and Higham §9.5.
- [ ] Benchmarked against UMFPACK on the CellularSheaves HJB local problems
      (`docs/scripts/hjb_benchmarks.jl`, 4-D boxes of 7–9 points per side):
      factor time, refactor time and solve time.
- [ ] Every existing CliqueTrees test passes, and the `ChordalLDLt` and
      `ChordalCholesky` code paths are untouched.

## Out of Scope

- Pivoting of any kind: threshold, static perturbation or delayed pivots. A
  later issue can add static pivot perturbation through the existing
  regularization interface.
- Tree-level parallelism (independent subtrees on separate tasks with
  per-task update stacks). This would benefit Cholesky and LDLᵀ equally and
  belongs in its own issue.
- Unsymmetric orderings (column ordering of ``A^TA`` as in COLAMD) and
  unsymmetric elimination DAGs: this issue is the symmetric-pattern
  multifrontal LU only.
- Switching CellularSheaves from UMFPACK. That is a follow-up once a
  CliqueTrees release ships `ChordalLU`.
