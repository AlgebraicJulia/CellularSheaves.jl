# 016 — CI compile time of the Mumblebee IPM under `--check-bounds=yes`

Status: diagnosed and measured, fixes not started (deferred). Recorded so the
numbers and the plan are not lost.

## Mathematical Background

None; this is a build and CI issue. The code involved is Mumblebee's
block-sparse matrix–vector product (`BlockSparseArrays/src/blas/gemv.jl`),
which the IPM calls at every Newton step on the KKT system.

## Codebase Orientation

| Where | Why |
|---|---|
| `test/ControlSheaves/PredictiveConsensus.jl` | First tests that call the IPM. 43 tests: about 1.5 min in a normal session, about 13 min under CI. |
| `.github/workflows/julia_ci.yml` | Calls the shared `AlgebraicJulia/.github` workflow, which runs `julia-actions/julia-runtest@v1` with its defaults (`check_bounds: yes`, `coverage: true`) and exposes no `check_bounds` input. |
| Mumblebee.jl `src/BlockSparseArrays/src/blas/gemv.jl` | `@generated` kernels unrolled over block sizes up to `GEMV_MAXN = 8` and tiles of `GEMV_FWD_TILE = 16`, selected by `if nrow == 1 … elseif nrow == 8` chains generated at compile time. |
| Mumblebee.jl `src/Mumblebee.jl` | No precompile workload, so every fresh process compiles the IPM. |

## Diagnosis (measured 2026-10-09, Julia 1.12.6, Mumblebee@9d308da)

All times are the first bounded `PredictedTrajectorySweeps` solve on a 3×3 grid
in a fresh process; the second solve takes 0.07 s plain, 0.14 s with bounds checks.

Compile time per method (`--trace-compile-timing`):

| | plain | `--check-bounds=yes` |
|---|---|---|
| `Mumblebee.IPM.solvepredictor!` | 16 s | 458 s |
| `Mumblebee.IPM.reinit!` | 11 s | 231 s |
| CellularSheaves code | 22 s | 27 s |
| whole first solve | 77 s | 793 s |

- Coverage is not the cause: `--code-coverage=@<pkg>` 98 s vs 98 s plain
  (`--code-coverage=user` would be 281 s, but CI does not use it).
- Cause: the static `if` chains make one GEMV call reach every unrolled kernel
  (about 64 per transpose mode), all compiled on first use whatever block
  sizes the problem has. `--check-bounds=yes` ignores `@inbounds`, so every
  unrolled load gets a bounds check and the IR LLVM must optimise grows about 25×.

Experiments with modified copies of Mumblebee:

| Variant | first solve, plain | first solve, bounds | 64-agent run time (central / sweeps, 1 thread) |
|---|---|---|---|
| current (8, 16) | 77 s | 793 s | 1.62 s / 1.61 s |
| `GEMV_MAXN = 4`, tile 8 | 28 s | 73 s | 1.56 s / 1.36 s |
| `GEMV_MAXN = 2`, tile 4 | 20 s | 35 s | 1.53 s / 1.63 s |
| current + PrecompileTools workload | 4.7 s | 6.3 s | — |

The workload variant moves the cost into a one-time precompile of Mumblebee:
75 s plain, 726 s with bounds checks. It is cached by `julia-actions/cache`
because Mumblebee is pinned. A workload in CellularSheaves would not help CI:
our own package image is rebuilt on every source change.

## Requested Implementation

In order of preference; 1 and 4 together are the recommended upstream PR.

1. **Compile only the block sizes used (Mumblebee).** Replace the generated
   `if nrow == m` / `if ncol == n` chains in `gemv.jl` with run-time dispatch
   on `Val(nrow)` / `Val(ncol)` behind a function barrier, so only sizes that
   occur are compiled (about 9 kernels for our 6/4/3 stalks instead of 64).
2. **Bounds-check-independent kernels (Mumblebee).** Load through pointers
   (`GC.@preserve` + `SIMD.vload` / `unsafe_load`) inside the unrolled kernels
   so their code size does not depend on `--check-bounds`. Larger change;
   pointer loads are unchecked by design.
3. **Lower the unroll limits (Mumblebee).** `GEMV_MAXN = 2` gives 35 s under CI
   flags with no run-time loss on our problems. Measured only on our
   problems; Mumblebee's own examples may rely on the unrolling.
4. **Precompile workload (Mumblebee).** `PrecompileTools` workload solving a
   small conic QP with free and second-order-cone stalks, re-solving after
   `reinit!`, and calling `frule!`. A tested version exists (draft from this
   investigation); combine with 1 or 3 to keep the one-time precompile cheap.
5. **CI fallback (`AlgebraicJulia/.github`).** Add a `check_bounds` input to the
   shared `julia_ci.yml` and pass `auto` here. Removes the 10× penalty for all
   packages, but gives up bounds checking in tests.

## Tests to Write

No new package tests. Verify by timing, in a fresh process:

```julia
# julia --check-bounds=yes --project=. -e '...'
t = @elapsed solve(ConsensusProblem(lq, x0), PredictedTrajectorySweeps(maxiter = 2))
@test t < 60          # currently 793 s
```

and by the CI test step on Ubuntu returning to about its `main` baseline
(9 min before PR #96).

## Verification Checklist

- [ ] First bounded solve under `--check-bounds=yes` below 60 s in a fresh process.
- [ ] 64-agent benchmark run time unchanged within noise (centralized and sweeps).
- [ ] Mumblebee's own examples' run time unchanged (if options 1–3 touch kernels).
- [ ] CI test step for PR #96 within a few minutes of the `main` baseline.
- [ ] PR #96 moved out of draft.

## Out of Scope

- Run-time optimisation of the IPM.
- The three Mumblebee correctness workarounds in `PredictiveConsensus`
  (`reinit!` keeping the barrier Hessian, `step_frac = 0.99` failures at μ′ > 0,
  NaN tangent for Δf = 0); report those upstream separately.
