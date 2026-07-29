# Per-iteration α oracle

Production-truth data on the *optimal augmentation α per IPM iteration*. At each iteration
`solve_logged` deep-copies the solver and takes one real `step!` at every α on a grid (half-decades
1e0–1e10), scores each `(state, ncraig)` — `state=0` iff every refinement status reached
force/floor, `ncraig` = summed base+refn CRAIG over the step — keeps the **lowest** α achieving the
best score, and advances. Real warm-started steps on full copies, so the cost is production truth
(no cold replay).

## Run

```bash
julia --project=examples examples/oracle/run_oracle.jl <key> <tol> <hsd|ipm>
# e.g.
julia --project=examples examples/oracle/run_oracle.jl e08 1e-8 hsd
```

Writes `<key>_<tol>_<solver>.csv` here — one row per `(iteration, α)`.

Keys are the example ids (`e01`…`e15`) plus `06`/`07`/`13` for e06/e07/e13 and `X03`/`X04`
adversarials — see `buildprob` in `run_oracle.jl`.

## Columns

- **score:** `state` (0 all-reached / 1 any stall-or-itmax), `ncraig`, `ipm_status`, `chosen`
  (marks the α picked that iteration).
- **candidate features** (what a controller could see): `mu`, `rho`, `force_tol` (=η·μ),
  `floor_tol`, `sigma2min` (decoupled Golub-Kahan estimate, `IPM.oracle_sigma2min`), and each
  solve's base/first-refine residuals `pres0/pres1`, `cres0/cres1`, and — HSD only — `wres0/wres1`
  (empty for IPM, which has no woodbury role). `res1 = NaN` means the base solve hit atol at pass 1
  (no second pass); it populates on α's where refinement runs ≥2 passes.
- **per-role counts/status:** `{p,c,w}base`, `{p,c,w}refn`, `{p,c,w}stat`.

## API

`solve_logged` returns the records in-memory (`Vector{NamedTuple}`); `write_oracle_csv(path, records)`
serializes them. Both live in `src/ipm/solver/oracle.jl`.
