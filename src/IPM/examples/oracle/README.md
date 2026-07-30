# Per-iteration α oracle

Production-truth data on the *optimal augmentation α per IPM iteration*. At each iteration
`solve_logged` deep-copies the solver and takes one real `step!` at every α on a grid (v3: half-decades
1e0–1e18, 37 candidates, uniform for all problems), scores each `(state, ncraig)` — `state=0` iff every
refinement status reached force/floor, `ncraig` = summed base+refn CRAIG over the step — keeps the
**lowest** α achieving the best score, and advances. Real warm-started steps on full copies, so the cost
is production truth (no cold replay).

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

### v3 additions (window ceiling instrumentation)

- **pass-1 residual, split by side:** `{p,c,w}_res0_dual`, `{p,c,w}_res0_prim` — the dual and primal
  components of the post-base-solve residual (same normalization the solver uses: dual `/(1+nc)`,
  primal `/(1+ng)`). `pres0`/`cres0`/`wres0` stay unchanged as their `max`. The window *ceiling* is
  dual-recovery error growing ∝ α — it lives entirely in the dual component and is invisible inside
  the max whenever the primal side dominates. `w_*` empty on IPM.
- **post-CRAIG constraint residual:** `r1_p`, `r1_c`, `r1_w` = ‖g−Bx‖₂ of the base solve *after* the
  CRAIG correction (recomputed, not from Krylov stats). With `r0_*` (pre-CRAIG) gives the exact
  contraction `rate = (r1/r0)^(1/niters)`.
- **refinement exit residual:** `p_res_exit`, `c_res_exit`, `w_res_exit` — the last `res` before the
  refinement loop returns, any status. Turns binary `state` into graded severity.
- **Cholesky factor diagonal extremes** (per candidate): `Ldiag_min`, `Ldiag_max` =
  min/max `|diag(F)|` — a κ(F) proxy from the factorization and a ρ-ladder early-warning.
- **IPM-level feasibility at step entry** (constant across candidates): `ipm_pres`, `ipm_dres`.
- **barrier-Hessian diagonal, pre-Q** (constant across candidates): `bar_hdiag_med`,
  `bar_hdiag_frac_mid` = median and mid-cluster fraction (`|log10 hⱼ|≤1`) of `hⱼ=dⱼ/pⱼ` over
  cone-constrained coords only (free-cone columns excluded). The degeneracy signature — computed
  *before* Q is folded into H, since a nonzero Q masks it. `hdiag_min/max` remain **post-Q** (the F
  the factorization sees).
- **meta:** `git_dirty`, `final_status`, `niter`; for X-instances, generator params + severity
  (`x_eps`/`x_spread`/`x_min_angle`/`x_colratio`/…).

Grid is 1e0–1e18 (37 half-decades) for all problems, to de-censor window ceilings that previously read
"≥ grid max".

## API

`solve_logged` returns the records in-memory (`Vector{NamedTuple}`); `write_oracle_csv(path, records)`
serializes them. Both live in `src/ipm/solver/oracle.jl`.
