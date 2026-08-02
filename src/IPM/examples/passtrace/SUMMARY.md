# passtrace — per-pass refinement instrumentation (E-series "targeted rerun")

git: `f0590af9` · 15 files × 2 runs = 30 CSV+meta pairs · sweep wall ~18 min (sequential)

Records the one term of the α cost-model that was never measurable: the Krylov (CRAIG) work done
*inside refinement passes*, and the per-pass entry residual it depends on. Every (problem, iteration, α)
row carries the oracle superset schema plus the new per-pass trace.

## New columns (per role p=predictor, c=corrector, w=woodbury[HSD only])
- `<role>_dres{1..6}` / `_pres{1..6}` / `_tres{1..6}` — the dual / primal / τ **entry residual** of refinement
  pass i (measured at the top of refinement iteration i, BEFORE that pass's CRAIG). `NaN` = the loop never
  reached pass i. `tres` is the HSD 3-row τ-component (NaN on the IPM 2-row system).
- `<role>_nkry{1..6}` — CRAIG iterations spent by pass i's inner solve (`solve_kkt!` return − 1). `−1` = pass
  i did not fire (loop broke at the force/floor/stall test before solving).
- `tau_stall`, `refine_cap` — the stall threshold and pass cap actually in force (run1 vs run2; never inferred).
- `normB2` — ordinary ‖B‖² (Frobenius). This is the OTHER sense of "σ²max"; the existing `sigma2max` is the
  **preconditioned** Ritz σ̂²max of B F⁻¹ Bᵀ, NOT ‖B‖². (§5 naming hazard resolved.)
- Base invocation (pass 0) is not in the trace: entry = `r0_{p,c,w}`, CRAIG = `pbase−1`, and its post-solve
  residual equals the pass-1 entry (`p_dres1`/`p_pres1`). `r1_{p,c,w}` = post-CRAIG base residual (equals r0
  when base CRAIG did nothing — that is the semantics, not a bug; §5 hazard resolved by leaving them alone).

## Two runs per file
- **run1** — production (`refine_stall=0.5`, `refine_itmax=10`). Baseline; matches current behaviour + new cols.
- **run2** — stall trigger STRICTLY inert (`refine_stall=Inf`, never fires) + `refine_itmax=20`. (Used Inf, not
  the spec's suggested 0.999, because 0.999 still stops a diverging pass early; Inf guarantees the spec's
  "some α hit the cap without converging".)

## Findings (predictor role; Σ ≡ per-pass refinement CRAIG = `prefn − ppass`)

1. **`a = log(entry_resid₁/ε) > 0` is tautological for pass-firing rows** (9 of 4972 exceptions, numerical).
   A pass fires *because* its entry residual exceeds ε=max(force_tol,floor_tol), so `a>0` by construction.
   `a>0` is therefore NOT a group discriminant — it holds in every group, both runs. The real discriminant
   is Σ>0 (whether the fired passes actually consume Krylov iters).

2. **HEADLINE — production refinement-CRAIG is ~0 across all 15 problems.** In run1, only **42 of 2478**
   pass-firing rows (1.7%) have Σ>0, with **max Σ=5**. The production stall (0.5) cuts refinement before CRAIG
   engages — this is empirically *why* the model's refinement-Krylov term was never testable. Under
   production settings the term is negligible.

3. **Σ>0 is strongly enriched in PRIMAL-driven passes, but not exclusively.** run1 cross-tab (driven × Σ>0):
   primal-driven 25/28 rows have Σ>0 (89%); dual-driven 17/2450 (0.7%). The plan's "all Σ>0 rows are
   primal-driven" is directionally right (huge enrichment) but **strictly false** — 17 dual-driven Σ>0 rows
   exist. Primal-driven rows are rare (62 total, 1.2%) and concentrated in `07_hsd` and `13_hsd`, matching
   the plan's ~1.5% expectation.

4. **The A/B grouping does NOT predict Σ-appearance → the central prediction is REFUTED.** With the stall
   inert (run2), the fraction of pass-firing rows developing Σ>0 overlaps fully between the "Sigma-should-
   appear" (A) and "control" (B) groups:

   | grp | file | run2 Σ_max | run2 frac Σ>0 |
   |---|---|---:|---:|
   | A | e02_hsd | 6000 | 0.145 |
   | A | e02_ipm | 7 | 0.018 |
   | A | e03_ipm | 1 | 0.008 |
   | A | e01_hsd | **0** | **0.000** |
   | B | e08_hsd | 217 | 0.112 |
   | B | e08_ipm | 14 | 0.058 |
   | B | X08b_hsd | 11 | 0.026 |
   | B | X08b_ipm | 14 | 0.013 |
   | B | X03_ipm | **0** | **0.000** |
   | C | 07_hsd | 29 | 0.151 |
   | C | 13_hsd | 15 | 0.102 |
   | C | SOS_P8_n9_s1_hsd | 26 | 0.106 |
   | D | 06_ipm | 7 | 0.021 |
   | D | SOS_P4_n13_s1_hsd | 23 | 0.113 |
   | D | X06b_hsd | **0** | **0.000** |

   Both A and B contain a member that never develops Σ (e01 / X03) and a member near ~0.11–0.14 (e02 / e08).
   Whatever theoretical quantity the A/B ends were meant to span, it does not gate whether Σ can be nonzero.

5. **Two divergence modes when the stall is removed.** (i) **CRAIG-engaging** — e02_hsd (Σ→6000, the CRAIG
   dimension cap), e08_hsd (217), 07/13/SOS: refinement passes run CRAIG hard into a diverging correction.
   (ii) **pure direct-solve** — e01_hsd, X03_ipm, X06b_hsd, e02_ipm: nkry=0 every pass yet the DUAL entry
   residual grows geometrically (e02_ipm α=3.2e15: p_dres 0.9→1.5→10→51→242→1290 over 6 of its 20 passes,
   primal at machine-zero throughout). The stall trigger was the single guard against both; run2's capped-row
   counts (up to 514 on 07_hsd) and the trajectory changes (SOS_P8: NEAR_OPTIMAL/17it → CONTINUE/15it;
   07_hsd: 17→22 it) show it is load-bearing.

## Solver changes (this branch, uncommitted)
- **Instrumentation (behaviour-preserving):** workspace-resident `PassTrace` buffers on `UzawaWorkspace`
  (kkt/uzawa.jl), filled in place by `refinekkt!` (ipm.jl) and `refinehsd!` (hsd.jl); read off the workspace
  by `_oracle_record` (oracle.jl). Zero-alloc, same pattern as the spectral harvest. No history / `step!` /
  return-signature changes. Nothing in solver logic reads the trace.
- **Robustness fix (harvest only):** `gk_spectral`/`ritz_spectral` guard against a non-finite Golub-Kahan
  bidiagonal and catch `LAPACKException` from the divide-and-conquer SVD (a diverging refinement direction at
  extreme α made it throw). Degrades to an all-NaN reading. The harvest is read-only diagnostics, so this
  cannot change convergence — but it was necessary: the harvest runs on every `solve_kkt!` including
  refinement passes, and run2's deep passes trip it.

## Reproduce
`examples/passtrace/run_passtrace.jl <key> <tol> <ipm|hsd>` (runs run1+run2). The 15-file sequential sweep:
`/Users/richardsamuelson/.claude/jobs/7f05e925/tmp/run_passtrace_all.sh`; analysis
`/Users/richardsamuelson/.claude/jobs/7f05e925/tmp/analyze_passtrace.py`.
