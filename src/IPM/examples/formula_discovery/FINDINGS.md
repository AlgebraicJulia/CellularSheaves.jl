# Distance-sensor formula discovery — findings log

Running record of the `distance_formula_discovery_spec.md` pipeline. Data: `oracle2/*.csv` +
`*.meta.json` → `dataset.pkl` (via `build_dataset.py`). Targets: `d_lo` = decades above the window
floor, `d_hi` = decades below the window ceiling. Evaluation: grouped 5-fold CV by problem
(`stage_a.py` etc.). Noise floor ≈ 0.2–0.3 decades (half-decade grid quantization) — "solved".

**Dataset:** 76,590 rows · **96 distinct problems** · 188 files (ipm 40k / hsd 36k). Zones: in 33k /
above 31k / below 12k. 6 adversarial families (X09/X04/X06/X08/SOS + fleet) × dials × seeds × 2
solvers. (Spec was written against 23 problems — fingerprint/overfit worries are far weaker here.)

---

## Stage A — forest reference + feature discovery  ✅

**Reference forest** (HistGradientBoosting depth 4, grouped 5-fold CV, MAE in decades):

| sensor | pooled MAE | operational slice | zones (below / in / above) |
|---|---|---|---|
| `d_lo` (floor)   | 0.61 | **0.50** (below+in) | 0.57 / 0.48 / 0.77 |
| `d_hi` (ceiling) | 0.62 | **0.62** (in+above) | 0.63 / 0.53 / 0.71 |

Above the ~0.25 noise floor → real learnable structure, not memorization. **Ceiling now as
well-predicted as the floor** (closes the spec's "ceiling is thinner" gap).

**Fingerprints don't help** — `core` (no `m/n/ng/alpha_anchor`) ≥ `full` (d_lo 0.610 vs 0.633; d_hi
0.620 vs 0.631). Sensors **transfer**; use core. (Payoff of 96 problems vs 23.)

**Cluster-permutation importance (core features):**
- **Both sensors dominated by the RESIDUAL family** (+5.1 d_lo / +4.1 d_hi): `r0_*, r1_*, pres*,
  cres*, *_res0_dual` all move together → floor/ceiling is fundamentally "residual crosses a
  tolerance." Stages C/D pick one parametrization from this cluster.
- **Floor (`d_lo`):** then `force_tol+mu` (+0.64), counter clusters `pbase` (+0.29), `cbase+ncraig`
  (+0.15). Classic floor law: residual vs forcing tolerance + refinement-count staircase.
- **Ceiling (`d_hi`): standout secondary is `Ldiag_min` (+1.10)** — far above `force_tol+mu` (+0.60)
  and `hdiag_min` (+0.33). `Ldiag_min` = Cholesky-factor min diagonal = the `κ(F)`/`λ_min(F)` proxy.
  **The ML independently rediscovered the λ_min-shrink ceiling mechanism** — and confirms the v3
  `Ldiag` instrumentation (added to *see* the ceiling) carries signal the residuals alone don't.

**Ingredient sets going into the formula stages:** floor = residual + `force_tol/μ` + counter
staircase; ceiling = residual + **`Ldiag_min`** + `force_tol/μ`.

---

## Stage B — capacity ladder + EBM shapes  ✅ (ladder)

**Capacity ladder** (grouped-CV MAE, core features):

| target | depth-4 | depth-2 | depth-1 (additive) | linear |
|---|---|---|---|---|
| `d_lo` (floor)   | 0.610 | 0.828 | 1.068 | 0.932 |
| `d_hi` (ceiling) | 0.620 | 0.808 | 1.103 | **2.847** |

- **Both sensors carry real interaction content** — depth-4 → depth-1 gap ≈ 0.46 for both. Not a pure
  sum of 1-D shapes; products/thresholds (residual-vs-tolerance crossing, counter×residual) matter.
- **Floor is near-linear-friendly** (linear 0.93 vs best 0.61) → a readable weighted-sum formula will
  capture most of `d_lo`.
- **Ceiling is strongly nonlinear** (linear **2.85** ≫ 0.62) → needs threshold/hinge + product terms,
  consistent with a roundoff *crossing* (`ε·κ(F)` vs `atol`), not a smooth sum. Expect hinges around
  `Ldiag_min`.

_(EBM shape functions appended below when the fit completes.)_

**EBM (full 76k, all 96 problems, weak settings interactions=0/outer_bags=4/max_rounds=3000):**
grouped-CV MAE **d_lo=0.729, d_hi=0.731** — additive-model level (between depth-1 stumps 1.07 and
depth-4 forest 0.61; ~0.11 interaction gap). Far better than diversity-starved subsamples (6 problems
→1.42; 14→1.09), confirming the bottleneck was **problem diversity, not rows**.

_Timing lesson (this cost a detour): EBM's `max_rounds` default=50000 causes super-linear blowup
(8s@2k→60s@4k); the fix is `max_rounds≈3000` → linear ~0.5s/1k, full 76k = **37s/fit**. `outer_bags`
(14→4) is a flat ~3.5× on top. `interactions` was 0 throughout — never the issue. Fit on numpy array ⇒
EBM has no feature names; index shapes by CORE position, not `term_names_`._

**Shape functions (x=feature value, y=±decades contributed):**
- `d_lo`: `alpha` linear ramp (+0.54/dec); `force_tol` ramp↑; `r0` ramp↓; **`pbase` saturating
  staircase** (0→−2.4, flat past ~17) — spec's counter staircase confirmed.
- `d_hi`: `alpha`↑; `force_tol` ramp↑; **`Ldiag_min` ramp↑** (−0.65 near-singular → +2.07 well-cond, the
  λ_min-shrink ceiling as a shape); `hdiag_min`↑; **`p_res0_dual` ramp↓** (+2.4→−1.7, the dual-side
  residual crossing — what v3's dual split was built for).

**Physical reading:** both sensors = α-geometry + tolerance ramp + residual ramp; **ceiling adds
conditioning (`Ldiag_min`/`hdiag_min`) + dual residual**, floor adds the **pbase staircase**. Knots for
Stage C: force_tol≈−5, Ldiag_min≈−2, pbase saturates≈17. Models pickled in `ebm_ckpt/full_{d_lo,d_hi}.pkl`.
