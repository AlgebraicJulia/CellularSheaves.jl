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

---

## Stage C — sparse SINDy-style formula  ✅

Library: survivors + physics composites (margin, logN, κ=Ldmax−Ldmin, ceilmargin=Ldmin+α) + hinges at
quantile AND EBM-indicated knots (force_tol≈−5, Ldiag_min≈−2) + counter indicators `[c≥k]`. Lasso *path*
→ pick α for target sparsity → OLS-refit → grouped-CV. _(Bug fixed: LassoCV given a one-shot generator cv
underregularized to 49 terms/MAE 2.05; switched to path + explicit sparsity.)_

**8-term formulas (grouped-CV MAE d_lo=1.17, d_hi=1.02):**
- `d_lo ≈ −0.29 + 0.70·logα − 0.94·[pbase≥3] − 0.36·[pbase≥2] + 0.27·(force_tol+4.6)₊ + 0.11·p_res0_dual − …`
  → α-geometry + **saturating pbase staircase** + tolerance hinge + residual. The floor law, as written.
- `d_hi ≈ +5.11 − 1.02·stall − 0.61·logα + 0.60·(Ldiag_min+4.4)₊ + 0.33·(force_tol+7.5)₊ + 0.26·(hdiag_min+7.1)₊ − 0.21·p_res0_dual + …`
  → **conditioning hinges (Ldiag_min, hdiag_min) + dual residual + tolerance** = the λ_min-shrink ceiling,
  exactly the EBM/Stage-A ingredients.

**Sparsity/accuracy frontier:** `d_lo` improves monotonically to ~0.92 @29 terms (≈ linear baseline).
**`d_hi` caps at ~0.99 (@12 terms) then OVERFITS** (collinear conditioning-hinges destabilize OLS: 1.79
@20). Both plateau **above** EBM (0.73) / forest (0.62): the additive main-effect form gets the physics
right but can't form **interactions** — which the ladder proved the ceiling needs (linear 2.85 → forest
0.62). → Stage D (PySR, with ×/÷) is the route to close the gap: discover `residual/tolerance`, `α·λ_min`.

---

## Stage D closeout — enumeration verdict + final formulas

**Reframe (A′):** A′ showed the missing structure is *known-type, unknown-instance* (interactions among
~15 named quantities), so the right tool is **fast linear enumeration**, not PySR genetic search. We ran
the enumeration menu; PySR was demoted and ultimately **skipped** (no compact term of any order to find).

### Gallery (grouped-5-fold CV MAE, decades; composite-augmented slate unless noted)

| model | d_lo | d_hi | notes |
|---|---|---|---|
| depth-4 forest        | **0.542** | **0.620** | accuracy ceiling; diffuse (100s of rules); pickled |
| additive EBM (GAM)    | 0.741 | 0.696 | smooth per-feature curves, interactions=0; ceiling past ≤0.75 target; pickled |
| EBM interactions=15   | 0.783 | 0.721 | *pairwise ceiling* (3-fold, handicapped) ≈ additive → pairwise adds ~nothing |
| RuleFit (3-way rules) | 0.760 | 0.957 | floor: 3-way helps (0.76<0.94); **ceiling: rules WORSE than smooth (0.96)** |
| product-SINDy (linear)| 0.939 | 0.935 | ~20 terms; linear products don't help |
| Stage-C 8-term        | 1.170 | 1.024 | the readable compact formula (stable) |
| dense-hinge distill   | 0.908* | 0.905* | UNSTABLE (collinearity): d_lo blows to 2.22@45t, d_hi ±31 coeffs — do not use |

### Verdict: the gap is not a compact term

Four independent tools agree there is **no compact interaction of any order** to discover:
- **H² < 0.07** on every candidate pair (incl. the physically-nominated `α×hdiag_min`, `α×Ldiag_min`).
- **i15 pairwise EBM ≈ additive EBM** (pairwise adds nothing beyond composites).
- **product-SINDy plateaus at 0.94** (no linear product helps).
- **RuleFit**: ceiling *worse* with rules (smooth wins); floor 3-way rules reach 0.76 ≈ EBM but the
  forest's 0.54 needs *hundreds* of rules → **irreducibly diffuse**, no compact 3-way term.

So: **ceiling = smooth additive function** (solved at the EBM level, 0.70, past the ≤0.75 target);
**floor = smooth mains + diffuse residual** (EBM/RuleFit ~0.75, forest 0.54 needs the full ensemble).
The compact-formula frontier is ~0.94–1.0 for both; EBM-parity is not *compactly/stably* distillable
(dense hinges are collinear → unstable).

### Physics finding — σ²min vindicated, as a coordinate not an interaction

`u = logα + logσ²min` (= log(α·σ²min), the Schur-side conditioning `κ_S ∝ 1/(α·σ²min)` from the
α-window theory) is a **top-5 feature** and drops the ceiling EBM 0.73→0.696. σ²min looked "useless"
the whole exploration because it's **marginally ~0 but multiplicatively essential** — every prior
analysis was marginal; the FAST/interaction lens is the first that could see products. But it enters
**additively via `u`**, not as an interaction (H²≈0). The `u` shape is monotone, sign-flips floor(↑)/
ceiling(↓), zero-crossing near **u≈5** (theory predicts u≈0 — a ~5-unit offset, likely the σ²min-
estimate / atol normalization; **open question below**).

### FINAL FORMULAS (stable, readable; grouped-CV MAE d_lo=1.17, d_hi=1.02)

```
d_lo ≈ −0.29 + 0.70·logα − 0.94·[pbase≥3] − 0.36·[pbase≥2]
       + 0.27·(log force_tol + 4.6)₊ + 0.11·log p_res0_dual − 0.09·margin
       [α-geometry + saturating pbase staircase + tolerance hinge + residual]

d_hi ≈ +5.11 − 0.61·logα + 0.60·(log Ldiag_min + 4.4)₊ + 0.26·(log hdiag_min + 7.1)₊
       + 0.33·(log force_tol + 7.5)₊ − 0.21·log p_res0_dual − 1.02·stall
       [conditioning hinges (λ_min-shrink) + dual residual + tolerance]
```
For ~0.94 add the composites (`u`, `tolratio`) as extra linear terms (product-SINDi n=12); for 0.70
use the pickled additive EBM; for 0.54–0.62 use the pickled forest. margin = log r0 − log max(tol).

### Reproducibility / artifacts (push these to git)

- `formula_discovery/*.py` — build_dataset, stage_a, stage_b, ebm_run, stage_c, product_sindy, i15_ebm,
  hstat, rulefit, distill, stage_ab_prime, stage_d0_pairs (pair-nomination), stage_d_resid (PySR, unused).
- `formula_discovery/dataset.pkl` — the 76,590×~55 training matrix (96 problems). Rebuild: `build_dataset.py`.
- `formula_discovery/ebm_ckpt/*.pkl` — additive-EBM full-fits (`full_{d_lo,d_hi}.pkl`), CV folds, and the
  depth-4 forests (`forest_{d_lo,d_hi}.pkl`, dict: model/features/target). Load: `pickle.load`.
- venv (not committable): `python -m venv v && v/bin/pip install pandas scikit-learn lightgbm interpret pysr`.
- Oracle CSVs the dataset is built from: `examples/oracle2/*.csv` (+ `.meta.json`).

### Open questions for a coworker

1. **The u≈5 offset.** Theory says the `u`-shape should bend at `u = logα+logσ²min ≈ 0` (where
   α·σ²min ~ O(1)). Measured zero-crossing ≈ 5. Derive the offset — is it `log(1/atol)`, the σ²min-
   estimate normalization (10-step Lanczos vs true), or the `κ_S` prefactor? A clean derivation would
   turn `u` into a physically-scaled sensor.
2. **The diffuse floor gap.** Forest 0.54 vs additive/RuleFit ~0.75 on d_lo, with H²<0.07 and no compact
   3-way term. Is it genuinely irreducible (many tiny effects), or is there a *change of coordinates*
   (like `u` was for σ²min) that makes it additive? The `hdiag_min×hdiag_max` pair (H²=0.069, the largest)
   is a lead — some conditioning combination.
3. **Would more oracle problems help?** 96 problems; the below/above zones are extrapolation-heavy. Does
   the diffuse gap shrink with more diverse ceilings (e.g. the X-family swept finer), or is it a genuine
   one-step-observability limit?
4. **X04 policy** (not chased here): X04-dial trajectories are information-limited at one-step range
   (measured: forest sensors 0.50 MAE still fail to control them). Formulas are for the sensible cells;
   X04 routes to the presolve track. Confirm/refute on the fleet-cost harness.
