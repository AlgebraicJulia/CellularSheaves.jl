# Distance-sensor formula discovery — findings log

Running record of the `distance_formula_discovery_spec.md` pipeline. Data: `oracle2/*.csv` +
`*.meta.json` → `dataset.pkl` (via `build_dataset.py`). Targets: `d_lo` = decades above the window
floor, `d_hi` = decades below the window ceiling. Evaluation: grouped 5-fold CV by problem
(`stage_a.py` etc.). Noise floor ≈ 0.2–0.3 decades (half-decade grid quantization) — "solved".

**Dataset:** 76,590 rows · **96 distinct problems** · 188 files (ipm 40k / hsd 36k). Zones: in 33k /
above 31k / below 12k. 6 adversarial families (X09/X04/X06/X08/SOS + fleet) × dials × seeds × 2
solvers. (Spec was written against 23 problems — fingerprint/overfit worries are far weaker here.)

---

# PHASE 1 — IPM / floor sensor (per-solver, canonical folds)  ✅ FROZEN

> **This is the current, authoritative result.** The Stage A–D sections below are earlier
> mixed-solver exploration (pre-pivot) — kept for history. Phase 1 rebuilds the floor (`d_lo`)
> sensor for the **IPM solver only**, on **leak-free canonical folds**, and freezes a legible formula.

## The evaluation is the artifact — `canonical_folds.json`

A week of this project taught us that *evaluation details silently decide everything*, so the folds are
committed, not recomputed. **Canonical family rule:** `family = problem.split('_')[0]`; if it matches
`X<digits><suffix>`, strip the suffix to the X-base (`X09bt`→`X09`, `X04b`→`X04`); else verbatim. This
collapses dial/twin variants so a family never leaks across folds via its near-twin → **21 families**
(was 31 with the naïve split). Balanced deterministic 5-fold map (no RNG) in `canonical_folds.json`.

- **X09 is 37% of all IPM rows** (14,911) — the merged ceiling-generation family. Whichever fold holds
  it dominates, so any 5-fold family-CV over 21 families is inherently high-variance (~±0.05).
- **Handshake (forest-5 on canonical folds):** in-sample **0.512**, family-CV **1.109**. (Earlier naïve
  31-family folds gave 0.936 — the twin-leakage was worth ~0.13 MAE. This is why we canonicalized.)
- **All Phase-1 numbers are quoted against these folds**, means with per-fold spread where it matters.

## Essential features — RFE redone honestly (the "8 vs 5" story)

The five `{log α, c_res0_dual, hdiag_max, bar_hdiag_med, margin}` were originally distilled from the
**mixed** dataset via a **problem-level 80/20** RFE — which is *leakier* than the 31-family folds (every
family scatters ~80/20 across train/test, so held-out problems always have near-twins in train). We
**redid the RFE on canonical folds** (`phase1_rfe_canonical.py`):

- **The greedy knee moved 5 → 8.** Honest curve is flat ~1.10–1.15 down to **n=8 (1.108)**, then breaks
  hard at n=7 (+0.19). Honest-8 set = the five **plus** `r0_c`, `tolratio`, and — notably — a counter,
  `cbase`. (The "no counters in this cell" note was itself a leaky-fold artifact.)
- **But the five are NOT refuted.** Greedy RFE is *path*-unstable in its tail: on canonical folds it drops
  `hdiag_max` too early and lands on a worse 5-set (1.308). Meanwhile **forest on our exact five = 1.109
  = forest-8 (1.108) = forest-25 (1.150)**. So the five are a fully valid essential set; there is a broad
  **5–8 feature plateau** (~1.11), and greedy just picks different members under different folds.

**Then we redid L0 with the honest eight** (`phase1_L0_eight.py`) to see if `r0_c`/`tolratio`/`cbase`
buy a better *formula*. **They do not** — the frozen stable-core is a **wash**:

| frozen stable-core (canonical CV) | terms | all-families | ex-X04 |
|---|---|---|---|
| 5-feature basis  | 7 | **1.204** | 0.944 |
| 8-feature basis  | 9 | 1.211 | **0.929** |
| off-the-dome champion (hand-built) | 7 | 1.320 | 1.048 |
| forest(5) reference | — | 1.109 | 0.809 |

The extra three are **near-redundant with the five, on BOTH models** — not merely formula-diffuse. Direct
measurement: **forest(5) = 1.109 vs forest(8) = 1.108** (the RFE's n=8 row is the forest on exactly the
five + `r0_c`/`tolratio`/`cbase`) — Δ 0.001, the three add nothing to the forest. And in the formula they
only reach microscopic terms (`cbase` → one 23%-active triple, coef **0.003**; `tolratio` → weather
triples; `r0_c` → a ≤0.013 touch-up). So they fail on both fronts: no forest-usable signal (5≈8) and no
formula-usable signal (5-core 1.204 ≈ 8-core 1.211). This is stronger than the classic "forest-essential
but formula-diffuse" case (where a feature helps the forest via interactions a linear form can't reach —
`σ²min` earlier was arguably that); these three don't even earn their keep in the forest. Per-fold they
*help* X09 extrapolation (1.23→0.99) but *hurt* X04 (2.6→3.6–4.1) — net wash. **8 features are no better
than 5 for either model**, so we freeze the simpler, more parsimonious five. (The RFE's "knee at 8" was a
greedy-elimination artifact — it dropped `hdiag_max` too early and propped the curve back up with these
near-equivalent substitutes; there is a broad 5–8 feature plateau at ~1.11.)

## How the formula is fit — L0 + the leakage correction

Basis: order-≤3 over the five, **centered at operating medians** (so main effects keep physical slopes
and interactions encode deviations — coefficients readable and fold-stable), hinges at **EBM bend
points** (`relu(f − knot)`, activity-% reported), pairwise + triple products. Selection: **greedy-forward
+ swap**, L0 (pick exactly *k* terms, no shrinkage — beats lasso's shrink-hedging tax, 1.58→).

**The correction that mattered:** an early single-support run scored **0.955** and "beat" the champion —
but that was inflated by *both* the easy 31-family folds *and* **support-selection leakage** (terms chosen
on inner folds drawn from all families, incl. test). The **nested** run (`phase1_nested.py`, select on
each fold's *training* families only) is honest, and it **overturned** that: per-fold greedy selection
generalizes to only **1.40** (loses to the champion's 1.32). Stability audit: only `la`, `mg` are
selected 5/5; every interaction term is ≤3/5 ("weather"). **The fix:** don't re-select per fold — *freeze
the terms that were stable across folds* (5/5 + 3/5), refit only coefficients. That fixed structure is the
deliverable and it honestly beats the champion.

## FROZEN FORMULA (IPM, floor; canonical CV MAE all **1.204** / ex-X04 **0.944**)

```
d_lo ≈ +5.511
       +0.529·la              (la = log₁₀α − 9.0)          [dominant: drop-Δ +0.33]
       −0.485·mg              (mg = margin + 3.86)         [dominant: drop-Δ +0.42]
       +0.025·relu(cd−10.52)·la                            [real:     drop-Δ +0.13]
       −0.034·relu(cd−10.52)·bm                            [real:     drop-Δ +0.09]
       −0.023·hx·bm           (hx = log hdiag_max − 1.63)  [minor:    drop-Δ +0.05]
       −0.003·cd·hx·mg        (cd = log c_res0_dual + 7.27)[near-dead: drop-Δ +0.035]
       −0.001·la·cd·bm        (bm = log bar_hdiag_med+0.87)[near-dead: drop-Δ +0.011]
```
`margin = log₁₀ max(r0_c, r0_p) − log₁₀ max(force_tol, floor_tol)`. All features centered at the medians
shown (medians: la 9.0, cd −7.27, hx 1.63, bm −0.87, mg −3.86). Output = signed decades above the window
floor (positive = safe; negative = below floor by |d|). Coefficients are all-data OLS; per-fold ranges are
tight (e.g. `mg` ∈ [−0.525, −0.470]). The last two terms are droppable (→ effectively a 5-term formula).

## Audits (all mandatory, on canonical folds)

- **By zone:** in-window **0.814**; **below 1.338 (bias +0.88 — optimistic near the floor)**; above 1.559.
  The formula is a good in-window sensor and degrades in the extrapolation zones — the "below" optimism
  (it says *safer than you are*) is the one caveat to flag for a controller.
- **By log α:** gently rising, 1.09 (low α) → 1.32 (high α) — worst where conditioning is hardest. Not
  flat, but a ~0.23 slope, explainable.
- **Per-family:** most families 0.24–1.04 (X09 itself only 1.038). **Two extrapolation holes dominate the
  pool: X04 = 2.67 and e15 = 3.49.** ex-X04 = 0.944; excluding e15 too would drop it further. X04 is the
  known adversary (routes to the presolve track); **e15 is a newly-surfaced hard family worth a look.**
- **Drop-test:** `mg`(+0.42) > `la`(+0.33) > `relu(cd−10.52)·la`(+0.13) > `relu(cd−10.52)·bm`(+0.09) >
  `hx·bm`(+0.05) > `cd·hx·mg`(+0.035) > `la·cd·bm`(+0.011). The physics: **α-geometry + residual-vs-
  tolerance margin** carry the formula; the `c_res0_dual`-hinge × {α, barrier-median} interactions add
  real signal; a 2-term tail is noise.

## Verdict + gap-to-forest

Decision rule (beat champion at ≤12 terms, tiebreak to fewer): **met** — 7-term formula 1.204 < champion
1.320. **Gap to forest: 0.095 (all) / 0.135 (ex-X04)** — within the 0.15 freeze band. That residual gap is
the honest **price of legibility**: the forest mines the diffuse `cbase`/`tolratio`/`r0_c` signal across
hundreds of rules; no compact term captures it (confirmed by the 8-feature wash). **Phase 1 frozen.**

**Artifacts:** `canonical_folds.json` · `ebm_ckpt/ipm_floor_forest5.pkl` (THE reference, provenance
inline) · `phase1_frozen_formula.json` (formula + audits) · scripts `canonical_folds.py`,
`make_canonical_forest.py`, `phase1_rfe_canonical.py`, `phase1_L0.py`/`phase1_L0v2.py` (single-support,
pre-canonical — superseded), `phase1_nested.py`, `phase1_L0_eight.py`, `phase1_freeze_probe.py`,
`phase1_freeze.py`, `rfe_pure_ipm.py` (leaky, superseded). Deprecated forests → `ebm_ckpt/*.DEPRECATED`.
**Next: Phase 2 (HSD/floor by transfer ladder).**

---

## Stage A — forest reference + feature discovery  ✅  _(earlier exploration; mixed-data, pre-pivot — superseded by Phase 1 above)_

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
