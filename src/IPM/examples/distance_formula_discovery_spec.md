# Distance-sensor formula discovery — self-serve spec

Goal: from the oracle CSVs alone, learn formulas for two sensors, with the machine—not the
analyst—selecting what matters ("maximally learned"):

- d_lo = log10(alpha) − log10(window lower edge)   (signed decades above the floor)
- d_hi = log10(window upper edge) − log10(alpha)   (signed decades below the ceiling)

Everything below is reconstructible from `src/IPM/examples/oracle2/*.csv` + `*.meta.json`.

## 1. Rows and labels

One training row per (iteration, candidate alpha) — every sweep row, not just chosen ones.
For each iteration: candidates with `state == 0` are successful; let mincraig = min ncraig
among them; the window = successful candidates with `ncraig <= mincraig + 1`; its edges are
the min and max alpha of that set (log10). Skip iterations with no successful candidate.
Labels: d_lo and d_hi as above, per candidate. Zone per row: below (alpha < lower edge),
in, above. Keep all zones in training.

Known label noise floor: the sweep grid is half-decade spaced, so edges are quantized;
~0.2–0.3 decades MAE is "solved". Do not chase below it.

## 2. Feature slate — hand the model everything

Include every per-row observable, log10-transformed if positive-continuous (guard with
+1e-300), raw if integer/flag:

- residuals: r0_c, r0_p, r1_c, pres0, cres0, p_res0_dual, p_res0_prim, c_res0_dual,
  c_res0_prim, c_res_exit, p_res_exit (and w-role variants where present)
- tolerances: force_tol, floor_tol — SEPARATELY (do not pre-combine; the max() composite
  measurably destroys information)
- structure: rho, hdiag_min, hdiag_max, bar_hdiag_med, bar_hdiag_frac_mid, Ldiag_min,
  Ldiag_max, mu, alpha
- counters (raw integers): cbase, pbase, cpass, ppass, ncraig (+ w-role where present)
- flags: any status containing REFINE_STALLED -> stall (0/1); is_hsd from filename
- meta (flag for audit, see §5): m, n, ng (log10(1+ng)), alpha_anchor

EXCLUDE (label leakage / oracle artifacts only): `chosen`, `score`, `state` used as a
feature for the *distance* targets is fine to include but note it near-encodes the zone;
anything derived from other candidates in the same sweep is forbidden (a deployed sensor
sees one row).

Transform note: tree models are invariant to monotone transforms (we verified: identical
MAE logged vs raw), so logging matters only for the linear/formula stages. Log everything
once and use one matrix everywhere.

## 3. Evaluation protocol (frozen before any fitting)

- Grouped cross-validation: folds split BY PROBLEM (prefix before first underscore);
  a problem never appears in both train and test. 5 grouped folds is fine; leave-one-
  problem-out if you have patience.
- Metric: MAE in decades, reported pooled AND sliced by zone (below/in/above). The
  operational slice is below+in for d_lo, in+above for d_hi — quote it, but report all.
- Every formula is re-scored under this protocol (refit per fold). Never compare a
  full-fit MAE against a CV MAE.
- Fix fold assignment once (seeded) and reuse for every method.

## 4. Pipeline

Stage A — forest reference + feature discovery.
Gradient-boosted trees (LightGBM / HistGradientBoosting, depth ~4) on the full slate.
Selection uses three instruments, not plain permutation ranking (validated on our data):
1. CLUSTER-LEVEL permutation: group features with |corr| > 0.8, permute whole clusters,
   score damage per cluster. Mandatory — individual permutation hides redundant families
   (we measured members scoring ~1% individually inside a family scoring +63%, and the
   dominant family at +247% that no individual number revealed). Keep every cluster with
   material damage; choose representatives later, or pass whole clusters onward.
2. Boruta-style shadows: append a permuted copy of every column, refit; a feature
   survives only if its permutation damage beats the best shadow's. This is the
   "all-relevant" set — it retains redundant family members, which Stages C/D want as
   alternative parametrizations.
3. Stability: repeat over 3 seeds; report the stable set and flag wobbling seats.
4. Optional but measured-best on our data: model-X knockoffs (use a proper package, e.g.
   knockpy; even a crude Gaussian implementation produced the best selected set we tested,
   with FDR control for free). Stability selection (subsampled-lasso frequencies) is a
   cheap competitive alternative on log-space features.
Do NOT use univariate filter screens (mutual information, mRMR) here: this problem's
signal lives in differences/ratios (residual vs tolerance), which have near-zero marginal
relevance — both filters dropped force_tol entirely and scored ~0.3 MAE worse than every
other selector. Measured contraindication, not taste.
Record CV MAE on the SELECTED set, not the full slate — selection improves forest
accuracy here (~0.07 MAE on our fleet: idle features are overfitting surface with only
23 problems), so the selected-set forest is the true accuracy ceiling for the gallery.

Stage B — capacity ladder + shapes.
Same data: depth-4 vs depth-2 vs depth-1 (additive) vs linear. Gaps tell you interaction
content vs 1-D nonlinearity content. Then fit a real EBM (`pip install interpret`,
ExplainableBoostingRegressor) on the survivors and READ the per-feature shape functions:
integer features give per-value tables (expect a graded, saturating staircase in the
counters), continuous features give curves whose bends locate your knots. EBM's
round-robin fitting shares load fairly across correlated features — trust its shapes over
one-shot boosted stumps.

Stage C — SINDy-style sparse selection.
Build a term library: every survivor; pairwise products of continuous survivors; hinge
terms (x−k)+ at a few quantile knots and at EBM-indicated bends; indicators [counter>=k]
for k=2..5 and ramps (counter−2)+; composites as ADDITIONAL menu items (margin =
log r0 − log max(tols), log(force_tol/floor_tol), log N = log r0 + log alpha) so raw and
composite parametrizations compete. Lasso across a small alpha grid, then OLS-refit on
the selected support and CV-score the support (lasso shrinkage biases coefficients).
Optional upgrade: L0 subset search (greedy/evolutionary over term subsets, OLS per
candidate, complexity penalty) — tends to pick ONE parametrization cleanly where lasso
hedges across redundant ones; we measured the lasso hedge worth ~0.05 MAE and the L0
formula 2x shorter.

Stage D — PySR (you have Julia; this is the step we could not run).
Feed the Stage-A survivors (logged), targets separately. Suggested config: binary
operators +, −, *, /; unary none to start (everything is already logged; add exp/log
later only if the Pareto front stalls); maxsize ~25; niterations 200+; use a validation
split by problem for model_selection, or export the Pareto front and CV-score each
expression yourself under §3 (recommended — PySR's internal selection is not grouped).
Two runs worth doing: survivors-only, and survivors-plus-raw-tolerances to see whether
it reinvents margin/tolerance-ratio on its own. Integer counters: pass them raw; if PySR
ignores them, add [c>=k] indicators as input columns (SR frameworks are poor at learning
thresholds from integers).

Stage E — the gallery.
Score every formula under §3 and tabulate: method, #terms, pooled MAE, zone MAEs.
Then rearrange all winners into common coordinates (substitute N = r0·alpha, margin =
r0/tol) and compare COEFFICIENTS: our runs found every method agreeing on the physics
(floor-law slope ≈ ½ in ingredient-split parametrizations, ≈ 1 when the alpha-geometry is
absorbed into a residual term; a tolerance-phase term; graded counter staircase; hsd
offset ≈ −0.5) with the MAE spread explained entirely by parametrization economics.
Your job is to confirm, refute, or extend that on the ceiling target especially — our
ceiling formula work is thinner than the floor.

## 5. Pitfalls (each one cost us a session; read before running)

- Fingerprints: m, n, ng, alpha_anchor, bar stats can act as problem-ID lookups with only
  23 problems. Protocol: run each stage with and without them; if a formula's accuracy
  leans on them, it will not transfer. Grouped CV catches some but not all of this.
- mu vs force_tol are collinear by construction (force_tol = min(0.1*mu/mu_1, 0.3));
  including both is fine for trees, poison for interpreting linear coefficients.
- Simpson/pooling: within-iteration, label and residual features are affine in the same
  variable (alpha), so pooled slopes mix within-iteration structure with cross-iteration
  intercept drift. A tight CI on a pooled slope is not evidence it is structural. Fit
  within-group or model intercepts before believing any slope.
- Zone slicing is mandatory: pooled MAE hides that each sensor is accurate on its own
  side of trouble and poor on the far side (which is operationally fine — check, don't
  average away).
- Selection optimism: choosing knots/terms/expressions on the same holdout you report is
  mild cheating; final numbers come from the frozen §3 protocol, selection from a
  different split or inner folds.
- In-sample vs CV: OLS coefficient displays can be full-fit; every MAE must be CV.
- Expression-search variance: rerun PySR / subset search with 2–3 seeds; report whether
  the winner is stable or a basin of equivalent forms.

## 6. Report back

The gallery table (both targets), the winning formulas as written expressions, EBM
counter tables and any shape that disagrees with: floor-law slope, saturating staircase,
tolerance-phase term, hsd offset ≈ −0.5. Disagreements are the interesting part.
