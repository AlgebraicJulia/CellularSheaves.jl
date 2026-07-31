# The four window-edge formulas

Self-contained reference. Everything here is computable from the CSVs in
`src/IPM/examples/oracle2/` alone. Each CSV row is one candidate α tried during one
solver iteration of one problem; `*_ipm.csv` and `*_hsd.csv` name the solver.

## 1. What the formulas estimate

At each solver iteration there is a window of usable α values: below its lower edge
("floor") the inner solve becomes too expensive; above its upper edge ("ceiling") the
arithmetic breaks down. The four formulas estimate, from quantities measurable at the
end of one iteration, the signed log-distance to each edge:

    d_lo = log10(alpha) - log10(floor)     (decades above the floor; negative = below it)
    d_hi = log10(ceiling) - log10(alpha)   (decades below the ceiling)

Ground-truth labels, from the CSVs: within one (file, iter) group, take rows with
`state == 0` (successful solves); let m = the minimum `ncraig` among them; the window is
the set of those rows with `ncraig <= m + 1`; floor and ceiling are the min and max
`alpha` of that set. Labels are quantized to the sweep's half-decade α grid, so ~0.25
decades of error is unavoidable noise.

## 2. Variable dictionary (formula symbol <- CSV column)

All logs base 10. Guard zeros: log10(max(value, 1e-300)).

    la   <- log10(alpha)
    ft   <- log10(force_tol)          per-iteration solve tolerance
    fl   <- log10(floor_tol)
    mg   <- max(log10(r0_c), log10(r0_p)) - max(ft, fl)        "margin"
    cd   <- log10(c_res0_dual)        corrector's dual residual, pre-refinement
    pd   <- log10(p_res0_dual)        predictor's dual residual, pre-refinement
    hx   <- log10(hdiag_max)          max diagonal of the scaled barrier Hessian
    bm   <- log10(bar_hdiag_med)      median diagonal of the barrier Hessian
    Ld   <- log10(Ldiag_min)          min diagonal of the Cholesky factor of F
    pp   <- ppass                     predictor refinement pass count (integer)
    rho_w<- (log10(r1_w) - log10(r0_w)) / max(wbase + craig_w, 1)   [hsd only]
    relu(x) = max(x, 0)

Centered variables (used by two formulas):
    la' = la - 9     ft' = ft + 5     pd' = pd + 7.5     Ld' = Ld + 4.5

## 3. The underlying law (one law, both solvers)

Define N = r0 * alpha (measured: r0 is proportional to 1/alpha within any sweep, slopes
-0.98/-0.99 per solver, so N is a property of the iterate, not of where you measure).
Both solvers' floors obey

    alpha_floor  ~  c * sqrt(N / eps)          (eps = force_tol)

Evidence: the sqrt (exponent ~0.5) appears in seven independent fits; within a solve,
both solvers' floors move with Delta(log N - log eps) at slope ~0.5 (ipm +0.66,
hsd +0.39). The solvers differ only in degree: hsd's renormalization damps N's
variance (within-solve 2.2 vs 3.2 decades; cross-problem 1.25 vs 1.60) and its
floor-coupling (corr 0.51 vs 0.82), and shifts the constant c by ~+0.46 decades.
Fitted as one template  d_lo = C + a*la + b*mg + g*relu(ft+6.42):

    solver    C       a(la)   b(mg)   g(ft)    CV
    ipm      -2.44    0.82    -0.30   0.04     1.59
    hsd      -3.73    0.99    -0.03   0.44     0.83

Same slots; hsd's margin coefficient is small, not absent. The deployment formulas
below are these rows sparsified plus per-cell correction terms. Ceilings share one
law outright (the ipm ceiling formula transfers to hsd verbatim).

## 4. The four deployment formulas

### 4.1 IPM floor ("the committee" — mean of two independent derivations)

    d_lo = -1.2
         + 0.5    * (la - mg)
         - 0.4    * relu(bm - 2)
         + 0.25   * relu(hx - 3)
         + 0.2    * relu(cd + 3.5)
         + 0.0125 * relu(cd + 10) * (la - 1.5*bm - 10)
         + 0.05   * relu(cd + 10) * relu(mg + 1)
         - 0.0125 * hx * bm

Reading: half the perch-adjusted gap (the law's sqrt), a curvature penalty past
bm = 2, a residual-scale credit past hx = 3, and a second altitude reading through
the dual residual (level + gated block + a term requiring dual and primal to agree).
Accuracy (fixed formula on the fleet): pooled 1.06; below 1.32 / in 0.78 / above 1.25.
Honest transfer estimate (coefficients refit excluding the tested family): ~1.31.

### 4.2 IPM ceiling (native 8-term; centered variables)

    d_hi = -1.15
         + 0.9   * relu(ft + 10)
         + 0.7   * Ld'
         - 0.4   * pd'
         - 0.4   * relu(pd + 10)
         - 0.07  * Ld' * la'
         + 0.05  * Ld' * ft'
         + 0.05  * ft' * la'
         + 0.025 * ft' * pd'

    with centered variables  la' = la - 9, ft' = ft + 5, pd' = pd + 7.5, Ld' = Ld + 4.5

Reading: conditioning (Ld) buys headroom; the dual residual is the altimeter
(altitude read through pd, self-calibrating per problem); tolerance sets the scale.
Accuracy: CV 0.93 exact / 0.91 fixed-artifact in this cute form — beats the reference
forest (1.08); zones 0.96 / 0.71 / 1.19. Applied to hsd rows (section 4.4) the cute
form scores 1.08. Note: this cell tolerates only gentle rounding (one decimal); coarse
rounding of the slopes costs ~0.4 because they multiply wide-range variables.

### 4.3 HSD floor

    d_lo = -3.7 + la + 0.5 * relu(ft + 6)

(Cute-rounded form; it slightly beats the raw fit — 0.766 vs 0.771 — and the raw
fit's five tiny interaction terms are dropped entirely, confirmed negligible.)

Reading: the same law with the margin term's pooled weight ~0 (see section 3); what
remains is the alpha-geometry (slope 1 = distance to a location that doesn't move with
the iterate at pooled level) and the tolerance shift at the law's half-exponent
(0.47). Accuracy: CV 0.77 (0.68 ex-X04); zones 0.88 / 0.64 / 0.87.
Frozen upgrade using the w-role's observed contraction rate (columns r0_w, r1_w,
wbase, craig_w):

    d_lo = -4.1 + la + 0.5*relu(ft + 6.4) - 0.4*rho_w + 0.02*rho_w*la
    (CV 0.79; earlier fit variants reached ~0.73 — treat 0.73-0.79 as the range) Dynamics note: this floor drifts within
a solve, tracking Delta(log N - log eps) at slope ~0.4; a one-row formula cannot use
that, but a controller carrying state across iterations can.

### 4.4 HSD ceiling

Use formula 4.2 unchanged on hsd rows. Measured: 1.11 vs hsd forest 0.83; a native
search found no stable compact formula that beats the transfer. The ceiling mechanism
is solver-shared; the residual gap is diffuse (no compact term captures it).

## 5. Caveats that travel with the formulas

- X04-family problems (near-dependent constraint rows) break the r0 ~ 1/alpha law
  the floor formulas rest on; all sensors degrade there (~2x). Quote ex-X04 numbers
  alongside pooled ones. The family's fix is structural (presolve), not sensing.
- Below the floor, both floor formulas run optimistic (report "safer than you are"
  by up to ~0.9 decades); a consumer should treat negative d_lo as a bound, not a value.
- Validity region: the fleet's data occupies a thin correlated manifold in feature
  space; outside it (unusual states), trust the spine terms (the section-3 template)
  and discount correction terms — they are fit only where data exists.
- Evaluation protocol behind every number: folds grouped by problem FAMILY (strip
  dial/seed suffixes, e.g. X04_e-7_s1 -> X04); problem-grouped folds leak sibling
  dials and inflate accuracy ~2x. Any re-evaluation must group by family.
  Exact fold assignment used for all quoted CV numbers (family: fold):
  06:3 07:2 13:3 SOS:2 X03:3 X04:1 X06:4 X08:3 X09:0 e01:2 e02:1 e03:4 e04:4
  e05:2 e08:4 e09:4 e10:3 e11:3 e12:4 e14:4 e15:2
- Reference forest recipe (for the "beats/trails the forest" claims):
  sklearn HistGradientBoostingRegressor(max_iter=250, max_depth=4,
  learning_rate=0.1, random_state=0), fit per solver-cell on that cell's essential
  features (4.1: la, cd, hx, bm, mg; 4.2: Ld, ft, pd; 4.3: la, ft, pd, Ld, pbase),
  scored under the fold assignment above.
