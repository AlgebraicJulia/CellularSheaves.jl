# Maximum-likelihood spline estimation: the first live epigraph lift

A companion to `nonnegative_spline.md`. Same sheaf — Bernstein segment stalks, `C^k`
continuity edges on a caterpillar spine, orthant shape cones — but the objective is a
**log-likelihood**, the corpus's first genuinely non-quadratic cost. Everything the
duality note proved about the epigraph/homogenizer lift (§5 there) has so far been
exercised only by the triangle verification; here the lift is the *builder*. This note
records only the deltas; the Bernstein machinery, shape cones, and exact/orthant
distinction are documented in `nonnegative_spline.md`.

The construction follows Papp & Alizadeh (JCGS 2014), whose shape-constrained spline
framework includes ML density estimation, and the Alizadeh–Papp line on estimating
Poisson arrival rates by maximum likelihood over nonnegative cubic splines — the same
two-variant (exact-SDP vs. linear-sufficient) structure this corpus uses.

---

## 1. Two likelihoods, one assembly

Observations `x_1..x_N` on `[0,1]`, spline `f = Σ P_i B_{i,n}` per piece as before.

**Density** (`mode = :density`). Maximize the log-likelihood of an i.i.d. sample:

```
    minimize    −Σ_j w_j log f(x_j)  +  ½ λr Σ_s (P^s)ᵀ G P^s
    subject to  C^k continuity,   ∫ f = 1,   P ≥ 0
```

The normalization is the same single `∫f = 1` pin edge the regression file already
ships; it is what makes the MLE well-posed (without it the likelihood is unbounded).

**Poisson intensity** (`mode = :intensity`). Events of an inhomogeneous Poisson
process with rate `λ(x)`; the log-likelihood is `Σ_j log λ(x_j) − ∫ λ`:

```
    minimize    ∫ λ  −  Σ_j w_j log λ(x_j)  +  ½ λr Σ_s (P^s)ᵀ G P^s
    subject to  C^k continuity,   P ≥ 0
```

The compensator `∫λ` is **linear** in the control points (`∫ piece = h/(n+1)·Σ P_i`),
so it rides in `c` on the segment stalks — no pin, no auxiliary anything. This is the
cleaner of the two: the inhomogeneity enters only through the wiring-edge rhs.

Everything else is shared, so the builder takes a `mode` flag and the two problems
differ by one pin edge versus one linear term.

---

## 2. The leaf: duality.md §5, executed

Each observation hangs one 3-dimensional cone leaf `(u_j, y_j, t_j)` off the segment
containing it, wired by a single 2-row edge:

```
    φ(x_j)ᵀ P^s − u_j = 0        (the spline value at the observation)
    y_j = 1                      (the homogenizer pin)
    leaf cone:  ExponentialCone = { (x,y,z) : x,y > 0,  y·log(x/y) ≥ z }
    objective:  −w_j on t_j
```

With `y = 1` the cone says `t_j ≤ log u_j`, and minimizing `−Σ w_j t_j` drives the
epigraph tight: `t_j = log f(x_j)` at the optimum. Read against the duality note:
`t_j` is Move 1's epigraph scalar, `y_j = 1` is Move 2's homogenizer slice `λ = 1`,
and both pins fold into the co-boundary through `g` — the rhs `(0, 1)` on the wiring
edge is "the apex broadcasts the constant," verbatim. The quadratic roughness stays
in `Q` (the rule of thumb from the MAP note: lift only the cost whose curvature the
solver cannot already see), so one problem carries `Q` + orthant + exponential cones
simultaneously — the widest cone mix in the corpus so far.

**Coordinate-order warning.** The library's exponential cone is *log-form*, bound
**last**: `(x, y, z)` with `y·log(x/y) ≥ z`, hence the leaf stalk order `(u, y, t)`.
MOSEK/MOI's cone is exp-form, bound **first**: `(x, y, z)` with `y·e^{x/y} ≤ z`, so
the same constraint in the JuMP reference reads `[t_j, 1, u_j]` — the reverse. The
two files keep both conventions explicit; confuse them and the model is silently
wrong. (cvxpy's `ExpCone(t, 1, u)` matches the MOI order.)

**Feasibility is by construction.** The constant spline `f ≡ κ` (all `P_i = κ`,
`κ = 1` for density) satisfies continuity and the pin exactly, sits in the interior
of every orthant, and with `t_j = log κ − 1` is strictly interior to every leaf cone
— a Slater point exactly analogous to the MAP note's uniform-independent belief.
Verified below to machine precision.

---

## 3. Where the cones certify what

The exp leaves certify positivity only **at the observations**; the orthant on the
segment stalk certifies `f ≥ 0` **everywhere** (and remains the shape prior — a
density or a rate must be nonnegative between data points too). So the orthant is
not redundant with the leaves; drop it and the estimate can dip negative in gaps.
The exact Karlin–Studden variant composes: `include` the exact file and replace the
orthant by certificate vertices, leaving the likelihood leaves untouched — the two
constructions touch disjoint parts of the graph.

Sheaf inventory (deltas over the regression instance):

| sheaf datum | MLE object |
|---|---|
| leaf vertex, stalk ℝ³ | `(u_j, y_j, t_j)`, cone `ExponentialCone` |
| wiring edge (2 rows) | `φᵀP − u = 0` and `y = 1`, rhs `(0,1)` |
| leaf objective | `c = (0, 0, −w_j)` |
| (density) pin edge | `∫f = 1` |
| (intensity) segment `c` | `+ s·h/(n+1)` per coefficient (the compensator, normalized §3.1) |

Leaves are degree-1 and eliminate with zero fill; the spine stays block-tridiagonal.
The cost of the likelihood is `N` extra 3-dim vertices and `2N` equality rows —
problem size scales with the *data*, not just the mesh. The `binned` helper caps
this: coalesce observations into bins, one leaf per occupied bin, weight `w_j` = the
count (the weight simply scales `c_t`). Binning is an approximation (evaluation
moves to bin midpoints), quantified below.

### 3.1 Scale the intensity model — a field note

The builders emit the intensity model in **normalized** form: substitute
`λ = s·f` with `s = Σⱼ wⱼ` and solve for `f`. This is a change of variables, not
of model — `F_raw = F_scaled − (Σw)·log s` exactly, verified below — and it
moves every exponential-cone argument to O(1): the fitted `f` has `∫f ≈ 1` (the
Poisson MLE fits total mass to the event count), where the raw `u_j = λ(x_j)`
grows with the sample and sits next to `y = 1` pins and `10⁻³`-scale compensator
coefficients.

Why this earns a section: the **raw** instances drove Mosek (JuMP-bridged; the
bridge is not implicated) into slow-progress stalls at benchmark sizes, while
the sheaf IPM converged. Badly scaled exponential cones are the canonical stall
mechanism for nonsymmetric-cone interior-point implementations, so the
normalized model is the benchmark-fair form for every solver. If a stall
survives normalization, first check `min_j f(x_j)` — events landing in
near-zero-intensity regions push their leaves into the cone's corner; if it
persists after that too, it has graduated from scaling artifact to a genuine
robustness differential and belongs in the benchmark table.

---

## 4. Verification

`mle_spline_oracle.py` (NumPy/CVXPY/Clarabel, standalone) solves the same instance
two ways: **native** — the un-lifted problem via `cp.log` — and **assembled** — the
`(Q, B, c, g, K)` system built index-for-index as `build_mle_problem` builds it,
exp-cone leaves and all. Agreement of the two validates the lift and every index in
the builder. `M = 6` quartic pieces, `C²`, `ρ = 2` roughness, `N ≈ 150` observations
from the bimodal ground truth:

```
[density]  N = 150 observations
  native (cvxpy log)  optimum : -73.18929832
  assembled lift      optimum : -73.18929318     |Δ| = 5.1e-06
  sections agree  ‖P − P′‖∞   : 6.6e-05
  epigraphs tight max|t−log f|: 1.6e-07
  primal feasibility ‖Bp − g‖ : 2.6e-07
  Slater point residual       : 1.6e-15
  min coefficient / ∫f        : -1.4e-08 / 1.000000
  binned (64) optimum         : -73.69689033     |Δ raw| = 5.1e-01

[intensity]  N = 156 observations   (normalized: λ = 156·f)
  native (cvxpy log)  optimum : 84.46984840
  assembled lift      optimum : 84.46985069     |Δ| = 2.3e-06
  sections agree  ‖P − P′‖∞   : 1.3e-05
  epigraphs tight max|t−log f|: 1.9e-08
  primal feasibility ‖Bp − g‖ : 2.7e-08
  Slater point residual       : 1.7e-15
  min coefficient / ∫f        : 3.5e-10 / 0.9535
  binned (64) optimum         : 84.40043822     |Δ raw| = 6.9e-02
  scaled − W·log s  vs  raw   : -703.307686 vs -703.307687   |Δ| = 3.7e-07
  exp-cone args f(x_j)        : [0.215, 2.717]   (raw λ(x_j): [33.5, 423.8])
```

Native and lifted optima meet; epigraphs are tight (`t_j = log f(x_j)`, the §5 lift
closing); the density integrates to one; the intensity's fitted mass
`s·∫f ≈ 148.8` estimates the expected event count, with the normalized and raw
optima agreeing through the `−W·log s` identity (and the section agreement
improving a thousandfold over the raw form — the conditioning gain is real for
Clarabel too, not just Mosek); and 64-bin coalescing (150 → ≤64 leaves) shifts
the optimum by well under 1% while preserving the structure.

---

## 5. Code

- **`mle_spline.jl`** — the builder in the corpus idiom: `include`s
  `nonnegative_spline_lp.jl` for the Bernstein/jet numerics, adds
  `generate_mle_instance` (rejection sampling for `:density`, Poisson thinning for
  `:intensity`), `binned`, `build_mle_problem`, a JuMP reference (`MOI` exp-cone
  order!), `mle_demo`, and a Mosek benchmark harness over raw vs. binned cases.
- **`mle_spline_oracle.py`** — the numerical ground truth quoted above.

Written against the CellularSheaves.IPM PR-67 API; not executed here (no Julia
toolchain). The first Julia run is the real shakedown for the `ExponentialCone`
iteration path — reuse `mle_settings()` as a starting point and expect the exp-cone
central path to want more iterations than the LP instances.

---

## References

- F. Alizadeh, D. Papp. *Estimating arrival rate of nonhomogeneous Poisson processes
  with semidefinite programming.* Annals of Operations Research (2013). *(ML intensity
  estimation over nonnegative splines; SDP-exact and linear-sufficient variants)*
- C. Kooperberg, C. Stone. *A study of logspline density estimation.* CSDA (1991).
  *(the log-linear alternative: model log f, positivity free, normalization inside the
  concave likelihood — a useful contrast to the direct nonnegative-spline MLE here)*
- D. Papp, F. Alizadeh. *Shape-constrained estimation using nonnegative splines.*
  Journal of Computational and Graphical Statistics 23(1), 211–231 (2014).
- R. T. Rockafellar. *Convex Analysis.* Princeton (1970). *(perspective/homogenization;
  the lift of duality.md §5)*
