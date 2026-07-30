# obstacle_option.jl — obstacle problems / American options

The example where Q stops being regularization and **is** the problem.
The objective is the PDE energy (dense-banded per-patch stiffness),
the cells are overlapping subdomain patches — Schwarz domain
decomposition made monolithic — and the orthant cell carries the best
cone semantics in the corpus: w = u − ψ is the option's
**early-exercise premium**, identically zero on the exercise region,
with the free boundary appearing exactly where the cone constraint
activates. No relaxation gap anywhere: after three gap species and
three holonomies, this is the pure-Q, pure-glue flagship.

## 1. Problem and symmetrization

Perpetual American put in log-price x = log S: the Black–Scholes
variational inequality Lu ≥ 0, u ≥ ψ = (K−eˣ)⁺, complementarity, with
L = −½σ²∂ₓₓ − μ∂ₓ + r, μ = r − σ²/2. L is not self-adjoint; the
integrating factor ρ = e^{2μx/σ²} fixes that:
ρLu = −(pu′)′ + rρu with p = ½σ²ρ, and since ρ > 0 the pointwise
complementarity is unchanged. The VI is then the QP

$$\min\ \tfrac12 w^\top A w + (A\psi)^\top w,\qquad w \ge 0,$$

with A = weighted stiffness + lumped mass — **positive definite
natively** (the mass term is the floor), so ε = 0: the first builder
in the corpus with no epsilon anywhere. The obstacle rides entirely in
shifts (c = Aψ and the boundary pins).

## 2. The sheaf

| datum | content |
|---|---|
| patch cells w_p | `PositiveCone`, the premium on P overlapping windows |
| Q_p | element-split stiffness (each interval weighted 1/#owners — the `tensor_spline` splitting) + node-split mass; Σ embed(Q_p) = A exactly |
| c_p | Q_p·ψ\|_p — the obstacle, patch-locally |
| agreement | shared node values on overlaps (the FD-nodal shadow of the spline corpus's jet matching) |
| pins | u(x_lo) = ψ, u(x_hi) = **Merton's value there** — the truncation value at S = 4 is 0.0498, not 0; pinning to 0 costs 5×10⁻² in ∞-norm (measured the hard way) |
| dual | the cone dual on premium cells is the **obstacle reaction measure** μ = Au + c ≥ 0, supported exactly on the contact set |

**Schwarz, precisely.** Lions (1988), Badea (SINUM 1991; 2006 for
nonquadratic energies), and Tai's linear-convergence theory analyze
*iterative* subdomain sweeps whose rate depends on the overlap width.
This builder is the *monolithic* QP rewritten on patches — exact by
construction, and measured overlap-invariant (values agree to
3.5×10⁻¹⁰ across overlaps 1–8). The relationship is: their iteration
converges to our formulation; our formulation quantifies what "the
Schwarz limit" is, and the overlap knob that governs their *rate*
governs nothing here.

## 3. Verification (obstacle_option_oracle.py; K=1, r=0.05, σ=0.3)

```
Merton: S* = γK/(1+γ), γ = 2r/σ²  →  S* = 0.526316 (ODE residual 2.1e-5)
n=801, 4 patches, overlap 4:
  patched ≡ monolithic  : value 2.5e-10, ‖Δu‖∞ 5.0e-7
  ‖u − PSOR‖∞           : 1.1e-5     (Cryer 1971)
  ‖u − Merton‖∞         : 1.2e-6
  free boundary          : S = 0.5300 vs S* = 0.5263  (one grid cell)
  contact measure        : ≥ +1.9e-4 on contact, ≤ 2.1e-7 off
American put T=1 (implicit-Euler chain, 200 obstacle QPs, n=601):
  0.098677  vs  CRR(2000) 0.098694   |Δ| = 1.7e-5
Game option δ=0.02: lower contact 514 nodes, upper 1 (both live)
Minimal surface vs Dirichlet (same bump): ‖Δu‖∞ = 2.8e-7   — identical
H¹ of the patch chain: 0
```

Rigor chain for the finance anchors: Merton's closed form (verified
against the ODE before anything trusted it), Cryer's PSOR as the
independent iterative reference, the CRR tree for finite horizon, and
Jaillet–Lamberton–Lapeyre for the American-option-as-VI foundations.

## 4. Modes

**`:game`** — two-sided obstacle ψ ≤ u ≤ ψ + δ: Israeli/game options
(Kifer 2000). Cells become `CofreeCone` u-patches plus two
`PositiveCone` corridor slack cells each — the spline corridor
pattern — and the Schwarz theory for exactly this problem is
Zeng–Zhou (SINUM **35**:600, 1998). Both contact sets are live at
δ = 0.02.

**`:minimal_surface`** — the classical obstacle problem
min ∫√(1+u′²): one 3-dim `SecondOrderCone` leaf per element
(t_e ≥ ‖(1, Δu/h)‖, sums-of-norms à la Lobo–Vandenberghe–Boyd–Lebret),
linear objective. The honest finding, measured then derived: **in 1D
the contact set cannot distinguish energies** — any ∫φ(u′) gives
straight lines off contact plus the same tangency condition, so the
minimal-surface solution *equals* the Dirichlet one (2.8×10⁻⁷). The
free-boundary difference between energies is a ≥2D phenomenon (out of
scope). The mode's role is therefore a rare kind of test: SOC leaves
and native Q are forced to produce the *same* function through
entirely different cone machinery, and do.

## 5. Time stepping as factorization reuse

The finite-horizon put is a chain of obstacle QPs with
G = M + Δt·A **fixed** across all 200 steps — only c changes
(c = Gψ − Mu^{prev}). For the IPM this is the cleanest
factorization-reuse showcase in the corpus: same Q, same B, same
cones, new linear term, two hundred times.

## 6. Solver notes

The genuinely new stress here is not the cone inventory (orthant,
plus Cofree in `:game`, plus SOC leaves in `:minimal_surface`) but
the **active-set geometry**: a macroscopic fraction of the premium
coordinates sits *on* the cone boundary at the optimum (the exercise
region — 514 of 801 nodes), and strict complementarity degenerates
exactly at the free boundary. Every prior orthant example had sparse
contact. Watch iteration counts vs the spline files; the free
boundary crossing patch seams (printed by `perpetual_put_demo`) is
the sheaf-specific variant of the stress. H¹ = 0 — full-rank regime.

First runs: `perpetual_put_demo()` (one solve, value + boundary +
seams), `american_put_demo()` (200 solves — the timing story),
`game_option_demo()`, `minimal_surface_demo()`.

## References

Merton (1973). Cryer, SIAM J. Control **9**:385 (1971). Jaillet,
Lamberton, Lapeyre, Acta Appl. Math. **21**:263 (1990). Kifer, Finance
Stoch. **4**:443 (2000). Kinderlehrer & Stampacchia. Lions, 1st DD
Symposium (1988). Badea, SINUM **28**:179 (1991); SINUM **44**:449
(2006). Tai (& Tseng, Xu): rate analysis for Schwarz on obstacle
problems. Zeng & Zhou, SINUM **35**:600 (1998). Lobo, Vandenberghe,
Boyd, Lebret, LAA **284**:193 (1998).

## Files

* `obstacle_option.jl` — standalone: `:dirichlet`/`:game` builders with
  time stepping, minimal-surface SOC builder, demos, JuMP reference.
* `obstacle_option_oracle.py` — self-contained; all numbers above.
