# IPM Examples

This directory contains benchmark examples demonstrating the sheaf-structured interior-point solver. Each example is a directory with two files:

- `def.jl` — definitions (structs, builders, gate tests, sweep function)
- `run.jl` — includes def.jl, calls run(), sample results in comments

Usage:
```bash
julia --project e01/run.jl              # Clarabel baseline
julia --project e01/run.jl --mosek      # + Mosek
julia --project e01/run.jl --quick      # quick (smaller) sweep
```

## Overview

| ID | Problem | Habitat | Dial |
|----|---------|---------|------|
| **E01** | Two-asset American basket put, Merton jumps (implicit-Euler chain) | dense-NOC QP | grid size nx |
| **E02** | Quantum marginal compatibility (XXZ ring) | non-chordal SDP | ring length N |
| **E03** | Tensor-spline shape-constrained regression | dense-QP | patch count M |
| **E04** | SOS-certified shape-constrained regression | dense SDP | intervals P |
| **E05** | Galerkin-compressed nonlocal TV | dense-map SOCP | signal length N |
| **E06** | Sqrt-loss spectral regression | dense-map SOCP chain | chain length P |
| **E07** | Poisson-TV photon-limited deconvolution | dense-map exp+SOC+NOC | signal length N |
| **E08** | Optimal execution, 3/2-power impact | dense-Q pow chain | horizon T |
| **E09** | Robust Kalman smoothing, burst-corrupted sensors | dense-Q + SOC chain | horizon T |
| **E10** | Compositional Lyapunov on a K×K torus | dense-Q dense-SDP, sum-decomposition edges | torus size K |
| **E11** | Multichannel FIR equalizer bank (MINT) | dense-Q + SOC, Cholesky ties | channels M |
| **E12** | Distributed sound-zone control, cooperative rendering | dense-Q + SOC, bounded degree | zones Z |
| **E12** | Min-fuel planetary soft landing with tracking | dense-Q SOC chain (LCvx) | grid N |
| **E13** | Gyre-aware robust flow recovery (punctured grid) | full-Hodge dense-map + dense-Q + SOC | grid size K |

## E01: Two-asset American basket put under Merton jump-diffusion

Finite-horizon American basket put, ψ = max(K − (S1+S2)/2, 0), priced by an
implicit-Euler chain: each timestep is a partial integro-differential LCP whose
matrix is dense because the Merton jump integral couples every grid node to
every other (d'Halluin–Forsyth–Vetzal-style PIDE pricing). A 2D integrating
factor makes drift+diffusion exactly self-adjoint, so each step is a genuine
QP; the small non-symmetric jump residue and any jump coupling between nodes
never co-resident in a patch are moved to the linear term and converged by
Picard sweeps (contraction ∝ λΔt; ~4 sweeps/step), so the fixed point solves
the full LCP exactly. The spatial domain is tiled into overlapping patches,
each carrying a nonnegative-orthant block over the early-exercise premium
v = u − ψ ≥ 0; Q is split into per-patch SPD blocks by a PSD element
decomposition (cell/edge/node pieces — the tensor_spline discipline), giving
dense-NOC habitat where general-purpose solvers factor a filled-in system
while this solver factors a small dense block per patch. Gate tests validate
kernel mass, the exact split identity, per-block SPD via Cholesky, European
prices against an in-file exact-simulation Monte Carlo oracle, and the IPM
chain against a Clarabel chain node-for-node. Full derivation and a
self-contained numpy validation in tools/b1_oracle.py.

**Dial:** Grid size nx (controlling patch dimension, hence block size);
benchmark rows time the whole 10-step chain.

## E02: Quantum marginal compatibility (XXZ ring)

Ground-state energy bounds from overlapping local density matrices on an XXZ spin ring. Each ℓ-site window carries a PSD block; consecutive windows overlap by ℓ-1 sites and must agree via partial-trace consistency. Unlike chordal moment hierarchies, partial-trace coupling is dense and non-decomposable — this benchmark isolates the non-chordal argument, showing that speedups persist where chordal decomposition gains nothing. Speedups: 1.6-1.8x vs Mosek, 4-5x vs Clarabel.

**Dial:** Ring length N at fixed ℓ=4.

## E03: Tensor-spline shape-constrained regression

Nonparametric regression with shape constraints (nonnegativity) using tensor-product Bézier surfaces. The domain is tiled into MxM patches, each carrying a control-net coefficient block. Per-patch Grams (data fit + thin-plate roughness penalty) are genuinely dense — the rare QP family where the dense-block condition holds through Q itself, not just cone scaling. Edges impose C^k smoothness via derivative-matching (jet) operators. Speedups: 1.9-2.5x vs Mosek, 2-2.5x vs Clarabel.

**Dial:** Patch count M (MxM grid).

## E04: SOS-certified shape-constrained regression

E03's reference formulation (Papp–Alizadeh 2014) done exactly: where E03 enforces
nonnegativity through a sufficient-condition control net, E04 certifies it — the
fitted piecewise polynomial (odd degree n per interval, C^k at knots) is
nonnegative everywhere by a Markov–Lukács sum-of-squares representation, which
for univariate polynomials is exact, and the oracle verifies the sampled-grid
relaxation converges to this optimum from below. Each interval is a vertex
whose stalk is one merged PSD Gram block of size n+1 (the two Markov–Lukács
Grams on the diagonal; the free cross block is exact by projection and the
optimizer zeroes it). The per-vertex Q is the data Gram plus roughness pulled
back through the Gram→coefficient map — dense-SDP habitat where the dense-block
condition holds through Q *and* the cone simultaneously. Edges are the knots,
carrying C^k jet-matching rows (row-normalized). Local Chebyshev bases keep
the pulled-back Grams well-conditioned at high degree (effective cond 4e2 vs
2.9e15 monomial at n=25, per Papp–Yildiz 2019). Gate tests validate the
representation rank and a Slater point, that the constraint is live
(unconstrained fit dips negative, SOS fit certified ≥ 0), the IPM against a
Clarabel solve of the same conic program, and the objective against a
sampled-relaxation lower bound (sandwich). Self-contained numpy validation in
tools/spline_sos_oracle.py.

**Dial:** Number of intervals P at fixed degree n=13 (block size 14); benchmark
rows time a single solve. This is the favorable regime for the sheaf solver:
many small identical blocks coupled by a thin path graph.

## E05: Galerkin-compressed nonlocal TV

Gilboa–Osher isotropic nonlocal TV, discretized in the dense-map habitat:
u is a piecewise polynomial — non-overlapping tiles, m = 32 local Chebyshev
coefficients each, C² jets at interfaces — instead of a nodal vector. Every
map is now dense **by nature**: point evaluation of a spectral representation
touches all m coefficients, so z-rows, jets, and the per-tile data Grams E'E
are honestly dense and the bloat factor is ~1 by construction, not by decree.
The compression is bought with a genuinely smaller problem (u-DOF ÷4), every
solver receives the same compressed conic program, and z-rows touch at most
two tile stalks (when a neighbor crosses a boundary). Measured statistical bonus: the truncated
basis filters fine-scale noise the TV term leaves behind, so the compressed
estimator *beats* the nodal one (MSE 0.0275 vs 0.0344 vs 0.0501 best local
TV). Oracle rigor: nodal reference with subgradient-selection KKT certificate;
monotone compression family with the Chebyshev resolution threshold appearing
exactly at ~π modes/wavelength; a completeness certificate (full basis, no
jets, reproduces the nodal optimum to +1.3e-6); sheaf assembly vs independent
monolithic solve to 5e-6; and a conditioning bound (equispaced Chebyshev
evaluation requires m ≤ 48). Self-contained numpy validation in
tools/galerkin_nltv_oracle.py.

**Dial:** Signal length N at fixed Tsz = 128, m = 32, K = 12; benchmark rows
time a single solve.

## E06: Square-root-loss spectral regression (whitened SOC chain)

Robust function fitting under block-corrupted noise: P Chebyshev patches with
C² jets, minimizing Σ_p ‖E_p θ_p − y_p‖₂ — norms, not squares. The sqrt loss
bounds each patch's influence (the gradient of ‖r‖ is a unit vector), so
corrupted patches cannot drag the global fit: the oracle's escalation curve
shows clean-region RMSE 0.0217 → 0.0233 for sqrt-loss as corruption grows
×1 → ×128 while least squares blows up 0.0216 → 0.1893, and robustness costs
~0.5% on clean data. The formulation demonstrates a general recipe — the
**whitening trick** — for keeping epigraph SOCPs in the dense-block habitat:
naively t ≥ ‖Eθ − y‖ needs a residual variable tied to θ by identity blocks
(internally sparse maps); instead, with E'E = LL',
b = L⁻¹E'y, c = y'y − ‖b‖², one has ‖Eθ−y‖² = ‖L'θ−b‖² + c, hence
(t, L'θ−b, √c) ∈ SOC(n+3), and eliminating θ = L⁻ᵀ(w+b) leaves one SOC stalk
(t, w, s) per patch with jets as dense maps through L⁻ᵀ (nonzero rhs) and a
one-row pin s = √c — no free stalks, no identity blocks, bloat ~1. Oracle
rigor: whitening identity to 5.9e-16; whitened conic vs direct norm-sum solve
to ‖Δf‖∞ = 4.9e-5; the escalation/saturation gates above; a stationarity
certificate via jet multipliers (1.2e-5). Self-contained numpy validation in
tools/sqrtloss_oracle.py.

**Dial:** Chain length P at fixed degree n = 24 (vertex count grows,
per-block work constant); benchmark rows time a single solve.

## E07: Poisson-TV photon-limited deconvolution (exp-cone habitat)

E05's Poisson twin, and the first example to exercise the exponential cone:
photon counts y ~ Poisson(τ·(Gu) + b) under a Gaussian PSF G with dark
background b > 0, reconstructed by exact penalized MLE — Poisson negative
log-likelihood (Shepp–Vardi territory) plus Gilboa–Osher nonlocal TV — over
E05's piecewise-Chebyshev tiles with C^k jets. Each positive-count bin
carries one 3-d exponential-cone stalk (t ≤ log μ, pinned second slot,
μ tied to ≤ 2 tiles by a dense PSF∘eval row); zero-count bins contribute
only the linear term Σμ and fold into c — no cone. TV stays SOC; u ≥ 0 is
sampled at tile nodes into per-tile NOC stalks. Three cone families in one
program, all small blocks on a thin path graph, and every map dense by
nature. The scene (dark plateaus + compact source + textured band) makes
the likelihood matter: at median count 3 the exact Poisson fit HALVES the
MSE of both Gaussianized workflows (0.0092 vs 0.0179 weighted-Gauss-TV vs
0.0178 no-deblur Anscombe-TV), and positivity is load-bearing in the
NODAL formulation — sub-graph-resolution spike dipoles are TV-invisible
and blur-null, so dropping the NOC sends min u to −220 and MSE to 98.4
(the "nearly black object" role of positivity;
Donoho–Johnstone–Hoch–Stern 1992). In the compressed space these
Nyquist-scale dipoles are unrepresentable (E05's ~π modes/wavelength
threshold), so basis truncation and positivity act as redundant
safeguards against the same null space — measured MSE ratio 1.0 when
dropping the NOC at m = 32. Oracle rigor (tools/poisson_tv_oracle.py, self-contained
numpy+cvxpy): MLEM certifies the conic MLE optimum one-sidedly to 3.5e-6;
Clarabel-vs-SCS cross-solver agreement 4.4e-8; subgradient-selection KKT
residual 2.0e-7; Gaussian-limit deformation into the SOC problem at
measured rate τ^−0.38; completeness (full basis, no jets) to −2.7e-11
excess; and the compression bonus recurs — m = 32 of Tsz = 128 (u-DOF ÷4)
beats nodal, MSE 0.0083 vs 0.0092, the truncated basis filtering
fine-scale Poisson noise exactly as it filtered Gaussian noise in E05.

**Dial:** Signal length N at fixed Tsz = 128, m = 32, K = 12, τ = 15;
benchmark rows time a single solve. Note the exponential-cone convention:
this IPM uses the log form (x₁, x₂, x₃), x₂ log(x₁/x₂) ≥ x₃; MOI/Clarabel
use the exp form, so baselines reverse coordinates.


## E08: Optimal execution with 3/2-power market impact (power cone)

The first power-cone example, and the return of the dense-Q habitat:
liquidate a portfolio over T periods under quadratic risk through a dense
factor-model covariance (Q_tt = 2γΣ per period — density in Q itself, as
in E03/E04, where E05–E07 carry it in the maps) and temporary impact
Σ η_i|v_{t,i}|^{3/2} — the empirically established impact power law
(Almgren–Thum–Hauptmann–Li 2005), used verbatim in Boyd et al.'s
multi-period trading framework (2017). Each asset-period pair is one
PowerCone(2/3) stalk (s, 1, v), s ≥ |v|^{3/2}; v-tie edges carry ±I_n
between adjacent holding stalks on a path graph over time. Unlike E07's exp
cone, this IPM's power cone matches the MOI convention, so baselines need
no coordinate reversal. Oracle rigor (tools/almgren_pow_oracle.py): the
suite's second ANALYTIC oracle — the continuous single-asset problem
conserves E = (η/2)|ẋ|^{3/2} − κx², reducing the BVP to a scalar bisection
(boundary layer removed exactly by x = √(E/κ)·sinh u), and the discrete
conic solution converges to it at measured rate Δt^1.95 in trajectory,
O(Δt) in cost; analytic limits (η×1e4 → TWAP to 3.9e-4, η→0 → immediate
liquidation to 1.6e-6, base regime mid at dev 0.673); exact decoupling
under diagonal Σ (9.2e-7) with factor coupling moving the solution by
0.375 — the dense Q is load-bearing, and the optimum transiently shorts a
hedge asset; three solver families agree (conic vs L-BFGS 5.4e-7, Clarabel
vs SCS 6.2e-6); C¹ stationarity residual 3.0e-7; monotone risk/impact
frontier over γ ∈ [0.25, 4].

**Dial:** Horizon T at fixed n = 8 (vertex count grows, per-block work
constant, cf. E06); n dials block size instead if preferred. Benchmark rows
time a single solve.


## E09: Robust Kalman smoothing under burst-corrupted measurements

E06's dynamic sibling, and the first occupant of the dense-Q + SOC cell of
the habitat matrix (E08 is dense-Q + pow; E05–E07 carry their density in the
maps with Q ≡ 0; here both mechanisms are live at once). Fixed-interval
smoothing of a linear-Gaussian state-space model — x_{t+1} = A x_t + w_t,
y_t = C x_t + v_t — with E06's square-root loss on the whitened measurement
residuals: prior and process terms stay quadratic, measurements contribute
λ‖Lv⁻¹(y_t − C x_t)‖ (norms, not squares). The sqrt loss bounds each
timestep's pull on the trajectory at λ, so measurement BURSTS — contiguous
windows where every sensor is corrupted (telemetry dropouts, EMI) — cannot
drag the estimate: the chain coasts across the outage on dynamics alone.
With the quadratic loss the same program IS the Rauch–Tung–Striebel
smoother (1965), which is what makes this family an analytic oracle; the
robust-loss lineage is Aravkin–Burke–Ljung–Lozano–Pillonetto's generalized
Kalman smoothing (Automatica 2017).

Everything is dense by nature: A = exp(Ac·Δt) for a damped spring-mass
chain (the matrix exponential is dense even though Ac is tridiagonal), W
from the Van Loan (1978) block exponential, C dense sensing rows, and
V = σ²(FF′ + I) a dense factor model (common-mode sensor drift, cf. E08's
factor Σ), so the whitened tie map M = Lv⁻¹C is dense too. Measured
densities at 1e-8: A 0.986, W 1.0, V 1.0, M 1.0.

Sheaf structure. Vertices t = 1..T are overlapping pair stalks
z_t = (x_{t−1}, x_t) ∈ R^{2n} (CofreeCone) carrying the dense process Gram
[A′ΩA, −A′Ω; −ΩA, Ω] with Ω = W⁻¹ (+ the prior P0⁻¹ on the first block of
z_1) — E01's overlapping-patch discipline for coupled quadratics — joined
by n-row agreement edges. Each timestep hangs one SOC stalk (u, w, s) of
dim m+2 off the chain: m dense tie rows w_t = M x_t − b_t (nonzero rhs,
E06's pattern) and a one-row pin s_t = δ. δ = 0 is the flagship sqrt loss;
δ > 0 with objective weight δ is the pseudo-Huber deformation
ψ_δ(r) = δ(‖(r, δ)‖ − δ) → ‖r‖²/2, whose limit is exactly RTS. Unlike E06
no Gram whitening is needed: the residual dimension m is small and the tie
is dense by nature already.

Gate tests validate the dense-by-nature densities, the two in-file
classical references against each other (RTS recursion vs banded
information-form solve), the exact assembly identity, the pseudo-Huber
deformation onto RTS, the robustness escalation, exact support recovery of
the corrupted bursts, and the IPM against a Clarabel solve of the same
conic program. Full validation in tools/robust_smoother_oracle.py
(self-contained numpy + scipy + cvxpy, fully executed; figures at n = 12,
m = 4, T = 200, λ = 1.5, three 8-step bursts):

* ANALYTIC ORACLE (the suite's third, after E01's MC and E08's ODE): three
  independent algorithms agree — RTS recursion vs information form
  2.0e-14; conic quadratic MAP vs RTS 4.6e-15.
* FORMULATION CERTIFICATE: the pair-stalk + agreement + SOC-tie assembly
  (exactly what e09/def.jl builds) equals the native norm-sum formulation to
  ‖Δx‖∞ = 7.4e-7; the pair-Gram split identity holds at random points to
  4.7e-16.
* DEFORMATION: x(δ) → RTS at measured rate δ^−1.98 (predicted −2).
* ESCALATION (corruption ×1 → ×128): robust clean-region state RMSE
  saturates 0.268 → 0.339 while RTS blows up 0.261 → 2.158; overall RMSE
  at ×128 is 0.554 vs 5.644 (10.2×); the clean-data price is +2.7%.
* SUPPORT RECOVERY + INFLUENCE CAP: the top-24 whitened residuals are
  exactly the 24 corrupted steps (margin ×2.8); the quadratic loss lets
  corrupted steps pull with ‖r‖ up to 95.9 vs the robust cap λ = 1.5
  (64×).
* KKT: subgradient-selection stationarity residual 1.35e-4 (70
  interpolation-active steps with ‖w_t‖ = 0, E05's discipline).
* CLASSICAL CERTIFICATE (DARE, the D-optimal oracle species — tools/doptimal_oracle.py): the influence of a
  boundary data perturbation on interior states decays at e^−0.2881/step
  vs the closed-loop ρ((I − KC)A) = 0.7519 predicted by the discrete
  Riccati equation — slope ratio 1.010.
* CROSS-SOLVER: Clarabel vs SCS ‖Δx‖∞ = 2.1e-6.

**Dial:** Horizon T at fixed n = 12, m = 4 (vertex count grows, per-block
work constant, cf. E06/E08); benchmark rows time a single solve.
DOF = 30T + 6, N1 = 17T − 7, median block 24.


## E10: Compositional Lyapunov certification on a K×K torus (`e10/`)

The suite's control-theory pillar, designed at the conjunction of the two
measured win-conditions — **dense Q** and **small balanced blocks** — on
the 2D topology where the fat-primal quantum-marginal experiment lost. A
K×K torus of subsystems with dense internal dynamics A_loc and dense
nearest-neighbor couplings is certified stable compositionally:
min Σ_v [tr P_v + (γ/2)‖W^{1/2}P_v W^{1/2}‖²_F] subject to
A′P + PA + I + S = 0 with P = diag(P_v) ⪰ 0 and S Agler-split into one
2d×2d PSD stalk per torus edge. W is the local controllability Gramian
(A_loc W + W A_loc′ + I = 0) — a canonical, dynamics-derived, fully dense
weight, so the Q blocks are dense by nature, and dense Q is free for the
sheaf solver (b³ per block regardless) while taxing baselines through
KKT fill or conic reformulation. Block accounting: vertex stalks svec 10
(P_v) and 36 (S_e) against row blocks of 10–16 — b_v ≈ b_e, no thin side
for a constraint-side factorizer to monopolize, with d a dial that
fattens both sides together. This is the **falsifiable test of the
thin-side account**: predicted IPM-competitive where the 136-vs-10
quantum torus was predicted (and measured) not to be.

The γ-term is not decoration: it makes the program strictly convex in P,
so the optimal certificate is **unique** (Tikhonov selection on the
linear SDP's optimal face) and cross-solver gates compare **states**, not
just values — a first among the suite's SDP examples. γ → 0 recovers the
pure linear certificate at measured rate γ^1.00.

Two new species. **Sum-decomposition edges** (vs. agreement, partial
traces, jets): the Agler split of the banded slack into clique-supported
PSD pieces, which on the torus's non-chordal cyclic band pattern is a
strict restriction — that *is* a compositional certificate, and the
conservatism is measured, not hidden. **A spectral (2D-DFT) oracle**: the
monolithic comparator (P block-diag, one big LMI) reduces exactly by
group averaging to per-wavenumber d-dim Hermitian LMIs at any K, so the
conservatism gap curve is computable everywhere: per-cell at γ = 0 it
runs 8.11% (K=3) → 3.36% → 2.90% → 3.10% (K=8), converging to its
continuum limit — E02's non-chordal argument relocated to the constraint
side, with a convergence curve attached. And by the private-slot theorem
(every row group pins a fresh S_e block), **coker B = 0 by design** —
verified by numerical rank at K = 3, 4 — making E10 a redundancy
*control*: same 2D plaquette topology, zero redundancy by design (the
clean counterpart to a coboundary-constrained program with nontrivial coker).

Validated by `tools/torus_lyapunov_oracle.py` (self-contained
numpy + scipy + cvxpy; fully executed; d = 4, ε = 0.15, γ = 0.1):
Gramian residual 8.9e-16; coker = 0 exact (B 378×738 / 672×1312, aspect
0.512, full row rank); Agler-split feasibility v̄ = 1.50529852; symmetry
oracle |v(K) − K²v̄| = 1.7e-12 / 3.2e-12 with state uniqueness 3.3e-9 /
6.2e-9 at K = 3 / 4 (K = 2 degenerate — wraparound doubles couplings —
documented); spectral vs. direct monolithic 5.1e-8; **dynamical
certificate** — V = x′Px monotone decreasing along ẋ = Ax with
V(t) ≤ V(0)e^{−t/λmax(P)}, max ratio 1.000; **instability
localization** — destabilizing one subsystem makes the program
infeasible and the minimal per-vertex relaxation concentrates exactly
there (t_broken = 1.036, elsewhere 0.0, healthy control all zeros);
deformation γ^1.00; formulation certificate 3.5e-11 (values) / 1.8e-8
(states); Clarabel vs SCS on states 1.1e-8.

First-run benchmarks: K dial 1.6–2.1× vs Clarabel, 0.79–0.96× vs
Mosek — the thin-side prediction held (the fat-primal torus at matched
DOF sat at 0.36× with Clarabel DNF). The d dial confirmed the
per-nonzero tax against Clarabel (1.70× → 3.08× as d = 4 → 8) but not
against Mosek (0.97× → 0.85×): Mosek's SDP Schur assembly exploits the
low rank (≤ ~4) of every constraint matrix here — Clarabel pays per
nonzero, Mosek per rank, the sheaf IPM per block. The d = 8 row is an
accidental control: svec-136 S stalks, the exact block size that sank
the fat-primal torus, hold 0.85× here because the rows are balanced
(aspect 0.53 vs 0.148) — balance, not block size, governs.

Dials: torus size K (DOF = 82K²) and subsystem dim d; γ, ε as physics
knobs. B has full row rank with a wide aspect — the clean structural
opposite of a coboundary-constrained (rank-deficient) program on the same graph.

Files: `e10/`, `tools/torus_lyapunov_oracle.py`.


## E11: Multichannel FIR equalizer bank design (`e11/`)

Short multichannel FIR equalizer-bank design under colored noise — regularized
MINT (Miyoshi–Kaneda 1988) at D = 1. M sensors observe a source through
reverberant FIR channels h_m; design one equalizer f_m per channel so the
composite response c = Σ_m h_m ∗ f_m approximates a pure delay, while
minimizing amplified ambient-noise power under an AR(2) colored-noise field.

Dense-Q + SOC habitat at larger block sizes (med_blk ≈ 65, vs E09's 6). The
noise autocorrelation Gram Θ is dense by nature (the noise is colored), giving
dense Q_m = ς_m² Θ + λI on each filter stalk. The constraint is a whitened
SOC ball: with G = T'T = LL' the joint convolution Gram, the reconstruction
budget ‖Tf − g‖ ≤ ε transforms to ‖L'f − b‖ ≤ εx (E06's whitening trick,
relocated from objective to constraint). The **exact slice split** partitions
the whitened residual by Cholesky row blocks — the gluing maps ARE the
Cholesky blocks of the joint Gram — giving zero conservatism but O(M²)
coupling mass. The **budget split** (--budget) uses per-channel whitening and
a triangle-inequality certificate, giving O(M) coupling mass at the price of
measured conservatism (~1.2x at small εx).

New edge species: Cholesky/convolution ties (dense triangular blocks from the
reverberant autocorrelations) and a hub SOC enforcing the ball constraint.
Gate tests include a native TRS (trust-region subproblem) oracle via secular
equation (Gay 1981; Moré–Sorensen 1983).

First-run benchmarks: IPM wins vs both Clarabel (3.4–6.8×) and Mosek
(1.5–3.2×) at M = 4–16. But the M-dial slope (DOF^2.29 vs Mosek's 1.72)
shows the M² coupling mass — the triangular slice hyperedges are visible
in the scaling.

Dials: channel count M (DOF = M(2Lf+2)+1) and block size Lf (--ldial).

Files: `e11/`, `tools/mint_equalizer_oracle.py`.


## E12: Distributed sound zones — cooperative multizone rendering (`e12/`)

Z zones on a path, each with a 4-speaker array and 3 control points;
reverberant FIR acoustics (Lh = 96) with inter-zone gain γ^{|Δz|}. Every
array within radius ρ = 1 of a program cooperates in rendering it
(pressure matching, Betlehem et al. IEEE SPM 2015), giving stalks f_{z',p}
(dim 64) with dense radiated-power Grams (sinc radiation coupling ×
heterogeneous AR(2) program spectra, Morse–Ingard) and one reproduction
ball per (zone, program): reference-field bright targets, dark leakage
caps, signal-scale radii at the global-LS floors. Per-array acoustic-output
budgets (SPL caps, same Gram family) connect the programs at each array.
Every ball gets E11's exact whitening + slice split LOCALLY: slice SOCs of
dim 65, hub SOCs of dim ≤ 5 — bounded arity, coupling mass O(Z). E12 is
E11 tiled, built to answer E11's M-dial reading (DOF^2.29 under O(M²)
global coupling).

Three structural observations recorded during design (oracle-measured):
(i) one-array-per-program renders the gluing exactly block-diagonal
(tri-density 1/(2ρ+1) = 0.338) — cooperative rendering restores 1.000;
(ii) without the budgets the problem decomposes exactly by program (Z
components, machine-checked) — and the cure is connectivity, not reach:
the physics is screening (hop ratio 0.21 at γ 0.25 vs 0.31 at 0.4,
the E09-DARE species); (iii) a signal-power budget (I_S ⊗ Toeplitz)
would be block-diagonal across speakers, tri-density 1/S — avoided by
capping acoustic rather than electrical power.

Oracle (`tools/soundzone_oracle.py`, 10 checks, fully executed): endpoint
η^0.98; assembly == native 1.8e-7 / 7.3e-11; TRS (Z = 1 reduction) with
Moré–Sorensen residuals ~1e-17; γ = 0 separates per array to 1.2e-9;
truncation reach γ^1.79 / γ^2.80 (increment 1.01 — one power of
γ per zone; the base exceeds naive γ^1 because the dark balls
shadow farther leakage); response-domain contrast +2.3..+4.3 dB per zone
(mean +3.3) over non-cooperative at budget utilization 0.83–0.94;
cross-solver 4.2e-8.

Benchmark (Z-dial 4 → 32, med blk 65). PRE-REGISTERED slope ~1.1–1.3
against E11's 2.29: measured IPM DOF^1.17 (HSD 1.19), Clarabel 1.07,
Mosek 1.55. IPM wins 2.14× vs Clarabel and 5.62× vs Mosek at Z = 32.
Honest caveat: Clarabel's slope is 0.10 SHALLOWER than ours — the win
over Clarabel on this dial is constant-factor, not exponent (extrapolated
crossover ~1e7+ DOF, far outside the regime; stated, not hidden). Mosek
mirrors E11 exactly: it loses the count dial (2.21× → 2.83× → 5.62×) as
decisively as it won E11's size dial.

Ablations: --joint (per-ball monolithic SOCs, dim ~194 vs 65): slope
nearly unchanged (1.22) but the win compresses 2.14× → 1.33× — the
balanced-blocks lesson operating at the constant. --nobudget
(decomposability control): IPM indifferent (1.16); only Mosek's presolve
shows decomposition-awareness (1.55 → 1.29); Clarabel does not (1.10).
--ldial (Lf 16 → 48 at Z = 8): all superlinear as predicted, but ratios
NOT stable — Clarabel steepest (2.77, IPM win grows to 4.02×), Mosek
shallowest (1.90); see the count-vs-size regime rule below.

Files: `e12/`, `tools/soundzone_oracle.py`.

---

## E13: Minimum-fuel soft landing with tracking (`e13/`)

3-DoF powered descent, fixed mass, exact-ZOH discretization — the convex
(SOC) lossless-convexification relaxation of Açıkmeşe–Ploen minimum-fuel
powered-descent guidance (JGCD 2007; the G-FOLD line) — plus a **dense
quadratic tracking term** `(γ/2) Σ (x_t − x̄_t)ᵀ W (x_t − x̄_t)`. `W` is the
inverse Van-Loan dispersion of accelerometer noise through the double
integrator: a dense 6×6 with real position–velocity cross terms, which forces
positions and velocities into one stalk — a genuinely dense-Q construction, the
habitat cell the pure-fuel landing (Q ≡ 0) lacks. Sheaf structure follows
E09's chain discipline: ℝ⁶ state stalks on the time path with a dense 6×6 Q
block each, SOC(4) thrust stalks on the ZOH dynamics edges, NOC(2) thrust
bounds, SOC(3) glideslope stalks.

**The dense-Q ranking inverts.** The quadratic folds cheaply into the sheaf
solver's block Hessian but is a reformulation tax for the conic solvers, so
HSD leads and the gap grows with N: at N = 240 → 960, HSD 28 → 117 ms,
Clarabel 33 → 158 ms (1.2–1.35×), Mosek-dual 54 → 455 ms (1.9 → 3.9×; slopes
HSD DOF^1.02, Clarabel 1.13, Mosek-dual 1.54). Mosek is benchmarked on its
dualized form (its best form here — dualizing recovers ~1.4× over the primal
QP but not the lead). This is the mirror image of the Q ≡ 0 case, where Mosek
would win — dense-Q is where "the sheaf IPM per block" pays off.

**Gated construction.** With `x̄` := the pure-fuel optimum, `x̄` minimizes both
terms simultaneously, so `x(γ) ≡ x̄` for all γ — an exact structural invariance
(matched at oracle accuracy 2.3e-5 m); γ → 0 recovers pure min-fuel at rate
γ^1.00; the fuel-vs-tracking frontier is monotone; and feasibility is
objective-independent (the min-landing-time boundary is unmoved by tracking).
**Caveat (measured, gated):** under tracking the LCvx relaxation is not tight
(‖u‖ < ρ1 at some nodes), so the benchmark object is honestly the relaxation.

**A near-degenerate face, and an open stall.** The min-fuel optimum is a
near-degenerate face: the strict-complementarity margin over ρ1-active thrust
arcs collapses as N⁻²·⁰⁰ (measured cross-solver at 1e-11 tolerance, reproduced
to the digit); γ lifts it ×8–×29. The plain affine IPM **stalls at N ≥ 240**
(a limit cycle at that face) while HSD solves every N. The margin law is real
but does **not** govern the stall — lifting the margin ×29 does not cure it,
and every formulation-level candidate (degeneracy, corner-smoothing, encoding,
augmentation) is falsified by direct experiment. A metrology note carries
through: boundary quantities (active sets, margins) must be measured at tight
(1e-11) tolerance, not benchmark-tolerance iterates that sit ~1e-4 off the
boundary. HSD is the solver at scale here.

Files: `e13/`, `tools/soft_landing_oracle.py`.

---

## E14: Gyre-aware robust flow recovery on a punctured grid (`e14/`)

Drifter-style flow measurements `g_e ∈ ℝ^T` on the edges of a K×K planar grid
with **islands** (removed faces → `dim H¹ = #islands`). Recover the Hodge
decomposition of the flow — potential `x` (vertices), vorticity `z` (faces), and
**gyre circulations `h`** (one block per island, the estimand) — while localizing
burst-corrupted edges:

```
min  Σ_e ½ u_eᵀ W_e u_e + λ Σ_e t_e + ridge(x,z,h)
s.t. (δ⁰x)_e·Φᵀ + (δ¹ᵀz)_e·Φᵀ + (ηh)_e·Φᵀ + u_e + v_e = g_e,   ‖v_e‖₂ ≤ t_e
```

`Φ` is an orthonormal smooth temporal basis, `W_e` are dense per-edge GP
precisions, and `η` is the orthonormal harmonic basis of the punctured complex.
This is the suite's **first full-Hodge design**: the constraint operator carries
`δ⁰`, `δ¹ᵀ`, *and* a harmonic stalk as first-class column blocks, so the normal
operator contains the complete degree-1 Hodge Laplacian plus harmonic completion
— structurally consistent for every `g` (no `ker Bᵀ` pathology). The gyre
coefficients are the **estimand**, not an obstruction, and an island-blind
ablation measures exactly what mismodeling `H¹` costs. Smooth-basis restriction
is the identifiability mechanism: white bursts lie outside `span(Φ)` at any
magnitude, so burst/clean separation is structural (measured margin 14.1×).

**A dense-Q + dense-column win.** Both sheaf solvers beat the conic baselines at
every K: at K = 8 → 16, IPM 11.5 → 45.8 ms, HSD 9.5 → 46.2 ms, Clarabel 13.1 →
64.7 ms, Mosek-dual 30.8 → 93.5 ms (slopes IPM DOF^1.04, HSD^1.18, Clarabel^1.20;
Mosek benchmarked dual-only). The dense GP precisions `W_e` fold into the block
Hessian for free but tax the conic solvers, and — crucially — the harmonic
stalk enters as a dense **column** of `B`, which a min-degree ordering eliminates
last, so `BᵀB` stays sparse. (Its transpose — gyres *prescribed* as constraints
— would put `ηᵀ` as a dense **row**, giving `BᵀB = L₁ + ηηᵀ` dense and poisoning
the factorization: dense columns are cheap, dense rows are not.)

**Gated construction.** `‖δ¹δ⁰‖ = 0` exactly; `dim H¹ = #islands`; the three
column blocks are Hodge-orthogonal (≤ 7e-16); LS gyre recovery sits at the ~2.5σ
noise floor; the island-blind ablation is a dose-response curve; group-lasso
localization has a 14.1× screening margin (burst support recovery is
draw-dependent — the shipped draw finds 5/6 vs the oracle's 6/6 — while the
structure and objective gates are exact).

Files: `e14/`, `tools/ocean_gyre_oracle.py`.

---

## Suite-level lessons

**Coupling mass sets the slope only under bounded topology.** E11 + E12
now demonstrate this as a three-point ladder within one problem family:
O(M²) global coupling → DOF^2.29 (E11 flagship); O(M) coupling through
an arity-M star hub → DOF^1.97 and a 6× LOSS to Clarabel (E11 --budget:
the hub SOC has dim M+2, a single (M+2)³ vertex cost, slope → 3
asymptotically); O(M) coupling through a binary hypot tree → DOF^0.99
and a 1.77× win (E11 --tree); O(Z) bounded-degree tiling → DOF^1.17
(E12). The habitat condition is therefore: dense-after-scaling blocks AND
bounded vertex arity AND bounded max block dimension — med_blk alone is
not sufficient (the star keeps med_blk at 65 while the max block grows
unboundedly; E10's imbalance lesson in topological form).

**The hypot tree is a reusable formulation species.** Any arity-k hub
s = ‖(s_1, .., s_k)‖ splits exactly into a depth-log(k) binary tree of
dim-3 combiners s_parent = ‖(s_left, s_right)‖ (nested norms; zero
conservatism). Cost: O(k) extra 3-dim SOCs. Use it whenever a budget or
ball aggregates more than O(1) channels.

**Count-dial vs size-dial regime rule (three independent sightings).**
Mosek is the strong opponent on BLOCK-SIZE growth (E11 M-dial 1.72; E12
--ldial 1.90, shallowest) and the weak one on BLOCK-COUNT growth (E12
Z-dial 1.55, steepest). Clarabel is the reverse (Z-dial 1.07, best;
--ldial 2.77, worst). The sheaf IPM's advantage is maximal on count dials
with bounded blocks (wins the exponent vs Mosek, the constant vs
Clarabel) and on size dials vs Clarabel; its exposed flank is size dials
vs Mosek (E11 flagship crossover at M ~ 16).

**Practical recommendation for E11-type problems now flips with scale:**
the flagship exact-slice split (zero conservatism) up to M ~ 16; the
--tree budget split (conservatism 1.04 measured) beyond — slope 0.99,
36.5 ms at M = 32 vs the flagship's 244.9 ms.

---

## Files

- `utils.jl` — Shared utilities (command-line parsing, timing, baselines)
- `e01/` — American basket put chain under Merton jumps (dense NOC)
- `e02/` — Quantum marginal (non-chordal SDP)
- `e03/` — Tensor-spline regression (dense QP)
- `e04/` — SOS-certified regression (dense SDP)
- `e05/` — Galerkin-compressed nonlocal TV (dense-map SOCP)
- `e06/` — Sqrt-loss spectral regression (whitened SOC chain)
- `e07/` — Poisson-TV photon-limited deconvolution (exp + SOC + NOC)
- `e08/` — Optimal execution with 3/2-power impact (dense-Q pow chain)
- `e09/` — Robust Kalman smoothing, burst corruption (dense-Q + SOC chain)
- `e10/` — Compositional Lyapunov on a K×K torus (dense-Q dense-SDP)
- `e11/` — Multichannel FIR equalizer bank (dense-Q + SOC, Cholesky ties)
- `e12/` — Distributed sound-zone control (dense-Q + SOC, bounded degree)
- `e13/` — Min-fuel soft landing with tracking (dense-Q SOC chain, LCvx)
- `e14/` — Gyre-aware robust flow recovery (full-Hodge: dense-map + dense-Q + SOC)

## Reference documentation

The `ref/` directory contains detailed mathematical background:

- `specialist_problems.md` — E01 (jump obstacle)
- `obstacle_option.md` — E01 foundation (1D obstacle/American option framework)
- `quantum_marginal.md` — E02 (XXZ spin chains, partial-trace coupling)
- `tensor_spline.md` — E03 (tensor-product Bézier surfaces, C^k jets)

## Project setup

```bash
cd src/IPM/examples
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```
