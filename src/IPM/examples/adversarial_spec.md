# Adversarial problem spec — α-controller stress campaign

Target directory: `examples/adversarial/`, one file per instance (`X03.jl`, `X04.jl`, …), each
exporting a `build_*` function returning an `IPMProblem` (via `IPMProblem(Q, B, c, g, K)` /
`blocksparse` / `IPM.allocblockdiag`), deterministic under a `seed` kwarg. Driver keys go in
`buildprob` in `run_oracle.jl` (X03/X04 already stubbed there with the signatures used below).

## Design constraint: immunity to the equilibrator

The Ruiz pass (`src/scaling/scaling.jl`) applies an unconstrained per-row scaling E and a
**single positive scalar per cone-vertex block** D = blkdiag(t_v·I), driving the ∞-norm of every
row of B and of every vertex block max(‖B_{·v}‖∞, ‖Q_v‖∞) to ≈ 1 (10 sweeps, tol 1e-3). HSD adds
one scalar on [c; g]. Therefore an instance is only adversarial if its pathology is invariant
under (E, block-constant D). Three invariants survive:

- **I1 — within-block spread:** the ratio of column norms *inside* one vertex block. D is
  block-constant; E can shuffle it between rows but with enough rows per column the spread is
  pinned (see X03 construction, which pins it).
- **I2 — row angles:** E·B·D preserves the angle structure of B's rows; near-parallel rows stay
  near-parallel at unit norm. σ_min(B̂) stays ~ ε.
- **I3 — Q/B conflict inside a vertex:** vnrm = max(B-block, Q-block), so a vertex with
  ‖Q_v‖ = 10^{2s} and ‖B_{·v}‖ = 1 gets t_v ≈ 10^{-s} (Q scales as t², B as t), leaving
  ‖B̂_{·v}‖ ≈ 10^{-s}: effectively-tiny constraint columns the equilibrator *created*.

Every instance below states which invariant it exercises and why the immunity argument holds.
Each also ships a manufactured optimum so correctness is checkable independently of the solver:
pick the intended primal p\* (interior or boundary as designed), set g = B p\*, pick dual
d\* ∈ K\* respecting the intended complementarity pattern and y\* freely, and set
c = Q p\* − Bᵀ y\* − d\*. Gate test per instance: solve at 1e-8 and assert
‖p − p\*‖∞ / (1+‖p\*‖∞) ≤ 1e-5 (looser where degeneracy makes p\* non-unique — then check
objective value and residuals instead).

Sizes: keep n in the low hundreds so a full oracle sweep (21 α × ~15 iters × deepcopy) stays
cheap. All randomness through `StableRNGs` or `Random.seed!(seed)` with a fixed default.

---

## X03 — `build_narrow(; npos=12, nfree=12, spread=8.0, rows_per_col=3, seed=1)` [I1]

**Idea.** One large PositiveCone block and one large CofreeCone block whose B-columns have
geometrically graded norms spanning `spread` decades — spread *inside* each block, where the
block-scalar cannot reach.

**Construction.** Two vertex blocks: v₁ = PositiveCone, dim npos; v₂ = CofreeCone, dim nfree.
m = rows_per_col·(npos+nfree) rows. For column j of block v, draw `rows_per_col` random rows and
place entries ±10^{−spread·(j−1)/(dim−1)} · (1 + 0.1·randn). Add one extra dense-ish anchor row
per block touching every column at that column's scale (keeps rows from being individually
rescalable into hiding the spread: every row contains one O(largest-in-row) entry, so E is pinned
near 1 after D normalizes block maxes — state this in a comment and assert post-equilibration:
run `equilibrate!` on a copy in the gate test and check the within-block column-norm ratio still
≥ 10^{spread−1}). Q = small ridge (1e-6·I per block) so the anchor's nQ is negligible and the
instance isolates B. Manufactured p\* strictly interior (pos part ≈ 1, free part ≈ 1).

**Predicted effect.** σ_min(B̂) ~ 10^{−spread} ⇒ the AL numerator constant N₀ inflates roughly
like σ_min^{−2}·(activity of the small columns); the plateau's lower edge rises by O(2·spread)
decades. At spread = 8 the plateau should be pushed toward or past the 1e10 grid ceiling —
extend the oracle grid to 1e14 for X03 runs.

**What it decides.** Whether N₀-universality (0.6–11 on the equilibrated e-set) is a property of
the equilibrator or of the AL structure. If the anchor·10 rule fails here while the r₀ controller
tracks the moving edge, that is the controller's headline case. Sweep spread ∈ {2, 4, 6, 8}.

## X04 — `build_degenerate(; n=16, degen=8, rowscale=8.0, eps=1e-6, seed=1)` [I2, ρ-ladder]

**Idea.** Near-dependent constraint rows. Ruiz normalizes row norms but preserves angles.

**Construction.** Base: a well-conditioned random B₀ (m₀ = n/2 rows) over a PositiveCone block of
dim n (plus a small CofreeCone block, dim 4, to involve both cone types). Append `degen` extra
rows: row i' = row i + eps·(random unit row), then multiply each appended row by
10^{uniform(−rowscale, rowscale)}. The `rowscale` part is deliberately *not* immune — the
equilibrator should remove it exactly; the gate test asserts post-equilibration row norms ≈ 1
while the principal angle between row i and row i' stays ≤ ~eps. This separates "equilibrator
fixes it" (rowscale) from the residual pathology (angle), which is the point of the knob.
g = B p\* with interior p\* keeps the duplicated constraints consistent.

**Predicted effect.** σ_min(B̂) ~ eps ⇒ BᵀB nearly singular. Large α makes F ≈ BᵀB: expect the
Cholesky ρ-shift ladder (never fired on the e-set — ρ stuck at rgmin throughout) to finally
escalate at the top of the grid, the upper wall to descend, and the plateau to be squeezed from
*both* ends. Watch the `rho` column: this is also the first dataset where the controller
interacts with ρ. Sweep eps ∈ {1e-4, 1e-6, 1e-8}.

## X05 — `build_qdominant(; nblk=8, dim=8, s=6.0, seed=1)` [I3]

**Idea.** Make the equilibrator itself manufacture tiny B-columns via the Q-vs-B vertex conflict,
and simultaneously stress the anchor formula raug·(nH+nQ)/nB², whose numerator this inflates.

**Construction.** nblk PositiveCone blocks of dimension dim. Half the blocks: ‖Q_v‖ ≈ 10^{2s}
(well-conditioned SPD, e.g. 10^{2s}·(I + 0.1·rand sym)), B entries O(1). Other half: Q_v = ridge,
B entries O(1). Constraint rows couple one loud and one quiet block each (so E cannot fix the
imbalance the vertex scalars create without breaking the quiet block's rows — assert
post-equilibration ‖B̂ columns of loud blocks‖ ≤ 10^{−s+1}). Manufactured interior optimum.

**Predicted effect.** Same σ_min mechanism as X03 but sourced from Q, plus a direct hit on the
anchor: nQ ~ 10^{2s} inflates alpha_anchor by ~10^{2s} while the plateau moves by ~10^{s·(…)} —
the anchor and the plateau move by *different* amounts, so the ×10 rule should misplace the
constant here even at s = 4. If it doesn't, the anchor formula is more robust than the analysis
suggests — either outcome is informative. Sweep s ∈ {2, 4, 6}.

## X06 — `build_corner_soc(; nsoc=6, dim=10, seed=1)` [cone-generality control]

**Idea.** Not a scaling trick — a control for cone-generality. The e-set oracle runs are
PositiveCone-dominated; the q ≈ 1, N ∝ μ^{0.55} law should be re-measured where the NT scaling
has intrinsic within-block spread. SOC instances whose optimum sits at the cone "corner"
(x₀ = ‖x̄‖ with x̄ ≈ 0) maximize the NT Hessian's internal spread, which no diagonal scaling
touches by design.

**Construction.** nsoc SecondOrderCone blocks; manufacture the optimum so half the blocks are at
the corner (choose d\* on the boundary ray accordingly) and half strictly interior. B random
well-conditioned, Q modest SPD. (If SDP cost is acceptable, a variant with two small
SemidefiniteCone blocks, rank-deficient at the optimum, is the same experiment for SDP.)

**Predicted effect.** If the μ^{−1/2} law and N₀ range survive here, the controller design is
cone-generic; if n (the N(μ) exponent) shifts, the feedforward constant is cone-dependent and
the r₀ controller's self-calibration becomes load-bearing rather than optional.

## X07 — protocol, not a problem: tolerance and warm-start travel

Run X03(spread=4) and two e-set problems (e01, e08) at feas_tol = gap_tol ∈ {1e-8, 1e-10, 1e-12}
with the extended grid, to measure the predicted two-sided plateau squeeze (lower edge ∝
tol^{−0.45}, upper wall descending ∝ tol) and locate the collapse tolerance directly. Separately,
exercise `reinit!`: solve, perturb g by 10%, re-solve warm — the second solve starts at small μ₁
where a fixed anchor is furthest from the instantaneous window; log the same oracle sweep on the
warm solve. This is the production-relevant case (chained QPs — e01/e08 are chains already) where
the controller can win on *speed*, not just robustness.

---

## Instrumentation and acceptance

Same oracle CSV schema as the current fleet (r0_\*, mu_next, step, hdiag, ρ, CRAIG status), plus:
extend `DEFAULT_ALPHA_GRID` to 1e14 for X03/X05 runs, and record in the meta the two
post-equilibration diagnostics each gate test computes (within-block column-norm ratio for
X03/X05; min row angle for X04) so the realized severity of each instance is in the dataset.

Per-instance acceptance questions, in decreasing order of importance:

1. Does q stay ≈ 1? (If not, the r₀ controller's exponent needs the measured q — record it.)
2. Where does the plateau intersection sit relative to (a) alpha_anchor·10 and (b) the r₀
   controller's trajectory at θ = 0.1? Success for the controller = stays in-window wherever any
   window exists; failure of the static rule is expected and fine.
3. Does the ρ ladder fire (X04), and does the α that triggers it match the predicted upper wall?
4. Does trajectory equivalence (mu_next spread ≤ ~2× within-window) still hold outside the floor
   regime? If it breaks mid-solve on any instance, flag it — the offline simulator's validity
   boundary is itself a result.

A null result — plateaus stay ≥ 3 decades wide and the ×10 anchor survives everything above —
is a publishable-quality answer too: it would say the equilibrator + anchor make the controller
unnecessary for this solver, and the spec is designed so that conclusion would be earned, not
assumed.


---

# ADDENDUM — numpy pre-validation results (read before implementing)

A faithful numpy port of the block-Ruiz equilibrator plus a miniature of the Uzawa solve
(Cholesky of F = H/α + BᵀB with the ρ ladder, CG-on-Schur as CRAIG stand-in, the refinement
loop with identical force/floor/stall logic) was run against small versions of every matrix
instance, with H from a manufactured central-path surrogate at μ ∈ {1e-2 … 1e-8}. Findings
that CHANGE the spec above:

**X04 is promoted to the primary adversarial.** Row angles are the one pathology the
equilibrator provably cannot touch, and the emulation confirms a clean severity ladder:
the plateau collapses (no feasible α on a 1e0–1e14 grid) at μ_collapse ≈ 1e2·eps²
(eps=1e-3 → collapse near μ=1e-7; 1e-4 → ~1e-5; 1e-5 → ~1e-3). Recommended production
settings: eps ∈ {1e-3, 1e-4} for "solvable at 1e-8 but fixed-anchor-hostile" (at eps=1e-3,
μ=1e-6 the feasible window was [3e12, 1e14] — far above any static anchor), and eps=1e-5
as a deliberate floor-regime/termination stress. Record μ_collapse per run and test the
eps² law.

**X03 is retained but demoted, with a mandatory redesign.** Two prototype failures to avoid:
(a) graded columns under the FREE cone make F near-singular at every α (tiny B and tiny
H=ridge on the same coordinates) — the problem becomes unsolvable, not controller-hostile;
grade under the PositiveCone only, where the barrier supplies H. (b) Pinning every row with
the SAME full-scale column makes rows near-parallel — accidental X04 contamination
(σmin ≈ 1e-17 at spread=2). Correct design: a band of ~m/2 full-scale columns, each row
pinned by a distinct one plus a second random top-band entry, graded tail after the band.
With that design the immunity holds exactly (post-equilibration spread = 10^spread) but the
measured plateau shift vs a spread=0 twin is only 1–2 decades: unconstrained row scaling
launders most within-block spread, because a row either carries a large entry (and doesn't
need its small columns) or is boosted wholesale by E. Expectation to test in production:
small effect; a null result here is confirmatory, not a failure of the run.

**X05 (qdominant) is DROPPED.** Doubly defeated in emulation: the quiet blocks keep every
row O(1) and satisfy the constraints without the loud columns (plateau came out wider than
control), and the intended anchor attack is void because alpha_anchor is computed from
post-equilibration norms, which are ≈ 1 by construction. Do not implement.

**Optional replacement X05′ — termination stress, clearly labeled as such:** the free-cone
version of X03 (graded free columns, ridge-only Q) produces instances where NO α reaches
force/floor at any μ — useful for exercising floor_patience/near-status logic and the
controller's behavior when no valid action exists, but it measures the termination design,
not α selection. Implement only if that codepath needs coverage.

**Calibration caveat for all predictions:** the surrogate exaggerates the lower edge's
μ-drift (≈ μ⁻¹ vs the production-measured μ⁻⁰·⁵⁵) because unit-random RHS lacks the real
solver's residual contraction. Use the prototype's numbers as differentials between an
instance and its own benign twin, never as absolute plateau positions. Every instance
should therefore ship with a matched benign-twin builder (spread=0 / eps=0) and the oracle
should run both.

X06 (SOC corners) and X07 (tolerance/warm-start travel) are unchanged — neither is testable
in the surrogate (X06 needs real NT scaling for non-symmetric-cone spread; X07 needs the
real trajectory) and both remain justified on the production data alone.

Prototype code (equilibrator port, uzawa miniature, instance builders) is in python and can
be handed over alongside this spec as a reference implementation for the gate tests.
