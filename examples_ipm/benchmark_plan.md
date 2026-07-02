# IPM Benchmark Suite: Analysis, Additions, and Diagnostics

Based on `examples_ipm/results.md` (2026-07-02 run), the benchmark sources, and the
solver internals (`src/IPM/kkt/uzawa.jl`, `src/IPM/ipm/`). The guiding principle:
**keep the losing benchmarks — they map the regime boundary — and add benchmarks that
live where the design actually pays off.**

---

## 1. Performance profile (what the results say)

The KKT solve is Uzawa: block-sparse multifrontal Cholesky over dense per-vertex
blocks, then CG in the coupling space, with `raug` as the augmentation parameter.
That architecture predicts every row in `results.md`:

| The solver wins when… | Evidence |
|---|---|
| Per-vertex blocks are **dense** and moderate-sized (BLAS-3 pays; generic sparse solvers face fill-in) | map_lp SparseMAP 3–7x; moment/CPTP SDPs up to 28x on complete graphs; tensor_spline regress 1.5–4x |
| **H¹ ≪ DOF** — the CG dual space is small | quantum_marginal (H¹=1): 2–7x; wasserstein_gauss dual: 2–6x |
| `raug` can be cranked to 1e5–1e6 (coupling well-conditioned ⇒ few CG iterations) | every winning row has raug ≥ 1e5 |
| The competitor must vectorize/dualize block structure the sheaf keeps native | graph-structured SDPs vs Mosek |

| The solver loses when… | Evidence |
|---|---|
| Blocks are **internally sparse** (tridiagonal stiffness, diagonal Q) — dense treatment does O(n³) work on O(n) nonzeros | obstacle 0.02–0.75x; "lifted diagonal Q" wasserstein rows 0.16–0.66x |
| **H¹ ≈ DOF** — CG space is the whole problem, and raug must drop to 1e-2 | obstacle (H¹ = 399–1599); game option 0.02x |
| Many **tiny cones/vertices** — per-leaf overhead dominates | num_sheaf meshes 0.16–0.5x; mle_spline with 1000 exp leaves 0.05x (58 binned leaves recover to 0.5–0.85x) |
| **Sequential re-solves** with no warm start | American put chain 0.23–0.45x vs OSQP warm-starting trivially |

The `raug` column is the single best predictor of the P/IPM ratio in the whole file.

One caution on the win column: the "15x" on `nonnegative_spline_exact` monotone is a
**Mosek-primal pathology** (735 ms primal vs 208 ms dual; vs Clarabel the same row is
~1.0x). Don't headline it — headline map_lp, the graph SDPs, quantum marginals, and
tensor_spline, which win against *both* reference solvers.

---

## 2. Benchmarks to ADD

### 2.1 Lévy / nonlocal obstacle problems (`obstacle/`) — highest priority *(builder draft: Appendix A)*

The current obstacle benchmark is the anti-pattern by construction: Q is a
tridiagonal Black–Scholes stiffness matrix, so OSQP's sparse LDLᵀ eats it directly.
Under a jump model the operator becomes **integro-differential — dense per patch** —
and the patch-dense-block design is exactly the right discretization. Same finance
story, large literature, and it *strengthens* the Schwarz-patch motivation already in
the file header.

Candidates, in increasing implementation effort:

1. **Merton jump-diffusion American put.** L = BS part + λ∫(u(x+y) − u(x))ν(dy)
   with Gaussian jumps. The jump integral is a dense Toeplitz-like block per patch.
   Closed-form-adjacent oracles exist (Kou/Merton series solutions; or just PSOR on
   the dense system, matching the existing `obstacle_option_oracle.py` pattern).
2. **CGMY / tempered-stable put.** Same structure, heavier literature cachet.
3. **Fractional Laplacian obstacle problem** ((−Δ)^s, s ∈ (0,1)). The cleanest
   "dense Q is the problem" instance; Caffarelli–Silvestre gives analytic structure
   for validation.

Keep the local BS rows as the honest stress case, with one sentence in the README:
*sparse-Q problems are the regime where general sparse solvers win; the nonlocal
rows are the regime the block structure targets.*

### 2.2 Warm-started time-stepping chain (`obstacle/`)

The 200-QP American put chain is currently cold-started per step against OSQP's
free warm starts. Two additions, cheap first:

1. **Amortize setup:** reuse the symbolic analysis, `UzawaWorkspace`, and scalings
   across the chain (one sparsity pattern, 200 numeric solves). Right now the chain
   presumably pays `make_kkt`'s `symbolic()` 200 times.
2. **Interior warm start:** initialize step k+1 from step k's (x, s, y), pushed back
   into the interior (shift s and z by a multiple of e sized by the previous μ, or a
   Mehrotra-style starting point centered on the previous solution). Consecutive
   implicit-Euler QPs differ by O(Δt), so even a crude shift should cut iterations
   by half or more.

Report the chain three ways: cold, setup-reused, warm-started. That's an honest and
interesting result regardless of who wins. If the chain also moves to the Merton
model (§2.1), the IPM may win it outright.

### 2.3 Binned MLE as the headline config (`spline/mle_spline.jl`)

The data already supports this: leaves scale with **distinct sites**, and the
N=1000/58-leaf rows run 5–10x faster than the 1000-leaf rows. Binning with weights
is statistically standard (Papp & Alizadeh do exactly this) and the `binned` support
already exists. Make binned the primary rows; keep one unbinned row explicitly
labeled as the per-leaf-overhead stress test. Optionally add a quadrature-based
likelihood variant (fixed leaf count independent of N) as the asymptotic story.

### 2.4 Vector-stalk NUM (`num/`)

`num_sheaf` is a semantics showcase (dual = TCP, Kelly's identity), not a speed
benchmark — scalar stalks give the block machinery nothing to bite on. Two options,
not mutually exclusive:

1. Label the existing file a correctness/interpretability demo in the README and
   drop it from speed summary tables.
2. Add a **multipath / MIMO variant**: per-source rate *vectors* over k paths, or
   per-node power-allocation SOC blocks (Chiang et al.'s "layering as optimization
   decomposition" has ready-made formulations). This puts dense per-vertex blocks
   into the same Kelly narrative.

### 2.5 Dense-Gram Wasserstein rows (`map/wasserstein_graph.jl`)

The "lifted" rows (quadreg, laplace) have diagonal Q — the known anti-pattern. Add
non-lifted counterparts with dense per-vertex Grams at the same DOF so the table
shows the lifted/dense contrast explicitly rather than just losing. (See also
diagnostic D6: `smooth bary` loses at raug 1e5, which does **not** fit the pattern
and needs an explanation before reworking this file.)

### 2.6 Scaling sweeps in the winning regimes

A fixed 5x at one size is weak evidence; a growing gap is strong evidence. Add size
sweeps (log-spaced, ~4–5 points each) and plot time vs size for IPM and the best
reference solver:

- `map_lp` grids: 20×20 → 200×200 (both SparseMAP and dense-channel modes)
- `quantum_marginal`: N = 6 → 24 at fixed locality ℓ ∈ {2, 3}
- `map_sdp` moment SDP: cycle/grid families, growing |V|
- `tensor_spline` regress: 4×4 → 32×32

The fill-in argument predicts the gap **grows** for the graph SDPs — if measured,
that's the strongest claim in the whole suite, and it comes from adding rows, not
removing any.

---

## 3. Reporting adjustments (no benchmark changes)

1. **Split the summary into "target regime" and "stress tests."** The current
   Summary already lists both — formalize it, and state the regime rule explicitly
   (dense blocks, small H¹/DOF, high raug feasible). Reviewers reward a solver that
   knows its own boundary.
2. Add **H¹/DOF** and a **block-size summary** (min/median/max stalk dim) as columns
   or per-table metadata — these are the covariates the whole story hinges on.
3. Report **min** alongside mean over runs (`bench_min_ms` already exists;
   `run_comparison` currently averages). Min is the standard for solver timing since
   it filters GC/OS noise.
4. Don't lead with the `nonnegative_spline_exact` 15x (Mosek-primal artifact, §1).

---

## 4. Diagnostics to run

Each item names the open question it answers. D1–D3 are prerequisites — they affect
the validity of numbers already in `results.md`.

### D1. Tolerance parity (validity — run first) *(patch draft: Appendix B)*

**Question: are the OSQP comparisons apples-to-apples?**
The IPM runs at `feas_tol = gap_tol = 1e-8` everywhere, but no example sets
`eps_abs`/`eps_rel` on OSQP, so OSQP runs at its defaults of **1e-3** — five orders
looser. Mosek and Clarabel default near 1e-8, so those columns are fairer.

- Re-run every OSQP table at `eps_abs = eps_rel = 1e-8` (and once at 1e-5 for a
  middle point).
- Record each solver's *achieved* primal/dual residuals and gap next to its time
  (the IPM already exposes `pres`/`dres` via `IPMHistory`; JuMP exposes solver
  attributes).
- Expectation: obstacle and the lifted-Q rows shrink from 10–50x to something much
  smaller; some 0.8x rows may flip. This may be the single highest-leverage change
  to the results file.

### D2. What does `ipm_solve()` actually time?

**Question: is IPM timing penalized by assembly/symbolic work the JuMP side
excludes?** In `run_comparison`, `jump_build(optimizer)` sits *outside* `@elapsed`
(model build excluded), while the `ipm_solve` closure is timed whole. Audit each
example's closure: if it includes sheaf assembly or `make_kkt` symbolic analysis,
split into `init` (build + symbolic) and `solve!` (numeric) and report both. This
matters most for the small-DOF rows (opf, num) where setup is a large fraction.

### D3. Per-solve time breakdown

**Question: in the losing rows, where does the time go?** Instrument (or profile)
one solve per benchmark into: symbolic analysis / numeric factorization
(`cholesky!`) / Uzawa CG (`cg!`) / iterative refinement (`refine_kkt!`) / cone
operations / everything else.

- Hypothesis A (obstacle): CG + refinement dominate because H¹ ≈ DOF and raug=1e-2.
- Hypothesis B (num, small opf): fixed overhead (symbolic, allocation, dispatch)
  dominates because the problems are tiny.
- The fix differs completely between A and B, so measure before optimizing.

### D4. Surface inner iteration counts in the tables

**Question: does the raug/conditioning story hold quantitatively?**
`refine_kkt!` already returns `kkt_iters` and `IPMHistory` records `npred`/`ncorr`.
Add three columns to every results table: IPM outer iterations, total CG iterations,
mean CG iterations per KKT solve. Then scatter total-CG-per-DOF against P/IPM across
all existing rows — if the regime story is right, this one plot explains the whole
suite.

### D5. raug sweep

**Question: how sharp is the raug cliff, and are the per-example raug choices
optimal?** For ~6 representative rows (one winner and one loser per group), sweep
raug over 1e-2…1e7 (log grid) and plot solve time + CG iterations vs raug. Answers:
(a) whether losers lose at *every* raug (structural) or only at the chosen one
(tuning); (b) whether an adaptive raug heuristic (e.g., from ‖B‖² / block spectral
estimates — `UzawaWorkspace` already computes `nrm = ‖B‖²`) could replace the
hand-picked per-example values, which currently look like results of manual tuning.

### D6. The `smooth bary` anomaly

**Question: why does `smooth bary` (raug=1e5, DOF 1220, |V|=4) lose 2–5x when the
profile says it should win?** It's the one row that violates the pattern. Pull its
`IPMHistory` and CG counts: candidates are (a) many IPM iterations (bad centering
from the exp leaves), (b) a few huge dense blocks where |V|=4 means the multifrontal
tree is trivial and the factorization is effectively one dense Cholesky, (c) CG
stalling despite high raug. Whichever it is, it refines the regime rule — worth
knowing before writing §2.5's new rows.

### D7. Warm-start ablation for the chain

**Question: how much of the 0.23–0.45x chain loss is warm starting vs per-solve
speed?** Run the 200-step chain in the three modes of §2.2 and, separately, force
OSQP cold (`warm_start = false` per step). If cold-OSQP ≈ IPM per-step, the entire
gap is warm starting and §2.2 item 2 is the fix.

### D8. Block-density ablation on obstacle

**Question: quantify the sparse-block penalty directly.** On the perpetual put,
time one numeric factorization + one CG solve of the current dense-block treatment
against a sparse Cholesky of the same A + raug·BᵀB (e.g., via `SparseArrays`/CHOLMOD)
at the same raug. The ratio is the exact price of densifying tridiagonal blocks, and
it calibrates how much §2.1's dense-operator variants should recover.

### D9. Solution agreement check

**Question: are all solvers solving the same problem to the same answer?** For every
row, assert |obj_IPM − obj_ref| / (1 + |obj_ref|) below a threshold tied to the
tolerances, and record it. Cheap insurance; also catches dualization sign/scale bugs
before they contaminate timing conclusions. The per-example `*_oracle.py` ground
truths already exist for correctness — this extends the check to the benchmark loop
itself.

---

## 5. Suggested order of work

1. **D1 + D2** (validity of existing numbers) → regenerate `results.md`.
2. **D3 + D4** (breakdown + iteration columns) → confirms/refutes the regime story
   with data already flowing through the solver.
3. **§2.1 Merton obstacle + §2.2 warm starts** (the two biggest narrative wins),
   with **D7/D8** as their controls.
4. **§2.3 binned MLE, §2.6 scaling sweeps** (cheap, high payoff).
5. **D5, D6** (raug heuristic + anomaly) as time permits.
6. **§2.4, §2.5** last — they depend on what D6 and D8 reveal.

---

## Appendix A — Draft: jump-diffusion obstacle builder (§2.1)

Written against the PR-67 API as used in `examples_ipm/obstacle/obstacle_option.jl`.
Not executed; the oracle plan is at the end.

### A.1 Model choice: assemble the Dirichlet form, not the PIDE

Discretizing the Merton PIDE directly gives a **nonsymmetric** jump matrix (the
integrating factor ρ symmetrizes the local part but not a general jump kernel), and
a QP needs symmetric PSD A. Instead, assemble the **Dirichlet form of the
(censored) jump diffusion** directly:

    a(u,u) = ∫ p u′² + r∫ ρu²                       (local part — patch_grams, verbatim)
           + (λ_J/2) ∬_{Ω×Ω} √(ρ(x)ρ(x′)) φ_s(x−x′) (u(x)−u(x′))² dx dx′   (jump part)

The jump part is symmetric PSD **by construction** (it is a weighted graph
Laplacian with dense Gaussian weights) — no symmetrization fudge, and the VI
minimizing a(u,u) over u ≥ ψ is the standard Dirichlet-form formulation of the
obstacle problem for jump processes (Bensoussan–Lions).

Honesty note for the README: this prices a *symmetrized* jump model with censored
jumps, not the risk-neutral Merton parameterization verbatim. Two clean framings:
(a) present it as a nonlocal obstacle benchmark with finance flavor and validate
against its own dense-PSOR oracle; or (b) use the exactly-symmetrizable Merton
family — a Gaussian jump density N(m, s²) conjugates to a symmetric kernel iff
m/s² = μ/σ² (μ = r − σ²/2) — if exact PIDE semantics are wanted. Either way the
solver-relevant structure (dense PD blocks) is identical.

### A.2 Structure: kernel bandwidth vs patch overlap

Truncate the Gaussian kernel at radius R = ⌈nσ·s/h⌉ cells (nσ = 5 ⇒ error
~1e-6·λ_J). Interacting pairs (i,j) satisfy |i−j| ≤ R, and since patches are
contiguous windows, **every interacting pair shares a patch iff olap ≥ R** — the
natural compatibility condition, mirroring the existing element-split rule one
level up (pairs instead of intervals). The jump scale s is then the benchmark's
**density dial**: per-patch blocks range from banded (small s) to fully dense
(R ≥ patch width), while OSQP's KKT sees O(n·R) nonzeros.

### A.3 Code

Drop-in extension of `patch_grams` (same return signature), plus a one-keyword
refactor of the builder.

```julia
# ---- instance -------------------------------------------------------------

struct JumpObstacleInstance
    base::ObstacleInstance
    λJ::Float64        # jump intensity
    s::Float64         # jump-size std dev (log-price units)
    nσ::Float64        # kernel truncation radius, in units of s
end

jump_obstacle_instance(; λJ = 0.5, s = 0.15, nσ = 5.0, kwargs...) =
    JumpObstacleInstance(obstacle_instance(; kwargs...), λJ, s, nσ)

# ---- per-patch Grams: local (verbatim) + jump Dirichlet form ---------------

"""Per-patch dense Grams for the symmetrized jump-diffusion obstacle problem.
Jump part assembled pairwise with 1/#owning-patch splitting (the `eown`
pattern lifted from intervals to pairs), so Σ_p embed(Q_p) = A_full exactly —
same invariant patch_grams already guarantees for the local part.
Requires olap ≥ kernel radius (asserted): every interacting pair must share
a patch, or cross-patch couplings would be orphaned."""
function jump_patch_grams(jinst::JumpObstacleInstance; dt = nothing)
    inst = jinst.base
    pat, Qs, Ms = patch_grams(inst; dt)            # local part, unchanged
    x = grid(inst); n = inst.n; h = x[2] - x[1]
    μ = inst.r - inst.σ^2 / 2
    ρ = exp.(2μ .* x ./ inst.σ^2)
    R = ceil(Int, jinst.nσ * jinst.s / h)
    @assert inst.olap ≥ R "olap ($(inst.olap)) < kernel radius ($R): raise olap or shrink s"
    φ(y) = exp(-y^2 / (2 * jinst.s^2)) / (√(2π) * jinst.s)
    # patches are contiguous windows ⇒ patches owning BOTH i < j are
    # exactly those with a ≤ i and b ≥ j
    npown(i, j) = count(((a, b),) -> a ≤ i && b ≥ j, pat)
    scale = dt === nothing ? 1.0 : dt              # time stepping: Q ↦ M + dt·A
    for (p, (a, b)) in enumerate(pat)
        Qp = Qs[p]
        for i in a:b, j in (i + 1):min(b, i + R)
            w = scale * jinst.λJ * √(ρ[i] * ρ[j]) * φ(x[j] - x[i]) * h^2 / npown(i, j)
            ii = i - a + 1; jj = j - a + 1
            Qp[ii, ii] += w; Qp[jj, jj] += w
            Qp[ii, jj] -= w; Qp[jj, ii] -= w
        end
    end
    return pat, Qs, Ms
end
```

Wiring: give `build_obstacle_problem` a `grams` keyword defaulting to the current
behavior — a two-line diff:

```julia
function build_obstacle_problem(inst::ObstacleInstance; dt = nothing,
                                uprev = nothing, grams = patch_grams)
    ...
    pat, Qs, Ms = grams(inst; dt)          # was: patch_grams(inst; dt)
```

```julia
# call site
jinst = jump_obstacle_instance(n = 801, P = 4, olap = 60, λJ = 0.5, s = 0.15)
prob, info = build_obstacle_problem(jinst.base;
                 grams = (inst; dt = nothing) -> jump_patch_grams(jinst; dt))
```

The c-vector logic is untouched: `c[p] = Qs[p] * ψ[a:b]` still assembles to A·ψ
globally because each Q_p only couples within-patch indices and Σ embed(Q_p) = A.
The time-stepping chain and the game mode inherit the jump Grams for free through
the same keyword.

### A.4 Validation & oracle plan

1. **Splitting exactness:** assert ‖Σ_p embed(Q_p) − A_full‖∞ ≤ 1e-12, A_full
   assembled globally without splitting (extends the existing ≤ 3.5e-10 check).
2. **λ_J → 0 continuity:** solution → Merton closed form as λ_J → 0 (reuse
   `merton_perpetual`); report ‖u_λ − u_Merton‖∞ at λ_J ∈ {1e-3, 1e-2}.
3. **Dense PSOR oracle:** extend `obstacle_option_oracle.py` with the same pairwise
   assembly (NumPy, dense) + projected SOR — the Cryer check, verbatim, on dense A.
4. **Monotonicity smoke test:** put value nondecreasing in λ_J at fixed spot.
5. **Contact-measure check:** cone dual ≥ 0 on contact, ~0 off it, as in the
   existing file.

### A.5 Suggested benchmark grid

| Sweep | Values | What it shows |
|---|---|---|
| jump width s (n=801, P=4) | 0.05, 0.10, 0.20, 0.40 | density dial: IPM flat, OSQP degrades with nnz |
| patches P (s=0.15) | 1, 4, 8 | P=1 is one dense Cholesky vs OSQP-dense; P>1 tests the coupling overhead |
| time-stepping chain (200 steps) | cold / setup-reuse / warm (§2.2) | interacts with D7 |

Report the per-patch density fraction as a column, and keep the λ_J = 0 rows —
they are the existing (honestly losing) local benchmark, now serving as the
regime-boundary control within the same table.

---

## Appendix B — Draft: tolerance-parity patch to `benchmark_utils.jl` (D1)

Three changes: a `--tol` CLI flag, tolerance-pinned optimizer factories, and
achieved-quality reporting in `run_comparison`. Plus one required cleanup outside
the file.

### B.1 CLI flag

```julia
# parse_benchmark_args(): add alongside nruns/nwarmup
    tol = 1e-8
    ...
    elseif startswith(arg, "--tol=")
        tol = parse(Float64, split(arg, "=")[2])
    ...
    return (open = open && !mosek, mosek = mosek,
            nruns = nruns, nwarmup = nwarmup, tol = tol)
```

Match the IPM side in each example: `IPMSettings(feas_tol = opts.tol,
gap_tol = opts.tol, ...)` so one flag moves all solvers together (useful for a
1e-5 midpoint run).

### B.2 Optimizer factories with pinned tolerances

```julia
using JuMP: optimizer_with_attributes

"""Optimizers with termination tolerances pinned to `opts.tol` so all
solvers chase the same accuracy. Note OSQP's DEFAULTS are eps_abs =
eps_rel = 1e-3 — five orders looser than the IPM's 1e-8 — so every
existing OSQP column was an unequal-tolerance comparison."""
function get_optimizers(opts; lp_only::Bool = false)
    tol = opts.tol
    primal = if opts.mosek
        optimizer_with_attributes(Mosek.Optimizer,
            "MSK_DPAR_INTPNT_TOL_PFEAS"      => tol,   # LP/QP path
            "MSK_DPAR_INTPNT_TOL_DFEAS"      => tol,
            "MSK_DPAR_INTPNT_TOL_REL_GAP"    => tol,
            "MSK_DPAR_INTPNT_CO_TOL_PFEAS"   => tol,   # conic path
            "MSK_DPAR_INTPNT_CO_TOL_DFEAS"   => tol,
            "MSK_DPAR_INTPNT_CO_TOL_REL_GAP" => tol)
    elseif lp_only
        optimizer_with_attributes(OSQP.Optimizer,
            "eps_abs"  => tol,
            "eps_rel"  => tol,
            "max_iter" => 200_000,   # default 4_000 will not reach 1e-8
            "polish"   => true)      # active-set polish; usually required at 1e-8
    else
        optimizer_with_attributes(Clarabel.Optimizer,
            "tol_feas"    => tol,
            "tol_gap_abs" => tol,
            "tol_gap_rel" => tol)
    end
    return primal, Dualization.dual_optimizer(primal)
end
```

Notes:
- `Dualization.dual_optimizer` accepts an `OptimizerWithAttributes` factory, so
  the dual inherits the same tolerances.
- OSQP `polish` adds a factorization to each solve — report obstacle once with
  and once without, since at 1e-8 unpolished OSQP frequently stalls at
  `ITERATION_LIMIT` (which must then be recorded as a non-solve, not a time).
- Mosek/Clarabel defaults are already ≈1e-8; pinning them is mostly about making
  the `--tol=1e-5` sweep meaningful.
- **Required cleanup:** several examples construct optimizers inline instead of
  via `get_optimizers` (e.g. `map_lp.jl:725-731`, `map_sdp.jl:841-847`). Route
  those through `get_optimizers(opts; lp_only = ...)` or the pinning silently
  won't apply.

### B.3 Achieved-quality reporting in `run_comparison`

Record what each solver *achieved*, not just how long it took; warn on
objective disagreement (diagnostic D9). Also switch to reporting min over runs
alongside the mean (§3.3).

```julia
function run_comparison(label::String, ipm_solve::Function, jump_build::Function;
                        optimizer, dual_optimizer = nothing,
                        nwarmup::Int = 2, nruns::Int = 5, verbose::Bool = true,
                        ipm_objective::Union{Nothing, Function} = nothing,
                        agree_tol::Float64 = 1e-5)
    # IPM
    res = nothing
    for _ in 1:nwarmup; res = ipm_solve(); end
    t = Float64[]
    for _ in 1:nruns
        push!(t, @elapsed (res = ipm_solve()))
    end
    t_ipm, tmin_ipm = 1000mean(t), 1000minimum(t)
    ipm_pres = res.history.pres[end]; ipm_dres = res.history.dres[end]
    @assert res.status in (OPTIMAL, NEAR_OPTIMAL) "IPM: $(res.status) on $label"

    # reference solver (primal shown; dual identical)
    t_primal = Inf; obj_p = NaN; status_p = nothing
    if optimizer !== nothing
        m = jump_build(optimizer)
        for _ in 1:nwarmup; optimize!(m); end
        tp = Float64[]
        for _ in 1:nruns
            m = jump_build(optimizer)
            push!(tp, @elapsed optimize!(m))
        end
        t_primal = 1000mean(tp)
        status_p = termination_status(m)
        obj_p    = objective_value(m)
        status_p == MOI.OPTIMAL ||
            @warn "reference solver did not reach OPTIMAL" label status_p
    end

    # D9: objective agreement (ipm_objective(res) = ½pᵀQp + cᵀp; the examples
    # hold prob, so pass e.g. res -> obj(prob, res) from the call site)
    if ipm_objective !== nothing && !isnan(obj_p)
        rel = abs(ipm_objective(res) - obj_p) / (1 + abs(obj_p))
        rel > agree_tol && @warn "objective disagreement" label rel
    end
    ...
end
```

And per results table, add columns: IPM `pres`/`dres` at termination, reference
solver termination status, and the D4 iteration counts (outer iters, total CG).
A row where the reference solver reports `ITERATION_LIMIT` or `ALMOST_OPTIMAL`
gets footnoted rather than counted as a loss/win.

### B.4 What to expect

The obstacle OSQP columns (currently 10–50x against) and the lifted-Q wasserstein
rows should compress substantially at equal 1e-8; some 0.7–0.9x rows elsewhere may
flip either way. Regenerate `results.md` in full after this patch — it supersedes
the 2026-07-02 numbers — and keep the old file as `results_2026-07-02_uneq_tol.md`
for the record.
