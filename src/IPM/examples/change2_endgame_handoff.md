# Handoff: Change 2 (HSD bordered) endgame cost regression

**Commit:** `acd98174` on `rhs/forcing`. **Scope:** HSD only (IPM untouched).

## TL;DR

Change 2 made the Woodbury column a bare backsolve, then drives it per-consumer to
targets `ptgt = (ptol−rb)/|dt|`, `ytgt = ytol/|dt|`. On hard problems (07, e08) the
**corrector's `|dt|` is O(1) in the endgame**, which collapses those targets to
machine precision and forces expensive column solves every step. This is a real,
structural solver-cost regression, measured at the **oracle best-α** (i.e. not a
controller artifact). Change 1 is clean (identical best-α cost); the live-controller
blowups (07: 3004→42223) are a *separate* controller-input problem, deferred.

## Best-α cost (oracle, fix_alpha, ground truth — controller excluded)

| best-α total CRAIG | 07 | e08 |
|---|---|---|
| pre-Change-1 (`5ce6646a`) | 82 | 47 |
| Change 1 (`e344d906`) | 82 | 68 |
| Change 2 (`acd98174`) | **140** | **303** |

Change 1 is cost-neutral on the solver. Change 2 raises the endgame per-step cost.

## Mechanism (verified from a traced best-α replay of e08)

`dt = Δτ = num/S` is the Schur-lift coefficient for the τ (border) row; the full
direction is `Δp = Δp0 + dt·Δp2`, so the column residual enters the assembled rows as
`dt·r2` — hence the column must be solved to `≈ tol/|dt|`.

**Why `|dt|` is O(1) for the corrector but ~1e-4 for the predictor:** it is `num`, not
`S`. Per-iter (e08, corrector = the driver; predictor stays loose all the way):

| iter | α | κ | ptol | ytol | num | S | \|dt\| | ytgt | ncraig |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 3.2e1 | 1.0 | 3.0e-1 | 9.0e-1 | 46.5 | 68.9 | 0.67 | 1.3 | 3 |
| 3 | 3.2e2 | 1.7e-2 | 4.1e-3 | 1.2e-2 | 7.4e-3 | 6.7e-1 | 0.011 | 1.1 | 3 |
| 9 | 3.2e6 | 2.4e-9 | 1.4e-9 | 4.1e-9 | -3.4e-9 | 2.4e-7 | 0.014 | 2.9e-7 | 3 |
| 10 | 3.16 | 3.6e-10 | 2.2e-10 | 6.4e-10 | -5.2e-10 | 5.2e-6 | 1.0e-4 | 6.4e-6 | **122** |
| 11 | 1e7 | 3.5e-10 | 2.1e-10 | 6.2e-10 | 2.8e-8 | 3.3e-8 | **0.86** | 7.2e-10 | 22 |
| 12 | 3.2e7 | 3.5e-10 | 2.1e-10 | 6.2e-10 | 6.8e-8 | 2.4e-8 | **2.87** | 2.1e-10 | 22 |
| 13 | 3.2e6 | 3.5e-10 | 2.1e-10 | 6.2e-10 | 5.4e-8 | 3.6e-8 | **1.52** | 4.1e-10 | 22 |
| 14 | 1e7 | 3.5e-10 | 2.1e-10 | 6.1e-10 | 2.6e-8 | 3.0e-8 | **0.88** | 7.0e-10 | 22 |
| 15 | 3.2e8 | 2.9e-10 | 1.8e-10 | 5.5e-10 | 1.1e-7 | 3.0e-8 | **3.80** | 1.4e-10 | 22 |
| 16 | 1e9 | 2.5e-10 | 1.9e-10 | 5.7e-10 | -8.9e-6 | 2.1e-8 | **417.7** | 1.4e-12 | 22 |
| 17 | 1e8 | 2.7e-10 | 1.7e-10 | 5.1e-10 | 4.0e-8 | 2.7e-8 | **1.50** | 3.4e-10 | 22 |

(predictor `num` at iters 11–17: ~1e-11 to 1e-12, i.e. `dt_pred ≈ 1e-3–1e-4.)

Reading the columns:
- **`S` collapses ~9 orders** (68.9 → ~2e-8) as the iterate converges: `S = δ + Δp2ᵀHΔp2 − 2QpᵀΔp2/τ`, `δ = pᵀQp/τ² + κ/τ`; κ falls 10 orders (1.0 → 2.5e-10), and `S` plateaus at ~2–3e-8.
- **`num_corr` does NOT vanish** — it plateaus at ~1e-8…1e-6, while `num_pred` → ~1e-11. The gap is the **Mehrotra second-order term** `(σμ − Δτa·Δκa)/τ` in the corrector's τ-RHS, which the 2-row base solve never touches (it doesn't see the τ row).
- So `dt_corr = num_corr/S ≈ (1e-8)/(2e-8) ≈ O(1)`, and `ytgt = ytol/|dt| ≈ ytol ≈ 1e-10` — the column must match the Newton tolerance itself → ~22 CRAIG/step. `dt_pred = (1e-11)/(2e-8) ≈ 1e-3` → loose → cheap.

**Not caused by:** S collapsing alone (predictor sees the same S and stays cheap), or `num` growing unbounded (it's O(1e-8), plateau not blowup; the single 417 spike at iter 16 has `ytgt < floor`, so the loop breaks there). It's the **ratio**: corrector `num` held up by the Mehrotra term, over a collapsed `S`.

**Iter-10 anomaly (`ncraig=122`):** oracle chose α=3.16 (tiny). At tiny α the augmented system `F=(1/α)H+BᵀB` is differently conditioned and the *base* solve is expensive there — unrelated to the column loop (`dt≈1e-4`, targets loose).

## Contrast with pre-Change-2 code

Old code solved the column **once to a fixed working tol** and let `refinehsd!` (3-row
IR) absorb the `dt·r2` injection — bounded cost regardless of `dt`. Change 2 replaced
that with "solve the column exactly as tight as `dt` demands," which inverts in the
endgame: the premise ("column over-solved on fixed data") holds early (small dt) but
fails late (dt≈O(1) ⇒ column driven to machine tol).

## Suggested fix direction (not yet implemented)

Break the `1/|dt|` coupling: cap how tight the column is ever driven (never below the
base tolerance) and let `refinehsd!` absorb the residual when `dt` is large — i.e.
recover the old bounded-cost behavior for large `dt` while keeping the loose-early
benefit. Needs measuring at best-α across 07/e08/e04 to confirm it removes the endgame
cost without breaking correctness.

## Reproduction

- Best-α replay: read `chosen` α per iter from `oracle/<key>_1e-10_hsd_c2.csv`, `init`
  with `fix_alpha=true`, set `s.α[]` per iter, call `IPM.step!`. Instrument the entry
  of the column loop in `src/kkt/bordered.jl:solvekkt!(bw,…)` to print
  `α, τ, κ, ptol, ytol, rb, num, S, dt, ptgt, ytgt`.
- Oracle sweeps used: `07_1e-10_hsd_c2.csv`, `e08_1e-10_hsd_c2.csv` (+ `_c1`, `_pre1`).

## Related, out of scope here

- **Live-controller blowup** (07: 3004→42223): Change 2 feeds the α-window the loose
  bare-backsolve residual (`wfres`), steering the controller to bad α. Deferred — plan
  is oracle-clean the solver first, fix the controller last.
- **Uninitialized-`y` bug** (already fixed in `acd98174`): the `craig=false` backsolve
  skipped `craig!`, which is the only writer of `y`; the dual recovery read
  uninitialized `y` → per-process nondeterministic blowups. Fixed by zeroing `y`.
