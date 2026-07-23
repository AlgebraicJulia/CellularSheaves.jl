# IPM solver — project ledger

Clean-room rebuild of the HSD conic IPM solver on `16c57d2` ("comp done"), per
`examples/final_spec.md`. Every mechanism proposed across the diagnostic → woodbury
→ raug → oracle → fix3 → ship → finish-line campaigns is mapped below to its final
disposition. The campaign evidence lives on branch `rhs/examples` (reports:
`DIAGNOSTIC_REPORT.md`, `RAUG_REPORT.md`, `ORACLE_REPORT.md`, `FIX3_REPORT.md`,
plus the `*_spec.md` design records).

## Shipped (C1–C6, this branch)

| # | Change | Mechanism | Evidence |
|---|---|---|---|
| C1 | Un-gated refinement | `refinehsd!` loses the `dynam_tol` floor exit (4–9 decades pessimistic); terminates on stall or `refine_itmax=10`; `floor_tol` only as atol backstop | FIX3_REPORT; P5(V10) floor 5.96e-9→~7.6e-10 @raug=1e2, 3.50e-12 @raug=1 |
| C2 | Mehrotra poor-affine override | `αa < α_bad=1e-2` ⇒ drop 2nd-order term + centering `σ_safe·μ` (σ_safe=0.5). No σ-clamp/floor. `nsafe` history column | B-1 grid; fires 0× absent the collapse signature |
| C3 | Predictor enforcement | Predictor refined to `force_tol` on the 3-row residual, same `refinehsd!` path, unconditional/self-gating, `η_pred=1`, no cap | fix3 exp-3; raug=10→OPTIMAL, raug=1→OPTIMAL (needs C2, mandatory order) |
| C5 | Honest reporting (revised — addendum) | `HSDResult.{mu,pres,dres,pobj,dobj}` in USER frame (pres/dres recomputed by un-scaling the residual *vector* — the stored embedding norm is wrong by up to 5.5× under equilibration); history gains embedding `gap = rτ−κ` column; `floored` dropped (accept-exit re-founds NEAR_OPTIMAL); infeasibility carve-out (certificates in obj slots) | DIAGNOSTIC I3 [CORR2]; addendum; classification 5/5 (status+gap, no flag) |
| C6 | G-form-native Uzawa | `G = βA + BᵀB`, β=1/α; the F-form does not exist. ρI ladder correctly scaled (‖G‖≈nB²) | B-4; exact-dual + corrupted-ỹc; resolves P6_base raug=1e6/1e10 init-failure to machine residual |

## Dropped by rule (measured harmful or inert)

| Mechanism | Why | Evidence |
|---|---|---|
| **Demand cap** | No waste to remove at tol=1e-8 (Δpasses=0); *actively harmful* at tol=1e-10 (tolerance-miss) | B-2 |
| σ-clamp | Byte-inert on the reproducer; the override is the active ingredient | B-1 (σmax 0.9≡0.99) |
| σ-floor (σ_min) | Steers — destabilizes raug=100 (spurious fires, 7.7× regression) | B-1 |
| `reltol` correction tol | Regressed refinement contraction | fix3 |
| `rowfloor` as behavior | Byte-for-byte redundant once ungate removes the floor exit | fix3 |
| Stall-gated predictor | Superseded — mechanism adjudicated (B) bad-region; safeguard+invariant supersedes | ship D-STOP → B-1 |
| Accept/reject veto | Never fires on its target (3-row-invisible corruption) | fix3 |
| Best-so-far / patience stall | Backfire — refinement 3-row res is not a proxy for outer step quality | fix3 |
| Woodbury guard/retry, `S_schur` division | No failure carried the signature; `R` unusable as raw hazard signal | woodbury V4/V2 |

## Deferred with evidence intact (back to the bench, not abandoned)

**Down-only raug backoff (would-be C8).** Core validated: V1 (backoff lands
OPTIMAL from the NEAR_OPTIMAL baseline), V2 (undershoot → `floored` accept-and-
record), the `raug_new = raug·(demand/stalled)` formula (no overshoot), and the
warm-start safety mechanism (safe by endgame `nH`-inflation timing: `α_eff = 788`
at nominal raug=3e-4). **Deferred at integration**: the trigger semantics are
demonstrated *incomplete* by the A2 gate — the bare `refstat==REFINE_STALLED`
signal fires on transient mid-solve stalls (P4 diverges: raug 1e2→3e-4 in one
premature hop), and the single-backoff-per-solve rung cannot reach demand from
raug=1e6 in one hop (P5 @ raug=1e6 stays at its 3e-7 floor). B-5's 3-iteration
pres-plateau qualifier — which `final_spec` dropped as "a trajectory heuristic" —
was load-bearing; removing it is the measured defect. **Suite-2 chartered** with
fixtures attached: P4(V20)@raug=1e2 (premature-trigger reproducer), P5(V20)@raug=1e6
(single-hop-reach reproducer), plus the V1/V2 regressions. Re-integrates with a
plateau-qualified trigger and a decision on multi-hop vs accept-and-record.

## Open (carried to the queue)

- **Machine-zero / STALLED-label edge.** A point driven to machine residual is
  labelled STALLED by the status enum. C5's user-frame `pres`/`dres` (~1e-14) report
  the truth, so a `status + gap` reading calls it solved — but the *status enum* is
  wrong. Per the C5 addendum this is now **assertable, not narrated**: the
  row-correspondence invariant (every exit path pushes a terminal row = the returned
  point; debug assertion ships) makes the anomaly's fingerprint — a classification
  reading state that does not describe the returned point — impossible, and its audit
  is chartered to resolve-or-file this. A status-classifier fix, not a solver fix.
- **Second P6 component:** raug=1e8 fails both forms (P6_base floor ~6.8e-5) — G-form
  resolves 1e6 and 1e10 but not 1e8. Non-monotone, undiagnosed.
- **Equilibration rewrite (maintainer-owned):** per-coordinate scaling on
  orthant/free blocks removes the intra-block span defect (P3_span2/5 → OPTIMAL at
  defaults). Detector shipped in B-6 (on `rhs/examples`); the fix is parked. Validate
  the intra-block ratio on the *input* `prob.B/Q` (Ruiz row-scaling masks it post-eq).
- **Schur-complement preconditioner:** metric/hard-fixture/warning-label pre-attached.
- **`comp/` extended-precision wiring:** gated on a sub-`F_eval` requirement + the
  BigFloat oracle tests; `comp/` kernels remain in `16c57d2`, unwired.

## A2 — outcome equivalence (rebuild C1–C6 vs experimental prototype, V=20,d=50, no backoff)

Reference: `rhs/examples` mirror config (ungate+override+enforce+G-form). Statuses
identical / floors within 2× — or rebuild strictly better. **No cell regresses.**

| fixture | cell | rebuild | reference | verdict |
|---|---|---|---|---|
| P5 | 1e2/1e-8 | OPTIMAL 5.2e-9 | OPTIMAL 5.2e-9 | identical |
| P5 | 1e2/1e-10 | NEAR 1.7e-10 | NEAR 1.66e-10 | identical (2×) |
| P5 | 1e6/1e-8 | NEAR 3.1e-7 | NEAR 3.12e-7 | identical |
| P5 | 1e6/1e-10 | NUME 3.1e-7 | NUME 3.12e-7 | identical |
| P4 | 1e2/{1e-8,1e-10} | OPTIMAL/NEAR ~1.1e-10 | OPTIMAL/NEAR ~2.1e-10 | identical (2×) |
| P4 | 1e6/{1e-8,1e-10} | OPTIMAL/NEAR 9.7e-10 | OPTIMAL/NEAR 2.9e-9 | identical status, rebuild better (3×) |
| P1_k16 | 1e2/1e-8 | **OPTIMAL 3.4e-13** | NEAR 1.35e-8 | rebuild better (machine) |
| P1_k16 | 1e2/1e-10, 1e6/* | STALLED ~1e-14 (machine) | NEAR/NUME ~1e-7 | rebuild better; STALLED = the machine-zero label edge |
| P3s2 | 1e2/{1e-8,1e-10} | NUME 3.9e-2 | NUME 1.3e-1 | identical (genuine span failure — expected) |
| P3s2 | 1e6/{1e-8,1e-10} | **OPTIMAL 9.7e-{10,12}** | NEAR/NUME 3.15e-7 | rebuild better |

The rebuild is equivalent (P5, P4) or strictly better (P1_k16, P3s2@1e6, reaching
machine/OPTIMAL where the toggle-laden prototype floored); the only non-identical
statuses that are *not* improvements are the machine-zero label edge (open above).
The backoff's absence leaves the high-raug cells at their natural raug floors —
the equivalence held without it, which is why C8 could be cleanly deferred.
