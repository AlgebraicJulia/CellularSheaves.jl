# snapshots3 — shifted-regime captures

git: `940e61bcec`  ·  captured 10/10  ·  wall 118.3s
Probed ρ-transitions (last-unshifted → first-shifted log α):
- 07_hsd : 9.9 → 10.0
- 07_ipm : 10.8 → 10.9
- e15_ipm : 11.8 → 11.9

One directory per (problem, iteration, α). A/B are α-independent (identical across a problem's three captures); vectors/meta are per-α.

| problem | solver | it | log α | ρ_engaged | ρ_value | A (n×n) | B (m×n) | fid p | fid c | fid w |
|---|---|---|---|---|---|---|---|---|---|---|
| 07 | ipm | 4 | 10.7 | false | 0.00e+00 | 1452×1452 | 1096×1452 | 8.6e-03 | 1.7e-02 | — |
| 07 | ipm | 4 | 10.8 | false | 0.00e+00 | 1452×1452 | 1096×1452 | 9.5e-06 | 2.7e-05 | — |
| 07 | ipm | 4 | 10.9 | true | 1.00e-09 | 1452×1452 | 1096×1452 | 2.6e-08 | 2.5e-09 | — |
| 07 | ipm | 4 | 12.0 | true | 1.00e-09 | 1452×1452 | 1096×1452 | 7.4e-09 | 3.1e-08 | — |
| 07 | hsd | 3 | 9.9 | false | 0.00e+00 | 1452×1452 | 1096×1452 | 1.1e-05 | 3.7e-04 | 8.6e-04 |
| 07 | hsd | 3 | 10.0 | true | 1.00e-09 | 1452×1452 | 1096×1452 | 5.7e-09 | 3.0e-10 | 1.8e-10 |
| 07 | hsd | 3 | 12.0 | true | 1.00e-09 | 1452×1452 | 1096×1452 | 2.0e-08 | 2.6e-08 | 2.7e-09 |
| e15 | ipm | 3 | 11.8 | false | 0.00e+00 | 4128×4128 | 3168×4128 | 2.1e-06 | 3.9e-05 | — |
| e15 | ipm | 3 | 11.9 | true | 1.00e-09 | 4128×4128 | 3168×4128 | 3.7e-08 | 5.4e-08 | — |
| e15 | ipm | 3 | 13.0 | true | 1.00e-09 | 4128×4128 | 3168×4128 | 3.7e-08 | 5.3e-08 | — |

**Fidelity: 9/10 captures pass ≤ 1e-2** (ρ-shift applied where ρ_engaged). The exception is `07_ipm_it4_a10.7` corrector = 1.7%: α=10^10.7 is the spec's *guessed* "last-before-shift" alpha but is two rungs below the true transition (10.8→10.9), where the **unshifted F is at cond ≈ 1.1e16 ≈ 1/eps** — the edge of numerical singularity that makes the ladder engage just above. The offline generic `\` loses ~1.7% on the corrector there; the solver's chordal factorization (which got `r0_c` fine at ρ=0) is more accurate. The captured inputs are correct — a reconstruction-conditioning artifact at the singularity boundary, not a capture error. The **tight bracket 10.8→10.9 is clean** (9.5e-6 → 2.6e-8).

## Transition bracketing (test 3)

- 07 ipm it4: ρ_engaged false@10.7 → true@12.0; tight bracket 10.8(false)→10.9(true) ✓
- 07 hsd it3: ρ_engaged false@9.9 → true@12.0; tight bracket 9.9(false)→10.0(true) ✓
- e15 ipm it3: ρ_engaged false@11.8 → true@13.0; tight bracket 11.8(false)→11.9(true) ✓

## Notes
- rho_engaged = (the applied shift row.ρ > 0), now that init_uzw! reports the ACTUAL shift (0 when the unshifted factorization succeeds, the ladder rung otherwise) rather than rgmin. Fidelity rebuilds F with that same ρ (0 unshifted, applied shift when engaged), confirming ρ is recorded correctly (acceptance test 1). All shifts landed on the first rung, ρ=rgmin=1e-9.
- Solver change (reporting only, behavior-preserving): init_uzw! now returns ρ_applied instead of the ladder's lower bound rgmin. ρ is report-only (no solver logic reads it), so this cannot change convergence. NOTE: existing sweeps (fine2/, oracle*/, snapshots/, snapshots2/) were generated before this fix and still carry the old rgmin-in-the-ρ-column convention.
- Extra meta fields per spec: rho_engaged, rho_value, norm_dp, norm_dy, r0_{p,c,w}, {p,c,w}_res0_dual, nc, ng.
- Capture point / A-B-vector conventions identical to snapshots/ and snapshots2/.
