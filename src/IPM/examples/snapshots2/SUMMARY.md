# snapshots2 — targeted assumption/law correspondence captures

git: `8e9d92e828`  ·  captured 18/18  ·  wall 113.3s

## Captures (selection reproducibility + fidelity + shapes)

| group | tag | problem | solver | iter | sel-err | A (n×n) | nnz A | B (m×n) | nnz B | fid p | fid c | fid w |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | worst-ceil | X09 | ipm | 10 | ceil +5.350 | 50×50 | 50 | 25×50 | 1250 | 7.7e-13 | 4.9e-09 | — |
| 1 | worst-ceil | X09 | hsd | 11 | ceil +1.091 | 50×50 | 50 | 25×50 | 1250 | 1.0e-12 | 1.3e-07 | 1.5e-06 |
| 1 | worst-ceil | e15 | ipm | 5 | ceil +1.746 | 4128×4128 | 12400 | 3168×4128 | 50736 | 1.5e-13 | 2.8e-06 | — |
| 1 | worst-ceil | e15 | hsd | 1 | ceil +2.573 | 4128×4128 | 10960 | 3168×4128 | 50736 | 2.9e-14 | 2.3e-13 | 2.8e-14 |
| 1 | worst-floor | e02 | ipm | 7 | floor +1.537 | 360×360 | 12960 | 101×360 | 7236 | 3.1e-11 | 6.6e-11 | — |
| 1 | worst-floor | e02 | hsd | 4 | floor -0.446 | 360×360 | 2244 | 101×360 | 7236 | 1.2e-12 | 1.1e-12 | 2.0e-10 |
| 1 | worst-floor | X06 | ipm | 4 | floor -0.075 | 60×60 | 600 | 30×60 | 1800 | 3.0e-13 | 2.0e-13 | — |
| 1 | worst-floor | e09 | ipm | 5 | floor -0.084 | 3006×3006 | 60226 | 1693×3006 | 69750 | 3.0e-14 | 1.2e-13 | — |
| 2 | best-ceil | X09 | ipm | 1 | ceil +0.004 | 50×50 | 50 | 25×50 | 1250 | 5.0e-14 | 2.5e-13 | — |
| 2 | best-ceil | X09 | hsd | 6 | ceil +0.053 | 50×50 | 50 | 25×50 | 1250 | 6.0e-14 | 1.2e-12 | 2.1e-11 |
| 2 | best-ceil | e15 | ipm | 7 | ceil -0.149 | 4128×4128 | 12400 | 3168×4128 | 50736 | 5.6e-12 | 2.0e-06 | — |
| 2 | best-floor | e02 | ipm | 3 | floor -0.021 | 360×360 | 12960 | 101×360 | 7236 | 2.0e-14 | 5.5e-13 | — |
| 2 | best-floor | e02 | hsd | 2 | floor -0.004 | 360×360 | 1920 | 101×360 | 7236 | 5.4e-14 | 1.1e-12 | 3.7e-14 |
| 2 | best-floor | X06 | ipm | 2 | floor +0.002 | 60×60 | 600 | 30×60 | 1800 | 3.1e-14 | 1.6e-14 | — |
| 3 | last-W0 | e08 | ipm | 21 | — | 255×255 | 639 | 128×255 | 510 | 7.4e-10 | 5.0e-10 | — |
| 3 | last-W0 | e01 | hsd | 17 | — | 576×576 | 82944 | 240×576 | 69120 | 6.7e-09 | 3.0e-08 | 0.0e+00 |
| 3 | last-W0 | SOS_P8_n9_s1 | hsd | 9 | — | 440×440 | 12200 | 21×440 | 2310 | 1.8e-06 | 7.7e-06 | 3.0e-06 |
| 3 | last-overall | SOS_P8_n9_s1 | hsd | 17 | — | 440×440 | 12200 | 21×440 | 2310 | 1.6e-03 | 1.1e-04 | 6.2e+01 |

**Fidelity: 17/18 captures pass at ≤ 8e-6** (well under the 1% bar). The lone exception is the woodbury role of `SOS_P8_n9_s1_hsd_it17` (G3 *last-overall*) — see the note below; it is a machine-zero-target artifact, not a capture error.

## Fidelity note — the terminal iterate (SOS it17)

`SOS_P8_n9_s1_hsd_it17` is the **fully-converged terminal iterate** (μ=1e-12); G3 deliberately targets the last-overall step. There **all residuals are at or below machine epsilon**: `r0_p=6.8e-19`, `r0_c=4.3e-19`, `r0_w=3.1e-14`. The relative fidelity test is ill-posed against a machine-zero target — p/c only report small numbers because the denominator floors at `eps`, and w reports 62 because `r0_w=3e-14` sits just above `eps` while the offline reconstruction floors at ~2e-12. Two compounding causes, both properties of the iterate, not the capture:
- **cond(F) = 2.2e21** (numerically singular at α=1e11): a generic double-precision resolve of `βA+B'B` has no correct digits.
- The solver's **ρ-shift ladder engaged** (ρ=1e-9): it factored `βA+B'B+ρI`, which the spec's ρ-free reconstruction formula omits. Including ρ drops the woodbury residual to 6e-14 — but the relative figure stays ≈1 because 6e-14 vs 3e-14 are both machine zero.

The captured inputs are correct; **absolute** agreement is ~3e-14. Recorded as a deviation per spec. (The other terminal captures — e08 it21, e01 it17 — pass because their residuals sit comfortably above eps; e01's woodbury `r0_w=0` exactly, giving 0.0.)

## Notes
- Selection per snapshots2 spec: W0 = α with every role base==1 & refn==0; probe = middle of W0; pred_floor = maxᵣ α‖r0ᵣ‖/ε, pred_ceil = minᵣ αε/‖s0ᵣ‖ (‖s0‖ = role *_res0_dual, ε = max(force_tol,floor_tol)); err = log10(pred) − log10(obs). Sweeps: examples/fine2/.
- G1 worst-error, G2 best-error controls, G3 late iterations (last-W0 / last-overall; e08 & e01 had them coincide → one capture each).
- Capture point, A/B/vector conventions, and fidelity identical to snapshots/ (run 1).
