# fine1 — fine-grid floor sweep

git: `ff903355`  ·  tol: 1.0e-10  ·  coarse: half-decade (37 pts)  ·  fine: 0.1-decade lattice

| problem | solver | iters | swept | fine½w | bracket-frac | det-maxΔ | consistency-exceptions |
|---|---|---|---|---|---|---|---|
| e01 | ipm | 13 | 13 | ±1.0 | 1.00 | 0.0e+00 | none |
| e15 | ipm | 7 | 4 | ±1.5 | 0.25 | 0.0e+00 | none |
| X08b | ipm | 9 | 5 | ±1.0 | 1.00 | 0.0e+00 | none |
| e01 | hsd | 17 | 17 | ±1.0 | 1.00 | 0.0e+00 | none |
| 06 | ipm | 9 | 9 | ±1.0 | 0.89 | 0.0e+00 | none |

## Skipped fine passes

- e15_1.0e-10_ipm iter 1: floor at coarse-grid minimum (not bracketed)
- e15_1.0e-10_ipm iter 6: floor at coarse-grid minimum (not bracketed)
- e15_1.0e-10_ipm iter 7: floor at coarse-grid minimum (not bracketed)
- X08b_1.0e-10_ipm iter 6: floor at coarse-grid minimum (not bracketed)
- X08b_1.0e-10_ipm iter 7: floor at coarse-grid minimum (not bracketed)
- X08b_1.0e-10_ipm iter 8: floor at coarse-grid minimum (not bracketed)
- X08b_1.0e-10_ipm iter 9: floor at coarse-grid minimum (not bracketed)

## Acceptance tests

1. **Bracketing** — fraction of swept iterations with ≥2 fine rows at each of plateau (nmin), nmin+1, and nmin+2, so both transitions and the full min+1 tread are interior to the range (see table). Problems below 2/3 were auto-widened to ±1.5 decades and re-run.
2. **Determinism** — max relative disagreement between α_floor re-solved on the frozen iterate and its coarse row, over (ncraig, pbase, r0_p). Must be ≤ 1e-12 (see det-maxΔ).
3. **Consistency** — iterations where ncraig rises in the interior of the fine α-lattice (edges allowed). Reported, not fixed — floor-staircase data.

## Structural notes (data, not analysis)

Max in-window predictor pbase (>1 ⇒ the plateau does genuine base-solve work) and per-iteration bracketing shortfalls (which tread level had <2 fine rows):

- e01_1.0e-10_ipm: max in-window pbase = 3; all 3 treads bracketed
- e15_1.0e-10_ipm: max in-window pbase = 2; shortfalls it2:<2@nmin+1, it4:<2@nmin+2, it5:<2@nmin+2
- X08b_1.0e-10_ipm: max in-window pbase = 2; all 3 treads bracketed
- e01_1.0e-10_hsd: max in-window pbase = 3; all 3 treads bracketed
- 06_1.0e-10_ipm: max in-window pbase = 4; shortfalls it3:<2@nmin

wall-clock: 69.5s
