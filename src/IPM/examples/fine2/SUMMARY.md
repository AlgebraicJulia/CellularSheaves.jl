# fine2 — two-edge fine-grid floor+ceiling sweep

git: `6ae3d512`  ·  42 runs (21 keys × ipm+hsd, ex-X04)  ·  Σwall 2600s  ·  fine 0.1-decade, adaptive width ±1→3.5

**Determinism (test 2): max Δ = 0.0e+00** across all runs (must be ≤1e-12).  
**Bracketing (test 1): mean floor = 0.86, mean ceil = 0.53.**

| problem | solver | iters | swept | floor-brk | ceil-brk | det-maxΔ | max-w | consistency (edge-aware) | wall |
|---|---|---|---|---|---|---|---|---|---|
| 06 | hsd | 9 | 7 | 0.71 | 0.71 | 0e+00 | ±3.5 | 2:ceil,3:ceil,4:ceil,5:ceil,6:ceil,7:ceil | 54s |
| 06 | ipm | 9 | 9 | 1.00 | 0.44 | 0e+00 | ±3.5 | 5:ceil,6:ceil,8:ceil,9:ceil | 46s |
| 07 | hsd | 17 | 14 | 0.71 | 0.43 | 0e+00 | ±3.5 | 1:ceil,2:floor,3:floor,3:ceil,4:floor,4:ceil,5:ceil,6:ceil,7:ceil,8:ceil,9:ceil,10:ceil,11:ceil,12:ceil,13:ceil,14:ceil | 73s |
| 07 | ipm | 8 | 8 | 1.00 | 0.00 | 0e+00 | ±3.5 | — | 55s |
| 13 | hsd | 17 | 17 | 0.82 | 0.94 | 0e+00 | ±3.5 | 1:ceil,3:floor,3:ceil,4:ceil,5:ceil,6:ceil,7:ceil,9:ceil,10:ceil,11:ceil,12:ceil,13:ceil,14:ceil,16:ceil | 54s |
| 13 | ipm | 19 | 17 | 0.88 | 0.29 | 0e+00 | ±3.5 | 1:ceil,4:ceil,5:ceil,6:ceil,7:ceil,8:ceil,9:ceil,10:ceil,11:ceil,13:ceil,14:ceil | 48s |
| SOS_P4_n13_s1 | hsd | 16 | 12 | 1.00 | 0.75 | 0e+00 | ±3.5 | 2:ceil,5:ceil,6:ceil,7:ceil,8:ceil,9:ceil,10:ceil,11:ceil,12:floor,12:ceil | 60s |
| SOS_P4_n13_s1 | ipm | 12 | 12 | 1.00 | 0.67 | 0e+00 | ±3.5 | 1:ceil,2:ceil,3:ceil,5:ceil,6:ceil,7:ceil,8:ceil,9:ceil,10:ceil,11:ceil,12:ceil | 52s |
| SOS_P8_n9_s1 | hsd | 17 | 14 | 1.00 | 0.86 | 0e+00 | ±3.5 | 7:ceil,9:ceil,10:ceil,11:ceil,12:ceil,13:ceil,14:floor,14:ceil | 59s |
| SOS_P8_n9_s1 | ipm | 12 | 12 | 1.00 | 0.42 | 0e+00 | ±3.5 | 1:ceil,3:ceil,4:ceil,5:ceil,6:ceil,7:ceil,8:ceil,9:ceil,10:ceil,11:ceil,12:ceil | 51s |
| X03 | hsd | 7 | 7 | 0.86 | 1.00 | 0e+00 | ±3.5 | 1:ceil,2:ceil,3:ceil,6:ceil | 33s |
| X03 | ipm | 6 | 6 | 1.00 | 0.50 | 0e+00 | ±3.5 | 1:ceil,2:ceil,4:ceil,5:ceil,6:ceil | 25s |
| X06b | hsd | 9 | 9 | 0.89 | 0.89 | 0e+00 | ±3.5 | 2:ceil,4:ceil,5:ceil,6:ceil | 52s |
| X06b | ipm | 8 | 8 | 1.00 | 0.12 | 0e+00 | ±3.5 | 1:ceil,2:ceil,4:ceil,5:ceil,6:ceil,7:ceil,8:ceil | 45s |
| X08b | hsd | 8 | 8 | 0.62 | 0.25 | 0e+00 | ±3.5 | 1:floor,2:floor,2:ceil,3:ceil,5:ceil,6:ceil,7:ceil,8:ceil | 35s |
| X08b | ipm | 9 | 9 | 1.00 | 0.11 | 0e+00 | ±3.5 | 5:ceil,6:ceil | 26s |
| X09 | hsd | 11 | 11 | 0.73 | 0.45 | 0e+00 | ±3.5 | 1:ceil,2:ceil,4:floor,4:ceil,5:ceil,7:ceil,8:ceil,9:ceil,11:ceil | 53s |
| X09 | ipm | 12 | 12 | 0.75 | 0.08 | 0e+00 | ±3.5 | 1:ceil,3:ceil,6:ceil,7:ceil,8:ceil,9:ceil | 46s |
| e01 | hsd | 17 | 17 | 1.00 | 0.41 | 0e+00 | ±3.5 | 1:ceil,2:ceil,3:ceil,4:ceil,5:ceil,7:ceil,10:ceil,13:ceil,14:ceil,15:ceil,16:ceil,17:ceil | 71s |
| e01 | ipm | 13 | 13 | 1.00 | 0.23 | 0e+00 | ±3.5 | 1:ceil,2:ceil,3:ceil,4:ceil,5:ceil,6:ceil,7:ceil,8:ceil,11:ceil,13:ceil | 59s |
| e02 | hsd | 7 | 6 | 1.00 | 0.67 | 0e+00 | ±3.5 | 1:ceil,2:ceil,6:floor | 57s |
| e02 | ipm | 7 | 7 | 1.00 | 0.29 | 0e+00 | ±3.5 | 1:ceil,2:ceil,3:ceil,4:ceil,5:ceil,6:ceil,7:ceil | 50s |
| e03 | hsd | 11 | 11 | 0.91 | 0.91 | 0e+00 | ±3.5 | 1:ceil,3:ceil,4:ceil,5:ceil,6:ceil,8:ceil,9:ceil,10:ceil,11:floor | 53s |
| e03 | ipm | 10 | 10 | 0.90 | 0.10 | 0e+00 | ±3.5 | 1:ceil,2:floor,2:ceil,4:ceil,5:ceil,6:ceil,7:ceil,8:ceil,9:floor,9:ceil,10:ceil | 47s |
| e04 | hsd | 16 | 13 | 1.00 | 0.85 | 0e+00 | ±3.5 | 1:ceil,3:ceil,4:ceil,5:ceil,6:ceil,7:ceil,8:ceil,9:ceil,10:ceil,11:ceil,12:floor,12:ceil,13:ceil | 64s |
| e04 | ipm | 11 | 11 | 1.00 | 0.55 | 0e+00 | ±3.5 | 1:ceil,2:ceil,6:ceil,9:ceil,10:ceil,11:ceil | 56s |
| e05 | hsd | 11 | 11 | 0.82 | 0.82 | 0e+00 | ±3.5 | 1:floor,2:floor,5:ceil,7:floor,9:ceil,10:ceil,11:floor,11:ceil | 114s |
| e05 | ipm | 10 | 10 | 1.00 | 0.30 | 0e+00 | ±3.5 | 7:ceil,8:ceil,10:ceil | 93s |
| e08 | hsd | 13 | 9 | 0.89 | 0.89 | 0e+00 | ±3.5 | 2:ceil,3:ceil,4:ceil,6:ceil,8:ceil,9:ceil | 54s |
| e08 | ipm | 21 | 21 | 0.86 | 0.57 | 0e+00 | ±3.5 | 3:ceil,4:ceil,5:ceil,8:ceil,10:ceil,11:ceil,13:floor,14:ceil,15:ceil,16:ceil,18:ceil,19:floor,19:ceil,20:ceil,21:ceil | 48s |
| e09 | hsd | 10 | 10 | 0.80 | 0.80 | 0e+00 | ±3.5 | 1:ceil,2:ceil,3:ceil,4:ceil,7:ceil,8:ceil,10:ceil | 68s |
| e09 | ipm | 17 | 17 | 0.82 | 0.12 | 0e+00 | ±3.5 | 1:ceil,2:ceil,4:ceil,5:ceil,10:ceil,11:ceil,13:ceil,14:ceil,15:ceil,16:ceil,17:ceil | 71s |
| e10 | hsd | 11 | 8 | 0.50 | 1.00 | 0e+00 | ±3.5 | 4:ceil | 64s |
| e10 | ipm | 9 | 9 | 0.67 | 0.22 | 0e+00 | ±3.5 | 1:ceil,3:ceil,5:ceil,6:ceil,8:ceil,9:ceil | 62s |
| e11 | hsd | 14 | 11 | 0.82 | 0.91 | 0e+00 | ±3.5 | 1:ceil,2:ceil,3:ceil,4:ceil,8:ceil,11:floor,11:ceil | 88s |
| e11 | ipm | 16 | 15 | 0.87 | 0.60 | 0e+00 | ±3.5 | 1:ceil,4:floor,5:ceil,6:ceil,8:ceil,11:ceil | 76s |
| e12 | hsd | 13 | 12 | 0.58 | 1.00 | 0e+00 | ±3.5 | 3:ceil,4:ceil,7:ceil,8:ceil,10:ceil,12:floor,12:ceil | 136s |
| e12 | ipm | 14 | 14 | 0.71 | 0.36 | 0e+00 | ±3.5 | 1:ceil,2:ceil,3:ceil,4:ceil,5:ceil,8:ceil,10:ceil,11:ceil,13:ceil,14:floor,14:ceil | 147s |
| e14 | hsd | 7 | 6 | 0.50 | 0.83 | 0e+00 | ±3.5 | 1:floor,1:ceil,3:ceil,4:ceil,5:ceil,6:floor,6:ceil | 63s |
| e14 | ipm | 9 | 9 | 0.89 | 0.22 | 0e+00 | ±3.5 | 1:ceil,2:ceil,3:ceil,4:ceil,5:ceil,6:ceil,7:ceil,8:ceil,9:floor,9:ceil | 57s |
| e15 | hsd | 7 | 7 | 1.00 | 0.57 | 0e+00 | ±3.5 | 2:ceil,3:ceil,4:ceil,5:ceil,6:ceil,7:ceil | 69s |
| e15 | ipm | 7 | 7 | 0.43 | 0.00 | 0e+00 | ±3.5 | 3:ceil,4:ceil,5:ceil,6:ceil,7:ceil | 68s |

## Acceptance

1. **Bracketing** — fraction of swept iterations with ≥2 fine rows at each of nmin, nmin+1, nmin+2 on that edge's lattice, after per-edge adaptive widening (±1→1.5→2.5→3.5). Floor and ceiling reported separately.
   - floor < 2/3: X08b_hsd, e10_hsd, e12_hsd, e14_hsd, e15_ipm
   - ceil  < 2/3: 06_ipm, 07_hsd, 07_ipm, 13_ipm, SOS_P8_n9_s1_ipm, X03_ipm, X06b_ipm, X08b_hsd, X08b_ipm, X09_hsd, X09_ipm, e01_hsd, e01_ipm, e02_ipm, e03_ipm, e04_ipm, e05_ipm, e08_ipm, e09_ipm, e10_ipm, e11_ipm, e12_ipm, e14_ipm, e15_hsd, e15_ipm  (ceiling rises can be steep — steep single-jump rises can't show 3 treads; this is structural, reported not fixed)
2. **Determinism** — edge-α re-solved on the frozen iterate vs its coarse row, max over (ncraig, pbase, r0_p). Column det-maxΔ.
3. **Consistency (edge-aware)** — floor lattice flagged for interior ncraig rises; ceiling lattice (rising side) flagged for interior drops. Reported as `iter:edge`, not fixed — staircase data.

## Notes
- Schema = oracle4 columns + `grid` (coarse|fine) + `edge` (coarse|floor|ceil); superset, loads unchanged.
- Both window edges labelled directly (no de-insetting); fit c and c2 on true edges.
- X04 excluded (designated pathology, unsuitable for the label baseline).
