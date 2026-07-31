# snapshots — solve-state capture

git: `92a4855edb`  ·  captured 7/7 snapshots  ·  wall 64.4s

X03 matched pair: ipm it2 (μ=0.1428) vs hsd it2 (μ=0.1403), Δlog10μ = 0.008 (≤ 0.5 ✓)

## Fidelity (test 1) + shapes (test 3)

| snapshot | A (n×n) | nnz A | B (m×n) | nnz B | fid p | fid c | fid w | status |
|---|---|---|---|---|---|---|---|---|
| e08_ipm_it9 | 255×255 | 639 | 128×255 | 510 | 2.8e-13 | 2.3e-12 | — | ok |
| e01_hsd_it8 | 576×576 | 82944 | 240×576 | 69120 | 6.3e-13 | 7.6e-12 | 0.0e+00 | ok |
| e15_ipm_it6 | 4128×4128 | 12400 | 3168×4128 | 50736 | 2.2e-13 | 8.3e-04 | — | ok |
| SOS_P8_n9_s1_hsd_it11 | 440×440 | 12200 | 21×440 | 2310 | 8.9e-08 | 2.8e-07 | 7.7e-08 | ok |
| X03_ipm_it2 | 24×24 | 24 | 12×24 | 57 | 1.7e-15 | 4.6e-14 | — | ok |
| X03_hsd_it2 | 24×24 | 24 | 12×24 | 57 | 1.4e-15 | 2.7e-14 | 1.7e-14 | ok |
| X04_ipm_it3 | 20×20 | 20 | 16×20 | 320 | 7.3e-06 | 2.6e-05 | — | ok |

**Fidelity: max relative Δ = 8.25e-04** (must be ≤ 1e-2).

## Test 2 (loadability)

Check per snapshot dir: `python -c "import scipy.io,pandas,json; scipy.io.mmread('A.mtx'); scipy.io.mmread('B.mtx'); pandas.read_csv('vectors.csv'); json.load(open('meta.json'))"`.

**All 7 load cleanly** (scipy 1.9 mmread for A/B, pandas read_csv for vectors, json.load for meta). Vector columns: IPM dirs `fp,gp,y0p,fc,gc,y0c` (6); HSD dirs add `fw,gw,y0w` (9). Directory sizes 16 KB – 3.2 MB, all far below the 200 MB skip threshold (no snapshot skipped).

## Notes
- A = full symmetric first block (H completed from its stored lower triangle — the operator the solver factors); no rescaling/extra symmetrization. B as the solver holds it.
- vectors.csv columns are NaN/empty-padded to the longest (fp,y0 length n; g length m); slice by n_rows_A / n_rows_B in meta.
- Fidelity reloads A,B from the written .mtx and reconstructs with the captured rhs; r0 targets are the logged pre-CRAIG residual norms.
