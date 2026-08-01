# fine2 — two-edge fine-grid floor+ceiling sweep

git: `9273d5b1`  ·  42 runs (21 keys × ipm+hsd, ex-X04)  ·  Σwall 2600s  ·  fine 0.1-decade, per-edge adaptive width ±1→3.5

## Key finding — the two window edges are different in kind

The **floor is a unit staircase** (CRAIG work rises one iteration per ~tread as α drops) — the fine lattice resolves it as nmin / nmin+1 / nmin+2. The **ceiling is a cliff**: at high α the augmentation F=(1/α)H+BᵀB collapses toward the singular BᵀB, so ncraig *jumps* (typically 2→4→6, skipping nmin+1) and the refinement state flips 0→1 (arithmetic breakdown) at the same α. So the two edges need different bracketing criteria, and the ceiling label is *sharper* than the floor (a cliff brackets to <0.1 decade trivially).

**Bracketing** — floor (3 unit treads): mean **0.86**.  ceiling (cliff: ≥2 plateau + ≥2 breakdown rows, transition interior): mean **0.99**  (42/42 ≥ 2/3).
**Determinism** — max Δ = **0.0e+00** across all runs (edge-α re-solved on the frozen iterate vs its coarse row; must be ≤1e-12).

| problem | solver | iters | swept | floor-brk | ceil-brk | det-maxΔ | max-w | consistency (edge-aware) | wall |
|---|---|---|---|---|---|---|---|---|---|
| 06 | hsd | 9 | 7 | 0.71 | 1.00 | 0e+00 | ±3.5 | 2:ceil,3:ceil,4:ceil,5:ceil,6:ceil,7:ceil | 54s |
| 06 | ipm | 9 | 9 | 1.00 | 1.00 | 0e+00 | ±3.5 | 5:ceil,6:ceil,8:ceil,9:ceil | 46s |
| 07 | hsd | 17 | 14 | 0.71 | 1.00 | 0e+00 | ±3.5 | 1:ceil,2:floor,3:floor,3:ceil,4:floor,4:ceil,5:ceil,6:ceil,7:ceil,8:ceil,9:ceil,10:ceil,11:ceil,12:ceil,13:ceil,14:ceil | 73s |
| 07 | ipm | 8 | 8 | 1.00 | 1.00 | 0e+00 | ±3.5 | — | 55s |
| 13 | hsd | 17 | 17 | 0.82 | 1.00 | 0e+00 | ±3.5 | 1:ceil,3:floor,3:ceil,4:ceil,5:ceil,6:ceil,7:ceil,9:ceil,10:ceil,11:ceil,12:ceil,13:ceil,14:ceil,16:ceil | 54s |
| 13 | ipm | 19 | 17 | 0.88 | 1.00 | 0e+00 | ±3.5 | 1:ceil,4:ceil,5:ceil,6:ceil,7:ceil,8:ceil,9:ceil,10:ceil,11:ceil,13:ceil,14:ceil | 48s |
| SOS_P4_n13_s1 | hsd | 16 | 12 | 1.00 | 1.00 | 0e+00 | ±3.5 | 2:ceil,5:ceil,6:ceil,7:ceil,8:ceil,9:ceil,10:ceil,11:ceil,12:floor,12:ceil | 60s |
| SOS_P4_n13_s1 | ipm | 12 | 12 | 1.00 | 1.00 | 0e+00 | ±3.5 | 1:ceil,2:ceil,3:ceil,5:ceil,6:ceil,7:ceil,8:ceil,9:ceil,10:ceil,11:ceil,12:ceil | 52s |
| SOS_P8_n9_s1 | hsd | 17 | 14 | 1.00 | 1.00 | 0e+00 | ±3.5 | 7:ceil,9:ceil,10:ceil,11:ceil,12:ceil,13:ceil,14:floor,14:ceil | 59s |
| SOS_P8_n9_s1 | ipm | 12 | 12 | 1.00 | 1.00 | 0e+00 | ±3.5 | 1:ceil,3:ceil,4:ceil,5:ceil,6:ceil,7:ceil,8:ceil,9:ceil,10:ceil,11:ceil,12:ceil | 51s |
| X03 | hsd | 7 | 7 | 0.86 | 1.00 | 0e+00 | ±3.5 | 1:ceil,2:ceil,3:ceil,6:ceil | 33s |
| X03 | ipm | 6 | 6 | 1.00 | 1.00 | 0e+00 | ±3.5 | 1:ceil,2:ceil,4:ceil,5:ceil,6:ceil | 25s |
| X06b | hsd | 9 | 9 | 0.89 | 1.00 | 0e+00 | ±3.5 | 2:ceil,4:ceil,5:ceil,6:ceil | 52s |
| X06b | ipm | 8 | 8 | 1.00 | 1.00 | 0e+00 | ±3.5 | 1:ceil,2:ceil,4:ceil,5:ceil,6:ceil,7:ceil,8:ceil | 45s |
| X08b | hsd | 8 | 8 | 0.62 | 1.00 | 0e+00 | ±3.5 | 1:floor,2:floor,2:ceil,3:ceil,5:ceil,6:ceil,7:ceil,8:ceil | 35s |
| X08b | ipm | 9 | 9 | 1.00 | 1.00 | 0e+00 | ±3.5 | 5:ceil,6:ceil | 26s |
| X09 | hsd | 11 | 11 | 0.73 | 1.00 | 0e+00 | ±3.5 | 1:ceil,2:ceil,4:floor,4:ceil,5:ceil,7:ceil,8:ceil,9:ceil,11:ceil | 53s |
| X09 | ipm | 12 | 12 | 0.75 | 1.00 | 0e+00 | ±3.5 | 1:ceil,3:ceil,6:ceil,7:ceil,8:ceil,9:ceil | 46s |
| e01 | hsd | 17 | 17 | 1.00 | 0.94 | 0e+00 | ±3.5 | 1:ceil,2:ceil,3:ceil,4:ceil,5:ceil,7:ceil,10:ceil,13:ceil,14:ceil,15:ceil,16:ceil,17:ceil | 71s |
| e01 | ipm | 13 | 13 | 1.00 | 1.00 | 0e+00 | ±3.5 | 1:ceil,2:ceil,3:ceil,4:ceil,5:ceil,6:ceil,7:ceil,8:ceil,11:ceil,13:ceil | 59s |
| e02 | hsd | 7 | 6 | 1.00 | 1.00 | 0e+00 | ±3.5 | 1:ceil,2:ceil,6:floor | 57s |
| e02 | ipm | 7 | 7 | 1.00 | 1.00 | 0e+00 | ±3.5 | 1:ceil,2:ceil,3:ceil,4:ceil,5:ceil,6:ceil,7:ceil | 50s |
| e03 | hsd | 11 | 11 | 0.91 | 1.00 | 0e+00 | ±3.5 | 1:ceil,3:ceil,4:ceil,5:ceil,6:ceil,8:ceil,9:ceil,10:ceil,11:floor | 53s |
| e03 | ipm | 10 | 10 | 0.90 | 1.00 | 0e+00 | ±3.5 | 1:ceil,2:floor,2:ceil,4:ceil,5:ceil,6:ceil,7:ceil,8:ceil,9:floor,9:ceil,10:ceil | 47s |
| e04 | hsd | 16 | 13 | 1.00 | 0.92 | 0e+00 | ±3.5 | 1:ceil,3:ceil,4:ceil,5:ceil,6:ceil,7:ceil,8:ceil,9:ceil,10:ceil,11:ceil,12:floor,12:ceil,13:ceil | 64s |
| e04 | ipm | 11 | 11 | 1.00 | 1.00 | 0e+00 | ±3.5 | 1:ceil,2:ceil,6:ceil,9:ceil,10:ceil,11:ceil | 56s |
| e05 | hsd | 11 | 11 | 0.82 | 1.00 | 0e+00 | ±3.5 | 1:floor,2:floor,5:ceil,7:floor,9:ceil,10:ceil,11:floor,11:ceil | 114s |
| e05 | ipm | 10 | 10 | 1.00 | 1.00 | 0e+00 | ±3.5 | 7:ceil,8:ceil,10:ceil | 93s |
| e08 | hsd | 13 | 9 | 0.89 | 0.89 | 0e+00 | ±3.5 | 2:ceil,3:ceil,4:ceil,6:ceil,8:ceil,9:ceil | 54s |
| e08 | ipm | 21 | 21 | 0.86 | 1.00 | 0e+00 | ±3.5 | 3:ceil,4:ceil,5:ceil,8:ceil,10:ceil,11:ceil,13:floor,14:ceil,15:ceil,16:ceil,18:ceil,19:floor,19:ceil,20:ceil,21:ceil | 48s |
| e09 | hsd | 10 | 10 | 0.80 | 1.00 | 0e+00 | ±3.5 | 1:ceil,2:ceil,3:ceil,4:ceil,7:ceil,8:ceil,10:ceil | 68s |
| e09 | ipm | 17 | 17 | 0.82 | 1.00 | 0e+00 | ±3.5 | 1:ceil,2:ceil,4:ceil,5:ceil,10:ceil,11:ceil,13:ceil,14:ceil,15:ceil,16:ceil,17:ceil | 71s |
| e10 | hsd | 11 | 8 | 0.50 | 1.00 | 0e+00 | ±3.5 | 4:ceil | 64s |
| e10 | ipm | 9 | 9 | 0.67 | 1.00 | 0e+00 | ±3.5 | 1:ceil,3:ceil,5:ceil,6:ceil,8:ceil,9:ceil | 62s |
| e11 | hsd | 14 | 11 | 0.82 | 1.00 | 0e+00 | ±3.5 | 1:ceil,2:ceil,3:ceil,4:ceil,8:ceil,11:floor,11:ceil | 88s |
| e11 | ipm | 16 | 15 | 0.87 | 1.00 | 0e+00 | ±3.5 | 1:ceil,4:floor,5:ceil,6:ceil,8:ceil,11:ceil | 76s |
| e12 | hsd | 13 | 12 | 0.58 | 1.00 | 0e+00 | ±3.5 | 3:ceil,4:ceil,7:ceil,8:ceil,10:ceil,12:floor,12:ceil | 136s |
| e12 | ipm | 14 | 14 | 0.71 | 1.00 | 0e+00 | ±3.5 | 1:ceil,2:ceil,3:ceil,4:ceil,5:ceil,8:ceil,10:ceil,11:ceil,13:ceil,14:floor,14:ceil | 147s |
| e14 | hsd | 7 | 6 | 0.50 | 1.00 | 0e+00 | ±3.5 | 1:floor,1:ceil,3:ceil,4:ceil,5:ceil,6:floor,6:ceil | 63s |
| e14 | ipm | 9 | 9 | 0.89 | 1.00 | 0e+00 | ±3.5 | 1:ceil,2:ceil,3:ceil,4:ceil,5:ceil,6:ceil,7:ceil,8:ceil,9:floor,9:ceil | 57s |
| e15 | hsd | 7 | 7 | 1.00 | 1.00 | 0e+00 | ±3.5 | 2:ceil,3:ceil,4:ceil,5:ceil,6:ceil,7:ceil | 69s |
| e15 | ipm | 7 | 7 | 0.43 | 1.00 | 0e+00 | ±3.5 | 3:ceil,4:ceil,5:ceil,6:ceil,7:ceil | 68s |

## Acceptance

1. **Bracketing** — floor: ≥2 fine rows at each of nmin, nmin+1, nmin+2. ceiling (cliff): ≥2 plateau rows and ≥2 rows past the cliff (ncraig≥nmin+2 or state=1), the transition interior. Per-edge adaptive widening ±1→1.5→2.5→3.5.
   - floor < 2/3: X08b_hsd(0.62), e10_hsd(0.50), e12_hsd(0.58), e14_hsd(0.50), e15_ipm(0.43)
   - ceil  < 2/3: none
2. **Determinism** — column det-maxΔ; 0 everywhere ⇒ the iterate is held exactly fixed across each sweep.
3. **Consistency (edge-aware)** — floor lattice flagged for interior ncraig rises; ceiling lattice flagged for interior drops. `iter:edge`, reported not fixed — the ceiling's near-cliff jumps and the floor's staircase noise both surface here.

## Notes
- Schema = oracle4 columns + `grid` (coarse|fine) + `edge` (coarse|floor|ceil); superset, loads unchanged.
- Both edges labelled directly at 0.1-decade (no de-insetting); fit c and c2 on true edges. The ceiling label = the state 0→1 / ncraig-jump boundary, not an ncraig tread.
- X04 excluded (designated pathology).
- `max-w = ±3.5` on every run is driven by the ceiling: the cliff has no nmin+1 tread, so the *floor's* 3-tread widening loop is separate — floor typically brackets by ±1.0; the ±3.5 reflects the ceiling loop exhausting its schedule (harmless — the cliff was already bracketed).
