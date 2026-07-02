# opf_sheaf.jl — AC optimal power flow as a ladder of covers

The industrial face of the corpus's obstruction story. The power
network G is the *problem graph*; every constraint and cost of AC-OPF
is linear in the entries of W = VVᴴ on G's pattern (rank W = 1 is the
dropped nonconvexity); and each relaxation in the standard OPF ladder
is a **choice of cover** of that pattern. The sheaf base graph is the
cover's nerve — one problem graph, three sheaves.

## 1. The three rungs

| rung | cover | cells | tight when | literature |
|---|---|---|---|---|
| 1 SOCP (Jabr) | lines | 4-dim Lorentz per line | radial networks | Jabr 2006; Sojoudi–Lavaei; Low's tutorial |
| 2 +cycles | lines + cycles | + real-embedded Hermitian W_C | (measured: closes all our meshes) | Kocuk–Dey–Sun, Oper. Res. 2016 |
| 3 chordal/full | cliques of a triangulation | PSD clique blocks | always (= full SDP, GJSW) | Lavaei–Low; Molzahn–Hiskens |

**Rung 1 cells and maps.** Per line e = (i,j), i<j: Jabr coordinates
u_i = |V_i|², c_e + j·s_e = V_i·conj(V_j); the 2×2 PSD condition
c² + s² ≤ u_iu_j becomes the Lorentz membership of
z_e = (u_i+u_j, u_i−u_j, 2c_e, 2s_e) — bound first, the library's
convention, and the corpus's first *load-bearing* `SecondOrderCone`
cells. Restriction maps with classical semantics: magnitude extraction
[½, ±½, 0, 0]·z_e = u_i ("|V_i|² as this line sees it") into shared bus
cells; 2-row injection hyperedges per bus assembling
S_i = Σ_k conj(Y_ik)W_ik in real coordinates (verified against complex
injections to 1.4×10⁻¹⁴); loads ride in g; generator cells (p, q) carry
the **native quadratic** ½c₂p² + c₁p — the flagship native-Q claim,
finally in industrial costume. Boxes are orthant slack leaves. The
base graph on line/bus cells is the *subdivision graph* of G plus
injection hyperedges and pendants.

**Rung 2 cells.** One Hermitian W_C per cycle, real-embedded as
M = [[X, −Z],[Z, X]] ⪰ 0 — the corpus's first complex embedding
(verified: same spectrum sign; assembled svec rows reproduce the
complex model's value *exactly*, 1.334255). Structure rows (A = D
block, Z antisymmetric) are degree-1 edges on the cell; agreement rows
tie X_ii to bus cells and X_ij, Z_ij to member lines' c, s. This is
Kocuk–Dey–Sun's cycle strengthening, sheaf-native: *enlarge the cover
along the cycle* — the exact move `wasserstein_gauss.jl`'s C4 repair
prototyped.

## 2. The dropped constraint is angular holonomy

Any rank-1 W satisfies Σ_C ±θ_e = 0 around every cycle, θ_e =
atan2(s_e, c_e) — voltage angles must be globally consistent. The SOCP
rung keeps each line's 2×2 physics and drops exactly this. Measured on
the seed-0 triangle (loads Pd = (0, .9, .4), Qd = (0, .35, −.25),
generator at bus 1, v ∈ [0.9, 1.1]):

```
  rung 1 SOCP     : 1.333280   holonomy defect −0.034636
  rung 2 +cycle   : 1.334255   holonomy defect −0.000000
  rung 3 full SDP : 1.334255   holonomy defect −0.000000
```

The cycle cell closes the gap exactly, and the defect vanishes with it.
On the C4 seed sweep the defect **tracks the gap seed by seed**
(defect 0.078 ↔ gap 2.3×10⁻³; defect 0.004 ↔ gap 1.3×10⁻⁵) — the
violated holonomy is the gap's dial, not just its symptom. Radial
star: gap −5×10⁻⁹ (tree-exactness). Quadratic cost (c₂ = 2): gap
3.6×10⁻³, sheaf assembly ≡ direct model to 2×10⁻⁸ in both cost modes.

This is the corpus's third holonomy — composed Monge maps
(`wasserstein_graph`), channel cycles (`quantum_marginal`), voltage
angles here — and the only one with an industrial literature measuring
it (Kocuk–Dey–Sun's cycle cuts reach 99.96% of the SDP value on IEEE
cases; our tiny meshes happen to close entirely).

## 3. Solver notes

The fullest cone mix in the corpus under one coboundary:
`SecondOrderCone` (lines) + `CofreeCone` (buses, generators) +
`PositiveCone` (box slacks) + `SemidefiniteCone` (rungs ≥ 2), with
native Q on the generator cell. H¹ of the assembled rung-1 triangle:
**0** (22×27, full rank) — boxes and injections leave no gauge, so OPF
joins `wasserstein_gauss` in the full-rank regime. Climbing a rung is
one kwarg (`cycles = [[1,2,3]]`): Q, c, and every rung-1 cell survive
untouched — the cover is the API.

First Julia runs: `opf_triangle_demo()` (three solves, the whole story:
values and defects above), `opf_radial_demo()`, `opf_c4_demo(seed = 1)`
(quadratic cost, biggest measured defect). Warning for the first run:
the SOC cells put the bound coordinate FIRST — if the triangle demo's
rung-1 value differs from 1.333280 by more than ~1e-5, suspect the
Lorentz ordering before anything else.

## 4. What is deliberately out of scope

True AC-OPF (rank-1, nonconvex) — this file is the relaxation ladder;
recovering feasible dispatches from near-rank-1 W is its own
literature. Shunts, transformers, thermal line limits — all linear in
W, all addable as rows without new machinery. Multi-generator cost
curves — the gen cell generalizes mechanically. Larger IEEE cases —
the builders are size-generic; instance parsers are not corpus
business.

## References

Jabr, IEEE TPS **21**:1000 (2006). Lavaei & Low, IEEE TPS **27**:92
(2012). Low, *Convex relaxation of OPF* (tutorial, 2014). Sojoudi &
Lavaei (SOCP = SDP on radial). Kocuk, Dey, Sun, Oper. Res. **64**:1177
(2016). Molzahn & Hiskens, FnT Electric Energy Systems (2019). GJSW
LAA **58**:109 for the chordal rung.

## Files

* `opf_sheaf.jl` — standalone: rung-1 builder with `cycles` kwarg for
  rungs 2–3, holonomy diagnostics, demos, JuMP reference.
* `opf_sheaf_oracle.py` — self-contained; all numbers above, including
  the real-embedding row verification.
