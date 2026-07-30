# T-junctions: subdivision restriction maps, and refinement as a partially filled hole

A companion to `tensor_spline.md`. Two-level adaptive refinement: selected coarse
cells of the patch grid are replaced by their four half-size children, producing
T-junctions — and the corpus's first restriction maps that are genuinely
**asymmetric across an edge**, with classical CAD semantics. Everything else
(Bernstein/thin-plate numerics, the lossy dials, the shapes) is documented in the
tensor note; this records the T-machinery and what it does to H¹.

---

## 1. The T-edge: de Casteljau as a restriction map

At a T-junction a coarse patch abuts two fine children along one line. Each child
gets one edge, saying: *the coarse boundary jet, restricted to the child's
sub-segment and reparametrized, equals the child's boundary jet.*
Restriction-to-half in Bernstein coordinates is the de Casteljau subdivision matrix
(splitting at `t = ½`), so a vertical T-edge reads

```
    kron(S_half, jetₓ(T = h)) · P_coarse  =  kron(I, jetₓ(T = h/2)) · P_fine
```

Two details carry the correctness. **Physical jets**: coarse and fine cells have
different widths, so C^k must be imposed in physical coordinates —
`jet_operator`'s hitherto-decorative `T` argument finally earns its keep. And a
**bookkeeping trap**: coarse–coarse lines must keep *one* full edge; splitting them
into two subdivided half-edges imposes the identical constraint set with
`(k+1)(n+1)` redundant rows per line, silently inflating H¹.

Against the duality note's setup, this is the first builder where the two endpoint
maps `A_ve` of an edge are structurally different objects (subdivision-composed jet
vs plain jet) — the "general restriction map" slot that only the MAP note's Markov
kernels had touched. One worry that does *not* transfer from the T-spline
literature: linear-independence / analysis-suitability. No reduced global basis is
ever constructed here; raw per-patch polynomials are constrained by exact jet
agreement, which is well-posed on any mesh configuration.

The subdivision matrices are verified in the oracle against the reparametrization
identity (`B(t/2)`, `B((1+t)/2)`) and the jet identities that make the geometry
work: `jet(T/2)·S_L = jet(T)` at the shared endpoint (corner jets transport
*exactly* through subdivision — the reason mixed corners behave like plain ones),
and the two halves stitch at the midpoint. All at machine precision.

---

## 2. H¹: the junction rule, and its n ≤ 2k anomaly

Measured first, stated second — this sheaf has now earned that discipline twice.

**The junction rule.** For every tested configuration with `n > 2k`:

```
    dim ker δᵀ = (k+1)² × #{interior fine-mesh vertices with all four
                            quadrants covered and ≥ 3 distinct patches meeting}
```

Fine corners (4 fine patches), T-points (coarse + two children), and mixed coarse
corners (3 coarse + 1 child) all contribute the *same* `(k+1)²` block — because
subdivision preserves endpoint jets exactly, the corner mixed-jet block is
over-determined around a mixed cycle precisely as around a plain one. The rule
subsumes the tensor-grid corner rule (a uniform corner is the 4-distinct-patch
case), and 2-patch mid-edge vertices contribute nothing. The "all quadrants
covered" clause is what excludes the hole-corner vertices of the lossy tensor
configurations (3 patches, but no closing cycle); on the full rectangles here it is
automatic.

**The ring regime.** At `n ≤ 2k`, a refined region enclosed by a ring of coarse
cells picks up extra dimensions — the tensor note's hole term, returned in a new
guise. But refinement *partially fills* the hole: the extras are bounded by, and
generally below, the empty-hole value `(2k−n+1)₊²` per coarse ring not already
spanned by pure-coarse corner cycles:

| configuration | measured extra | empty-hole bound |
|---|---|---|
| 3×3, center refined, n=4 k=2 | **+1** | +1 |
| 3×3, center refined, n=6 k=3 | **+1** | +1 |
| 3×3, center refined, n=3 k=2 (cubic C²!) | **+3** | +4 |
| 4×4, central 2×2 block refined, n=4 k=2 | **+0** | +1 |
| 5×5, two separated cells refined, n=4 k=2 | **+2** | +2 |

The gradation is the story: an *empty* hole (tensor note) contributes the full
`(2k−n+1)₊²`; a hole *filled with four children* contributes 1 of 1 at `n = 2k`
but only 3 of 4 at `(n,k) = (3,2)`; a hole filled with a *finer* 2×2 block of
refined cells contributes nothing. Refinement interpolates between hole and solid,
and H¹ measures how much of the coarse-level transport the refined interior can
absorb. Note that `(3,2)` is the classical cubic-C² family, so the delicate regime
is not exotic. The exact `n ≤ 2k` law is open — a sharp combinatorial question the
oracle adjudicates per instance; `h1_predict_tj` returns the junction rule, exact
for `n > 2k` and a lower bound otherwise.

For the solver, nothing changes qualitatively from the tensor note: the redundancy
is consistent (`g = 0`), the dual is an H₁-coset, and the Uzawa CG pins its
basepoint; T-meshes just add mixed-cycle members to the coset.

---

## 3. Verification

```
subdivision self-test: reparam 1.7e-16  endpoint-jets 0.0  midpoint-stitch 0.0

[2x2 refine (1,1), n=4 k=2, regress nonneg]
  native / assembled optimum : -19.93673799 / -19.93673799   |Δ| = 1.4e-14
  sections ‖ΔP‖∞ / ‖Bp‖      : 2.0e-12 / 1.7e-13
  T-line jump max|∂ₓ^r f|Δ   : 2.1e-14      (r = 0..k, physical, dense sampling)
  H1 measured / predicted    : 36 / 36

junction rule (n > 2k or ring-free)          — 9 configurations, all exact,
  including mixed patterns, k = 1..2, n = 4..5, diagonal and corner refinement
ring regime (n ≤ 2k)                          — table of §2, measured

adaptivity, FIXED data budget (N = 2400):
  coarse   4x4   dof =  400   fit-mse = 1.34e-04     (underfits the bumps)
  adaptive 4x4   dof =  700   fit-mse = 1.12e-04     (refined under the bumps)
  fine     8x8   dof = 1600   fit-mse = 2.12e-04     (overfits: ~37 pts/patch)
```

The T-line jump check is the one an index bug cannot survive: the physical
`∂ₓ^r f`, `r = 0..k`, sampled densely on both sides of a coarse–fine line, agrees
to 1e-14. And the adaptivity table is the reason T-junctions exist, with the honest
statistical reading: at a fixed data budget the adaptive mesh beats *both* uniform
meshes — the coarse one underfits, the fine one overfits, and refinement spends
coefficients only where the truth has structure.

---

## 4. Code

- **`tjunction_spline.jl`** — `subdiv_half`, `generate_tjunction_instance`
  (`refine = [(a,b)]` or a `refined` mask, `refine_by_truth` for the adaptive
  demo), a shared `tj_edges` enumerator consumed by both `build_tjunction_problem`
  and the JuMP reference (so the two models use byte-identical maps),
  `h1_predict_tj`, `tjunction_demo`, and `adaptive_demo`. Scope: `:regress`,
  `:nonneg`/`:free`, full rectangle, uniform `k`, two levels. Monotone leaves, the
  `:intensity` likelihood, and the tensor note's lossy dials compose per-patch and
  per-edge without touching the T-machinery; multi-level refinement is bookkeeping
  plus a 2:1-balance decision, deferred until a use case demands it.
- **`tjunction_spline_oracle.py`** — the ground truth above, including
  `coarse_ring_term` for the ring-regime bound.

Written against the CellularSheaves.IPM PR-67 API; not executed here.

---

## References

- G. Farin. *Curves and Surfaces for CAGD.* 5th ed. (2002). *(de Casteljau
  subdivision, tensor patches)*
- D. Forsey, R. Bartels. *Hierarchical B-spline refinement.* SIGGRAPH (1988).
  *(the origin of two-level local refinement)*
- C. Giannelli, B. Jüttler, H. Speleers. *THB-splines: the truncated basis for
  hierarchical splines.* CAGD 29 (2012). *(hierarchical meshes in analysis; the
  basis-level issues this construction sidesteps)*
- T. Sederberg, J. Zheng, A. Bakenov, A. Nasri. *T-splines and T-NURCCs.*
  SIGGRAPH (2003).
