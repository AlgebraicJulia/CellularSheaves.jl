# Tensor-product surfaces on lossy grids: the sheaf whose H¹ reads the restriction maps

A companion to `nonnegative_spline.md`, one dimension up. Same discipline — cost on the
cells, agreement coupling, cones on the stalks — but the base graph is now a **grid of
patches**, possibly lossy, and two genuinely new things appear: a sheaf H¹ that is
large, structured, and *not* a function of the base topology alone; and the first real
test of the solver's dual-coset path. The Bernstein machinery, orthant sufficiency, and
the 1D exact variant are documented in the spline note; this records the 2D deltas.

Following Papp & Alizadeh (JCGS 2014) and the Alizadeh–Papp multi-dimensional
arrival-rate line (Ann. Oper. Res. 2013), which covers exactly the bivariate
nonnegative-spline ML estimation the `:intensity` mode implements.

---

## 1. The instance family: grids, and every useful way to break them

The domain `[0,1]²` is tiled by `Mx × My` cells; each **active** cell carries a
bidegree-`(n,n)` Bézier patch, stalk `ℝ^{(n+1)²}` (control net, x-index fastest).
Adjacent active patches are joined by continuity edges matching transversal boundary
jets — `kron(I, jetₓ)` across vertical shared lines, `kron(jet_y, I)` across
horizontal ones. Everything is a Kronecker product of the 1D numerics: evaluation
`kron(φ_y, φ_x)`, differences `kron(I, Δ)`, and the thin-plate roughness
`∫ f²ₓₓ + 2f²ₓᵧ + f²ᵧᵧ = kron(G₀,G₂)/h⁴ₓ + 2·kron(G₁,G₁)/h²ₓh²ᵧ + kron(G₂,G₀)/h⁴ᵧ`.

The semantics demand only *per-edge* jet agreement, not a grid — so the builder takes
an active-cell mask and a per-edge smoothness `k_e`, and one code path covers:

| lossiness | `k_e` / mask | meaning |
|---|---|---|
| none | uniform `k` | C^k surface on a rectangle (CAD/spline standard) |
| missing cells | mask | polyomino domains, holes (irregular-domain fitting) |
| missing edges | `k_e = −1` | **faults**: the surface may jump (geological faults, discontinuity-preserving regression) |
| thinned edges | `0 ≤ k_e < k` | **creases**: continuous ridge/feature lines (terrain, CAD) |

Objectives: `:regress` is least-squares + thin-plate (pure `Q` + orthant); `:intensity`
is a spatial-Poisson log-likelihood reusing the `mle_spline.jl` leaf construction
verbatim — one `(u, y, t)` exponential-cone leaf per event, `φ = kron` of Bernstein
evals, compensator `s·∫f` linear in `c` — emitted in the normalized form
`λ = s·f` of `mle_spline.md` §3.1 (the raw form stalled Mosek in the field).

---

## 2. H¹: per-cycle, but not topological

Redundant-but-consistent rows of `B` are exactly `ker δᵀ` — the sheaf's H¹, the
dual-coset dimension of `duality.md` §7. On these instances it decomposes per cycle of
the base graph, with **cycle-dependent weights** (all verified in the oracle):

- **Corner 4-cycle** (4 active cells, 4 live segments): contribution
  `(min(kx_lo, kx_hi) + 1) · (min(ky_lo, ky_hi) + 1)` — the corner mixed-jet block,
  reachable two ways around the cycle. Uniform order: `(k+1)²`.
- **Unit hole** (8-ring, uniform order): contribution `max(2k − n + 1, 0)²`. Transport
  around the ring must pass *through* patches side-to-side, and a bidegree-`n` patch's
  left and right jet functionals overlap in only `(2k−n+1)₊` indices per coordinate.
- **Trees, faulted grids, L-shapes**: zero.

Two consequences. First, in the practical regime `n > 2k` holes are *free* — carve the
domain arbitrarily and only genuine 4-patch corners create redundancy; at the boundary
case `n = 2k` a hole contributes exactly `1`. Second, and this is the point worth
stating as a slogan: **this sheaf's H¹ is not determined by the topology of the base
graph.** Two homologous cycles contribute different amounts depending on the stalk
degree and the restriction maps — the constant-sheaf intuition
`H₁ = cyclerank × dim` fails quantitatively (a full-grid corner gives `(k+1)²` where
the constant sheaf would give the stalk dimension; a hole gives `(2k−n+1)₊²` where the
constant sheaf would again see one full cycle). The corridor's sheaf was tidy but
decorative and the MAP note's cohomology lived in the *positivity*; here, for the
first time, the **linear** cohomology itself genuinely reads the restriction maps.

Sanity anchor: on the full grid the corner rule gives `(Mx−1)(My−1)(k+1)²`, which
matches the dimension count
`rows − [Mx·My·(n+1)² − dim(tensor spline space)]` exactly.

---

## 3. What this asks of the solver

The redundancy is consistent (`g = 0` on continuity edges — exact trivially), so per
§7 of the duality note the primal is unaffected and the dual `y` is a coset mod
`ker δᵀ`. Reading the Uzawa path (`kkt/uzawa.jl`): the factorized operator
`F = A + αBᵀB` stays SPD regardless of `B`'s row rank (A carries `Q ⪰ εI` plus the
cone scaling); the CG right-hand side `r = g − BF⁻¹(f + αBᵀg)` lies in
`range(B) = range(S)` whenever `g` does; and on the range the augmentation clusters
the Schur spectrum near `1/α`. So CG should converge fast and return (approximately)
the **minimum-norm** `y` — which *is* the "solver pins the dual coset through its
basepoint" of §7, now load-bearing. The one failure mode theory can't exclude is
floating-point stagnation of CG on a singular system near the `√eps` tolerance:
that is what the default 2×2 demo (`H¹ = 9`) probes, and the lossy dials give a
bisection tool — faults/holes reduce H¹ continuously to zero, the full grid is the
maximal stress case. Run the probe before anything larger.

---

## 4. Shape constraints in 2D — and the missing exact variant

- **`:nonneg`** — orthant on the stalk, native, still only *sufficient* (convex-hull
  property tensors), still tightened by degree elevation and subdivision.
- **`:monotone_x` / `:monotone_y`** — difference leaves `kron(I,Δ)P − leaf = 0`,
  `leaf ≥ 0`; with `C⁰`-or-better continuity everywhere, per-patch monotonicity glues
  to global monotonicity (not across faults, by design).
- **No exact variant exists.** Karlin–Studden is strictly univariate. Fixed-degree
  certificates for nonnegative *bivariate* polynomials on a box are
  Schmüdgen/Putinar-type **relaxations** whose certificate degree is unbounded in
  general — there is no 2D companion to `nonnegative_spline_exact.jl`, and the note
  says so rather than gesture at one. The orthant is the honest model.
- **Convexity (deferred)** changes character: it is a pointwise PSD-Hessian
  condition. A clean sufficient recipe exists in-idiom — degree-elevate the three
  Hessian entries to a common bidegree and require each 2×2 Bernstein coefficient
  matrix PSD, i.e. one `SecondOrderCone(3)` leaf per coefficient, `(n+1)²` leaves per
  patch — but it is conservative and bulky, so it is documented here and left out of
  the builder.

---

## 5. Verification

`tensor_spline_oracle.py` mirrors `build_tensor_problem` index-for-index (vertex
order: patches, monotone leaves, observation leaves; edge order: shape, observation,
continuity E-then-N) and solves each instance natively and assembled:

```
[probe 2x2 n=4 k=2 regress nonneg]            (the rank-deficiency probe)
  native / assembled optimum : -12.41466998 / -12.41466998   |Δ| = 8.9e-15
  sections ‖ΔP‖∞ / ‖Bp−g‖    : 0.0 / 1.0e-12
  H1 measured / predicted    : 9 / 9
  coset: dim=9   ‖Bᵀν‖∞ = 8.8e-15   |⟨g,ν⟩| = 0.0

[donut 3x3 n=4 k=2]      optima agree to 0.0;  H1 = 1 / 1   (hole term (2k−n+1)² = 1)
[monotone_x 2x2]         |Δ| = 5.7e-07;        H1 = 9 / 9
[intensity L-shape 2x2]  |Δ| = 7.8e-07;  epigraphs tight 8.5e-09;  H1 = 0 / 0
                         normalized λ = 118·f;  scaled − W·log s = raw, |Δ| = 5.0e-06
                         exp-cone args f ∈ [0.34, 5.9]   (raw λ ∈ [40, 696])

H1 formula sweep (rows − rank B):
  full 3x3, n=4 k=2  (4 corners × 9)            36 / 36
  donut 3x3, n=5 k=2 (hole term 0: n > 2k)       0 /  0
  fault 3x3, n=4 k=2 (corner killed: 3 × 9)     27 / 27
  crease k_e=0, 3x3  ((0+1)(2+1) + 3×9)         30 / 30
  crease k_e=1, 3x3  ((1+1)(2+1) + 3×9)         33 / 33
  4x4 two diagonal holes (rings intersect)      18 / 18
```

The per-cycle formula — corner min-rule, hole `(2k−n+1)₊²`, faults zero — holds on
every configuration, including mixed-order creases and interacting holes; the null
vectors of `Bᵀ` pair to zero with `g` (the dual coset is well-defined, `duality.md`
§8's `⟨g,n⟩ = 0` line); and native/assembled optima coincide across all four modes.

---

## 6. Code, and the T-junction extension

- **`tensor_spline.jl`** — the builder: `generate_tensor_instance` (mask + `faults` +
  `creases` kwargs, regression sampling or Poisson thinning on the active region),
  `build_tensor_problem`, `h1_predict`, `tensor_demo` (defaults to the 2×2 probe),
  a JuMP/Mosek reference sharing `patch_QC`, and a benchmark sweep whose large grids
  are the first genuinely 2D test of the multifrontal factorization (real separators,
  frontal blocks `~ grid-width × (n+1)²`, versus the corpus's tridiagonal spines).
- **`tensor_spline_oracle.py`** — the ground truth quoted above.

Deferred, deliberately: **T-junctions** (a coarse patch abutting two fine ones). The
restriction map becomes jet ∘ de Casteljau subdivision matrix — non-square,
non-identity restrictions with classical CAD semantics (THB/LR-spline territory) —
about twenty lines, but it changes the character of the example, so it waits until
the uniform version has survived a real solve. Written against the
CellularSheaves.IPM PR-67 API; not executed here (no Julia toolchain).

---

## References

- F. Alizadeh, D. Papp. *Estimating arrival rate of nonhomogeneous Poisson processes
  with semidefinite programming.* Annals of Operations Research (2013).
  *(multi-dimensional ML intensity over nonnegative splines)*
- G. Farin. *Curves and Surfaces for CAGD.* 5th ed., Morgan Kaufmann (2002).
  *(tensor-product patches, subdivision, creases/feature lines)*
- J. Hansen, R. Ghrist. *Toward a spectral theory of cellular sheaves.* J. Appl.
  Comput. Topology 3 (2019). *(sheaf cohomology on graphs)*
- D. Papp, F. Alizadeh. *Shape-constrained estimation using nonnegative splines.*
  JCGS 23(1) (2014).
- K. Schmüdgen. *The K-moment problem for compact semialgebraic sets.* Math. Ann. 289
  (1991). *(why no fixed-degree exact 2D certificate)*
