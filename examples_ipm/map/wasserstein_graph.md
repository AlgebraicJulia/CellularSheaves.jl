# wasserstein_graph.jl — graph-structured optimal transport as a sheaf

The pairwise (local-consistency) relaxation of multimarginal optimal
transport with a graph-structured cost — and the liberation of the
Wasserstein barycenter from its star. One assembly, three readings:

* **all nodes pinned** — the pairwise relaxation of graphical MMOT
  (Haasler–Singh–Zhang–Karlsson–Chen, IEEE Trans. Inf. Theory
  **67**:4647, 2021);
* **star, center free** — the Wasserstein barycenter (Agueh–Carlier),
  whose tractability becomes a *corollary*: star ⊂ tree ⊂ tight;
* **general graph, some pins** — Wasserstein propagation (Solomon et
  al., ICML 2014), where the couplings are the honest variables and
  the model is exact on every topology.

This is `map_lp.jl` wearing transport clothes, and the file leans into
the correspondence rather than hiding it: the minimal gap witness below
*is* the frustrated triangle, reread.

## 1. The problem and the sheaf

Cost graph $G = (V, E)$, node supports $x_k$ of size $n_k$, pairwise
costs $C_e$, marginals $\nu_k$ pinned on a subset. The pairwise model:

$$\min_{\pi, \mu}\ \sum_{e=(i,j)} \langle C_e, \pi_e\rangle
\quad\text{s.t.}\quad
\mathrm{rowmarg}(\pi_e) = \mu_i,\ \ \mathrm{colmarg}(\pi_e) = \mu_j,\ \
\mu_k = \nu_k\ (k\ \text{pinned}),\ \ \pi, \mu \ge 0.$$

| sheaf datum | transport object |
|---|---|
| coupling cell $\pi_e$, one per cost-graph edge | `PositiveCone`, dim $n_i n_j$, column-major vec; carries $c = \mathrm{vec}(C_e)$ |
| marginal cell $\mu_k$, **free nodes only** | `PositiveCone`; pinned marginals ride in $g$ (the corridor's "$g = \delta$ of a reference" move) |
| two base-edges per coupling | $\mathrm{rowmarg} = \mathrm{kron}(\mathbf 1^\top, I)$, $\mathrm{colmarg} = \mathrm{kron}(I, \mathbf 1^\top)$ |
| normalization | implied: a pin's mass propagates through every coupling; ≥1 pin per component asserted |
| exact MMOT reference | the **one-cell degenerate sheaf**: the joint tensor as a single `PositiveCone` vertex with $K$ marginalization pin edges — the transport analogue of `map_sdp`'s single-clique Shor block |

## 2. Tightness and the gap: the theorems

**Trees are exact.** Haasler–Ringh–Chen–Karlsson (SIAM J. Control
Optim. **59**:2428, 2021, §4): the unregularized MMOT with cost
decoupling over a tree equals the sum-of-bimarginal-costs problem —
i.e. this sheaf. Measured: path $K{=}4$, random marginals, gap
$1.5\times10^{-7}$ (LP tolerance). The same paper exhibits the cycle
counterexample (pairwise transports that fail path-consistency around a
cycle), which is our gap in the primary source.

**Cycles must gap.** Altschuler–Boix-Adserà (Discrete Optim. 2021):
MMOT with pairwise-decomposable cost is NP-hard in general —
tractability requires low treewidth (Fan–Haasler–Karlsson–Chen,
AISTATS 2022). So the pairwise relaxation is *necessarily* loose on
cyclic cost graphs, by the same $P \ne NP$ logic that makes MAP-LP
loose. The obstruction is again a failure of positive local sections to
glue: locally consistent couplings that are the two-dimensional
marginals of **no** joint law.

**The gap is holonomy** (measured). Uniform marginals on a $K$-cycle,
quadratic costs, with `nanti` edges wanting $y = 1{-}x$ and the rest
$y = x$. Each coupling can independently be its optimal Monge map (all
marginals stay uniform — locally consistent, pairwise value 0), but a
joint law exists iff the *composed* Monge map around the cycle is the
identity, i.e. iff `nanti` is even:

```
  binary agreement triangle : pairwise 0.000000  exact 1.000000   (map_lp: 0 vs 1, verbatim)
  K=3 n=8  nanti=1 / 2      : gap 0.326531 / 0.000000
  K=4 n=8  nanti=1 / 2      : gap 0.255102 / 0.000000
  K=5 n=6  nanti=1 / 2      : gap 0.240000 / 0.000000
```

Cycle *length* parity is irrelevant; the parity that matters is the
number of orientation-reversing edges — the composed-map holonomy, the
transport sibling of the synchronization-holonomy remark in
`local_consistency_sheaf.md` §2.

## 3. Barycenter mode

Star with pinned leaves, free center, edge costs $\lambda_k \cdot
|x - y|^2$. Verified: McCann's midpoint (barycenter of $\delta_{.25}$,
$\delta_{.75}$ at $\lambda = (\tfrac12,\tfrac12)$ on a grid containing
$0.5$) puts mass 1.000000 at the midpoint; a 3-leaf mixture agrees with
POT's `ot.lp.barycenter` in **objective** to $1.3\times10^{-8}$ while
the centers differ by $1.4\times10^{-2}$ — the barycenter LP is
degenerate, and the $\varepsilon I$ block picks the min-norm optimum
(a feature: SparseMAP's projection reading, unique and differentiable).

## 4. Entropic dial: a second gap, on trees

`build_wg_entropic_problem` adds $(1/\eta)\sum \pi\log\pi$ per coupling
via one `ExponentialCone` leaf per coupling coordinate: leaf $(x,y,z)$
with $x = 1$ pinned, $y$ wired to $\pi_{ab}$, objective $-1/\eta$ on
$z$; at the optimum $z = \pi\log(1/\pi)$ (epigraph tight) — the
`mle_spline.jl` pattern with the roles of $x$ and $y$ swapped. This
gives Sinkhorn-comparable problems, and exposes the *regularization*
gap orthogonal to the topological one: on a tree, where the
unregularized problems coincide, the entropic pairwise and entropic
multimarginal problems differ (HRCK's smoothing observation). Measured
on a pinned-ends path ($K{=}3$, $n{=}15$):

```
  η=10: middle-marginal entropy  pairwise 2.5083   multimarginal 2.2390
  η=40:                          pairwise 1.9808   multimarginal 1.6471
```

— pairwise regularization is strictly more diffusive, at every $\eta$,
converging to the common LP value as $\eta \to \infty$. Two gaps, two
independent knobs: topology (η-independent) and regularization
(topology-independent).

## 5. H¹: one Kantorovich gauge per coupling

Measured first: path all-pinned 3, triangle all-pinned 3, star
center-free 2, 4-cycle one-node-free 3. So the deficiency is **per
coupling, not per cycle**: within one coupling, the row- and
column-marginal blocks both imply the same total mass — the classical
transportation-polytope redundancy, dual to the Kantorovich potential
gauge $(f, g) \mapsto (f + c,\, g - c)$. Free nodes glue the gauges:
a free $\mu_k$ column forces the incident gauge constants to balance,
removing one dimension per independent balance equation:

$$H^1 = |E| \;-\; \mathrm{rank}\big(\text{signed incidence restricted
to free nodes}\big),$$

implemented as `wg_h1_predict` and matching all four measurements. Note
what changed from `map_lp`: there, per-node normalization pins ate this
redundancy; here normalization is implied and the gauge survives — so
*every* instance of this file, even on trees, runs Uzawa on a
consistently rank-deficient $B$. The file is therefore also a
volume probe of the rank-deficiency tolerance, at LP cost.

## 5½. Proper block-diagonal Q (measured)

With Q = εI the QP is a projection and the stalks are quadratically
decoupled — all structure lives in B. This file also carries *genuine*
per-cell Grams, each with a literature home, so the vertex blocks of the
KKT system stop being scalar:

* **First-class η** — quadratically regularized OT (Blondel–Seguy–Rolet
  AISTATS 2018; Lorenz–Manns–Meyer). The solution is the Euclidean
  projection of −C/η onto the polytope: unique, differentiable, and
  **sparse** where entropic is dense. Measured (n = 25 single coupling):
  LP support 49 (= 2n−1), quad η = 0.05: 53, η = 0.5: 83 — against
  entropic η = 20: **519 of 625**.
* **Laplacian-regularized OT** (Courty–Flamary–Tuia–Rakotomamonjy,
  TPAMI 2017) — "transported neighbors stay neighbors". The Kronecker
  Gram $Q_e = 2\eta[\alpha\,(X_tX_t^\top \otimes L_s) + (1{-}\alpha)
  (L_t \otimes X_sX_s^\top)]$ matches POT's implemented `emd_laplace`
  objective exactly (Gram self-test $1.7\times10^{-13}$, Q ⪰ 0 to
  $8\times10^{-14}$), and our QP optimum beats POT's conditional-gradient
  value by $2.6\times10^{-5}$ on the same instance — the IPM as the
  *exact* solver for a problem POT solves approximately.
* **Evidence Grams** (inverse-OT / matching econometrics: Dupuy–Galichon;
  Stuart–Wolfram): observed plan reads $y = \Phi\,\mathrm{vec}(\hat\pi)
  + \text{noise}$ enter as $Q_e \mathrel{+}= \lambda\Phi^\top\Phi$,
  $c \mathrel{-}= \lambda\Phi^\top y$ — the spline files' data-Gram
  pattern on couplings. Measured: $\|\pi - \hat\pi\|_1$ drops
  1.95 → 1.47 when the evidence is switched on.
* **H¹ smoothness prior on free marginals** (`wg_smooth_barycenter`):
  center roughness $\Sigma(\Delta\mu)^2$ falls 0.250 → 0.100 → 0.013 as
  $\lambda_{\text{smooth}}$ sweeps 0 → 0.02 → 0.2, support spreading
  4 → 8 → 11 of 21 — the transport cost and the prior visibly trading.

All four touch one cell's coordinates each, so Q stays block-diagonal
per vertex — the normal form intact, the blocks finally dense. (Side
effect worth watching in Julia: Uzawa's $F = A + \alpha B^\top B$ now
has dense orthant-cell A-blocks, a path nothing else in the corpus
exercises.)

## 6. Verification summary (wasserstein_graph_oracle.py, CLARABEL)

Marginalization maps $2.2\times10^{-16}$; sheaf assembly vs native
matrix-variable model $|\Delta| = 0$; tree gap $1.5\times10^{-7}$;
McCann midpoint exact; POT objective agreement $1.3\times10^{-8}$;
the holonomy gap table, entropy table, and structured-Q measurements above; H¹ 3/3/2/3 against the
gauge law. Exact MMOT references are brute-force joint LPs
($n^K \le 3.3\times10^4$ variables).

## 7. Field notes

* First Julia runs: `wg_gap_demo()` (three solves incl. the one-cell
  exact references; the binary triangle should print gap 1.000000),
  then `wg_barycenter_demo()`, then `wg_entropic_demo()` (first live
  mix of `PositiveCone` + `ExponentialCone` leaves outside the spline
  family; 450 leaves at the demo size).
* The one-cell exact builder is deliberately naive ($\prod n_k$
  variables) — it is a *reference*, not a method; keep $K \le 4$,
  $n \le 10$. The junction-tree route to exact MMOT at scale is the
  HRCK algorithm, out of scope.
* Propagation mode needs no new code: any graph, any pin subset. On
  cyclic graphs its couplings are honest variables, so there is no
  relaxation gap to worry about — only the (intended) modeling choice
  that consistency is pairwise.
* The $\varepsilon$-shift correction for reported optima mirrors the
  quantum file: exact $\varepsilon$-optimum $\le$ LP optimum
  $+\ \varepsilon\cdot(\text{cells})/2$; demos print the linear part.

## References

Haasler, Singh, Zhang, Karlsson, Chen, IEEE Trans. Inf. Theory
**67**:4647 (2021). Haasler, Ringh, Chen, Karlsson, SIAM J. Control
Optim. **59**:2428 (2021). Altschuler & Boix-Adserà, Discrete Optim.
**42**:100669 (2021). Fan, Haasler, Karlsson, Chen, AISTATS 2022.
Agueh & Carlier, SIAM J. Math. Anal. **43**:904 (2011). Solomon,
Rustamov, Guibas, Butscher, ICML 2014. Niculae & Martins (SparseMAP)
for the $\varepsilon$/projection reading.

## Files

* `wasserstein_graph.jl` — standalone builder: pairwise, one-cell exact
  MMOT, entropic; demos; JuMP reference.
* `wasserstein_graph_oracle.py` — self-contained (POT optional, used
  only as a cross-check); all numbers above.
