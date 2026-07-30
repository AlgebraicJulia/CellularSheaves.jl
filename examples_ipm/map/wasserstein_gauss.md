# wasserstein_gauss.jl — Bures–Wasserstein transport as a sheaf of covariances

`wasserstein_graph.jl` after the map_lp → map_sdp cone swap: orthant
coupling cells become PSD **joint-covariance blocks**, marginalization
becomes **principal-block compression**, and the graphical gap becomes
**PSD completability** — the correlation-completion story and the
transport story fuse into one file.

## 1. The cells and why they are exact

Coupling per cost-graph edge: $J_e = \begin{pmatrix}\Sigma_i & K\\
K^\top & \Sigma_j\end{pmatrix} \succeq 0$, cost
$\langle C_e, J_e\rangle$ with $C_e = \begin{psmallmatrix}I & \mp I\\
\mp I & I\end{psmallmatrix} = \mathbb E\|x \mp y\|^2$ — linear. No
support grid: cell size is the marginal's *dimension*. For quadratic
pairwise costs the objective reads only second moments, and every
consistent second-moment profile is realized by a Gaussian
(Givens–Shortt; Olkin–Pukelsheim), so:

* 2-marginal: the SDP optimum **is** Bures² — verified against the
  closed form to $3.6\times10^{-7}$;
* barycenter (star, free center PSD cell): agrees with the
  Álvarez-Esteban fixed point to $6\times10^{-5}$ in $\Sigma$
  (Agueh–Carlier; Bhatia–Jain–Lim);
* graphical MMOT: the **exact** problem is one pattern-covariance SDP
  (`build_gw_exact_problem`, a single cell) — poly size at every $K$,
  the only multimarginal setting where the whole ladder can be
  instrumented exactly.

Equivalently this is the **degree-1 moment relaxation**: for
non-Gaussian marginals the same SDP evaluates to **Gelbrich's lower
bound** on the true $W_2^2$ (verified: degree-1 SDP = Gelbrich closed
form to $10^{-8}$), and the moment ladder climbs past it — measured on
bimodal marginals: Gelbrich $0.0216 \le$ degree-2 $0.0978 \le$ true
$W_2^2$ $0.1231$ (oracle §5; the degree-2 cell is a $6\times6$
pseudo-moment matrix with order-4 marginal pins).

## 2. The gap is sign holonomy; the repair is the cover

With `anti` edges (cost $\mathbb E\|x+y\|^2$) on a cycle of pinned
$N(0,1)$ marginals, each coupling is independently optimal at
$K = \pm1$ — locally consistent, pairwise value 0 — while a joint law
needs the sign product around the cycle to be $+1$:

```
  K=3: gap 3.000000 / 0 / 3.000000   (nanti = 1 / 2 / 3)
  K=4: gap 2.343146 / 0             (= 8 − 4√2 at nanti = 1)
```

Same parity law as the discrete-support transport holonomy (odd number
of orientation-reversing edges, any cycle length) — *not* the
binary-MRF odd-cycle-length law. Adding the **cycle cell** (the joint
covariance over the cycle's nodes, with compression agreements to each
member coupling) closes the C4 gap exactly: 0 → 2.343146 = exact.
Refine the cover, restore exactness — GJSW operationalized, and the
concrete rung-2 prototype for the OPF ladder.

## 3. The fusion: completion mode

`tg` targets switch the objective to $\tfrac12\sum_e\|J_e - T_e\|^2$
(Q = I on couplings, transport cost off) — nearest pattern covariance.
On the C4 family $(s,s,s,-s)$: the edge-cover sheaf is **blind**
(distance $\approx 10^{-9}$ for all $|s| \le 1$; edge blocks only share
pinned diagonals), while the exact projection turns on at the
Barrett–Johnson–Tarazaga arccos threshold $s^* = 1/\sqrt2$:
$4\times10^{-8}$ / $6.7\times10^{-4}$ / $3.45\times10^{-2}$ at
$s = 0.70/0.72/0.80$ — the corpus's only obstruction with an analytic
switch-on point, now living inside the transport file it always
belonged to.

## 3½. Proper Q (measured)

Two dense Gram families with founding-literature pedigree:

* **Higham weights** on completion: $\tfrac12\|W^{1/2}(J_e - T_e)
  W^{1/2}\|_F^2$ — the *weighted* nearest correlation matrix (Higham
  2002; Qi–Sun). In svec coordinates the Gram is the **symmetric
  Kronecker** $W \otimes_s W$, built column-wise as
  $\mathrm{svec}(W\,\mathrm{smat}(e_k)\,W)$ (`gw_symkron`; self-test
  $1.8\times10^{-15}$). Measured invariance worth keeping: the weighted
  C4 distances rescale ($6.7\times10^{-4} \to 1.4\times10^{-3}$,
  $3.45\times10^{-2} \to 7.5\times10^{-2}$) but the $1/\sqrt2$
  threshold does not move — W changes the projection geometry, not the
  feasible cone.
* **Shrinkage priors** on free covariance cells (Ledoit–Wolf flavor):
  $\tfrac12\lambda\|\Sigma_k - P\|^2$ via $Q \mathrel{+}= \lambda I$,
  $c = -\lambda\,\mathrm{svec}(P)$. On the barycenter center with
  $P = (\mathrm{tr}\,\Sigma_{fp}/d)I$: $\lambda = 0$ recovers the fixed
  point to 4 decimals, $\lambda = 5$ sits at distance 0.06 from the
  prior — the transport term and the prior trading exactly as the
  smooth-barycenter demo did in the discrete file.

Both stay block-diagonal per vertex; both give the `SemidefiniteCone`
cells dense A-blocks in Uzawa's $F$ — the PSD counterpart of the
Laplacian-Gram probe in `wasserstein_graph.jl`.

## 4. H¹ = 0, and why that is interesting

Measured: 0 for the 2-marginal, the all-pinned triangle, and the
center-free star. Mechanism: the two compression blocks of a coupling
constrain *disjoint* sub-blocks of $J_e$ — unlike discrete couplings,
whose row- and column-marginal blocks share the total-mass functional.
The Kantorovich mass gauge has **no Gaussian analog** (couplings are
not normalized objects), so where `wasserstein_graph.jl` runs Uzawa
rank-deficient on every instance ($H^1 = |E| -$ free balance), this
file runs full-rank everywhere. Two transport sheaves, opposite
deficiency regimes — a clean differential probe for the solver.

## 5. Field notes

* First runs: `gw_bures_demo()` (closed-form check in one solve), then
  `gw_barycenter_demo()`, `gw_gap_demo()` (includes the one-cell exact
  reference and the cycle-cell repair), `gw_completion_demo()`.
* Means separate: $\|m_i - m_j\|^2$ is a tiny standalone QP with closed
  forms; cells here are centered covariances (standard practice).
* Entropic Gaussian OT has closed forms (Janati–Muzellec–Peyré–Cuturi,
  NeurIPS 2020) but needs log-det terms outside the current cone
  language — cited, not built.
* Gromov–Wasserstein is quadratic in $\pi$ but **indefinite** —
  excluded by the normal form, not by taste; its SDP lifts are a
  different project.
* The degree-2 moment rung lives in the oracle only (`moment_bound`);
  building moment cells in Julia is mechanical (Hankel equalities are
  ordinary rows) but deferred until a use case demands non-quadratic
  costs.

## References

Givens & Shortt (1984); Olkin & Pukelsheim (1982); Gelbrich (1990).
Agueh & Carlier, SIAM J. Math. Anal. **43**:904; Álvarez-Esteban,
del Barrio, Cuesta-Albertos, Matrán, J. Math. Anal. Appl. **441**:576
(2016); Bhatia, Jain, Lim, Expo. Math. **37**:165 (2019). Grone,
Johnson, Sá, Wolkowicz, LAA **58**:109 (1984); Barrett, Johnson,
Tarazaga, LAA **192**:3 (1993). Janati, Muzellec, Peyré, Cuturi,
NeurIPS 2020. Alfonsi, Coyaud, Ehrlacher, Lombardi, Math. Comp.
**90**:689 (2021); Khoo, Lin, Lindsey, Ying (moment MMOT).

## Files

* `wasserstein_gauss.jl` — standalone: edge-cover builder with cycle
  cells and completion mode, one-cell exact MMOT, demos, JuMP reference.
* `wasserstein_gauss_oracle.py` — self-contained; all numbers above.
