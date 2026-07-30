# tjunction_sdp.jl — SDP shape certificates on the T-junction mesh

The 2D sibling of `nonnegative_spline_exact.jl`, and the corpus's first
mixed **PSD + subdivision** problem. Two shapes, selected by
`inst.shape`:

* `:nonneg_sdp` — nonnegativity of the fitted surface, by a
  degree-truncated Schmüdgen certificate on the box;
* `:convex_sdp` — convexity, by an SOS-matrix certificate for the
  Hessian. This one has no LP-sufficient alternative worth having: PSD-ness
  of a polynomial matrix is natively semidefinite.

Everything T-junction — subdivision restriction maps, physical jets,
coarse-line bookkeeping, adaptive refinement — is inherited from
`tjunction_spline.jl` untouched. Only the per-vertex data changes.

## 1. The certificates, and what they are not

**Honesty first.** In one variable, nonnegativity on an interval *is* an
SDP (Lukács/Karlin–Studden) — that is why the 1D file is called `_exact`.
In two variables it is not: nonnegative ≠ sum of squares (Hilbert 1888;
Motzkin's $u^4v^2+u^2v^4-3u^2v^2+1$), and the Schmüdgen/Putinar
representations that do exist on the box require unbounded multiplier
degree. What we implement is the *degree-truncated* Schmüdgen
preordering: sufficient, never necessary. The compensation is that it is
strictly tighter than the coefficient-nonnegativity LP already in the
corpus — and the demo measures the gap.

**Nonnegativity** ($n = 2m$ even). Per patch, in local coordinates,

$$f = \sigma_0 + u(1{-}u)\,\sigma_1 + v(1{-}v)\,\sigma_2 +
      u(1{-}u)v(1{-}v)\,\sigma_3,\qquad
  \sigma_a = z_a^\top G_a z_a,\ G_a \succeq 0,$$

with $z_a$ the tensor Bernstein bases of bidegree $(m,m)$, $(m{-}1,m)$,
$(m,m{-}1)$, $(m{-}1,m{-}1)$ — so each summand lands exactly in bidegree
$(n,n)$. The Gram→coefficient maps $\Lambda_a$ come from two univariate
identities,

$$b_i^p b_{i'}^p = \tfrac{\binom{p}{i}\binom{p}{i'}}{\binom{2p}{i+i'}}
  b_{i+i'}^{2p}, \qquad
  u(1{-}u)\, b_s^{2p} = \tfrac{\binom{2p}{s}}{\binom{2p+2}{s+1}}
  b_{s+1}^{2p+2},$$

which compose to the single weight
$w = \binom{p}{i}\binom{p}{i'}/\binom{n}{i+i'+e_u}\cdot(\text{same in }v)$
targeting coefficient $(i{+}i'{+}e_u,\ j{+}j'{+}e_v)$.

**LP ⊆ SDP, strictly.** For even $n$, split $u^I(1-u)^{n-I}$ by parity of
$I$: $I$ even gives a perfect square, $I$ odd gives $u(1{-}u)$ times a
square — only the multipliers $1$ and $u(1{-}u)$ ever arise, i.e. exactly
our four tensor blocks. So every coefficientwise-nonnegative $P$ has a
*diagonal*-Gram certificate, and the LP cone embeds. Strictness is cheap:
$(u-\tfrac12)^2$ has Bernstein coefficients $(\tfrac14,-\tfrac14,\tfrac14)$
— a square that no amount of elevation to $n=4$ brings into the LP cone.

**Convexity** ($n = 2m$ even, $k \ge 1$). The local Hessian
$H_{\mathrm{loc}} = \begin{psmallmatrix} f_{uu} & f_{uv}\\ f_{uv} & f_{vv}
\end{psmallmatrix}$ gets the matrix version of the same certificate
(Scherer–Hol):

$$H_{\mathrm{loc}} = S_0 + u(1{-}u)S_1 + v(1{-}v)S_2 + u(1{-}u)v(1{-}v)S_3,
\qquad S_a = (I_2 \otimes z_a)^\top G_a (I_2 \otimes z_a),$$

with $G_a \succeq 0$ of size $2D_a$; the entry $(s,t)$ of $S_a$ reads off
the $(s,t)$ block of $G_a$, so the same $\Lambda$ machinery applies
blockwise. Entries of $H$ have mixed bidegrees $(n{-}2,n)$, $(n{-}1,n{-}1)$,
$(n,n{-}2)$; all are degree-elevated to a common $(n,n)$ by
$A_2 = E_{n-1}E_{n-2}\,\bigl(n(n{-}1)\Delta^2\bigr)$ and
$A_1 = E_{n-1}\,(n\Delta)$, giving coefficient maps
$C_{uu} = I \otimes A_2$, $C_{uv} = A_1 \otimes A_1$, $C_{vv} = A_2 \otimes I$.
Two globalization facts close the argument: the physical Hessian is the
congruence $D\,H_{\mathrm{loc}}\,D$ with $D = \mathrm{diag}(1/w_x, 1/w_y)$,
so local PSD suffices per patch at any refinement level; and a $C^1$
piecewise-convex function on a convex domain is convex, so $k \ge 1$
continuity across patch lines and T-lines (which the base file provides)
upgrades per-patch certificates to global convexity. The builder asserts
$k \ge 1$.

## 2. Sheaf shape: certificates are vertex furniture

Per patch: the coefficient vertex $P$ (now `CofreeCone` — the cone moved
to the Grams), four `SemidefiniteCone` Gram vertices ($\mathrm{svec}$
coordinates: column-major **lower** triangle, diagonal first per column,
off-diagonals $\times\sqrt2$ — the library's `svec!` in `cone/sdp.jl`),
and one certificate hyperedge

$$P_{\text{map}}\,P \;-\; \textstyle\sum_a \Lambda_a\,\mathrm{svec}(G_a) = 0,$$

with $P_{\text{map}} = I$ (nonneg) or the stacked
$[C_{uu}; C_{uv}; C_{vv}]$ (convex). A symmetric weight matrix $W$
becomes an svec functional by the rule diag → $W_{cc}$, off-diag →
$\sqrt2\,W_{rc}$; that one rule handles diagonal *and* off-diagonal Gram
blocks uniformly. All continuity/T-edges act on the $P$ vertices exactly
as in the base file.

**H¹ is invariant.** The products $z_\alpha z_\beta$ of the block-0 basis
span all of bidegree $(n,n)$, so $\Lambda_0$ is surjective and
$\Lambda_0^\top$ injective; a left-null vector of the coboundary must
therefore vanish on every certificate edge, and the Gram vertices (each
touched by exactly one edge) contribute nothing. So `h1_predict_tj` —
junction rule plus measured ring terms — applies verbatim. Verified below,
including in the ring regime.

## 3. Verification (tjunction_sdp_oracle.py, CLARABEL)

Trust order: certificate maps against direct $z^\top G z$ evaluation →
Hessian maps against lower-degree differentiation and finite differences →
svec-sheaf assembly against an *independent* native model (PSD matrix
variables, raw ordered-pair index bookkeeping, no svec anywhere) →
semantics on dense grids → H¹ ranks.

```
certificate maps: scalar 1.11e-16   matrix-block 3.33e-16
hessian maps: elevation-vs-direct 2.66e-15   vs finite-diff 1.39e-07

[nonneg certificate, 2x2 refine (1,1), n=4 k=2, same data as LP]
  sheaf-svec vs native optimum : -21.32320150 / -21.32318811   |Δ| = 1.34e-05
  sections ‖ΔP‖∞ / ‖Bp‖        : 2.37e-05 / 3.12e-13
  sandwich  free ≤ SDP ≤ LP    : -21.496342 ≤ -21.323201 ≤ -19.936738   (LP−SDP gap = 1.39e+00)
  min Bernstein coeff  LP/SDP  : 2.33e-11 / -9.05e-01   (SDP < 0 allowed)
  min f on dense grid  SDP     : 3.81e-07   (certified ≥ 0)
  min f on dense grid  free    : -1.27e+00   (unconstrained)

[convex certificate (SOS-matrix Hessian), 2x2 refine (1,1), n=4 k=2]
  sheaf-svec vs native optimum : -12.61048770 / -12.61027220   |Δ| = 2.15e-04
  sections ‖ΔP‖∞ / ‖Bp‖        : 2.97e-04 / 1.49e-11
  min eig H (phys, dense grid) : free -1.214e+03   certified 1.218e-07
  fit-mse vs convex truth      : free 3.997e-03   certified 3.652e-04

H1 invariance (certificate edges must add no left-null vectors):
  2x2 refine (1,1) n=4 k=2           base 36   nonneg 36   convex 36
  3x3 refine ctr   n=4 k=2 (ring)    base 82   nonneg 82   convex 82
```

Three headlines worth keeping:

* **The gap is real and large.** On identical data the LP forgoes 1.39
  of objective relative to the SDP, and the SDP optimum carries a
  Bernstein coefficient of **−0.905** while its surface bottoms out at
  $+3.8\times10^{-7}$ on a dense grid — strict cone containment, made
  tangible. (The unconstrained fit dips to −1.27, so the constraint is
  genuinely active.)
* **Convexity works where nothing else exists** — free fit's physical
  Hessian reaches eigenvalue $-1.2\times10^{3}$ near a noisy corner;
  certified fit: $+1.2\times10^{-7}$, across coarse *and* half-width fine
  patches through the congruence argument.
* **Shape as regularization.** With deliberately tiny roughness
  ($\lambda_r = 10^{-6}$, $\sigma = 0.06$, $N = 250$), the convexity prior
  cuts fit MSE against the convex truth by **10×** (4.0e−3 → 3.7e−4). The
  certificate is not just a safety guarantee; on data whose truth has the
  shape, it is the better estimator.

## 4. Field notes

* **Even degree only.** The four-block certificate needs $n = 2m$; the
  builder asserts it. Odd $n$ has a Lukács-style variant with multipliers
  $u$ and $(1{-}u)$ and bidegree-$\tfrac{n-1}{2}$ blocks — same
  machinery, different `cert_blocks`; not implemented.
* **ε picks the Gram.** Certificates are non-unique in $G$; the
  $\tfrac{\varepsilon}{2}\|\mathrm{svec}(G)\|^2 = \tfrac{\varepsilon}{2}\|G\|_F^2$
  regularization (mirroring the patch blocks' $\varepsilon I$, and needed
  by Uzawa's SPD requirement anyway) selects the minimum-Frobenius
  certificate and keeps the problem strictly convex.
* **Cost.** Per patch at $n=4$: nonneg adds Grams of sizes $9,6,6,4$
  (svec total 97 next to 25 coefficients); convex doubles the matrix
  sizes to $18,12,12,8$ (svec total 363). The certificate edge adds
  $(n{+}1)^2$ rows (nonneg) or $3(n{+}1)^2$ (convex) per patch.
* **Composability.** Nothing prevents attaching both certificates, or
  attaching them only on selected patches (e.g. convexity only where the
  data demands it) — certificates are per-vertex furniture and the
  T-machinery never sees them. Per-patch selectivity is a one-line edit
  to the builder loop.
* First Julia run: `nonneg_sdp_demo()` then `convex_sdp_demo()`; both
  print the H¹ prediction (36 for the default 2×2 config) which the IPM's
  dual coset behavior should reflect exactly as in the base file. The
  JuMP reference `build_jump_tjsdp` accepts any SDP-capable optimizer
  (Clarabel.jl; Mosek once the exp-cone scaling verdict is in).

## References

Schmüdgen (1991), *The K-moment problem for compact semi-algebraic sets*.
Scherer & Hol (2006), *Matrix sum-of-squares relaxations for robust
semi-definite programs*. Hilbert (1888) / Motzkin (1967) for
nonnegative ≠ SOS. Lukács representation as in Karlin–Studden; see the
1D notes in `nonnegative_spline.md`.

## Files

* `tjunction_sdp.jl` — includes `tjunction_spline.jl`; certificate maps,
  builder, demos, JuMP reference.
* `tjunction_sdp_oracle.py` — imports `tjunction_spline_oracle.py`,
  `tensor_spline_oracle.py`, `mle_spline_oracle.py` (ship the four
  together); all numbers above.
