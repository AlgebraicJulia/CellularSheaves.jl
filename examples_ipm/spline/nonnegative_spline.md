# Nonnegative and shape-constrained Bézier splines as a cellular sheaf

A self-contained description of a shape-constrained spline-estimation problem — fit a
piecewise-Bézier curve that is provably nonnegative, monotone, or convex over its whole
domain — and how to pose it as a constrained section of a cellular sheaf. Unlike a
safe-corridor problem, where the cone constrains an *affine image* of the control points and
therefore needs an auxiliary slack, here the cone constrains the control points *themselves*:
it lives directly on the vertex stalk. That makes this the clean "cone-on-vertex" companion to
the "cone-on-edge" corridor.

Two faithful formulations are given: an **orthant (LP)** relaxation that is genuinely
cone-on-vertex, and an **exact (SDP/SOCP)** characterization. The exact one specializes to
second-order cones for low-degree pieces; that special case is discussed explicitly.

The construction follows Papp & Alizadeh, *Shape-Constrained Estimation Using Nonnegative
Splines* (JCGS 2014), whose exact characterization rests on the classical Karlin–Studden
representation of nonnegative univariate polynomials.

---

## 1. The problem

We are given noisy observations `(x_j, y_j)`, `j = 1..N`, of an unknown smooth function `f`
on an interval `[x_lo, x_hi]`, together with a *shape* assumption:

- **nonnegativity**: `f(x) ≥ 0` for all `x` (densities, intensities, occupancy/rate fields);
- **monotonicity**: `f' ≥ 0` (growth curves, dose–response, cumulative distributions);
- **convexity / concavity**: `f'' ≥ 0` or `≤ 0` (production frontiers, utility, cost curves).

We estimate `f` by a polynomial spline of degree `n` with continuity `C^k`, choosing the
piece coefficients to minimize a loss

```
    minimize    data-fit(f)  +  λ · roughness(f)
    subject to  f is C^k across knots
                f satisfies the shape constraint everywhere on [x_lo, x_hi]
```

where `data-fit` is the residual sum of squares `Σ_j (f(x_j) − y_j)²` and `roughness` is a
min-derivative penalty `∫ (f^(ρ))²` for some order `ρ` (e.g. `∫(f'')²` for a min-curvature
prior). Both terms are quadratic in the coefficients, so this is a convex quadratic program
over a cone. Density estimation adds one linear constraint (`∫ f = 1`); see §7.

The point of interest here is the shape constraint, which is a *pointwise* condition holding
at infinitely many `x`, yet — via Bernstein control points — collapses to a finite conic
constraint on the coefficients.

---

## 2. Bernstein pieces (the basis that makes shape constraints cheap)

Each spline piece is a Bézier segment. On the local coordinate `t ∈ [0,1]` (each knot interval
is affinely mapped to `[0,1]`, which keeps every constraint independent of knot spacing), a
degree-`n` piece is

```
    f(t) = Σ_{i=0}^{n} P_i · B_{i,n}(t),      B_{i,n}(t) = C(n,i) t^i (1−t)^{n−i}
```

with `n+1` control points (Bernstein coefficients) `P_0..P_n`. Three facts drive everything:

1. **Partition of unity + convex hull.** The `B_{i,n}` are nonnegative and sum to 1, so
   `f(t)` is a convex combination of the `P_i`. Hence
   ```
        all P_i ≥ 0   ⟹   f(t) ≥ 0 on [0,1].
   ```
   This is *sufficient* for nonnegativity (not necessary — see §4).

2. **Derivative = differenced control points.** The derivative of a degree-`n` Bézier is a
   degree-`(n−1)` Bézier whose control points are `n·ΔP_i = n·(P_{i+1} − P_i)`. Iterating,
   the `m`-th derivative has control points proportional to the `m`-th differences `Δ^m P`.
   Therefore
   ```
        monotone (f' ≥ 0)  ⟸  ΔP_i  ≥ 0  for all i,
        convex   (f'' ≥ 0) ⟸  Δ²P_i ≥ 0  for all i.
   ```
   Shape constraints on `f` become nonnegativity constraints on `f`, `ΔP`, or `Δ²P` — i.e. the
   *same* orthant condition applied to a fixed linear image of the control points.

3. **Bernstein ↔ monomial is a fixed linear map.** Writing the piece in the monomial basis,
   `f(t) = Σ_ℓ c_ℓ t^ℓ`, the coefficients are `c = M·P` with
   ```
        M[ℓ, i] = C(n,i) · C(n−i, ℓ−i) · (−1)^(ℓ−i)     (for i ≤ ℓ, else 0).
   ```
   This is needed only for the *exact* formulation (§5).

Why Bézier rather than B-splines: at equal degree and knots the two have the same degrees of
freedom, but the cone of nonnegative-coefficient B-splines is a *strict subcone* of the
nonnegative-coefficient Bézier splines (Papp 2011, Thm 3.6). Bézier gives a better
nonnegativity model at no extra cost, so we use it.

---

## 3. The sheaf

The estimator is a section of a cellular sheaf whose base graph is a caterpillar: a spine of
segment vertices joined by continuity edges, with pendant cone-leaves for the shape constraint
(and, in the exact case, for the certificate).

**Vertices (stalks).**

- **Segment vertex** `v_s`, one per spline piece `s = 1..M`. Stalk `= ℝ^{n+1}` (the Bernstein
  control points of that piece). This vertex carries the objective: its quadratic block is the
  per-piece Gram `Q_s = λ·G_ρ + (data-fit Gram)`, where `G_ρ` is the min-derivative Bernstein
  Gram and the data-fit Gram accumulates `Σ φ(x_j) φ(x_j)^T` over the data points falling in
  that piece (`φ(x)` is the Bernstein evaluation vector). Because each data point lies in one
  piece, `Q` stays block-diagonal per vertex.
- **Shape / certificate vertices** — pendant leaves, contents depend on the variant (§4, §5).

**Edges (restriction maps).** Each edge carries a restriction map on each endpoint stalk; the
coboundary of the edge is the difference of the two, required to lie in the edge's cone.

- **Continuity** `(v_s, v_{s+1})`: the jet operators `R` (right/`t=1`) and `L` (left/`t=0`),
  each `(k+1)×(n+1)`, giving `R·P^s − L·P^{s+1} = 0` — i.e. value and first `k` derivatives
  match across the knot. Edge cone `= {0}`.
- **Boundary pins** (degree-1): fix an endpoint jet where required (e.g. `f(x_lo)`, or a
  known boundary value); the map is a row of `L` or `R`, the constant rides in the right-hand
  side. Density estimation's `∫ f = 1` is a single such linear edge.
- **Shape edges**: attach the shape cone (§4, §5).

**Objective.** `min ½ Σ_s (P^s)^T Q_s P^s + c^T P` subject to `B·P = g` and the per-vertex
cones, where `B` is the coboundary assembled from the maps above. This is exactly the solver's
normal form `min ½ pᵀQp + cᵀp s.t. Bp = g, p ∈ K`.

The graph is a tree of caterpillars: eliminating the leaves leaves a block-tridiagonal spine,
so the node factorization is `O(M)`.

---

## 4. Variant A — orthant (LP)

Impose the *sufficient* Bernstein-coefficient condition. This is the clean cone-on-vertex form.

**Nonnegativity.** The segment vertex's own cone *is* the orthant:

```
    cone(v_s) = ℝ₊^{n+1}          (PositiveCone on the stalk)
    objective Gram Q_s on the same vertex
```

No auxiliary vertex, no reification. The cone is a property of the control points themselves —
this is the sense in which the problem is "properly cones-on-vertices." Contrast the corridor,
where the constraint `A·P ≤ b` acts on a wide affine image and must be reified with a slack;
here nonnegativity is native to the Bernstein coordinates.

**Monotone / convex.** The constraint is `Δ^m P ≥ 0`, a linear image of the stalk, so it does
need a carrier. The clean reading is *the derivative curve is itself a nonnegative Bézier*:
attach a leaf vertex `d_s` holding the `m`-th-difference control points, with

```
    edge:   Δ^m · P^s − d_s = 0            (D^m is the m-th difference operator)
    cone(d_s) = ℝ₊^{n+1−m}                 (PositiveCone)
    cone(v_s) = free (CofreeCone), carrying Q_s
```

The leaf is degree-1 and eliminates with zero fill, so it is cheap — but note that, unlike
nonnegativity, monotone/convex are one linear map away from the stalk. (Nonnegativity is the
special case where that map is the identity.) Combine constraints by attaching more than one
leaf (e.g. monotone *and* convex).

**Cross-knot shape.** Within a piece the leaf handles the difference; at a knot the condition
must also hold *across* the join. The `C^k` continuity edges already tie the boundary jets, so
the cross-knot difference is determined by them — just make sure the difference sequence is
read on the joined control points, not reset per piece.

**Fidelity.** The orthant is a *relaxation*: it optimizes over a strict subcone of the truly
shape-constrained splines. Concretely, `(x−0.5)²` is nonnegative on `[0,1]` but its degree-2
Bernstein coefficients are `(0.25, −0.25, 0.25)` — the orthant rejects a valid nonnegative
function. Two mitigations:

- **Sieve density.** Papp & Alizadeh (Thm 1) show the nonnegative-Bernstein-coefficient
  splines are *dense* in the nonnegative continuous functions as the mesh shrinks (the
  Bernstein basis qualifies because it sums to 1). So any target is approximable; the gap is
  only for individual boundary functions at fixed resolution.
- **Degree elevation.** For a function *strictly* positive on the interval, there is a degree
  above which all Bernstein coefficients become positive, so elevating the degree drives the
  orthant toward the exact cone — and does so benignly, since Bernstein conditioning improves
  under elevation. (A function with an interior zero, like `(x−0.5)²`, always has a
  nonpositive coefficient, since the minimum coefficient is `≤ min f ≤ 0`; density still
  gives approximation, just not exact membership.)

Cones/constraints in this variant are all linear, so the model is an LP (a QP with an LP cone,
given the quadratic objective).

---

## 5. Variant B — exact (SDP / SOCP)

Impose the *exact* pointwise condition `f ≥ 0 on [0,1]` per piece via the Karlin–Studden
representation of nonnegative univariate polynomials. On the scaled interval `[0,1]` it reads
(monomial coefficients `c = M·P`, indices from 0):

**Odd degree `n = 2k+1`.** `f ≥ 0` on `[0,1]` iff there exist symmetric PSD matrices
`X, Y ⪰ 0` of size `(k+1)×(k+1)` with

```
    c_ℓ = Σ_{i+j=ℓ} Y_ij  +  Σ_{i+j=ℓ−1} (X_ij − Y_ij),      ℓ = 0..n.
```

**Even degree `n = 2k`.** `f ≥ 0` on `[0,1]` iff there exist `X ⪰ 0` of size `(k+1)×(k+1)` and
`Y ⪰ 0` of size `k×k` with

```
    c_ℓ = Σ_{i+j=ℓ} X_ij  +  Σ_{i+j=ℓ−1} Y_ij  −  Σ_{i+j=ℓ−2} Y_ij,      ℓ = 0..n.
```

(These are the `a=0, b=1` specializations of Karlin–Studden's Propositions 1–2; the general
`[a,b]` form carries `a` and `b` in the coefficients, but the scaled representation removes
them.) Monotone/convex use the identical machinery applied to the derivative piece's
coefficients.

**As a sheaf.** The certificate matrices are auxiliary cone vertices; the coefficient identity
is an equality edge:

```
    segment vertex v_s : P^s ∈ ℝ^{n+1}, free, carrying Q_s
    aux vertex X_s : PSD matrix, cone = SemidefiniteCone(k+1)
    aux vertex Y_s : PSD matrix, cone = SemidefiniteCone(k+1 or k)
    edge (v_s, X_s, Y_s):   M·P^s  −  KS(X_s, Y_s)  =  0        (n+1 scalar rows)
```

where `M` is Bernstein→monomial and `KS(·)` is the linear antidiagonal-sum map above. So the
exact variant *does* reify — it hangs a semidefinite cone off each piece — structurally the
same move as the corridor slack, but with a PSD cone instead of an orthant. This is why the
orthant, not the exact form, is the "clean cone-on-vertex" one; the exact form is the
semidefinite-on-auxiliary-vertex cousin.

### 5.1 The SOCP special case (degree ≤ 3)

The certificate blocks have size `⌈(n+1)/2⌉`. For **cubic or lower** pieces (`n ≤ 3`) that is
`2×2` (or `1×1`), and a `2×2` PSD matrix is a single Lorentz cone:

```
    [[a, b], [b, c]] ⪰ 0    ⟺    (a+c, a−c, 2b) ∈ Q₃      (second-order cone, bound first)
```

So for `n ≤ 3` each certificate block becomes a `SecondOrderCone(3)` rather than a
`SemidefiniteCone`, and the exact model is an **SOCP** — no semidefinite solver needed. This is
not a separate construction: it is the exact characterization with every block small enough to
Lorentz-reduce. For `n ≥ 4` the blocks are `3×3` or larger, no reduction exists, and the exact
model is a genuine **SDP**.

| degree | certificate block size | exact cone | model |
|--------|------------------------|------------|-------|
| `n = 1,2,3` | `2×2` (and `1×1`) | `SecondOrderCone(3)` | SOCP |
| `n ≥ 4` | `⌈(n+1)/2⌉ ≥ 3` | `SemidefiniteCone` | SDP |

For the min-crackle pairing (`n = 9`) the exact form is an SDP with `5×5` blocks; the SOCP face
appears only if a cubic instance is added to the suite.

---

## 6. The three faces at a glance

There are **two** formulations. SOCP is not a third — it is the low-degree face of the exact
one.

| | orthant (LP) | exact, `n ≤ 3` (SOCP) | exact, `n ≥ 4` (SDP) |
|---|---|---|---|
| condition | Bernstein coeff `≥ 0` | Karlin–Studden | Karlin–Studden |
| exactness | sufficient (strict subcone) | exact | exact |
| where the cone lives | **on the segment stalk** | on an aux vertex | on an aux vertex |
| cone type | `PositiveCone` | `SecondOrderCone(3)` per block | `SemidefiniteCone` per block |
| reification | none (nonneg) / 1 diff-leaf (mono/convex) | aux certificate vertex | aux certificate vertex |
| cost | cheapest | cheap | heaviest |
| tighten by | degree elevation | — (already exact) | — (already exact) |

---

## 7. Applications and instances

- **Nonnegative regression.** Data fit + roughness, `P ≥ 0`. Orthant is the natural model;
  exact if boundary-tight nonnegativity matters.
- **Monotone / convex regression.** Add the `ΔP ≥ 0` and/or `Δ²P ≥ 0` leaves. Multiple shape
  constraints just add leaves. (For a concave `f`, nonnegativity on `[a,b]` reduces to
  `f(a) ≥ 0, f(b) ≥ 0` — two linear pins — a nice simplification.)
- **Density estimation.** Nonnegativity plus one linear edge `∫ f = 1`. For a spline the
  integral is linear in the coefficients (`∫₀¹ B_{i,n} = 1/(n+1)`, so on a knot interval of
  width `h` the piece integral is `h/(n+1) · Σ_i P_i`), so `∫ f = 1` is a single equality
  edge summing all pieces. The "integrate to one" set is also bounded, which keeps the problem
  well-posed. A least-squares data term (fit to a histogram/pilot estimate) keeps the objective
  quadratic; a maximum-likelihood objective `−Σ log f(X_i)` is convex but not quadratic and
  needs an exponential-cone/`log` capable solver.

---

## 8. Relationship to the corridor problem

The corridor and the shape-constrained spline are the two ends of the cone-placement spectrum,
and the difference is exactly *invertibility of the map between the stalk and the constraint*:

- **Corridor**: the constraint `A·P ≤ b` uses a **wide** face matrix `A` (more faces than
  dimensions). It cannot be folded into the stalk, so it is reified with a slack — a
  *cone-on-edge* problem written in cone-on-vertex form.
- **Shape (nonneg)**: the constraint `P ≥ 0` is on the coordinates themselves — a **native**
  orthant on the stalk. Monotone/convex use `Δ^m P ≥ 0`, and `Δ^m` is invertible-with-anchor,
  so it is a fixed change of coordinates rather than a wide projection; the orthant still lives
  on (a linear image of) the stalk, reified only by a cheap degree-1 leaf.
- **Shape (exact)**: the Karlin–Studden certificate reintroduces an auxiliary cone vertex — so
  the *exact* shape problem is again a reified, cone-on-auxiliary-vertex problem, this time
  with a semidefinite (or, for cubics, second-order) cone.

So the orthant nonnegativity model is the cleanest cone-on-vertex instance in the whole suite,
and the exact model is its semidefinite reification — useful precisely because it exercises the
`SemidefiniteCone` (or, at low degree, `SecondOrderCone`) path directly on the certificate.

---

## 9. Implementation summary

Numerics reused from the Bézier corridor generator: the min-derivative Bernstein Gram, the
difference operator, the jet operators `L`/`R`, and Bernstein evaluation. Two builders:

- `nonnegative_spline_lp.jl` — orthant variant. `shape ∈ {:nonneg, :monotone, :convex}`.
  `PositiveCone` on the segment stalk for `:nonneg`; a `PositiveCone` difference-leaf for
  `:monotone`/`:convex`. Pure LP cone.
- `nonnegative_spline_exact.jl` — exact variant. Builds the Karlin–Studden certificate per
  piece and auto-selects the cone by degree: `SecondOrderCone(3)` blocks for `n ≤ 3` (SOCP),
  `SemidefiniteCone` blocks for `n ≥ 4` (SDP). One code path, one degree branch.

Both share: `C^k` continuity edges on the spine, a data-fit + roughness objective, and an
optional `integrate-to-one` edge for density estimation.

Two solver-API conventions to confirm against your cone implementations: `SecondOrderCone`
takes the bound as its **first** coordinate (`t ≥ ‖u‖`); `SemidefiniteCone` expects a
symmetric matrix in some vectorization (the exact file uses the standard `svec` with
off-diagonals scaled by `√2` so that `⟨A,B⟩ = svec(A)ᵀ svec(B)` — adjust if your cone expects
a different layout or scaling).
