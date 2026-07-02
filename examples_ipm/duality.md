# Duality for Sheaf Section Problems: the Conic Lift is Universal

*Two problems over a cellular sheaf, and one duality joining them. The **general** form puts a
separable convex cost on the cells and couples them by a co-boundary equality; this is
Rockafellar's monotropic program, and its dual is again a separable problem on the same sheaf —
**conjugate** costs evaluated at the **divergence** of an edge cochain. The **conic linear**
form — linear objective, cones on the cells — is the special case where each cell cost is
"linear plus a cone indicator," and its dual is the familiar "co-boundary-image lands in the dual
cone." The closing fact is that the second form is **universal**: every problem of the first kind
lifts to one of the second, one dimension up, and the lift commutes with duality — so the simple
conic program is the right normal form on **both** faces at once. This note supersedes the conic
lifting note, folding it into the primal–dual picture the solver actually runs.*

---

## 0. Setup and conventions

A cellular sheaf $\mathcal F$ on a graph: each vertex $v$ carries a stalk $\mathcal F(v)$, each
edge $e$ a stalk $\mathcal F(e)$, and each incidence a restriction map
$A_{ve} := \mathcal F_{v\lhd e}\colon \mathcal F(v)\to\mathcal F(e)$. Cochains live on
$C^0=\bigoplus_v\mathcal F(v)$ (vertices) and $C^1=\bigoplus_e\mathcal F(e)$ (edges); the
co-boundary is $\delta\colon C^0\to C^1$.

Orient each edge $e$ with a tail (**lower**) endpoint $t(e)$ and a head (**upper**) endpoint
$h(e)$, and collect per vertex the upper/lower stars

$$
\mathrm{adj}^+(v)=\{e:h(e)=v\},\qquad \mathrm{adj}^-(v)=\{e:t(e)=v\}.
$$

Orientation is a free choice; the restriction maps are intrinsic. The co-boundary and its adjoint
$\delta^{\mathsf T}=\partial$ are then sign-free — head minus tail on edges, in minus out at
vertices:

$$
(\delta x)_e = A_{h(e)\,e}\,x_{h(e)} - A_{t(e)\,e}\,x_{t(e)},
\qquad
(\delta^{\mathsf T}y)_v = \sum_{e\in\mathrm{adj}^+(v)} A_{ve}^{\mathsf T}y_e
                       \;-\sum_{e\in\mathrm{adj}^-(v)} A_{ve}^{\mathsf T}y_e .
$$

For the constant sheaf ($A_{ve}=I$) these are $x_j-x_i$ and the signed graph divergence.

The two problems, side by side:

$$
\text{(P)}\quad \min_{x\in C^0}\ \sum_i f_i(x_i)\ \ \text{s.t.}\ \ \delta x=g
\qquad\Longleftarrow\qquad
\text{(P}_{\text{LP}}\text{)}\quad \min_{x\in C^0}\ \langle c,x\rangle\ \ \text{s.t.}\ \ \delta x=g,\ x\in K .
$$

(P) is general: each $f_i\colon\mathcal F(i)\to\mathbb R\cup\{+\infty\}$ is closed proper convex,
so a cone constraint $x_i\in K_i$ is simply an indicator summand inside $f_i$. (P$_{\text{LP}}$)
is the special case $f_i=\langle c_i,\cdot\rangle+\iota_{K_i}$. The arrow is the lift (§5); it runs
the other way as a specialization.

Throughout, regularity (a relative-interior / Slater point) is assumed — guaranteed in this corpus
because $g$ is **exact by construction**, $g=\delta x_0$ for an apex/homogenizer, so the affine
slice $\{\delta x=g\}$ meets $\operatorname{ri}\operatorname{dom}\sum f_i$.

---

## 1. The general problem and its dual (monotropic)

Attach a multiplier $y\in C^1$ — one per edge equation — and split the Lagrangian over vertices:

$$
L(x,y)=\sum_i f_i(x_i)-\langle y,\delta x-g\rangle
=\langle g,y\rangle+\sum_i\big[f_i(x_i)-\langle (\delta^{\mathsf T}y)_i,\,x_i\rangle\big].
$$

Each vertex minimizes independently, and $\inf_{x_i}[f_i(x_i)-\langle w_i,x_i\rangle]=-f_i^*(w_i)$
is a Legendre conjugate. With $w:=\delta^{\mathsf T}y$ the divergence,

$$
\boxed{\ \text{(D)}\qquad \max_{y\in C^1}\ \ \langle g,y\rangle-\sum_i f_i^*\big((\delta^{\mathsf T}y)_i\big)\ }
$$

The dual is **the same shape as the primal, read on the cosheaf**: a separable sum over the *same
vertices*, of the *conjugate* costs $f_i^*$, evaluated at the *divergence* $\delta^{\mathsf T}y$ of
an edge cochain. The multiplier $y$ is the only object that lives on edges; the cones/costs stay on
the vertices. This is monotropic self-duality — cost-on-sections dualizes to
conjugate-cost-on-divergences, $\delta\leftrightarrow\partial$ — exactly the cosheaf with stalks
$\mathcal F(v)^*$ and boundary $\partial=\delta^{\mathsf T}$.

---

## 2. KKT for the general problem

Two conditions, regularity assumed:

$$
\boxed{\quad
\delta x=g
\qquad\text{and}\qquad
(\delta^{\mathsf T}y)_i\in\partial f_i(x_i)\ \ \forall i .
\quad}
$$

There is **no separate complementarity line** — it is absorbed into the subdifferential inclusion.
By Fenchel–Young, at each vertex

$$
(\delta^{\mathsf T}y)_i\in\partial f_i(x_i)
\iff
f_i(x_i)+f_i^*\big((\delta^{\mathsf T}y)_i\big)=\big\langle x_i,(\delta^{\mathsf T}y)_i\big\rangle
\iff
x_i\in\partial f_i^*\big((\delta^{\mathsf T}y)_i\big),
$$

and that per-vertex equality is what plays the role of complementarity. In the upper/lower
notation the inclusion reads, cell by cell,

$$
\sum_{e\in\mathrm{adj}^+(i)} A_{ie}^{\mathsf T}y_e-\sum_{e\in\mathrm{adj}^-(i)} A_{ie}^{\mathsf T}y_e
\ \in\ \partial f_i(x_i):
$$

the tensions $y_e$ on the star of $i$, pulled back by the adjoint restrictions, must form a
subgradient of the local cost. **Edges carry agreement** ($\delta x=g$, reading $x$ at two
endpoints); **vertices carry the subgradient balance** (reading $x_i$ and the $y_e$ on the star).

**Weak duality in one line.** With $\delta x=g$, $\langle g,y\rangle=\langle x,\delta^{\mathsf T}y\rangle=\sum_i\langle x_i,w_i\rangle$, so

$$
p-d=\sum_i f_i(x_i)+\sum_i f_i^*(w_i)-\langle g,y\rangle
=\sum_i\underbrace{\big[f_i(x_i)+f_i^*(w_i)-\langle x_i,w_i\rangle\big]}_{\ge\,0\ \text{(Fenchel–Young)}}\ \ge 0 .
$$

Zero gap $\iff$ every per-vertex Fenchel–Young gap vanishes $\iff$ the inclusions hold. Each
summand is a Bregman-type divergence on a stalk; in the conic case it degenerates to a bilinear
pairing (§3).

---

## 3. The conic linear program as the special case

Take $f_i=\langle c_i,\cdot\rangle+\iota_{K_i}$. The conjugate is a dual-cone indicator,
$f_i^*(w_i)=\iota_{K_i^*}(c_i-w_i)$, so (D) collapses to a feasibility-constrained linear program.
Writing the dual slack $s:=c-\delta^{\mathsf T}y$:

$$
\text{(P}_{\text{LP}}\text{)}\quad \min_x\ \langle c,x\rangle\ \ \text{s.t.}\ \ \delta x=g,\ x\in K
\qquad\Longleftrightarrow\qquad
\boxed{\ \text{(D}_{\text{LP}}\text{)}\quad \max_y\ \langle g,y\rangle\ \ \text{s.t.}\ \ s:=c-\delta^{\mathsf T}y\in K^*\ }
$$

The cones are **still on the vertices** — $s\in K^*=\prod_v K_v^*$ — and $y$ is **still free on the
edges**. Nothing migrated to the edges; this is the cosheaf reading of §1 with the conjugate
specialized to a cone indicator. (The genuinely edge-coned problem — cone-constrained flows,
$\partial f=g$ — is a *different* primal, out of class; see §7.)

**KKT** (the redundant Jordan condition removed):

$$
\boxed{\quad
\delta x=g,\ \ x\in K
\qquad
\delta^{\mathsf T}y+s=c,\ \ s\in K^*
\qquad
\langle x,s\rangle=0 .
\quad}
$$

Cell by cell, in the upper/lower notation:

$$
\begin{aligned}
&\textbf{agreement (edges):} && A_{je}x_j-A_{ie}x_i=g_e,\quad x_v\in K_v,\\
&\textbf{balance (vertices):} && s_v=c_v-\!\!\sum_{e\in\mathrm{adj}^+(v)}\!\!A_{ve}^{\mathsf T}y_e+\!\!\sum_{e\in\mathrm{adj}^-(v)}\!\!A_{ve}^{\mathsf T}y_e\ \in K_v^*,\\
&\textbf{complementarity:} && \langle x_v,s_v\rangle=0 .
\end{aligned}
$$

Scalar complementarity is the certificate: with $x_v\in K_v$, $s_v\in K_v^*$ each
$\langle x_v,s_v\rangle\ge0$, so $\langle x,s\rangle=0$ forces termwise complementarity. The
vector form $x_v\circ s_v=0$ is **equivalent** to $\langle x_v,s_v\rangle=0$ under cone membership
(for symmetric cones, $\langle x,s\rangle=\operatorname{tr}(x\circ s)$), hence redundant in the KKT;
and for the exponential/power stalks the lift produces there is **no** Jordan product, so the
scalar pairing is the only complementarity that exists. It earns its keep solely in the iteration:
the central path $x_v\circ s_v=\mu e_v$ is the vector relaxation that squares the Newton step, and
its trace $\langle x,s\rangle=\nu\mu\downarrow0$ is the gap monitor.

**Gap identity.** $p_{\text{LP}}-d_{\text{LP}}=\langle x,s\rangle$ on a shared $x$ — the §2 Bregman
sum degenerated to a single bilinear pairing.

### 3.1 The native quadratic is just one conic cost

A block-diagonal $f_i=\tfrac12 x_i^{\mathsf T}Q_i x_i+c_i^{\mathsf T}x_i$ (optionally $+\iota_{K_i}$)
is smooth, so $\partial f_i$ is the singleton $\{Q_ix_i+c_i\}$ and the inclusion sharpens to the
*equation* $\delta^{\mathsf T}y+s=Qx+c$. Conjugating out $x$ gives the Dorn dual
$\max_{x,y}-\tfrac12 x^{\mathsf T}Qx+\langle g,y\rangle$ s.t. $Qx+c-\delta^{\mathsf T}y\in K^*$ —
the quadratic reappears with a flipped sign (it is $-f_i^*$), and the cones stay on the vertices. In
the Newton system the only footprint of $Q$ is to thicken the vertex $(1,1)$ block from $-H$ (the
cone scaling) to $-(Q+H)$; both summands are vertex-block-diagonal, so the augmented KKT matrix
keeps the exact bipartite sparsity of $\delta$. So the QP is **not** a separate generalization — a
quadratic is a conic-representable cell cost (rotated second-order cone), hence one more instance of
the lift below.

---

## 4. Dictionary: the two faces of one sheaf

| object | primal (section) | dual (cosheaf) |
|---|---|---|
| variable support | vertices $x\in C^0$ | edges $y\in C^1$ |
| operator | co-boundary $\delta$ | boundary $\partial=\delta^{\mathsf T}$ |
| stalk data | $\mathcal F(v)$, cone $K_v$ | $\mathcal F(v)^*$, cone $K_v^*$ |
| per-cell cost | $f_i$ | conjugate $f_i^*$ at $(\delta^{\mathsf T}y)_i$ |
| feasibility | $\delta x=g$ (agreement) | $c-\delta^{\mathsf T}y\in K^*$ (balance) |
| coupling | $\langle x_i,s_i\rangle=0$ / $(\delta^{\mathsf T}y)_i\in\partial f_i(x_i)$ | (same, read the other way) |
| freedom | sections $\ker\delta=H^0$ | cycles $\ker\delta^{\mathsf T}=H_1$ |

A primal–dual interior-point method carries **both** columns simultaneously: $x$ the section, $y$
the edge tension, $s$ the dual vertex slack. The convergence certificate is the residual triple
$\|\delta x-g\|$, $\|\delta^{\mathsf T}y+s-(Qx+c)\|$, $\langle x,s\rangle$ — primal infeasibility,
dual infeasibility, gap — driven to zero together.

---

## 5. The lift: every general problem becomes a conic LP

Two stalk-local moves turn (P) into (P$_{\text{LP}}$).

**Move 1 — epigraph (linearize the objective).** One scalar $t_i$ per vertex:

$$
\min_{x,t}\ \sum_i t_i\quad\text{s.t.}\quad \delta x=g,\ \ (x_i,t_i)\in\operatorname{epi}f_i .
$$

**Move 2 — homogenize (epigraph $\to$ cone slice).** Every closed convex set is a slice of its
homogenization. For the epigraph,

$$
\operatorname{epi}f_i=\{(x_i,t_i):(x_i,t_i,1)\in K_i\},
\quad
K_i=\operatorname{cl}\{(x_i,t_i,\lambda):\lambda>0,\ \lambda f_i(x_i/\lambda)\le t_i\},
$$

with $\lambda f_i(x_i/\lambda)$ the **perspective** of $f_i$ — positively homogeneous, so $K_i$ is a
closed convex cone (Rockafellar 1970). Pin the homogenizer $\lambda_i=1$. Stacking
$z_i=(x_i,t_i,\lambda_i)$, the pins fold into an augmented co-boundary exactly as an apex broadcasts
a constant through $g$, and the result is

$$
\boxed{\ \min_z\ \langle \mathbf 1,t\rangle\ \ \text{s.t.}\ \ \tilde\delta z=\tilde g,\ \ z\in\textstyle\prod_i K_i\ }
$$

a conic **linear** program with the co-boundary as its equality block and a product cone on the
augmented vertex stalks — **the same type, one dimension up**.

Why it stays the same type:

- **Stalk-local.** Each $f_i$ contributes only its own $t_i,\lambda_i,K_i$; no new coupling crosses
  an edge, so $\tilde\delta$ is $\delta$ plus block-diagonal local rows. The benign dimension
  increase — the opposite of a global-affine accumulator blow-up.
- **Degree-1 costs are free.** Norms, gauges, support functions need no homogenizer; cone
  indicators need no lift at all (that is the completion family, $f_i=\iota_{K_i}$). A non-conic
  convex set needs exactly one homogenizing coordinate (the elliptope/correlation case).
- **Data rides in the linear objective** $\langle\mathbf 1,t\rangle$; the co-boundary still carries
  only consistency plus the $\lambda=1$ pins.

The standard conic library covers the costs that arise — affine/quadratic and norms (SOC), $\ell_p$
and geometric terms (power cone), entropy / log-sum-exp / log-det (exponential and PSD cones) — so
the lift reaches every separable convex cost on the cells (Ben-Tal–Nemirovski 2001;
Nesterov–Nemirovski 1994; Boyd–Vandenberghe 2004).

---

## 6. The lift commutes with duality — why the simple form is universal

The lift is not just a primal convenience; it is **dual-compatible**, and that is what makes the
conic LP the right normal form on both faces.

Specialize the §1 dual to the conic LP and you recover §3's $s\in K^*$. Go the other way and the
structure closes: the per-vertex dual cone $K_i^*$ of the lifted primal is itself the
**homogenization of $\operatorname{epi}f_i^*$**. So the conic dual of the *lifted primal* is the
*lift of the unlifted dual* — the homogenization of $f_i$ and of its conjugate $f_i^*$ are polar
cones, and conic complementarity $\langle z_i,s_i\rangle=0$ is precisely the Fenchel–Young inclusion
$(\delta^{\mathsf T}y)_i\in\partial f_i(x_i)$ written one dimension up. Lifting and dualizing
commute:

$$
\begin{array}{ccc}
\text{(P) general} & \xrightarrow{\ \text{lift}\ } & \text{(P}_{\text{LP}}\text{) conic} \\[2pt]
\big\downarrow{\scriptstyle\text{ dual}} & & \big\downarrow{\scriptstyle\text{ dual}} \\[2pt]
\text{(D) conjugate} & \xrightarrow{\ \text{lift}\ } & \text{(D}_{\text{LP}}\text{) }s\in K^* 
\end{array}
$$

The conic LP is therefore universal in the strong sense: its primal absorbs every separable convex
cost on the cells (Move 1–2), and its dual simultaneously *is* the lift of the monotropic dual —
$f_i^*$ at the divergence becomes "vertex residual in $K_i^*$." One factorization, one pair of cones
$(K_i,K_i^*)$ per stalk, serves both the section problem and its conjugate.

---

## 7. Geometry, cohomology, and the two freedoms

The primal and dual feasible directions are Hodge-orthogonal. The primal $x$ ranges over
$x_0+\ker\delta=x_0+H^0$ — sections up to global symmetries; the dual divergence $w=\delta^{\mathsf T}y$
ranges over $\operatorname{im}\delta^{\mathsf T}=(\ker\delta)^\perp=(H^0)^\perp$. Fenchel–Young
couples a point of the section space to a point of its orthocomplement: the splitting
$C^0=\ker\delta\oplus\operatorname{im}\delta^{\mathsf T}$, with the nonlinearity living pointwise on
the stalks.

This surfaces **two cohomological freedoms**, dual to each other:

- **Primal: $H^0=\ker\delta$.** The section is determined only up to global sections (the formation
  guide's rigid motions); the cost selects which translate.
- **Dual: $H_1=\ker\delta^{\mathsf T}$.** The multiplier $y$ enters (D) only through
  $w=\delta^{\mathsf T}y$ and through $\langle g,y\rangle$, so $y$ is determined only up to the
  **cycle space** $\ker\delta^{\mathsf T}$. Adding any divergence-free edge cochain (a circulation,
  a harmonic $1$-cycle) leaves the dual unchanged.

The dual coset is well-defined **because $g$ is exact**: $g=\delta x_0\in\operatorname{im}\delta$
and any cycle $n\in\ker\delta^{\mathsf T}=(\operatorname{im}\delta)^\perp$ give $\langle g,n\rangle=0$,
so $\langle g,y+n\rangle=\langle g,y\rangle$. The same exactness that makes the primal feasible
($[g]=0\in H^1$) makes the dual coset consistent. A solver pins the primal coset through the cost
and the dual coset through its basepoint/regularization, but the two indeterminacies are genuinely
$H^0$ and $H_1$ of one sheaf.

---

## 8. Verification

A triangle sheaf — three vertices, edges $e_1=(1\!\to\!2),e_2=(2\!\to\!3),e_3=(1\!\to\!3)$,
vertex and edge stalks $\mathbb R^2$, rotation restriction maps with trivial holonomy (so
$H^0=\mathbb R^2$). Mixed cell costs: quadratic $\tfrac12\|x_1-a_1\|^2$, robust
$\|x_2-b_2\|_1$, exponential $\sum e^{(x_3)_k}$ — exercising rotated-SOC, SOC/box, and exponential
cones. Right-hand side $g=\delta x_0$ (exact). Solved three ways: direct convex (P), the monotropic
dual (D) reconstructed from the equality multiplier $y$, and the epigraph lift.

```
(P) direct primal opt            : 2.22216709
(D) dual  <g,y> - sum f*(δᵀy)     : 2.22216708     gap 5.3e-09
    primal feasibility ‖δx − g‖   : 3.8e-11
Fenchel–Young gaps (v1,v2,v3)     : 2.8e-17, 4.7e-09, 6.8e-10   (KKT inclusions hold)
    incl v1  ‖δᵀy − ∇f1‖          : 3.0e-11        (quadratic: ∇f1 = x1 − a1)
    incl v3  ‖δᵀy − ∇f3‖          : 3.1e-05        (exponential: ∇f3 = exp(x3))
    dual feas v2  ‖δᵀy‖_∞ ≤ γ     : 1.000 ≤ 1.0    (ℓ1 dual cone, binding)
cohomology                        : dim H⁰ = 2,  dim H₁ = 2
    ⟨g, n⟩ = 0  for n ∈ ker δᵀ    : 2.3e-16        (exactness ⇒ dual coset well-defined)
    dual obj after y += 1.7·cycle : 2.22216708     (unchanged — H₁ freedom)
lift: epigraph linear-obj opt     : 2.22216708     |direct − lifted| = 4.1e-09
    epigraphs tight  t − f(x*)    : 1.5e-9, 2.8e-9, 9.9e-9
```

Primal and dual meet; the per-vertex Fenchel–Young gaps vanish (the KKT inclusions of §2); the
$\ell_1$ vertex sits on the boundary of its dual $\infty$-ball; the lift reproduces the optimum with
tight epigraphs; and the dual objective is invariant under adding a cycle, with $\langle g,n\rangle=0$
on $\ker\delta^{\mathsf T}$ confirming the $H_1$ coset freedom of §7.

---

## 9. Scope: the discipline, on both faces

The lift is *useful* exactly when the convex cost sits **on the cells** $x_i$ (the domain of
$\delta$), so that after lifting the cones live on the augmented cell stalks and the solver
factorizes a sparse up-Laplacian $\delta^{\mathsf T}W\delta$ (plus a vertex-block-diagonal $Q$).
The boundary of the class is one line — *is the cost on the cells, or on their disagreement?* — and
it reads the same on the dual side:

- **Cost on a difference $\delta x$** (a cone in the *codomain* of $\delta$ — "cones on edges").
  Needs a slack, $B=[\delta,-I]$, no longer a strict co-boundary. Eliminating the slack still gives
  $\delta^{\mathsf T}W\delta$ (fast), but it is outside the clean section semantics — and dually it
  is a cost on $y$ rather than a conjugate at $\delta^{\mathsf T}y$. Examples: total variation /
  fused lasso, isotonic regression, metric nearness, synchronization.
- **Conservation $\partial f=g$** (cone-constrained edge flows). Here variables live on edges and
  $B^{\mathsf T}B=\delta\delta^{\mathsf T}$ is the *edge* Laplacian, fill $\sim\sum_v\deg(v)^2$ —
  dense, no reformulation fixes it. This is the genuinely edge-coned primal that the §3 dual is
  *not*: its cone sits on edges no matter which side you dualize. Examples: min-cost flow, DC
  optimal power flow, spring/electrical networks. Strictly out.

In short: the cost and the cones belong on the cells; the co-boundary supplies only the equality
$\delta x=g$; and the conjugate cost and dual cones land back on the cells too. Cost on the cells,
agreement coupling, sparse $\delta^{\mathsf T}\delta$ — the class is self-dual, and the conic LP is
its normal form.

---

## References

- A. Ben-Tal, A. Nemirovski. *Lectures on Modern Convex Optimization.* SIAM (2001).
- S. Boyd, L. Vandenberghe. *Convex Optimization.* Cambridge University Press (2004). *(perspective, conjugacy, conic representability)*
- P. Goulart, Y. Chen. *Clarabel: An interior-point solver for conic programs with quadratic objectives.* arXiv:2405.12762 (2024).
- J. Hansen, R. Ghrist. *Toward a spectral theory of cellular sheaves.* Journal of Applied and Computational Topology **3** (2019), 315–358.
- Y. Nesterov, A. Nemirovski. *Interior-Point Polynomial Algorithms in Convex Programming.* SIAM (1994).
- R. T. Rockafellar. *Convex Analysis.* Princeton University Press (1970). *(perspective, homogenization, conjugate duality)*
- R. T. Rockafellar. *Network Flows and Monotropic Optimization.* Wiley (1984). *(separable convex duality over a subspace)*
