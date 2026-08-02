# The augmentation window: justification and estimate

This note concerns the solution of a single saddle-point (KKT) system by
augmentation and Krylov correction, inside a refinement loop, and the
choice of the augmentation parameter $\alpha$.

We identify the set of $\alpha$ for which the augmentation alone suffices —
for which the system is solved to tolerance by one backsolve, with no
Krylov iteration and no refinement pass — and show that this set is an
interval whose two endpoints are computable from quantities the solve
records, with no fitted constants. We then show that admitting Krylov
iterations extends this interval downward, and admitting refinement passes
extends it upward, in two graded families obeying the same crowding law; and
that since both kinds of work are paid in the same unit, the quantity to
minimize is their sum, of which the window is the special case where both
vanish. Part IV states the procedure that consumes all of this: from one
completed solve, which $\alpha$ should have been used, returned as an
interval with a budget and a quality on each endpoint.

Sections 1–25 are self-contained and purely mathematical: no measurements
are reported, and every constant is derived from a stated assumption. They
concern one invocation of the solver. Section 26 is an appendix of a
different kind — it concerns a *sequence* of solves along an interior-point
trajectory, and it does report measurements, because the questions it raises
are not settled by derivation; it is marked provisional throughout. A computation requiring several
solves at a common $\alpha$ must intersect the intervals obtained for each;
that intersection is not treated here.

## 1. The solver

We solve

$$\begin{bmatrix} A & -B^\top \\ B & 0 \end{bmatrix}
\begin{bmatrix} x \\ y \end{bmatrix} =
\begin{bmatrix} f \\ g \end{bmatrix}$$

with $A$ symmetric positive semidefinite and $B$ of full row rank on the
relevant subspace. The system is indefinite and is not factored directly.

The semidefiniteness of $A$ is doing quiet work and is worth marking: it makes
the augmented operator positive definite for every $\alpha > 0$, so the
parameter is constrained from below by cost alone. Everything below relies on
that, and nothing below treats an indefinite $(1,1)$ block.

It is an assumption, and not a mild one in the setting the note points at.
An indefinite $A$ imposes a definiteness threshold
$\gamma_{\min} = -\lambda_{\min}(A)/\lambda^*_{\min}(B^\top B) > 0$
[RC, Thm. 1, p. 7] below which no $\alpha$ is admissible at all — a hard
lower bound, unlike the cost floor of §5, and one that binds before anything
derived here. Interior-point methods reach that regime routinely: the
leading block may lose definiteness and the constraint Jacobian may lose rank
as the iterations proceed [RC, Sec. 4.1, p. 7], which is the situation
[GS] exists to analyse. So the hypothesis is not a formality that a wider
treatment would absorb; it excludes a regime the cited implementations spend
part of their time in, and everything below — both endpoints, both graded
families, the cost — is conditional on it.
For a chosen $\alpha > 0$ (write $\beta = 1/\alpha$) form the **augmented
operator**

$$F = \beta A + B^\top B + \rho I = \beta A_\rho + B^\top B,
\qquad A_\rho = A + \alpha\rho I,$$

symmetric positive definite, admitting a Cholesky factorization. Here
$\rho \ge 0$ is a regularization chosen so that the factorization succeeds;
$\rho = 0$ is the unregularized case, to which every statement below
reduces. Since $\lambda_{\min}(F) \ge \rho$ while the top of the spectrum is
unchanged, $\kappa(F) \le \sigma_{\max}^2/\rho$ independently of $\alpha$,
which is the term's purpose. The system solved is always the original one:
$\rho$ enters through the factorization alone, and its perturbation is left
to the refinement loop. Its consequences are collected in §12.
Multiplying the first block row by $\beta$ and adding $B^\top$ times the
second yields the equivalent single equation

$$F x = \beta f + B^\top(g + \beta y). \qquad (1)$$

Equation (1) has two unknowns and is solved iteratively in $y$.

1. **Initial backsolve.** Given a starting multiplier $y_0$, solve
   $F x_0 = \beta f + B^\top(g + \beta y_0)$ by one backsolve and form the
   **initial constraint residual**
   $$r_0 = g - B x_0 .$$
2. **Krylov correction.** A Golub–Kahan/CRAIG iteration, preconditioned by
   $F^{-1}$, solves $B \delta x = r_0$ for the minimum-$F$-norm correction,
   stopping when the constraint residual is within the inner tolerance
   $\varepsilon$. Its iterates satisfy the invariant
   $$F \delta x = B^\top \delta y \qquad \text{at every iterate,}$$
   where $(\delta x, \delta y)$ are the primal iterate and its multiplier.
   We write $n$ for the number of iterations performed. The letter is also
   the dimension of $A$ where the operator is described ($B$ is
   $p \times n$); the two never appear in one expression, and the
   participation ratio of §17, which ranges over the $n$ rows of the
   perturbation, is the dimension.
3. **Recovery.** Assemble $x = x_0 + \delta x$, let $r = g - Bx$, and set
   $$y = y_0 + \alpha(\delta y + r).$$

**The stopping test.** The Krylov iteration is *entered* only if its
stopping criterion fails at the initial point:

> **(T1)** The iteration tests $\|g - Bx\| \le \varepsilon$ before its
> first step; if the test passes, $n = 0$.

The tolerance is **absolute** throughout, and that is a restriction on the
solves this note addresses rather than a normalization. Under a relative
test, $\|r_k\|/\|r_0\| \le \varepsilon_{\mathrm{rel}}$, the entry check can
never pass — the initial relative residual is $1$ — so $n \ge 1$ always and
$W_0$ is empty by construction rather than by (15). The graded family
survives but reindexes, to
$\alpha_k = (c_k/\varepsilon_{\mathrm{rel}})^{1/k}$: no $k = 0$ member,
crowding as $1/k$ rather than $1/(k+1)$, and — decisively — no $\ell_p$, so
the multiplier gap that is Part I's entire measured content does not enter.
That is a different object, not a rescaling of this one, and it is not
treated here. Implementations do use relative tests: [RC, Sec. 5.1, p. 12]
stop CG at $10^{-12}$ in relative residual norm, with no absolute pre-check.
A reader applying the estimate to such a solve must supply the absolute
reading itself.

**The refinement loop.** The assembled pair is checked against the original
two-block system, and further passes of the same machinery are applied to
the residual if the check fails:

> **(T2)** The refinement test is applied to the residual components
> separately: the pass fires unless *every* component of the two-block
> residual is within $\varepsilon$. In particular no component may borrow
> slack from another.

(T1) and (T2) are premises about the algorithm, not approximations. They
are the source of both constants in what follows.

**Norms.** We take $\varepsilon$ to be measured in a fixed norm in both
tests. If the two tests use inequivalent norms, the equivalence constants
enter as fixed factors multiplying each endpoint below, shifting the
interval but not its structure.

### The operator every spectral quantity is read on

The solve is of the augmented system

$$\begin{pmatrix} A + \alpha B^\top B & B^\top \\ B & 0\end{pmatrix}
  \begin{pmatrix} x \\ y \end{pmatrix} =
  \begin{pmatrix} c + \alpha B^\top d \\ d \end{pmatrix},$$

whose solution is that of the unaugmented system for every $\alpha$; here
$B$ is $p \times n$ and $A$ is $n \times n$. Write

$$F = \beta A + B^\top B, \quad \beta = 1/\alpha, \qquad S = BF^{-1}B^\top,$$

so that the $(1,1)$ block is $\alpha F$ and $\sigma(S) \subset (0,1)$.

**Every spectral quantity in this note is read in the variable $\mu$ of $S$,
and none is invariant under rescaling that operator.** The first CG step
length $a_1 = (r_0^\top r_0)/(r_0^\top S r_0)$ gives the elasticity
$\theta = 1/a_1 = \langle\mu\rangle_w$, which is dimensionless and — since
$\sigma(S) \subset (0,1)$ — at most $1$. The levels are the running product

$$e_i = \frac{\alpha\sqrt{b_i}}{a_i}, \qquad c_k = \prod_{i \le k} e_i ,$$

in which each term carries one factor of $\alpha$, so $c_k$ carries
$\alpha^k$ as [GG] and (21) require. The $\alpha$
belongs inside the product rather than applied to it afterwards: the
constants of (21) are norms in the variable $\tau = \alpha(1-\mu)$,
which is affine in $\mu$ with slope of magnitude $\alpha$, so a monic norm of
degree $k$ there is $\alpha^k$ times the corresponding norm in $\mu$. $e_i$
performs that change of variable one degree at a time. Running the
iteration on $\beta S$ instead of $S$ moves every one of these by a power of
$\alpha$.

$S$ is $\alpha$ times the matrix conventionally called the Schur complement
of this system, $-B(A + \alpha B^\top B)^{-1}B^\top$, up to sign; where that
object is referred to below it is named in full.

§1 writes the unaugmented system as
$\begin{psmallmatrix} A & -B^\top \\ B & 0\end{psmallmatrix}$ with the
original right-hand side, where the block above carries $+B^\top$ and an
augmented one. The two differ by the sign of $y$ and by which form of the
right-hand side is displayed; $F$, $S$, $\mu$ and every quantity built on
them are identical, and every result below transfers unchanged. The
one place the convention is visible is the dual residual, which is
$f - Ax + B^\top y$ there and $c - Ax - \alpha B^\top(Bx-d) - B^\top y$
here; both are the first block row of the original system, and $\|s_0\|$
means that row's norm in either.

The augmentation itself is standard. Golub and Greif [GG, eq. (2.1), p. 2078]
obtain it by multiplying the constraint row by $BW$ and adding it to the
first, for a general weight $W$; the case $W = \gamma I$ treated there
[GG, eq. (2.7), p. 2080] is the present $\alpha$, and their $B$ is the
transpose of the $B$ used here. The restriction to a scalar weight is a real
one: on a test problem where the null space of $A$ is available,
$W = (B B^\top)^{-1}$ — which makes the added term an orthogonal projector —
attains a condition number of $25.5$ against $44.2$ for the best scalar
$\gamma$ [GGV, Ex. 4.1, p. 788], and a rank-two weight aligned with $A$'s
null vectors is both sparser and better clustered than any scalar choice
[GGV, Ex. 4.3, pp. 789–790]. Everything below is stated for $W = \alpha I$,
in which the parameter is one number and can be read off a completed solve;
a non-scalar weight has no single $\alpha$ to estimate and is out of scope. What is new below is not the parameter but
how it is chosen: the augmented-Lagrangian literature selects $\gamma$ from
an a priori analysis of the operator, beginning with norm balance
[GG, Sec. 2, p. 2079] and continuing through the preconditioning literature,
where the analysis of the augmented system yields $\gamma = O(\nu^{-1})$
[BO, p. 2102]. That the a priori route is not the end of the matter is
argued by its own authors: they judge the result of limited practical use,
report instead a scaling rule $\gamma \approx \|w\|$
[BO, Rem. 4.3, p. 2102], and for small viscosity reduce $\gamma$ by hand,
running $1, 1, 1, 0.1, 0.02$ across four decades of $\nu$
[BO, Table 6.2, p. 2109]. The pattern persists: the augmentation is the
basis of the HyKKT method for interior-point solvers on GPUs
[PS, Sec. 3.1, p. 12], where $\gamma$ is again fixed by experiment, set to
$10^7$ throughout a large benchmark [PS, Sec. 6.2, p. 30] after a sweep over
four decades. The judgement is stated outright by another group working the
same block: selecting $\gamma$ in practice is difficult and requires problem
heuristics or trial and error, and a better method for selecting it is left
to future work [RC, Sec. 4.1, p. 9; Sec. 8, p. 18]. The estimate here
addresses that gap directly: it
recovers the parameter from a completed solve rather than from an analysis
of the operator, and returns an interval rather than a point.

## 2. The window

**The cost.** The unit of work is the triangular-solve pair against the
Cholesky factor: one is consumed by each backsolve and one by each Krylov
iteration, while the multiplications by $B$ and the residual evaluations are
of lower order. A solve performing $j$ backsolves — an initial one and
$j-1$ refinement passes — and $n$ Krylov iterations in total therefore costs
$n + j$ solve-pairs, and the object of the analysis is the $\alpha$
minimizing that sum. Iterations and passes are the same currency; neither is
privileged, and in particular there is no reason to forbid refinement.

**Definition.** For $k \ge 0$ and $j \ge 1$ let

$$W_{k,j} = \{\alpha > 0 : n(\alpha) \le k \ \text{ and }\ j(\alpha) \le j\},$$

and call $W_{0,1}$ — the set of $\alpha$ at which one backsolve suffices,
with no iteration and no refinement pass — **the window**, written $W_0$.

These definitions are absolute, not relative: they do not refer to a minimum
attainable cost, so they require no prior knowledge of the cost curve.
Sections 3–7 determine $W_0$, whose endpoints are the two convergence tests
read literally. Sections 9 and 10 determine the two graded families, which
extend the window downward by iterations and upward by passes; §11 assembles
them into the cost. The window is not the set of usable $\alpha$: it is the
set at which both counts are simultaneously minimal, a stronger condition
than optimality, and it may be empty while the cheapest solve is cheap.

By (T1) and (T2), $\alpha \in W_0$ iff both entry conditions hold at the
initial backsolve — the first denying the iteration any work, the second
denying refinement any:

$$\underbrace{\|r_0(\alpha)\| \le \varepsilon}_{\text{primal}}
\qquad\text{and}\qquad
\underbrace{\|s_0(\alpha)\| \le \varepsilon}_{\text{dual}} \qquad (2)$$

where $s_0$ denotes the first-row (dual) component of the two-block
residual of the assembled pair. The asymmetry between them is structural
and is the subject of the two parts that follow: the Krylov iteration acts
on the primal row, so the primal residual is the quantity it controls and
would reduce if entered; the dual row is untouched by the iteration, and
its residual is limited by the arithmetic. The remaining components of the
two-block residual are treated in §7.

Throughout, $\ell$ denotes an $\alpha$-independent quantity determined by
the solve.

# Part I — the primal condition

## 3. An exact identity for the initial residual

**Proposition 1.** *Let $(x^*, \lambda^*)$ solve the two-block system
exactly. Then, in exact arithmetic,*

$$r_0 = -\frac{1}{\alpha} S (y_0 - \lambda^*),
\qquad S = B F^{-1} B^\top. \qquad (3)$$

*Proof.* The exact solution satisfies both block rows, hence satisfies (1)
with $y = \lambda^*$: $F x^* = \beta f + B^\top(g + \beta\lambda^*)$.
Subtracting the defining equation of $x_0$,

$$F(x_0 - x^*) = \beta B^\top (y_0 - \lambda^*)
\Longrightarrow
x_0 - x^* = \beta F^{-1} B^\top (y_0 - \lambda^*).$$

Apply $B$ and use $B x^* = g$:

$$r_0 = g - Bx_0 = B(x^* - x_0)
= -\beta B F^{-1} B^\top (y_0 - \lambda^*). \qquad\square$$

The operator $S$ is symmetric positive semidefinite with eigenvalues

$$\mu_i = \frac{\alpha\sigma_i^2}{1 + \alpha\sigma_i^2} \in (0,1),
\qquad (4)$$

$\sigma_i$ being the singular values of $B$ in the metric induced by the
splitting. They are identified exactly, and not only asymptotically, by a
Sherman–Morrison argument: writing $S_\bullet(\gamma)$ for the augmented
Schur complement of the unaugmented one $S_\bullet$,
$S_\bullet(\gamma) = S_\bullet(I + \gamma S_\bullet)^{-1}$, so its eigenvalues
are $\beta_i/(1+\gamma\beta_i)$ with $\beta_i$ those of $S_\bullet$
[GGV, Prop. 2.7 and Cor. 2.8, p. 785]. Since this note's $S$ is $\alpha$ times
that object, the comparison with (4) is term by term and gives

$$\sigma_i^2 = \beta_i = \lambda_i\big(B A^{-1} B^\top\big),
\qquad \tau_i = 1/\sigma_i^2 ,$$

for every $\alpha$. So the limiting measure of §13 is carried on the
*inverse* spectrum of the unaugmented Schur complement, and $\alpha_\infty$
and $c_1$ — its capacity and its standard deviation — are properties of that
single object.

Two riders. The expression names $A^{-1}$, which §1's hypothesis does not
supply: where $A$ is singular the object is defined by the perturbation
$A(\epsilon) = A + \epsilon B^\top B$ and the limit taken, as
[GGV, §2.2, p. 784] do, and the $\sigma_i$ above are that limit's. And the
same identification arrives independently as a large-$\alpha$ expansion,
$\gamma S = (I+C)^{-1}$ with $C = \gamma^{-1}(BA^{-1}B^\top)^{-1}$
[RC, Thm. 4, p. 10], which [BO] attribute to Benzi and Liu; the exact form
above is preferred here because the note reads the measure at finite
$\alpha$. Directions with $\alpha\sigma_i^2 \gg 1$ have $\mu_i \approx 1$
and are said to be **resolved** by the augmentation; those with
$\alpha\sigma_i^2 \ll 1$ have $\mu_i \approx \alpha\sigma_i^2$ and
constitute the **tail**.

Identity (3) supplies, with no assumption, the prefactor $1/\alpha$:
whatever $\|S(y_0-\lambda^*)\|$ does, doubling $\alpha$ halves the residual
through the explicit factor.

## 4. Self-similarity, and inverse proportionality

Write $\delta = y_0 - \lambda^*$ and expand (3) in the eigenbasis of $S$:

$$\|r_0(\alpha)\|^2 = \sum_i \Big(\frac{\mu_i}{\alpha}\Big)^2 \delta_i^2,
\qquad
\frac{\mu_i}{\alpha} = \frac{\sigma_i^2}{1+\alpha\sigma_i^2}
\approx
\begin{cases} 1/\alpha, & \alpha\sigma_i^2 \gg 1, \\
\sigma_i^2, & \alpha\sigma_i^2 \ll 1, \end{cases}$$

so that

$$\|r_0(\alpha)\|^2
= \frac{\|P\delta\|^2}{\alpha^2} + W(\alpha), \qquad (5)$$

with $P = P(\alpha)$ the orthogonal projector onto the resolved directions
and $W = \sum_{\alpha\sigma_i^2 < 1} \sigma_i^4 \delta_i^2$ the tail
contribution. Resolved directions contribute inversely with $\alpha$; tail
directions contribute constants; the partition boundary
$\alpha\sigma_i^2 = 1$ moves down the spectrum as $\alpha$ decreases.

The assumption comes in three tiers, ordered by strength. The weakest is
what the derivation consumes; each stronger one implies the ones below and
adds something the derivation does not need.

> **A1 (scaling, local) — used.** *Let $\alpha$ be the parameter at which
> the solve was performed and $\alpha_0$ the endpoint sought. On the
> interval between them,* $\|r_0\| \propto 1/\alpha$.

> **A2 (scaling, global).** *The same, for every $\alpha$ in the operating
> range.*

> **A3 (mechanism).** *The gap's mass lies on directions the augmentation
> resolves:* $\|P\delta\| \approx \|\delta\|$ *and* $W$ *negligible. Then
> $\|r_0\| = \|\delta\|/\alpha$, and the constant of proportionality is
> named.*

A3 $\Rightarrow$ A2 $\Rightarrow$ A1. Only A1 enters the derivation of the
floor endpoint: the law needs the residual to scale inversely with $\alpha$
between where it was read and where the endpoint lies, and nothing more. A2
is cleaner to state, being free of any reference to where the solve happened
to be probed. A3 is the spectral account — it explains *why* the scaling
holds, supplies the constant $\|\delta\|$, and is what one would try to
establish for a class of operators. It is also strictly more than is needed:
any mechanism producing inverse scaling would serve equally well.

Under A3,

$$\|r_0(\alpha)\| = \frac{\ell_p}{\alpha}, \qquad \ell_p = \|\delta\|,
\qquad (6)$$

an exact inverse proportionality with an $\alpha$-independent constant. The
content of A is quantified by the **elasticity**

$$\theta(\alpha) \equiv -\frac{d \ln \|r_0\|}{d \ln \alpha}
= \frac{\sum_i \mu_i w_i}{\sum_i w_i}
= \langle \mu \rangle_w, \qquad
w_i = \Big(\frac{\mu_i}{\alpha}\Big)^2 \delta_i^2, \qquad (7)$$

obtained by differentiating (4) — $d\ln(\mu_i/\alpha)/d\ln\alpha = -\mu_i$
— so that $\theta$ is the residual-weighted mean of the spectrum of $S$.
Since $\mu_i \in (0,1)$ we have $\theta \in (0,1)$ identically: **the
residual can never fall faster than inversely with $\alpha$**, and A3 is the
statement $\theta \approx 1$, i.e. that the residual's weight sits on
resolved directions. As $\alpha$ decreases the boundary in (4) sweeps into
populated spectrum, $\theta$ decays, and the proportionality (6)
overstates the rate at which the residual falls.

## 5. The primal endpoint

By (T1) the primal condition of (2) is $\|r_0(\alpha)\| \le \varepsilon$.

**That this condition is a half-line requires no assumption.** By (7) the
elasticity is $\theta = \langle\mu\rangle_w$, an average of the spectrum of
$S$ against non-negative weights; and $\sigma(S) \subset (0,1)$ by §1. So
$\theta > 0$ strictly, for every $\alpha$ and every non-zero $r_0$, and
$\|r_0\|$ is strictly decreasing. The set $\{\|r_0\| \le \varepsilon\}$ is
therefore an interval unbounded above, whatever the residual's rate of
decrease may be. This is the same distinction §8 draws about contiguity:
the *structure* of the condition is derived from Proposition 1 alone, and
assumption A is needed only to say *where* the endpoint falls, not that
there is one. Its endpoint is where the residual meets the tolerance:

$$\boxed{\ \alpha \ge \alpha_0 = \frac{\ell_p}{\varepsilon}
= \frac{\alpha\|r_0(\alpha)\|}{\varepsilon}\ } \qquad (8)$$

the second form expressing the endpoint through the recorded residual at
any $\alpha$ where A holds — A entering here and not above. (In finite
precision the computed $\|r_0\|$ stops decreasing once it reaches the
rounding of forming $g - Bx_0$; the monotonicity is a statement about the
exact quantity, as everywhere in Part I.) Equivalently, in ratio form,

$$\frac{\alpha}{\alpha_0} = \frac{\varepsilon}{\|r_0(\alpha)\|}. \qquad (9)$$

**Remark (no free constant).** The number 1 implicit in (8) — the endpoint
is where the residual equals *one* tolerance — is (T1) read literally, not
a calibration: the iteration is entered exactly when the initial residual
fails the test it is given. Under a variant that tested a multiple of
$\varepsilon$, that multiple would appear in (8) unchanged in form.

**Remark (one-sidedness of the A-error).** By (7), $\theta \le 1$, so below
$\alpha_0$ the true residual is never smaller than (6) predicts. The
estimate (8) therefore errs, when it errs, by understating how far below
the endpoint a small $\alpha$ lies; it does not misplace the endpoint
itself while A holds near it.

# Part II — the dual condition

## 6. The dual residual is a finite-precision quantity

For a pair assembled naively through $F$ — any $x$ solving (1) with its
$y$ — the first-row residual is *not* a finite-precision quantity: the
defining equation forces it to equal $\alpha B^\top(g - Bx)$ identically,
an exact-arithmetic coupling to the constraint residual. The recovery step
of §1 cancels precisely this coupling.

**Proposition 2 (the dual row).** *In exact arithmetic the assembled pair of
§1 — $x = x_0 + \delta x$, $y = y_0 + \alpha(\delta y + r)$ with
$r = g - Bx$ — satisfies*
$$f - Ax + B^\top y = \alpha\rho\, x,$$
*for every $n \ge 0$, including $n = 0$; in particular the first block row
holds identically when $\rho = 0$.*

*Proof.* The defining equation of $x_0$ gives
$\beta(Ax_0 - f - B^\top y_0) + \rho x_0 = B^\top(g - Bx_0)$, i.e.
$f - Ax_0 + B^\top y_0 = -\alpha B^\top r_0 + \alpha\rho x_0$. The
invariant $F\delta x = B^\top\delta y$ gives
$\beta A \delta x = B^\top \delta y - B^\top B \delta x - \rho\,\delta x$,
i.e. $A\delta x = \alpha B^\top(\delta y - B \delta x) - \alpha\rho\,\delta x$.
With $r = r_0 - B\delta x$,

$$f - Ax + B^\top y
= \alpha B^\top(r - r_0 + B\delta x) + \alpha\rho(x_0 + \delta x)
= \alpha\rho\, x. \qquad\square$$

The dual residual is therefore the sum of an exact term proportional to
$\rho$ and a rounding term. Both are proportional to $\alpha$ and to the
solution produced, so both are absorbed into the single level $\ell_d$ of
§7; the consequences of the split are taken up in §12.

At $n = 0$ the recovery $y = y_0 + \alpha r_0$ alone performs the
cancellation. Consequently the computed $s_0$ consists entirely of
accumulated rounding error: every exact-arithmetic contribution vanishes by
construction. This licenses treating it as a sensor of attainable accuracy.

## 7. Conditioning, attainability, and the dual endpoint

> **S (two-scale spectrum).** *$B$ has a nontrivial near-kernel on which $A$
> acts with scale $h > 0$, while the top of $F$'s spectrum is pinned at
> $\sigma_{\max}^2$ independently of $\alpha$.*

S belongs to the mechanism layer (§24, group 3), not to the derivation. It
accounts for the *cost* side of the picture — the conditioning of $F$, and
hence the graded family — but not for the dual residual, which by (10b)
depends on $\|F\|$ alone and needs no statement about the near-kernel.

On the near-kernel the only support is $\beta A$, so the corresponding
eigenvalues of $F$ scale as $h/\alpha$, and

$$\kappa(F) \approx \frac{\sigma_{\max}^2}{h/\alpha}
= \frac{\alpha \sigma_{\max}^2}{h} \propto \alpha. \qquad (10)$$

This is the mirror of (6): raising $\alpha$ improves the primal residual
and degrades the conditioning, each in direct proportion.

Three tiers again, in the same order — but here the mechanism is a short
derivation rather than an assumption, and it delivers the global form
outright.

> **C1 (scaling, local) — used.** *On the interval between the parameter at
> which the solve was performed and the endpoint $\alpha_c$ sought,*
> $\|s_0\| \propto \alpha$.

> **C2 (scaling, above the switch).** *$\|s_0\| \propto \alpha$ for every*
> $\alpha > \alpha_\star$, *where*
> $$\alpha_\star = \frac{\lambda_{\max}(A)}{\sigma_{\max}^2}. \qquad (10a)$$

C2 implies C1 whenever the transport interval lies above $\alpha_\star$, and
unlike a bare global assertion it is checkable: (10a) is computable from two
spectral quantities of the operator pair, with no reference to where the
solve was probed. Under a regularized factorization the switch moves down
and the proportionality holds over a correspondingly wider range, the
regularization contributing a term exactly linear in $\alpha$ where the
rounding term is flat.

**C3 (mechanism).** Both statements follow from Proposition 2 and a standard
backward-error argument, with no further assumption. Let $\hat x_0 = x_0 +
\Delta x_0$ be the computed backsolve and $\hat y = y_0 + \alpha(g - B\hat
x_0)$ the recovered multiplier. Expanding the dual row and using
Proposition 2 to kill the exact-arithmetic part,

$$f - A\hat x_0 + B^\top \hat y
= -\big(A + \alpha B^\top B\big)\Delta x_0
= -\alpha F \Delta x_0 .$$

The computed solve is exact for a perturbed operator,
$\Delta x_0 = -F^{-1}\Delta F\, x_0$ with
$\|\Delta F\| \le c_u u \|F\|$, so the inverse cancels:

$$\|s_0\| \approx \alpha\,\|\Delta F\, x_0\|
\;\lesssim\; c_u\, u\, \alpha\, \|F\|\, \|x_0\|
\;=\; c_u\, u \max\big(\lambda_{\max}(A),\; \alpha \sigma_{\max}^2\big)\,\|x_0\|,
\qquad (10b)$$

the last step using $\|F\| = \max(\lambda_{\max}(A)/\alpha,\,
\sigma_{\max}^2)$.

**The constant $c_u$.** Here and throughout, $c_u$ is the backward-error
constant of the *solve* — the Cholesky factorization together with the two
triangular solves — defined by the property that the computed $\hat x_0$ is
the exact solution of a relatively perturbed system, $\|\Delta F\| \le c_u u
\|F\|$. The bound is relative, so it is indifferent to whether $F$ or
$\alpha F$ is named as the factored operator: both sides scale together, and
(10b) is unchanged. Four things about it.

*It is not universal.* Standard backward-error analysis of a Cholesky solve
bounds $c_u$ by a modest polynomial in the dimension, growing linearly in the
worst case; the realized value is far smaller, and a probabilistic analysis
of the accumulated roundings replaces the worst-case growth by roughly its
square root. So $c_u$ is a dimension-dependent quantity, not a constant of
the arithmetic, and this note does not pin it.

*It is independent of the parameter, the tolerance and the right-hand side.*
It is a property of the factorization applied to $\alpha F$, so it does not
move with $\alpha$ except through whatever the factorization does — in
particular it does not carry the $\alpha$-dependence that the transports
turn on, which is why (10b)'s two branches are $\alpha$-branches and not
$c_u$-branches.

*The estimate never uses it.* This is the reason its value can be left open.
The ceiling is $\alpha_c = \alpha(\varepsilon/\|s_0\|)^{1/\varphi}$, which
reads $\|s_0\|$ from the solve rather than predicting it, so $c_u$ cancels
out of every quantity the algorithm of §21 returns. It appears only in
statements *about* $\|s_0\|$ — the switch (10a), the equivalence of §16, the
ceiling's flat branch — which explain the estimate rather than compute it.

*The one exception is a design guideline, not a return value.* The shift
threshold $\rho_\star = c_u u \sigma^2_{\max}$ of §12 is a statement about
where a regularization ladder's first rung should sit, and it does carry
$c_u$ linearly. A ladder designed against an underestimate of $c_u$ places
rungs below the level at which they do anything.

Equation (10b) is flat in $\alpha$ below
$\lambda_{\max}(A)/\sigma_{\max}^2$ and linear above it. That is C2, with
the switch (10a) identified, and C1 follows on any interval above it.

Three things are worth noting about what does *not* appear in (10b). There
is no condition number: the $F^{-1}$ introduced by the backsolve is
annihilated by the $F$ that the recovery step reintroduces, so the dual
residual is governed by $\|F\|$ alone. There is no near-kernel scale $h$;
the bottom of $F$'s spectrum is irrelevant to it. And there is no condition
on which directions the right-hand side excites — the perturbation acts on
$x_0$ directly, and the bound follows whatever $x_0$ happens to be.

The practical content of the switch is that below $\alpha_\star$ the term
$\beta A$ dominates $F$, the achievable accuracy does not move with the
parameter, and the computed dual residual sits at a constant level carrying
no information about where the ceiling lies. Above it $B^\top B$ dominates,
$\|F\|$ grows in proportion to $\alpha$, and the residual grows with it.
The ceiling law reads that growth; a reading taken below $\alpha_\star$
extrapolates from a flat signal and the law fails accordingly.

Only C1 enters the derivation of the endpoint: the ceiling law needs the
dual residual to scale linearly between the reading and the endpoint, and
the constant of proportionality cancels. C2 and C3 establish that scaling
and locate where it begins.

Because the computed quantity cannot fall below the evaluation's own
rounding level, what is available is

$$\|s_0\| = \max(\ell_d \alpha,\ \text{evaluation floor}), \qquad (11)$$

a proportional ramp above a constant floor.

By (T2) the dual condition of (2) is $\|s_0\| \le \varepsilon$, which by
(11) is a half-line bounded above:

$$\boxed{\ \alpha \le \alpha_c = \frac{\varepsilon}{\ell_d}
= \frac{\alpha\varepsilon}{\|s_0(\alpha)\|}\ } \qquad (12)$$

again with the second form expressing the endpoint through the recorded
residual. In ratio form,

$$\frac{\alpha_c}{\alpha} = \frac{\varepsilon}{\|s_0(\alpha)\|}. \qquad (13)$$

**Remark (no free constant; the role of (T2)).** As in (8), the number 1 is
a test read literally. Here it depends on (T2) specifically: because the
refinement test admits no borrowing between components, the dual component
must clear $\varepsilon$ on its own, and the endpoint is its crossing. Under
a test that combined components — a norm of the stacked residual, say — the
dual allowance would be reduced by the share the other components occupy,
and (12) would carry that share as a factor.

**Remark (which component binds).** The other components of the two-block
residual do not compete for the endpoint. The second row is the constraint
residual, already within $\varepsilon$ by the primal condition; any further
scalar components of the assembled system are likewise controlled by the
iteration or by construction. The dual row is the one the iteration does
not act on, and by Proposition 2 it is pure roundoff; hence it is the
component whose failure defines the upper endpoint.

**Remark (one-sidedness).** By (11) the computed $\|s_0\|$ never falls
below the ramp, so (12) never overstates $\alpha_c$: when the evaluation
floor dominates — small $\alpha$, far below the endpoint — the estimate is
a valid lower bound on the true endpoint rather than an equality.

# Part III — the window and its graded family

## 8. The window is an interval

**Theorem.** *Under A, B and C, and with (T1), (T2),*

$$W_0 = \Big[\ \frac{\ell_p}{\varepsilon},\ \frac{\varepsilon}{\ell_d}\ \Big],
\qquad (14)$$

*an interval, nonempty if and only if*

$$\ell_p \ell_d \le \varepsilon^2, \qquad (15)$$

*with width (as a ratio)*

$$\frac{\alpha_c}{\alpha_0} = \frac{\varepsilon^2}{\ell_p \ell_d}. \qquad (16)$$

*Proof.* By §5 the primal condition holds exactly on $[\alpha_0, \infty)$
and by §7 the dual condition on $(0, \alpha_c]$; $W_0$ is their
intersection, which is $[\alpha_0, \alpha_c]$ when $\alpha_0 \le \alpha_c$
and empty otherwise. Substituting (8) and (12) gives (15) and (16).
$\square$

Contiguity is thus derived, not assumed: it follows from the monotonicity
of the two residuals in $\alpha$, which in turn follows from (6) and (10).
Both quantities in (15) are properties of the solve — $\ell_p$ the size of
the multiplier gap the solve must close, $\ell_d$ the amplification the
factorization contributes per unit $\alpha$ — and both are measured by the
recorded residuals through (8) and (12).

**Precision exhaustion.** When (15) fails there is no $\alpha$ at which the
augmentation alone suffices — though by §11 the solve remains possible, and
usually cheap, at a cost of one step on one staircase or the other: at every $\alpha$ either the backsolve leaves
too large a constraint residual or the factorization cannot deliver a dual
residual within tolerance. At equality the two endpoints coincide, and at
that point both residuals equal $\varepsilon$. Since $\ell_p$ and $\ell_d$
are properties of the problem and the arithmetic while $\varepsilon$ is
imposed, (15) states the precision available to the tolerance demanded: the
window closes from both ends as $\varepsilon$ tightens, and the quantity
$\ell_p \ell_d / \varepsilon^2$ is the natural measure of how close the
configuration is to exhaustion.

**Leaving the window.** Neither endpoint is a boundary of correctness.
Below $\alpha_0$ the iteration is entered and converges, at a cost
quantified in §9. Above $\alpha_c$ refinement fires and converges too, at a
cost quantified in §10: a pass's own dual residual is proportional to the
correction it computes, and those contract, so the passes are not held at
the base solve's level. Leaving the window downward buys iterations; leaving
it upward buys passes; both are paid in the same currency, and §11 weighs
them. The asymmetry that does survive is narrower and lies further out — the
primal family terminates, while the dual family meets a wall at
$\alpha_{\mathrm{wall}}$ beyond which no number of passes suffices.

## 9. The graded family below: Krylov iterations

For $k \ge 1$, $W_k$ relaxes only the primal condition: the dual condition is
unaffected by the iteration count, so all $W_k$ share the endpoint
$\alpha_c$.

Since the correction is sought in the Krylov space generated by $S$ from
$r_0$, every iterate's constraint residual is a polynomial in $S$ applied to
$r_0$: $r_j = p_j(S) r_0$ with $p_j \in \Pi_j$, $p_j(0) = 1$. Write

$$R_k(\alpha) = \min_{j \le k} \frac{\|r_j\|}{\|r_0\|}
\ \ge\ R_k^{\mathrm{opt}}(\alpha) = \min_{\substack{p \in \Pi_k \\ p(0)=1}}
\frac{\|p(S) r_0\|}{\|r_0\|} \qquad (17)$$

for the reduction the iteration actually achieves within $k$ steps, measured
in the norm of the stopping test. The inequality is not an artifact of
bookkeeping: CRAIG selects its polynomial by minimizing the error in the
$S$-norm — equivalently the residual in the $S^{-1}$-norm — not the residual
in the norm the test uses, so $R_k$ exceeds the min–max value
$R_k^{\mathrm{opt}}$ in general. Under A3 the discrepancy is slight, since the
residual's weight then sits where $\mu \approx 1$ and the two norms agree
there to within $1 + O(1/\alpha)$; but $R_k^{\mathrm{opt}}$ is a *lower*
bound on the cost, so substituting it for $R_k$ overstates how far the window
extends. This error has the same sign as the transport error of §10 and as
the degradation of A below the point of measurement; the three compound
rather than cancel, as recorded at the end of this section.

By (T1) applied at each step, $n(\alpha) \le k$ iff
$R_k(\alpha)\|r_0(\alpha)\| \le \varepsilon$, so

$$W_k = [\alpha_k, \alpha_c], \qquad
R_k(\alpha_k)\, \|r_0(\alpha_k)\| = \varepsilon, \qquad (18)$$

and $W_0 \subseteq W_1 \subseteq W_2 \subseteq \cdots$, a nested family with
a common upper endpoint.

**The rate of descent.** Expanding in the eigenbasis of $S$, for any
admissible $p$,

$$\frac{\|p(S)r_0\|^2}{\|r_0\|^2} = \frac{\sum_i p(\mu_i)^2 w_i}{\sum_i w_i}
= \big\langle p(\mu)^2 \big\rangle_w ,$$

the same residual weights $w_i$ as in (7). Taking $p(t) = (1-t)^k$ and using
$1 - \mu_i = 1/(1+\alpha\sigma_i^2)$: on resolved directions this is
$O(1/(\alpha\sigma_i^2))$, while on tail directions it is $O(1)$ — so the
bound must be taken in the weighted mean, not the maximum, and it is
A3, which places the weight on resolved directions, that makes the tail's
$O(1)$ contribution negligible. Hence
$\langle p(\mu)^2\rangle_w = O(\alpha^{-2k})$, and since the method does at
least as well as any fixed admissible polynomial in its own norm,

$$R_k(\alpha) = \frac{c_k}{\alpha^{k}}\,(1 + o(1)), \qquad (19)$$

with $c_k$ determined by the residual's spectral measure and independent of
$\alpha$. Note the dimensions: $R_k$ is a ratio, so $c_k$ carries the units
of $\alpha^k$; it is not a pure number, and not slowly varying in $k$.

**The constants are orthogonal-polynomial norms.** Both $R_k$ and $c_k$ have
exact characterizations in terms of the orthogonal polynomials of the
residual's spectral measure. Let $\{\pi_j\}$ be the monic orthogonal
polynomials of the normalized measure $\hat w$ placing mass $w_i$ at
$\mu_i$. Minimizing $\langle p(\mu)^2\rangle_{\hat w}$ over $p \in \Pi_k$
subject to $p(0) = 1$ is, in the orthonormal basis
$\varphi_j = \pi_j/\|\pi_j\|$, the minimization of $|c|^2$ subject to
$\sum_j c_j \varphi_j(0) = 1$, whence

$$R_k^{\mathrm{opt}}(\alpha)^2
= \Big(\sum_{j=0}^{k} \varphi_j(0)^2\Big)^{-1} \qquad (20)$$

— the Christoffel function of the measure evaluated at the origin, exactly,
for every $\alpha$ and every $k$. For the asymptotic constant, substitute
$\mu = 1 - t$ with $t_i = 1/(1+\alpha\sigma_i^2) = \tau_i/\alpha$,
$\tau_i = \alpha/(1+\alpha\sigma_i^2) \to 1/\sigma_i^2$: the constraint
$p(0)=1$ becomes normalization at $t = 1$, and placing the $k$ roots at the
scaled nodes leaves a monic polynomial in $\tau$ divided by $\alpha^k$.
Hence

$$c_k = \|\pi_k^{\tau}\|_{\hat w}
= \prod_{j=1}^{k} e_j, \qquad (21)$$

the $L^2(\hat w)$-norm of the degree-$k$ monic orthogonal polynomial of the
**limiting measure** — nodes $\tau_i = 1/\sigma_i^2$, weights
$\propto \delta_i^2$ — written as the running product of that measure's
Jacobi off-diagonals $e_j$, which is the $e_j$ of §1: the squared
off-diagonal is the recurrence coefficient, and $e_j = \alpha\sqrt{b_j}/a_j$
expresses it in the CG quantities, where $b_j$ is the CG ratio and not the
recurrence coefficient. The two are related by
$b^{\mathrm{rec}}_j = e_j^2 = \alpha^2 b_j/a_j^2$ and are nowhere
interchanged. The units are those of $\tau^k$, i.e. of $\alpha^k$, as
required.

*The node variable, exactly.* At a finite $\alpha$ the nodes are
$\tau_i = \alpha/(1+\alpha\sigma_i^2) = \alpha(1-\mu_i)$, which is affine in
$\mu$ with slope $-\alpha$; the recursion works in $\mu$ and therefore
returns norms in exactly this variable, one factor of $\alpha$ per degree.
The limiting form $1/\sigma_i^2$ is its $\alpha\sigma_i^2 \gg 1$ limit, and
the two differ by the factor $1/\mu_i$. Inside the window $\mu$ clusters near
$1$ and the distinction is immaterial, but it is not immaterial at moderate
$\alpha$, where the standard deviations in the two variables can differ by
several-fold. What the recursion computes is the affine variable; $c_k$ as
used below is that quantity, and the identification with the limiting
measure is asymptotic. The limiting measure is the residual's own measure
with the augmentation divided out, and its Jacobi coefficients are precisely
what the Golub–Kahan recursion of §10 computes: the constants of the graded
family and the quadrature of §10 are one object, seen twice.

An explicitly chosen polynomial (the $(1-t)^k$ above, or a Chebyshev
polynomial on the cluster) yields an upper bound on $c_k$ in place of (21),
and by (22) an upper bound on $\alpha_k$: the resulting estimate of the
window is conservative, understating the extension each iteration buys.

Combining (19) with (6) in (18),

$$\boxed{\ \alpha_k = \Big(\frac{\ell_p c_k}{\varepsilon}\Big)^{1/(k+1)}\ }
\qquad (22)$$

**Termination, and the intermediate regime.** The limiting measure is
supported on the $N \le m$ distinct values among the $\tau_i$ carrying
weight. Two consequences follow from (21), and they are different in kind.

*Termination.* For $k \ge N$ the monic orthogonal polynomial vanishes
identically, so $c_k = 0$, $R_k = 0$, and

$$W_k = (0,\ \alpha_c] \qquad (k \ge N):$$

with enough iterations the correction system is solved exactly and the primal
condition ceases to bind at any $\alpha$. The graded family therefore does
not accumulate; it terminates.

$N$ is left abstract here, but it is not structureless. In the preconditioned
spectrum of the same augmentation the eigenvalue $1$ carries multiplicity
$n - m + \operatorname{nullity}(A)$ [GGV, Thm. 2.5, p. 783] — clustering that
*improves* as $A$ loses rank. Repeated $\sigma_i$ are repeated $\tau_i$, so
whatever raises that multiplicity lowers $N$ and moves termination earlier.
A rank-deficient $(1,1)$ block, which is the difficult case for the
factorization, is the easy case for the primal family.

*The intermediate regime.* For $k \ll N$ the discrete measure is
indistinguishable from its continuum envelope, and there the classical
asymptotics of orthogonal polynomials apply: for a measure on an infinite
compact set $E$, $\|\pi_k\|^{1/k} \to \operatorname{cap}(E)$, the
logarithmic capacity of the support. Writing

$$\alpha_\infty = \lim_{k} c_k^{1/k} = \operatorname{cap}(E),
\qquad
\operatorname{cap}\big([\tau_-,\tau_+]\big) = \frac{\tau_+ - \tau_-}{4}
= \frac{1}{4}\Big(\frac{1}{\sigma_-^2} - \frac{1}{\sigma_+^2}\Big)$$

for the interval case — a scale with the units of $\alpha$, determined by the
extent of the weighted spectrum — and substituting
$\log c_k \approx k \log \alpha_\infty$ into (22) with
$\Lambda = \log(\ell_p/\varepsilon) = \log \alpha_0$,

$$\log \alpha_k \approx \log\alpha_\infty
+ \frac{\Lambda - \log\alpha_\infty}{k+1},
\qquad
\log\alpha_k - \log\alpha_{k+1} \approx
\frac{\Lambda - \log\alpha_\infty}{(k+1)(k+2)}. \qquad (23)$$

Equivalently, in exponentiated form,

$$\alpha_k \approx \alpha_\infty
\Big(\frac{\alpha_0}{\alpha_\infty}\Big)^{1/(k+1)}: \qquad (24)$$

the levels interpolate geometrically between the endpoint $\alpha_0$ of the
window itself and the scale $\alpha_\infty$, the logarithmic distance to the
latter contracting as $1/(k+1)$. They are thus not equally spaced — each
iteration buys a smaller extension than the last, the increments falling as
$1/k^2$ — and the distance is measured from $\alpha_\infty$, not from the
origin. The value $1$ has no meaning here: $\alpha$ carries units, and a
problem-independent accumulation point would be dimensionally impossible.

**Consequence: the budget for a given $k$.** Comparing (24) with the fixed
upper endpoint, $W_k$ is nonempty once

$$k + 1 \ \ge\ \frac{\log(\alpha_0/\alpha_\infty)}
{\log(\alpha_c/\alpha_\infty)}, \qquad (25)$$

the number of iterations that must be budgeted when $W_0$ is empty. If
$\alpha_\infty \ge \alpha_c$ the right-hand side is not positive and (25)
carries no information: the intermediate regime never delivers a window, and
one is available only as $k$ approaches $N$ — that is, only by solving the
correction system essentially exactly. This is not an impossibility, but it
is a change of regime, from a few iterations to a number comparable with the
dimension. In finite precision the exact termination at $k = N$ is not
attained, so the intermediate description governs whatever budget is
practical.

**Validity, and the signs of the errors.** The asymptotic (19) rests on A3 —
the weight concentrating on resolved directions. Equations (23)–(25) rest in
addition on $k \ll N$, and degrade as $k$ approaches $N$, where the true
behaviour is the termination above; they degrade also for measures whose
support is not well described by its convex hull, since the capacity is then
not that of an interval.

The several approximations do not have a common sign, and it is worth
recording which way each points. Four understate the level $\alpha_k$,
placing it lower than it truly is — an optimistic estimate of how far the
window extends:

- *A1 below the point of measurement.* By (7), $\theta \le 1$, so
 as $\alpha$ falls the true residual exceeds $\ell_p/\alpha$; more work is
 needed than (6) supposes, and the true level lies above the computed one.
- *The norm of (17).* Substituting the min–max value $R_k^{\mathrm{opt}}$ for
 the achieved $R_k$ credits the method with a reduction it does not deliver
 in the tested norm.
- *Transport (§10).* Node lumping and undiscovered spectrum both understate
 $R_k$ at parameters below the one the iteration ran at.
- *A single calibrated scale in place of $c_k$ (§10).* Since $c_k^{1/k}$
 ascends from the standard deviation of the limiting measure towards its
 capacity, replacing $c_k$ by $c_1^{\,k}$ understates it.

These share a sign and therefore compound: they do not cancel, and an
estimate carrying several of them is optimistic by their product. Two others point
the opposite way. Bounding $c_k$ by an explicit polynomial in place of (21)
overstates it, hence overstates $\alpha_k$; and the $o(1)$ of (19) overstates
$\alpha_k$ as the levels approach $\alpha_\infty$, where the resolved-limit
substitution loses its accuracy. Since this last grows in precisely the
regime where the others grow, no net sign can be asserted near
$\alpha_\infty$. In the asymptotic regime, with $c_k$ taken from (21) rather
than bounded, what remains is the optimistic group, and the computed
$\alpha_k$ is properly read as a lower bound on the level.

## 10. The graded family above: refinement passes

Section 9 relaxed the primal condition by admitting Krylov iterations. The
dual condition relaxes too, by admitting refinement passes, and the
structure that emerges is the mirror image.

**A pass and its contraction.** Suppose $\alpha$ exceeds $\alpha_c$, so the
dual component fails and a pass fires. The pass solves the original system
with the recorded residual as right-hand side. Its primal component is
already within $\varepsilon$, so to leading order the right-hand side is
$(s_j, 0)$, and the second block row of (1) forces the correction to lie in
$\ker B$. Writing $Z$ for a basis of that space and projecting the first
row onto it,

$$Z^\top A_\rho Z \,(Z^\top x_j) = Z^\top s_j,
\qquad\text{hence}\qquad
\|x_j\| \le \frac{\|s_j\|}{\lambda_Z},
\qquad \lambda_Z = \lambda_{\min}\big(Z^\top A Z\big), \qquad (26)$$

with $A_\rho = A$ in the unregularized case treated here. By §7 that pass
leaves in its turn a dual residual $\ell_d$ times its own $\alpha$-scaled
solution — the level of §7 is proportional to $\|x\|$, so writing
$\rho_{\mathrm{eff}} = \ell_d/\|x\|$ for the perturbation per unit solution,
the pass leaves $\alpha\rho_{\mathrm{eff}}\|x_j\|$. Composing with (26), the
dual residual contracts by

$$\varrho_d = \frac{\alpha\,\rho_{\mathrm{eff}}}{\lambda_Z} \qquad (27)$$

per pass. Two features of (27) are worth naming. It is *independent of the
Krylov work* done inside a pass: the correction $x_j$ lies in $\ker B$,
about which the iteration's equation $B\,\delta x = r_0$ says nothing, so
iterations cannot accelerate it. And it *grows* with $\alpha$, where the
primal side's contraction improves with $\alpha$ — the two mechanisms pull
oppositely, which is the source of everything below.

**The levels.** The base solve leaves a dual residual $\ell_d\alpha$, and
$j$ backsolves reduce it to $\ell_d\alpha\,\varrho_d^{\,j-1}$. Requiring
this to meet $\varepsilon$ and solving for $\alpha$:

$$\boxed{\ \alpha_c^{(j)}
= \big(\alpha_c\big)^{1/j}\,\big(\alpha_{\mathrm{wall}}\big)^{1-1/j},
\qquad
\alpha_c = \frac{\varepsilon}{\ell_d},
\quad
\alpha_{\mathrm{wall}} = \frac{\lambda_Z}{\rho_{\mathrm{eff}}}
= \frac{\lambda_Z \|x\|}{\ell_d}\ } \qquad (28)$$

or, in logarithms,

$$\log\alpha_c^{(j)} = \log\alpha_{\mathrm{wall}}
+ \frac{\log\alpha_c - \log\alpha_{\mathrm{wall}}}{j}. \qquad (29)$$

This is (24) reflected. The family of upper endpoints interpolates
geometrically between the one-backsolve ceiling $\alpha_c$ and a scale
$\alpha_{\mathrm{wall}}$, the logarithmic distance contracting as $1/j$; the
increments therefore fall as $1/j^2$, exactly as the primal levels crowd
about $\alpha_\infty$. The span of the whole family is

$$\log\frac{\alpha_{\mathrm{wall}}}{\alpha_c}
= \log\frac{\lambda_Z\|x\|}{\varepsilon},$$

so refinement can buy as many decades of ceiling as there are between the
tolerance and the near-kernel scale of $A$ — and no more.

**The wall.** Here the mirror is imperfect, and in the direction that
matters. The primal family terminates: for $k$ large enough the correction
system is solved exactly and $W_k$ extends to all small $\alpha$. The dual
family does not. By (27) $\varrho_d \ge 1$ once
$\alpha \ge \alpha_{\mathrm{wall}}$: the passes stop contracting, and no
number of them reaches the tolerance. $\alpha_{\mathrm{wall}}$ is a genuine
limit, not an accumulation point of levels that could be passed with
patience, and it is the first hard bound on $\alpha$ that this note has
produced. Its position is set by the same three quantities as everything
else: the near-kernel scale of $A$, the perturbation per unit solution, and
nothing about the tolerance.

## 11. The cost, and where it is least

Assembling the two families gives the cost of (2) as a function of
$\alpha$. The cost is *additive* in the two counts rather than
multiplicative, and since everything below rests on that, it is worth the
two lines it takes.

**Why $n + j$ and not $n \times j$.** A refinement pass is a full solve of
the correction system, so in principle it carries its own Krylov count and
the total is $\sum_{\text{passes}}(1 + n_{\text{pass}})$. What collapses
that sum is the pass's *initialization*. A pass has no prior multiplier
estimate — it solves for a correction — so it runs with $y_0 = 0$, and
Proposition 1 then gives its initial constraint residual as
$r_0^{\text{pass}} = \alpha^{-1} S\lambda^{*,\text{pass}}$, whence

$$\ell_p^{\text{pass}} = \alpha\|r_0^{\text{pass}}\|
= \|S\lambda^{*,\text{pass}}\| \;\le\; \|\lambda^{*,\text{pass}}\| ,$$

the inequality because $\sigma(S) \subset (0,1)$. So the pass's workload is
bounded by the multiplier of the system it is correcting — and that
multiplier is not the previous residual itself but $K^{-1}$ applied to it.
The relevant block of $K^{-1}$ acts on $\ker B$, where it is controlled by
the reduced Hessian, so

$$\ell_p^{\text{pass}} \;\lesssim\; \frac{\text{previous residual}}{\lambda_Z} ,$$

and additivity holds to the extent that the residual contracts faster than
$1/\lambda_Z$ inflates. That is worth stating rather than eliding, because
it makes the failure mode a familiar one: $\lambda_Z$ is the same quantity
whose collapse degenerates $\alpha_{\mathrm{wall}}$ in §18 and shrinks the
span $M$ in §17. A near-singular reduced Hessian is therefore not three
difficulties but one, appearing in three places — and the estimate records
$\alpha_{\mathrm{wall}}$, so the one condition is visible in what it
returns. Only the
first backsolve carries the full $\ell_p$, inherited from the *initial*
multiplier gap $y_0 - \lambda^*$; every later pass inherits a gap that has
already contracted. Once $\ell_p^{\text{pass}}/\alpha$ falls below
$\varepsilon$ the entry test (T1) passes and the pass costs one backsolve
and nothing more, which is the regime after the first pass.

The claim is checkable rather than asserted: $\ell_p^{\text{pass}}$ is
recorded by each pass, and additivity holds exactly to the extent that the
sequence contracts. Where it does not — a pass whose correction system is
no easier than the original — the cost is the sum above and $n + j$
understates it. To leading order, therefore,

$$\mathrm{cost}(\alpha) = n(\alpha) + j(\alpha),
\qquad
n(\alpha) = \min\{k : \alpha_k \le \alpha\},
\qquad
j(\alpha) = \min\{j : \alpha \le \alpha_c^{(j)}\}, \qquad (30)$$

the sum of a descending staircase and an ascending one, whose treads are
given by (24) and (29). Both arms are logarithmic in $\alpha$ and both
crowd as $1/k^2$, so the sum is a broad, flat-bottomed bowl rather than a
sharp minimum: the cheapest $\alpha$ is not a point to be hit but a range,
bounded below where the primal staircase flattens and above where the dual
staircase starts to climb.

Three consequences follow, and they revise the reading of §8.

*The window is a sufficient condition, not a necessary one.* $W_0$ is
$\{\alpha : n = 0, j = 1\}$, the set where both staircases stand at their
minimum simultaneously. That is strictly stronger than minimizing their
sum. When (15) fails and $W_0$ is empty, the cost (30) is still finite and
usually small — one arm simply pays a step so the other need not.

*The exhaustion criterion (15) is not a limit of the method.* It states
when one backsolve ceases to suffice, not when the tolerance ceases to be
attainable. The latter is governed by $\alpha_{\mathrm{wall}}$ on one side
and by the exact termination of §9 on the other.

*The endpoints are demoted, not discarded.* $\alpha_0$ and $\alpha_c$ are
the points at which each staircase reaches its lowest tread — $n = 0$ and
$j = 1$ respectively — and so remain the natural targets when they are
simultaneously reachable. Where they are not, (30) is the object to
minimize and the endpoints are its two anchors.

## 12. Regularization

Collecting the consequences of $\rho > 0$. Everything follows from
$A_\rho = A + \alpha\rho I$ of §1 — the note's statements with $A$ replaced
by an operator that depends on $\alpha$.

**The dual level.** By Proposition 2 the dual residual carries an exact term
$\alpha\rho x$ alongside the rounding term of (10b), so

$$\ell_d = \big(\rho + c_u u \sigma_{\max}^2\big)\|x\|,
\qquad \rho_\star = c_u u \sigma_{\max}^2. \qquad (31)$$

$\rho$ enters exactly where machine epsilon does. Below $\rho_\star$ it is
invisible; above it, $\alpha_c$ falls as $1/\rho$. Three further remarks.
The ceiling law (12) is *unchanged in form* — it reads the recorded $s_0$
as before and returns a smaller number. The scaling C1 is *strengthened*:
$\alpha\rho\|x\|$ is exactly linear where the rounding term is flat below
$\alpha_\star$, so the switch (10a) moves down to $u\lambda_{\max}(A)/\rho$
and the proportionality holds over a wider range. And because
$\alpha\rho x$ is an exact vector identity rather than a bound, it holds in
whatever norm the test of (T2) employs, with no equivalence constant: if
that test reads $\|s_0\|_\bullet \le \varepsilon\nu$ for an
implementation-fixed scaling $\nu$, then

$$\alpha_c^{\rho} = \frac{\varepsilon\,\nu}{\rho\,\|x\|_\bullet}
\qquad (32)$$

is computable *before* $\alpha$ is chosen — the only endpoint in this note
that need not be measured. It omits the rounding term and so overstates the
endpoint; it is sound for $\rho \gg \rho_\star$.

**The primal level.** The identity of §3 acquires a bias:
$r_0 = -\alpha^{-1}S_\rho\delta + \rho BF^{-1}x^\*$, whose second term is
$\alpha$-independent in the limit. Hence $\|r_0\|$ does not fall to zero but
onto a floor $\rho b_\infty$, $b_\infty = \|B(B^\top B + \rho I)^{-1}x^\*\|$,
so

$$\alpha_0 = \frac{\ell_p}{\varepsilon - \rho b_\infty}, \qquad (33)$$

unsatisfiable at any $\alpha$ when $\rho b_\infty \ge \varepsilon$. A1 thus
fails at *large* $\alpha$ as well as small. The mechanism is the
saturation: with $A_\rho$ the spectrum of $S$ approaches
$\sigma_B^2/(\sigma_B^2+\rho)$ rather than 1, so directions with
$\sigma_B^2 \ll \rho$ are never resolved. Two ways to handle the mixed
reading: the bias is bounded a priori by
$\rho b_\infty \le \tfrac12\sqrt{\rho}\,\|x\|$, which if small against
$\varepsilon$ settles the matter; otherwise $\rho BF^{-1}x$ costs one
solve-pair and is subtracted from $r_0$ as a vector. Unlike the dual term,
this damage is repairable by iterations, since $B\,\delta x = r_0$ is solved
whatever $r_0$ is.

**The wall becomes a knee.** In (27) the denominator is
$\lambda_{\min}(Z^\top A_\rho Z) = \lambda_Z + \alpha\rho$: the perturbation
shifts the very operator that bounds the correction. Hence
$\varrho_d = \alpha\rho_{\mathrm{eff}}/(\lambda_Z + \alpha\rho) < 1$ for
every $\alpha$, and **refinement always converges under regularization**,
however slowly. The wall of §10 becomes a knee at $\lambda_Z/\rho$; the
graded ceiling has no accumulation point, and the usable condition is
$\alpha\rho \lesssim \lambda_Z$ — a condition on the *product*, so that the
large $\alpha$ the primal side favours is what makes a given $\rho$
expensive. In the cost of §11, $\rho$ raises the ascending arm without
touching the descending one, and moves the optimum leftward.

**Window.** Combining (31) and (33), $W_0$ is nonempty iff
$\ell_p\ell_d \le \varepsilon(\varepsilon - \rho b_\infty)$, degrading (15)
on both factors.

An instance is on record. On one power-flow sequence [RC, Fig. 4, p. 13]
report that a Ruiz-scaled problem required $\delta_1 \approx 1$ to make the
augmented block definite, and abandon the method there on the grounds that
the regularization would carry as much weight as the problem
[RC, Sec. 5.2, p. 13]. In the terms above that is both failures at once: a
$\rho$ of order one collapses $\alpha_c^{\rho}$ by (32), and $\rho b_\infty$
overruns $\varepsilon$ so that (33) is unsatisfiable at any $\alpha$. Their
four remaining sequences ran at $\rho = 0$.

**A shift on the other block is not treated.** Implementations also
regularize the Schur complement itself, solving with $S + \delta_2 I$ in place
of $S$ when $B$ loses row rank [RC, Alg. 1, line 17, p. 9]. That shift moves
every node of the measure on which §9 and §13 read $c_k$, $\alpha_\infty$ and
$\theta$, so it is not the $\rho$ of this section under another name and
nothing here transports across it.

## 13. The spectral content of the cost

The constants $c_k$ of (19), characterized in (21) as norms of the limiting
measure's orthogonal polynomials, are not extra data: the iteration that
produced $r_0$ computes a finite representation of that measure.

**The blind case.** The representation has $k$ nodes after $k$ steps, so a
solve *inside* the window computes none: $n = 0$ is precisely the condition
defining $W_0$, and a backsolve that passes the entry test leaves no
spectral record whatever. The levels $\alpha_k$ below such a solve are
therefore not observable from it, and $\|r_0\|$ — the one number it does
record — determines $\alpha_0$ and nothing else. What is available cheaply
is the first level, from two moments of the measure. Since the weights $w_i$
are the squared components of $r_0$ in the eigenbasis of $S$,
$\langle\mu^j\rangle_w = r_0^\top S^j r_0 / r_0^\top r_0$ for every $j$; in
particular

$$\langle \mu \rangle_w = \frac{r_0^\top (S r_0)}{r_0^\top r_0},
\qquad
\langle \mu^2 \rangle_w = \frac{\|S r_0\|^2}{r_0^\top r_0},$$

both obtained from a single application of $S$ — one triangular-solve pair
and two multiplications by $B$ —

$$R_1^{\mathrm{opt}} = \Big(1 -
\frac{\langle\mu\rangle_w^2}{\langle\mu^2\rangle_w}\Big)^{1/2}
= \frac{\mathrm{sd}_w(\mu)}{\langle\mu^2\rangle_w^{1/2}}, \qquad (34)$$

the one-step min–max reduction exactly, hence $c_1$ and $\alpha_1$ by
(19) and (22). The same computation returns the elasticity
$\theta = \langle\mu\rangle_w$ of (7), so one application of $S$ reports
both the validity of A at the current $\alpha$ and the first level below
the window. Deeper levels require correspondingly more applications, at one
node each; this is the price of the graded family in the regime where the
solve itself is free.

**Calibration in place of measurement.** If even one application of $S$ is
unwanted, the first level may be estimated from a constant calibrated once
over a family of problems. By (21) the constant $c_1 = \|\pi_1^\tau\|$ is
the standard deviation of the limiting measure — of $\tau = \alpha(1-\mu)$,
the limiting $1/\sigma^2$
under the weights $\delta_i^2$ — and is therefore a property of the operator
pair and the direction of the multiplier gap alone: it does not depend on
$\alpha$, and it does not depend on $\varepsilon$. Given a calibrated
$\hat c_1$, (22) supplies

$$\alpha_1 \approx \Big(\frac{\ell_p \hat c_1}{\varepsilon}\Big)^{1/2}
= \Big(\frac{\alpha \|r_0\| \hat c_1}{\varepsilon}\Big)^{1/2}
\qquad (35)$$

from the recorded residual alone. Three properties recommend this over
calibrating a constant of the endpoint itself. It is $\varepsilon$-free,
where a constant expressing the endpoint as a fixed margin ratio is not: such
a formulation posits $\alpha_1 = \ell_p q/\varepsilon$ for a fixed $q$,
which by (22) requires $q = (\varepsilon c_1/\ell_p)^{1/2}$ — not a
constant at all, but a quantity falling as $\varepsilon^{1/2}$ and varying
with the multiplier gap. It is *damped*: $c_k$ enters $\alpha_k$ through the
$(k+1)$-th root, so an error of a factor $\phi$ in the calibration displaces
$\alpha_1$ by $\phi^{1/2}$ and $\alpha_2$ by $\phi^{1/3}$, where a constant
entering linearly would transmit $\phi$ in full. And it prices the whole graded family, not
only its first level: $c_1$ is the measure's standard deviation and
$\alpha_\infty$ its capacity, which for interval-like support agree to a
factor between $1$ and $\sqrt2$, so $c_1$ may stand for $\alpha_\infty$
in (24), giving

$$\alpha_k \approx \hat c_1
\Big(\frac{\alpha_0}{\hat c_1}\Big)^{1/(k+1)},
\qquad \alpha_0 = \frac{\alpha\|r_0\|}{\varepsilon},
\qquad (36)$$

from the same single constant. Equation (36) is exact at $k = 0$, where it
returns $\alpha_0$, and at $k = 1$, where it returns
$(\hat c_1 \alpha_0)^{1/2}$ and so reproduces (35); for $k \ge 2$ its only
approximation is $c_k \approx c_1^{\,k}$, the intermediate-regime statement
of §9. That approximation is one-signed: $c_k^{1/k}$ ascends from the
standard deviation towards the capacity, so $c_k$ exceeds $c_1^{\,k}$ and
(36) places $\alpha_k$ below its true value, joining the optimistic group of
§9 rather than opposing it. The exponent contains the resulting error, which
does not accumulate with $k$: the widening gap between $c_k$ and
$c_1^{\,k}$ is taken to the power $1/(k+1)$.

The validity condition is checkable from the estimate itself. Equation (19)
requires $c_1/\alpha \ll 1$ at the level being located, i.e.
$\alpha_1 \gg c_1$; as $\alpha_1$ descends towards $c_1$ the asymptotic
fails, which is the approach to $\alpha_\infty$ of §9 seen from below. Since
$\hat c_1$ and $\alpha_1$ are both in hand, the ratio reports whether the
estimate is in its domain.

**The quadrature.** Over $k$ iterations the Golub–Kahan recursion assembles
the Lanczos tridiagonal $T_k$ of $S$ with starting vector $r_0/\|r_0\|$.
Its eigenpairs give the Ritz values $\mu^{(k)}_j$ — approximations to the
eigenvalues $\mu$ of $S$, and written in that variable for that reason —
and their weights,

$$(\mu^{(k)}_j,\ w_j), \qquad w_j = \|r_0\|^2 (e_1^\top v_j)^2,
\qquad j = 1,\dots,k, \qquad (37)$$

which are the nodes and weights of the $k$-point Gaussian quadrature of the
residual's spectral measure on $S$. Since (17) depends on $r_0$ only
through that measure, and since the Krylov space of dimension $k$ is
determined by its first $2k$ moments — exactly those the quadrature
integrates exactly — the pairs (37) reproduce the iteration's own residual
history for the first $k$ steps in exact arithmetic. The representation is
lossless at $\alpha$, not an approximation of it.

**Transport.** The dependence on $\alpha$ in (4) is a bijection between
$\mu_i$ and $\sigma_i^2$, so the nodes may be stripped of the parameter and
re-dressed at any $\alpha'$:

$$\hat\sigma_j^2 = \frac{\mu^{(k)}_j}{\alpha(1-\mu^{(k)}_j)},
\qquad
\hat\delta_j^2 = w_j\Big(\frac{\alpha}{\mu^{(k)}_j}\Big)^2, \qquad (38)$$

$$\hat\mu_j(\alpha') = \frac{\alpha'\hat\sigma_j^2}{1+\alpha'\hat\sigma_j^2},
\qquad
\hat w_j(\alpha') = \Big(\frac{\hat\mu_j}{\alpha'}\Big)^2 \hat\delta_j^2,
\qquad (39)$$

(39) being identity (3) applied per node. Evaluating (17) on the
transported pairs gives $R_k(\alpha')$, hence $c_k$ and, through (22), the
levels $\alpha_k$ — all from the single solve, with no assumption of
self-similarity spent on the transported directions, since (39) carries
each node's $\alpha$-dependence explicitly.

**What transport does and does not determine.** The representation (37) is
exact at $\alpha$ for the first $k$ steps. Transported, it is not: the map
(39) reweights *within* each node, but a node standing for a cluster of
distinct $\sigma^2$ transports as a single point, and the reweighting of
that cluster is only as accurate as the lumping. Moreover the measure is
known only where the Krylov space explored it — mass on directions the
starting vector did not excite, or that $k$ steps did not reach, is absent
from (37) altogether and cannot be recovered by any map. Both effects
concern spectrum at the low end, which is where descending $\alpha$ moves
weight; hence the transported estimate of $R_k(\alpha')$ for
$\alpha' \ll \alpha$ is systematically optimistic, and the levels (22)
computed from it are lower bounds on the true $\alpha_k$. This is a
statement about the representation, not about any particular
implementation: a $k$-node quadrature determines the measure's first $2k$
moments and nothing more.

## 14. What a solve determines

The endpoints and their families are expressed in four quantities:
$\ell_p$, $\ell_d$, the constants $c_k$, and $\lambda_Z$. They are not
equally available.

**Determined, at no cost.** $\ell_p = \alpha\|r_0\|$ and
$\ell_d = \|s_0\|/\alpha$ are read from residuals the algorithm forms
regardless — one as the Krylov iteration's right-hand side, before any work;
one as the refinement loop's own convergence test, after it. They fix
$\alpha_0$ and $\alpha_c$, the lowest tread of each staircase, exactly and
without calibration. The two differ in kind, however: by Proposition 1
$\ell_p$ is an exact quantity, while by Proposition 2 $\ell_d$ is a single
realization of a rounding process. Successive passes furnish further
readings of the latter — normalized by the correction each produced, since
the level is proportional to it — where the former admits no such
repetition. They are readings of *related* quantities and not repeated draws
of one, each pass solving a smaller system than the last; §17 records what
that costs, which is that they neither sharpen the ceiling nor calibrate its
scatter.

**Determined by the iteration, up to its own depth.** By (21) the constants
$c_k$ are products of the Jacobi coefficients of the limiting measure, and
by §13 the Golub–Kahan recursion assembles exactly those. Iteration $i$
produces both $a_i$ and $b_i$, hence the off-diagonal
$e_i = \alpha\sqrt{b_i}/a_i$, so a solve performing $k$ iterations yields
$c_1,\dots,c_k$ — one per iteration, not one fewer, the $k$-th belonging to
$T_{k+1}$ rather than to $T_k$. Since the solve already witnesses that $k$
suffices, the informative ones are $c_1,\dots,c_{k-1}$: **the levels lying
above where it ran**, which are the ones its Krylov space explored. It says
nothing about levels below, except through the transport of §13, which is
one-signed. A single iteration therefore fixes $\alpha_1$, and only
$n = 0$ leaves the floor family unknown.

**Not determined by the iteration at all.** $\lambda_Z$ is an eigenvalue on
$\ker B$, and the equation $B\,\delta x = r_0$ that the iteration solves is
precisely the one that says nothing about that subspace. No amount of Krylov
work reveals it. But $\alpha_{\mathrm{wall}}$ need not be assembled from
$\lambda_Z$: a single observed pass gives the ratio
$\varrho_d = \|s_{j+1}\|/\|s_j\|$ and hence, by (27),

$$\alpha_{\mathrm{wall}} = \frac{\alpha}{\varrho_d},$$

the whole upper family from one pass. Two qualifications attach. Inequality
(26) is a worst case, attained only when the residual aligns with the
smallest eigendirection of $Z^\top A Z$, so the resulting
$\alpha_{\mathrm{wall}}$ errs conservative. And the observed ratio equals
$\varrho_d$ only where the dual mechanism dominates the augmentation's own
contraction, which is to say at large $\alpha$; below that the ratio
measures the primal side instead and carries no ceiling information.

**The duality.** Each family is observable exactly where it binds. The floor
levels are revealed by a solve that grinds, the ceiling levels by a solve
that passes, and each is silent about the other. At the cost minimum of §11
a solve does neither — both $n$ and $j$ stand near their least values — so
**the cheapest solve is the one that says least about why it is cheapest.**
The blind case of §13 is the primal half of this statement; the requirement
of an observed pass is the dual half. Away from the minimum, in either
direction, the solve pays for information about the side it has moved
toward.

**A halt on insufficient contraction.** A rule that abandons refinement when
the residual fails to fall by a prescribed factor per pass is, by (27), a
test of $\alpha$ against a multiple of $\alpha_{\mathrm{wall}}$. The graded
family (28) is defined on $\alpha < \alpha_{\mathrm{wall}}$, where by
construction the passes contract and no such halt occurs; the two are
consistent, and a halt firing at $\alpha$ well below
$\alpha_{\mathrm{wall}}$ reflects a threshold stricter than the analysis
requires rather than a phenomenon the analysis omits.


# Part IV — the estimate

The sections to this point derive what the two residuals mean and how they
move with $\alpha$. What follows is the procedure that consumes them: given
one completed solve, which $\alpha$ should have been used.

The two endpoints are the same objects throughout, under one set of names:
$\alpha_0$ is the primal endpoint of §5, and is also the $k = 0$ member of
the graded family $\alpha_k$ of §9; $\alpha_c$ is the dual endpoint of §7,
and the $j = 1$ member of $\alpha_c^{(j)}$ of §10. Part IV states them in
the form the algorithm uses — with the transport exponents of §16 in place
of the fixed unit exponents assumed in Parts I and II, which are their
$\theta = \phi = 1$ case.

## 15. The answer is an interval

That the answer is a range rather than a point is not new. Golub and Greif
conclude that there appears to be a band of $\gamma$ over which at least two
of the three condition numbers they track, and possibly all three, sit near
their minima [GG, Sec. 5, p. 2090]; the observation is later given a theorem,
$\kappa(A + wbb^\top)$ being convex in $w$ with at most one minimum
[GGV, Thm. 3.4, p. 786], and an exact plateau in the rank-one case
[GGV, Prop. 3.1, p. 785; and §16 below]. What follows sharpens the
observation into a computable interval with endpoints and a cost.

The sharpening is not an analogue of the convexity result, and the difference
is worth marking, since it is where §20's warning comes from — though the
parallel runs closer than it first appears, and then breaks at a definite
step. Their object is a condition number, smooth and convex in the
parameter, so a unique minimum follows and the plateau is an artifact of the
max–min form. This note's smooth interpolant is *also* convex in
$\log\alpha$ — §20 shows both its terms are, with positive second
derivative on the interval between the poles — and so it too has a unique
minimum, available in closed form. What differs is that the object of
interest is not that interpolant but its quantization, and quantization
costs more than convexity: it costs unimodality. §20 gives the mechanism.
The flatness here is therefore the phenomenon and not a degeneracy of it,
and it is a flatness that may occupy several separated intervals rather than
one.

Cost is flat across the window. Every $\alpha$ at which the solve does no
Krylov work and fires no refinement pass costs exactly one backsolve, and
these form an interval; no point within it is preferable to any other on
cost grounds. The estimate is therefore that interval and not a point
inside it. Placement — where in the interval to sit, and what margin to
keep from its edges — is a separate question, decided by how the interval
moves between solves rather than by the cost within it, and it is out of
scope.

The interval is also *this* solve's and *this* tolerance's. A computation
performing several solves at a common $\alpha$ must intersect the intervals
obtained for each, as §1 notes and this note does not treat; and $\varepsilon$
enters as given, so how the interval moves when the tolerance tightens is a
property of the schedule rather than of the solve. Both are constraints on
what a caller may do with a return, not qualifications on the return itself.

**One component of the optimal set, chosen deliberately.** The
cost-minimizing $\alpha$ need not form a single interval: by §20 the two
staircases round separately, so the minimizing set is in general a union of
intervals. The components are not arbitrary. They are indexed by how the
same budget is split — one is $(k, j') = (2,1)$, the next $(1,2)$, the next
$(0,3)$, all costing $3$ — and since $\alpha_k$ falls with $k$ while
$\alpha_c^{(j)}$ rises with $j$, they are ordered in $\alpha$ by that split:
iteration-heavy at the bottom, pass-heavy at the top.

The estimate returns the **lowest**, by the tie-break of §21 step 7. That is
the safe choice and not merely a convention: by §17 an $\alpha$ below the
optimum costs Krylov iterations, which the next solve reports at once, while
one above risks a refinement stall, which is a failed step. Spending the
budget on iterations rather than passes buys that asymmetry for nothing,
the total being equal across components by construction.

The return is therefore sound and deliberate, but not exhaustive: an equally
cheap $\alpha$ exists above, separated by a band on which the cost is
higher, and no comparison against the returned interval can conclude that a
parameter outside it is worse.

Reporting an interval also degrades in the right way. Every quantity the
solve fails to determine costs the estimate coverage rather than
correctness: it shrinks the reported interval, or blanks an endpoint
outright, but never displaces it. A solve that reveals little therefore
returns a smaller or emptier answer rather than a confident wrong one, and
every point it does report still delivers the quoted cost. The direction
guard of §18 is the clearest instance — a level it cannot vouch for is raised
to the probe, narrowing the interval — and the containment is what the
return contract below requires.

**The return value** is a triple $(\mathrm{lo}, \mathrm{hi}, \mathrm{cost})$:
the interval of parameters, and the budget in solve-pairs that the interval
assumes. Either endpoint may instead be $\star$; $\mathrm{cost}$ never is.

    (lo, hi, cost)   with lo <= hi
                     an interval was found; solving anywhere in it costs
                     `cost` solve-pairs

    any entry = *    that quantity was not established by this solve

There is no third outcome. The estimate takes a **converged** solve, and a
converged solve is a witness (below) that some budget suffices, so no input
it accepts can be reported infeasible. An inverted return, $\mathrm{lo} >
\mathrm{hi}$, occurs only as the withheld half-answer of step 6, where one
endpoint is a bound and the comparison is undecided rather than decided
against.

**The completed solve is a witness, and bounds the answer.** The estimate
takes converged solves only. The input ran at $\alpha$, spent $n$ iterations
and $j$ backsolves, and reached the tolerance. That is a proof — not an estimate — that the one-point interval
$[\alpha,\alpha]$ is delivered at budget $n+j$. The estimate therefore never
returns a budget exceeding $n+j$, and has no infeasible outcome at all:
the witness is carried as a floor on the return and seeds
the search, so the worst available answer is an exact one-point interval at
a known cost rather than a refusal. Equality, $\mathrm{cost} = n+j$, is the
statement that given what this solve reveals the probe was already at a
local best.

Two consequences are worth naming because they are otherwise invisible.
Every level with $k \ge n$ clamps to $\alpha$ under the guard of §18, so the
search never selects $k > n$: the estimate can only ever recommend *fewer*
iterations than the probe spent. And since $(n, j)$ is feasible whenever the
levels are correct — $\alpha_n \le \alpha$ by the definition of $n(\alpha)$,
and $\alpha \le \alpha_c^{(j)}$ by the definition of $j(\alpha)$ — any
returned budget above $n+j$ is evidence of a defective level, not of a
costly problem. That is the direction limit of §14 in operational form.

Each endpoint carries, besides its value, one of four qualities:

    x                exact
    x^-              a lower bound: the truth is at least x
    x^+              an upper bound: the truth is at most x
    *                not established

The signs are not decoration. Every compromised quantity in this note fails
in a *known direction*, and a bound with a direction remains usable where an
unsigned warning would only be discardable:

    alpha_c, off-ramp -- hence below the switch,
      and so below the ceiling (step 2)             lower  ^-
    alpha_wall from a pass failing the screen (5)   lower  ^-
    a floor level raised by the guard (4)           upper  ^+
    a level extrapolated past the harvest (4)       upper  ^+

Two consequences follow. Comparisons are permitted only when they are
decided by the signs. $\mathrm{lo} \le \mathrm{hi}$ survives a
$\mathrm{hi}^-$, since the true ceiling is higher still, and equally a
$\mathrm{lo}^+$, since the true floor is lower still; in both cases the
reported interval is a subset of the true one. $\mathrm{lo} > \mathrm{hi}$
survives neither — which is exactly the asymmetry the algorithm turns on at
step 6, here a property of the values rather than a remark. And exhaustion,
being the claim that two known endpoints cannot meet, would require both
entries exact — and is in any case unavailable, the precondition of §21
ruling it out.

`cost` is neither signed nor nullable. It is a count rather than a
transported quantity, it is always a number, and it never exceeds the
witnessed $n+j$. Where an interval was named it is the budget that interval
assumes; where only a half-answer was, it falls back to the witness and
reads as *no interval was priced, but $n+j$ at the probe is established*.
There is no return on which it is $\star$: the precondition makes $n+j$
known whatever the endpoints do, so a solve that determines nothing about
where to go still reports what it cost to get where it is. Only the two
endpoints are nullable.

Three conventions make this unambiguous.

$\star$ means **not established**, never *unbounded*. A caller must not read
a $\star$ floor as a licence to descend: in the bias-floor regime of §16 the
floor is invalid rather than absent, and there is a real limit the estimate
simply cannot locate.

One asymmetry is worth stating because the algorithm turns on it: the
ceiling's meter can only cause a *false empty*, never a false nonempty
(§16), so a nonempty window is reported whatever the meter says, while an
apparently empty one is withheld when the meter is low.

$\mathrm{lo} > \mathrm{hi}$ **with both endpoints exact** cannot occur. The
ordinary empty-window case — $\alpha_0 > \alpha_c$, which is common — is not
returned; it is the internal trigger to spend budget, and what comes back is
the graded interval, or failing that the witness, either of which satisfies
$\mathrm{lo} \le \mathrm{hi}$. The only inverted return is step 6's withheld
$(\alpha_0, \alpha_c^-, n+j)$, in which the ceiling is a bound and the
comparison is undecided.

So an inversion always carries a $\star$ quality, and a caller reading the
values alone will misread it as a verdict on the configuration when it is a
statement about the reading. The qualities are the rule, and `cost` does not
help distinguish the cases: it is a number on every return.

The $\star$ is localized: a reading that invalidates the floor leaves the
ceiling intact, and conversely. A half-answer is returned rather than a
refusal, since one usable endpoint is strictly better than none:

    (*, alpha_c, n+j)        floor not established — the bias-floor regime
                             of §16, where alpha_0 is invalid rather than
                             merely bounded.  This is now the only route to
                             a blanked floor: an absent harvest still leaves
                             alpha_0, and a clamped level is returned as
                             alpha_k^+ rather than discarded (step 7)
    (alpha_k^+, alpha_c, k+1) floor level clamped by the guard: the true
                             level lies at or below it
    (alpha_0, alpha_c^-, n+j) ceiling off its ramp: understated, and an
                             apparent emptiness withheld (step 6)
    (alpha_0, *, n+j)        ceiling not established — an unshifted reading
                             cannot be transported across a shift engagement,
                             so step 2 has no value to transport and blanks
                             the endpoint. The budget is witnessed even so
    (*, *, n+j)              neither — e.g. a degenerate tolerance

## 16. Measurement

The solve records $\alpha$, the tolerance $\varepsilon$, its initial
constraint residual $\|r_0\|$, and the dual component $\|s_0\|$ of the
two-block residual evaluated after it. Both residuals are formed by tests
the algorithm performs regardless: the first is the Krylov iteration's
right-hand side, checked before any work; the second is what the refinement
loop examines to decide whether to run. Neither costs anything. From them,

    alpha_0 = alpha * ( ||r0|| / eps )             # floor:  no Krylov work
    alpha_c = alpha * ( eps / ||s0|| )^(1/phi)      # ceiling: no refinement pass

with one transport exponent and one meter,

    phi   =   d log||s0|| / d log alpha      # an input, default 1
    theta = - d log||r0|| / d log alpha      # measured, = 1/a_1  (Sec. 5);
                                             # read, not applied -- see below

recovering the plain forms at $\theta = \phi = 1$, and with $a_1$ the first
CG step length on $S = BF^{-1}B^\top$ (§1) — the identity $\theta = 1/a_1$
holds on that operator and on no rescaling of it. Both leading constants
are $1$: they are the two convergence tests read literally, not calibrated
quantities, and the norms are whichever the tests themselves use — the
transport is exact in any norm because the underlying relations are
identities on vectors, not bounds. What is not free is the pair of
exponents, and the two differ in what it costs to know them.

**The floor's exponent is measured, and is not applied.** $\theta \le 1$
always — it is $\langle\mu\rangle_w$ and $\sigma(S) \subset (0,1)$ by §1, so
the bound is a consequence of the normalization rather than a separate
assumption — so the residual falls *more slowly* than assumed. Under a
locally constant elasticity the transport would be exact at

    alpha_0 = alpha * ( ||r0|| / eps )^(1/theta)

and $\theta$ is free whenever the solve took any Krylov step, being the
reciprocal of the first CG step length on $S$ (§1, §19). Nevertheless **the
plain form is the default and the corrected form is not**, for a reason that
is a property of $\theta$ rather than of the arithmetic.

The exponent is not constant on the transport interval. $\theta$ rises
towards $1$ as $\alpha$ approaches the floor, since that is where the
augmentation resolves the directions carrying the residual's weight; so the
probe's $\theta$ is a *lower bound* on the interval's, and $1/\theta$ inflates
the entire logarithmic distance by a factor the interval does not warrant.
The damage is exactly quantified, and needs no threshold to state. The
correction displaces the endpoint by

$$\Big(\frac{1}{\theta} - 1\Big)\log\frac{\|r_0\|}{\varepsilon},$$

which is to be judged against the floor family's tread,
$\Lambda/((k+1)(k+2))$, by the dominance argument of §17. That comparison
settles the matter without a constant, and it settles it against applying
the correction at all. Below a tread, the displacement cannot move the
answer — the reported interval is the same either way, so the correction
buys nothing. Above a tread it does move the answer, and by the paragraph
above it moves it the wrong way, since $1/\theta$ over-inflates a distance
whose true exponent is larger than the probe's. There is therefore no band
in which the correction both matters and helps: it is irrelevant where it is
safe, and harmful where it is not.

The rule that follows is the note's rule elsewhere. $\hat\phi$ meters the
ceiling and does not correct it (below); §19 admits a bought $\theta$ as a
validity check and refuses it as an exponent; §21 declines to buy $c_1$
because a quantity adequate for judging is not thereby adequate for moving.
$\theta$ on the floor is the same case, and the same answer: **report it, do
not raise to it.** Use $1/\theta$ only as the upper arm of the bracket below,
and only where that arm's width is below a tread — which is the same
comparison again, and the only gate the note needs.

*The direction of the plain form's error, and why only one direction
arises.* Since $1/\theta \ge 1$, the exponent raises
$(\|r_0\|/\varepsilon)^{1/\theta}$ above $\|r_0\|/\varepsilon$ when that
ratio exceeds $1$ and lowers it when it does not: the plain form understates
the floor for a probe below the window and *overstates* it for one inside.
Only the first case arises. By (T1) the iteration is entered only when
$\|r_0\| > \varepsilon$, and $\theta = 1/a_1$ exists only if a step was
taken; so wherever the correction is available the ratio exceeds $1$, and
wherever the ratio does not, the correction is unavailable and the plain
form is used unmodified. The one way to reach the reversed case is to *buy*
$\theta$ at $n = 0$, which §19 forbids for this reason. The bias is one-signed on the branch the algorithm
can reach, and the closure is exact rather than probabilistic — it is (T1)
that supplies it, not an assumption about where probes fall.

**What this changes in the assumptions.** The floor no longer requires
$\|r_0\| \propto 1/\alpha$; it requires only that the *exponent* be
constant between the probe and the endpoint, and it measures that exponent
rather than positing it. The old form is the special case $\theta = 1$, so
this is a strict weakening — and a material one, since $\theta$ is not close
to $1$ wherever the tail is populated. The assumption is exchanged, not
removed: over a short transport a locally constant exponent is a mild
hypothesis, over a long one it is not, and $\theta$ traverses a wide range
across the operating interval. That is the same statement as the correction
helping most for a nearby probe and least for a distant one.

**A conditional bracket, with a gate.** Where $\theta$ increases with
$\alpha$ — nodes and weights both moving, the nodes dominating — the two forms
bracket the true floor when transporting upward: the plain form understates
it and the corrected form, using an exponent smaller than the interval's true
one, overstates it. The monotonicity is not a theorem — the weights can
redistribute enough to reverse the mean — so treat the pair as a usual
bracket to be reported, not as bounds to be relied on.

The gate is the tread, not a value of $\theta$. The bracket's width is the
displacement above, $(1/\theta - 1)\log(\|r_0\|/\varepsilon)$; report the pair
where that is below $\Lambda/((k+1)(k+2))$ and the plain value alone where it
is not, since a bracket wider than a tread reports nothing but its own
inapplicability. The rule is per-problem — it moves with $\Lambda$, with $k$,
and with how far the probe sits from the endpoint — and carries no fitted
number; $\theta$ is carried as the reading that enters it.

**The ceiling's exponent is an input.** No single solve measures $\phi$: it
is a slope in $\alpha$, and one reading gives a point. Two readings at
different $\alpha$ determine it, which the surrounding computation
ordinarily has — the same history that supplies the ramp check below. The
estimate therefore takes $\phi$ as a parameter, defaulting to $1$.

The default is not arbitrary, and it comes with a meter. By the conditioning
argument $\|s_0\| \approx c_u u\,\alpha\|F\|\,\|x_0\|$, so

$$\phi = 1 + \frac{d\log\|F\|}{d\log\alpha} + \frac{d\log\|x_0\|}{d\log\alpha},$$

whose last term vanishes as $x_0$ settles. It is worth saying why the
exponent is not $2$. Golub and Greif show $\kappa_2(K(\gamma))/\gamma^2 \to
\|B\|^2$ [GG, Cor. 2.4, p. 2080], so the conditioning of the whole system
grows quadratically, and a reader who takes $s_0$ to scale with it will
expect $\varphi = 2$. But backward stability bounds the *residual* by
$\|\Delta F\| \le c_u u\|F\|$, in which the norm appears and the
condition number does not; $\kappa$ enters the error instead. The solve
factors $\alpha F$, so $\|F\|$ is what the rounding rides on, and the
exponent is $1$. (Their Prop. 2.1 [GG, p. 2079], that the augmentation
perturbs $K^{-1}$ by exactly $-\operatorname{diag}(0, W)$, is the exact
statement about the $\alpha$-dependence of the dual block from which
Cor. 2.4 follows, and the same identity is what [BO, Lem. 4.1, p. 2101]
rests on, cited there to [GG, Prop. 2.1, p. 2079].) The linear growth is
observed independently: sweeping $\gamma$ over four decades on a large
interior-point instance, [PS, Sec. 6.1.2, p. 28] report the relative residual
norm rising linearly with the parameter while the Krylov count falls — the
ceiling and the floor, traded against one another. The same trade is visible
earlier and on both arms of a single sweep: over eight decades of $\gamma$ on
a power-flow sequence, [RC, Fig. 2(a), p. 12] show the CG count falling from
some $10^3$ to a mean of $9.4$, while [RC, Fig. 3, p. 12] show the backward
error and relative residual of the original system rising across the same
range, the authors concluding that only an intermediate band of four decades
is usable [RC, Sec. 5.1, p. 13]. That band is the flat-bottomed bowl of §11
observed rather than derived. And $\|F\| = \|\beta A + B^\top B\|$
interpolates between its two terms, giving

    phi_hat = T_B / (T_A + T_B),     T_A = beta * ||A||,  T_B = ||B^T B||

— one when $B^\top B$ dominates, zero when $\beta A$ does, and one half at
the switch $\alpha_\star = \lambda_{\max}(A)/\sigma^2_{\max}$. That switch is
the classical parameter choice: for positive semidefinite $A$ it is exactly
$\|A\|/\|B\|^2$, which Golub and Greif propose on the grounds that it
equalizes $\|A\|$ with the norm of the added term [GG, Sec. 2, p. 2079], and
which they find lies inside the good range in both of their numerical
examples [GG, Sec. 4.1, p. 2087; Sec. 4.2, p. 2089]. In the interior-point
setting the same choice acquires a rate: with $\Xi$ the duality measure,
[PS, Assumption 3(b), p. 20] require $\gamma = \Theta(\Xi^{-1})$ and identify
that requirement with Golub and Greif's norm balance
[PS, Remark 2, p. 13] — so the switch tracks the barrier parameter, rising by
a decade for each decade the tolerance falls.

**The switch is the top of a second window.** The norm balance is not merely
a good point but an endpoint of an interval, and the interval is derived
rather than reported. For $A$ positive semidefinite of nullity one and the
rank-one augmentation along its null vector,

$$\kappa_2\big(A + w u_n u_n^\top\big)
  = \frac{\max(\lambda_{\max}(A),\, w)}{\min(\lambda^+_{\min}(A),\, w)}$$

exactly — strictly decreasing below $\lambda^+_{\min}(A)$, **constant on**
$[\lambda^+_{\min}(A),\, \lambda_{\max}(A)]$, strictly increasing above
[GGV, Prop. 3.1, p. 785]. Carrying the same strategy to the $B^\top B$
augmentation, as [GGV, §3, p. 787] do in proposing $w = \|A\|/\|B\|^2$,
divides that plateau by $\sigma^2_{\max}$, and its upper endpoint is then
$\alpha_\star$ exactly. (The plateau is proved for the rank-one case; its
transfer to $B^\top B$ is by analogy there, not by theorem.)

So there are two windows over the same parameter, and they differ in
everything but shape. The conditioning plateau
$[\lambda^+_{\min}(A)/\sigma^2_{\max},\ \alpha_\star]$ is computable a
priori from the spectra, spans $\kappa$ of $A$ on its range, and is flat
because $\kappa(F)$ is; the window $W_0$ of §8 is measured from one solve,
spans $\varepsilon^2/\ell_p\ell_d$, and is flat because cost is. They share
the endpoint $\alpha_\star$, and the sharing is not coincidence: it is where
$\beta A$ and $B^\top B$ exchange dominance in $F$, which is simultaneously
where $\kappa(F)$ starts to climb and where $\|F\|$ — hence the dual
residual, hence the ceiling — starts to grow. The two analyses read the same
crossing through different tests. This note's contribution is not the
crossing but the second endpoint: the conditioning account gives no lower
edge tied to a tolerance, and no cost.

Here $\alpha_\star$ arrives from a third direction,
as the point at which the ceiling begins to bind; of which
$\hat\phi = \tfrac12$ is the smooth restatement. Both norms may be taken
from diagonals, $\max_i A_{ii}$ and $\max_i (B^\top B)_{ii}$, which are free:
the second is assembled with the operator and the first is the Hessian
diagonal.

**This is the ceiling's validity meter**, the counterpart of $\theta$ on the
floor, and it should be read before the ceiling is used. It is not, however,
a threshold on $\hat\phi$: as on the floor, what decides is the displacement
it induces measured against a tread. Nor is $\hat\phi$ to be used as an
exponent, and the reason is the mirror of §19's.

**Why $\hat\phi$ is read and not applied.** Since $\hat\phi \le 1$, the
exponent $1/\hat\phi \ge 1$; and at $j = 1$ the dual test passed, so
$\varepsilon/\|s_0\| \ge 1$ and raising it to that power *increases*
$\alpha_c$. Correcting therefore moves the ceiling up — the
anti-conservative direction, the one that risks a stall rather than an
iteration. It is the same argument as $\theta$ at $n = 0$ turned through the
other endpoint: in both cases the correction is available exactly on the
branch where the ratio lies on the wrong side of $1$, and in both the answer
is to report the meter and not to raise to it. The two refusals are one rule
stated twice, and the rule is that a quantity adequate for judging a reading
is not thereby adequate for moving it. Transporting with the default exponent
against a true exponent $\varphi$ displaces the ceiling by exactly

$$\log\alpha_c^{\text{used}} - \log\alpha_c^{\text{true}}
  = \Big(1 - \tfrac{1}{\varphi}\Big)\,L,
  \qquad L = \log\frac{\varepsilon}{\|s_0\|},$$

to be judged against the ceiling family's tread $M/(j(j+1))$: the reading is
sound when the displacement falls below it and marked when it does not.
Taking $\varphi \approx \hat\phi$, that is a per-problem test with no fitted
number, moving with $M$, with $j$, and with the distance $L$ being
transported. It replaces any threshold on $\hat\phi$ itself, which would be
a condition on the wrong quantity — a low $\hat\phi$ is harmless over a short
transport and ruinous over a long one, and only the product distinguishes
them.

which grows like $L/\varphi$ as $\varphi \to 0$: decades rather than tenths
whenever $L$ is itself order-decades.

**The magnitude is $\hat\phi$'s business; the direction is the probe's.**
Off the ramp $\varphi < 1$, so the first factor is negative and the sign of
the displacement is $-\operatorname{sign}(L)$: it is fixed by whether the
dual test passed at the probe, which is the same bit as whether a refinement
pass fired. A probe that sat *below* the ceiling ($j = 1$, $L \ge 0$)
understates it; one that sat *above* it ($j \ge 2$, $L < 0$) overstates it.
The determination is exact and needs no reference to $\hat\phi$.

**The switch never lies above the ceiling.** Only the first of those two
cases can arise, and the reason is a fact about the two scales worth stating
on its own. Below the switch $\|F\| \approx \beta\|A\|$, so
$\|s_0\| \lesssim c_u u \lambda_{\max}(A)\|x\|$ is flat in $\alpha$. This
is the lower half of the mechanism of §7 (C3, (10b)) — derived from
the exactness of the dual row and a backward-error argument, not assumed —
whose two halves are flat below the switch and linear above it. The dual
test therefore returns the same verdict at every $\alpha$ below the switch,
and so

$$\alpha_c < \alpha_\star \iff \varepsilon < c_u u \lambda_{\max}(A)\|x\|
  \iff \text{the dual test fails at every } \alpha .$$

The right-hand condition is not a low ceiling but the absence of one: no
$\alpha$ admits a refinement-free solve, and the tolerance is unattainable.
Only the forward direction is needed below, and it uses (10b) as the upper
bound it is; the converse additionally requires that the bound be attained,
which (C3) gives to within a factor of order one.
So in any problem with a ceiling at all, $\alpha_\star \le \alpha_c$. Two
consequences. An off-ramp reading has $\alpha < \alpha_\star \le \alpha_c$,
hence $j = 1$ and $L \ge 0$, so **the degradation is one-signed —
understatement — and the $j \ge 2$ branch above is unreachable**; the failure
is excessive conservatism rather than a ceiling that cannot be met. And the
ramp always overlaps the operating range, so $\hat\phi$ is a gate that
admits something rather than a screen that might exclude everything.

Prefer moving
the probe to correcting the exponent; a supplied $\phi$ repairs the
arithmetic but not the fact that the sensor is reading its own noise floor.

The single dimensionless summary is

    X = alpha_0 / alpha_c = l_p * l_d / eps^2

the exhaustion ratio. $X \le 1$ is exactly the condition that the window is
nonempty, so no separate test of regime is needed: the wide-window and
narrow-window cases are one formula returning a wide or narrow interval,
and $X > 1$ is the empty case of §18. The tension $X$ measures is the one
Golub and Greif describe qualitatively — $\gamma$ large enough to clear the
ill-conditioning of $A$, but not so large that the rank-deficient added term
dominates [GG, Sec. 2.3, p. 2083]. Benzi and Olshanskii state the upper side in
the same terms, that the $(1,1)$ block becomes increasingly ill conditioned
for large $\gamma$ because the added term has a nontrivial null space, so
the parameter is better kept moderate [BO, p. 2099]. $\alpha_0$ and
$\alpha_c$ are that statement given two endpoints.

### Regularization

If the factorization applied a shift $\rho > 0$, the recorded residuals are
not the quantities the transports assume, and the two ends are affected
differently.

*Ceiling — no change to the formula.* The shift contributes an exact term
$\alpha\rho x$ to the dual residual. It is proportional to $\alpha$, as the
rounding term is, so $\alpha_c = \alpha\varepsilon/\|s_0\|$ still transports
and returns a smaller and correct ceiling. The requirement is that the
reading and the endpoint lie on the same side of the shift's engagement: a
reading taken with no shift cannot be extrapolated to a regime where the
ladder has engaged, because the level itself jumps. When the shift is known
in advance, the ceiling may be predicted before $\alpha$ is chosen, from
$\rho$, the tolerance, the test's own scaling and an estimate of $\|x\|$ —
the one endpoint that need not be measured at all.

*Floor — the formula changes.* The shift contributes to $\|r_0\|$ a bias
that is constant in $\alpha$, and a constant cannot be transported along a
law that assumes inverse proportionality. Applying the uncorrected formula
to a reading in which the bias dominates over-predicts the floor by however
far the reading sits above the ramp.

Two remedies, in order of cost.

    screen (free):    if  sqrt(rho) * ||x|| / 2 << eps  the bias cannot
                      matter; use alpha_0 unmodified

    correct (one     bias <- rho * B F^{-1} x       # one solve-pair
    solve-pair):     r0   <- r0 - bias              # vector subtraction
                     alpha_0 <- alpha * (||r0||/eps)

A third check is free when the history holds two readings at different
$\alpha$: on the ramp the product $\alpha\|r_0\|$ is invariant, on the bias
floor $\|r_0\|$ itself is. Which invariant holds identifies the regime; how
well separated the two are in a given problem is not established here. A
floor-regime reading should be refused rather than transported.

That a shift is available at all is a property of the augmentation: for
positive semidefinite $A$ with $K$ nonsingular, every $\nu > 0$ however small
makes $A + \nu B^\top B$ positive definite [GG, p. 2084]. The alternative is
to remove the singularity by reducing the system rather than shifting it, by
deflating the $(1,1)$ block by exactly its nullity [GG, Sec. 3, pp. 2084–2086];
that route is not taken here, since it presumes the nullity is small and
known, and the estimate is meant to run on a solve that has already happened.

*The switch–ceiling ordering survives.* The shift contributes to $\|s_0\|$ a
term exactly linear in $\alpha$, so it joins the linear branch of (10b) and
leaves the flat branch alone: below the switch the level is still
$c_u u \lambda_{\max}(A)\|x\|$, independent of $\rho$. What moves is the
switch, down to $c_u u \lambda_{\max}(A)/(\rho + c_u u \sigma^2_{\max})$,
which is $\approx u\lambda_{\max}(A)/\rho$ once the shift dominates. The
equivalence of §16 therefore holds verbatim with the relocated switch — its
right-hand side is unchanged — and so do its consequences: an off-ramp
reading still understates, and $\alpha_c^+$ remains unreachable. A shift
widens the ramp rather than invalidating the argument.

The shift is invisible below $\rho_\star = c_u u \sigma^2_{\max}(B)$, where
it is dominated by the rounding term it adds to. A ladder whose first rung
sits far above $\rho_\star$ spends ceiling for nothing.

## 17. The two ends differ in kind

The floor rests on an exact identity and an assumption with a computable
meter. The ceiling rests on a residual that is *entirely rounding* — by the
exactness of the dual row it has no exact-arithmetic part at all — so a
single reading is one draw from a distribution rather than a measurement of
a constant.

$s_0$ is formed once, before any refinement pass has run, so the estimate has
one reading and no choice about it. Its scatter cannot be *calibrated* from
a single solve, and $\hat\phi$ sorts readings on the ramp from readings off
it without bounding the error of either.

The qualification matters because at $j \ge 2$ the passes do supply further
readings, as §14 says — and $j \ge 2$ is the regime in which the ceiling
binds. What they do not supply is *repeated draws*. Each pass solves a
smaller system, so its rounding is drawn from a different population;
normalizing by the correction each produced makes the readings comparable in
scale but not in kind. Two consequences follow, and they are the same
consequence twice. Aggregating them sharpens nothing — it introduces a bias
toward the unsafe side, which is why the single reading is taken above — and
their spread is not the dispersion either: it tracks it, but with a scale
factor not bounded here. So the passes give a correlated proxy of unknown
calibration, which is more than nothing and much less than a measurement.
That is the sense in which the scatter has no calibrated route from a single
solve.

What removes the difficulty is not precision but **dominance**. Cost is flat
within a tread, so an error in $\alpha_c$ is free unless it crosses a tread
boundary.

Both treads follow from the level laws. Consecutive ceiling levels sit at
$(1/j - 1/(j+1))M$ and consecutive floor levels at
$(1/(k+1) - 1/(k+2))\Lambda$, so with
$M = \log(\alpha_{\mathrm{wall}}/\alpha_c)$ and
$\Lambda = \log(\alpha_0/\alpha_\infty)$ the scatter is irrelevant while

$$j(j+1) \ll \frac{M}{\operatorname{sd}(\log\alpha_c)}, \qquad
  (k+1)(k+2) \ll \frac{\Lambda}{\operatorname{sd}(\log\alpha_c)} ,$$

that is while $j \lesssim \sqrt{M/\mathrm{sd}}$ and
$k \lesssim \sqrt{\Lambda/\mathrm{sd}}$. These are per-problem quantities,
not constants, and only half of each is recorded: $M$ and $\Lambda$ follow
from the poles and endpoints of §22, while $\operatorname{sd}$ has no route
from a single solve at all. So the conditions are not evaluated. They are
used in the direction that needs no number — both are quadratic in the
budget, so a dispersion merely *small* against the spans suffices at small
budget, and the closure asks nothing sharper. What follows supports that
smallness; it does not measure it.

**Two inputs, not one.** The closure needs the dispersion small *and* the
spans not degenerate, and these are different quantities. A perfectly spread
perturbation still fails the ceiling condition at $j = 2$ if $M$ is a
fraction of a decade. But $M = \log(\lambda_Z\|x\|/\varepsilon)$ collapses
exactly when $\lambda_Z$ does, which is the near-singular regime §18 already
names and $\alpha_{\mathrm{wall}}$ already records — so the second input is
one the estimate has in hand, and the closure's failure mode is a nearly
singular reduced Hessian rather than any misbehaviour of the arithmetic.

Two features make the rest robust. Both conditions are quadratic in the
budget, so the slack at small budget is large; and both spans widen as the
tolerance tightens — $\alpha_0$ rises and $\alpha_c$ falls at the same rate
in $\varepsilon$, so $\Lambda$ and $M$ grow together and the two thresholds
move in step rather than one improving at the other's expense. Since $j$
stays small (§18), the ceiling side has slack to spare. The floor side is the
tighter, $k$ running higher than $j$, and it is the one that can bind — but
only deep in the empty-window branch at large $k$, which is precisely where
§23's item on the graded branch's excess already applies, and where adjacent
levels are in any case indistinguishable and no particular $k$ should be
trusted. The two known imprecisions coincide rather than compound.

**The dispersion itself has an account**, though a qualitative one, and it is
not needed in any sharper form. By Proposition 2 the dual row is exactly
$s_0 = \Delta\,x$ for the backsolve's backward perturbation, so $\|s_0\|^2$
is a sum of squares over the perturbation's rows and its logarithm
concentrates. What governs that concentration is not $n$ but the
**participation ratio** of the row scales,

$$n_{\mathrm{eff}} = \frac{\big(\sum_i \eta_i^2\big)^2}{\sum_i \eta_i^4}
  \;\in\; [1, n], \qquad
  \operatorname{sd}(\log\|s_0\|) \;\simeq\;
  \tfrac{1}{2}\sqrt{\operatorname{Var}(Z^2)}\; n_{\mathrm{eff}}^{-1/2},$$

the componentwise form of the Cholesky backward error putting
$\eta_i \propto (|L||L^\top||x|)_i$ — a fresh letter, $\sigma$ being reserved
throughout for the generalized singular values of $(B, A^{1/2})$ and $s_0$
for the dual residual — and $Z_i$ the components normalized by
those scales. The estimate's safety at small budget therefore rests on
assumption (D) of §24: that $n_{\mathrm{eff}}$ is large, so that the
dispersion is not of order one.

*The exponent is a second moment, not a central limit.* The relative
dispersion of a sum of squares $\sum_i\eta_i^2 Z_i^2$ is
$\sqrt{\sum\eta_i^4}\big/\sum\eta_i^2$ times
$\sqrt{\operatorname{Var}(Z^2)}$, which is $n_{\mathrm{eff}}^{-1/2}$ times a
constant; halving it passes to the logarithm. Nothing in that computation
asks the $Z_i$ to be Gaussian, and the rate is therefore not hostage to a
concentration picture $\|s_0\|$ does not satisfy. What Gaussianity would
supply is the constant, $\sqrt{\operatorname{Var}(Z^2)}/2 = 1/\sqrt2$, and
that is the part not to be relied on: the $Z_i$ are neither independent nor
identically shaped, and rounding errors generated by one factorization are
correlated. Correlation deflates the effective count below
$n_{\mathrm{eff}}$, so the estimate above is optimistic and the sign is
known — the true dispersion is at or above what $n_{\mathrm{eff}}$ predicts,
never below. That is the note's usual currency, and the rate is used only
in that direction: as an order-of-magnitude lower bound on the dispersion,
against treads whose spacing is measured in decades.

Note what (D) is *not*: it is not a condition on
the conditioning of $A$ or of $F$. Ill-conditioning collapses
$n_{\mathrm{eff}}$ only if it fails to be matched by a compensating spread in
the solution — and where the two are aligned, as they are wherever a
complementarity relation ties $\|x\|$ to the operator's small directions, the
row scales stay comparable however badly conditioned the operator is. That
alignment is the same structure that makes ill-conditioning benign for the
accuracy of the computed step.

$n_{\mathrm{eff}}$ is cheaply computable — $L$ and $x$ are both in hand, and
forming $|L|(|L^\top||x|)$ costs two triangular multiplies, the work of a
backsolve rather than of nothing — so (D)
is checkable per problem should it ever be in doubt — and checkable
absolutely rather than relative to the dimension, since what enters the
concentration is $n_{\mathrm{eff}}$ itself and not its ratio to $n$. It is
therefore not merely evidence but an estimate of the dispersion, one-signed
in the safe direction, available from a solve that has already happened. The
estimate does not compute it and nothing it returns depends on it, but a
caller in doubt may evaluate the conditions above rather than reason about
them. What remains open is narrow, and is stated in §23.

**Asymmetry of consequence.** Which direction the error takes is not
symmetric. Leaving the window downward costs Krylov iterations, which the
next solve reports immediately. Leaving it upward fires refinement, and a
refinement that starts too far above the attainable level does not merely
cost a pass — it stalls, and a stalled pass is a failed step. The estimate
reports $\alpha_c$ as measured; a margin below it belongs to the placement
policy, which is out of scope here, but the direction the margin should take
is settled by this asymmetry and not by the estimator's scatter.

## 18. When the window is empty

The quantities whose conditioning varies with the parameter are three — the
system, the $(1,1)$ block, and the Schur complement — and they move in
different directions as it grows [GG, Sec. 2.3, p. 2083, with the behaviour
plotted at Figs. 4.1, p. 2087 and 4.3, p. 2090]. The floor and the ceiling
below are a two-sided, measurable form of that picture.

$X > 1$ means no $\alpha$ admits a zero-work solve. Two purchases recover an
interval, and they are paid in the same currency:

- **Krylov iterations lower the floor.** Admitting $k$ of them gives
  $\alpha_k = (\alpha_0 c_k)^{1/(k+1)}$.
- **Backsolves raise the ceiling.** Admitting $j$ of them — an initial one
  and $j-1$ refinement passes — gives
  $\alpha_c^{(j)} = \alpha_c^{1/j}\,\alpha_{\mathrm{wall}}^{\,1-1/j}$.

The estimate is the interval $[\alpha_k, \alpha_c^{(j)}]$ for the pair
$(k,j)$ minimizing $k + j$ among those that make it nonempty, together with
that budget. Reporting $k$ alone would treat a refinement pass as forbidden
rather than as a unit of the same price, which it is not.

In practice $j$ stays small — the useful ceiling purchase is one extra
backsolve, occasionally two, while $k$ does the work. The search should be
written as a joint one but will terminate at low $j$; an elaborate ceiling
climb will not fire.

### The floor levels

$c_k$ is the $L^2$ norm of the degree-$k$ monic orthogonal polynomial of the
limiting spectral measure, and carries the units of $\alpha^k$: not a pure
number, and not slowly varying in $k$.

It needs no separate computation. The Lanczos recursion underlying CG
supplies the step lengths $a_i$ and the ratios $b_i$, and by §1

$$c_k = \prod_{i \le k} e_i, \qquad e_i = \frac{\alpha\sqrt{b_i}}{a_i},$$

accumulates as a running product with the $\alpha$-scaling already inside it.
Record the coefficients rather than the Ritz pairs: the pairs are downstream
of the coefficients, recovering $c_k$ from them is a round trip, and the
coefficients arrive incrementally so the sequence can stop as soon as it
drops below the ceiling.

Three limits on what the harvest supplies.

- *Depth.* Iteration $i$ produces both $a_i$ and $b_i$, hence the Jacobi
  off-diagonal $\sqrt{b_i}/a_i$; a solve running $k$ iterations therefore
  determines $c_1,\dots,c_k$ — one per iteration, not one fewer. Since the
  solve already witnesses that $k$ suffices, the informative ones are
  $c_1,\dots,c_{k-1}$: the levels lying **above** where it ran, which are
  the ones its Krylov space explored. It determines none below. A single
  iteration is therefore enough to fix $\alpha_1$, and only $n = 0$ leaves
  the floor family unknown.
- *Direction, and the guard.* Transporting the harvest upward is exact;
  downward it is systematically optimistic, since mass on directions the
  iteration never excited is absent from the representation and no map
  recovers it. **Guard: a level below the probe is not reportable.** Where
  the levels are located by searching — solving $\hat n(\alpha') = k$ for
  $\alpha'$ — the search must be clamped at $\alpha' \ge \alpha$, returning
  $\alpha$ itself if $\hat n(\alpha) \le k$ already; an unclamped search
  descends into the optimistic region, finds a spuriously low $\alpha'$ at
  which it believes $k$ iterations suffice, and recommends it. Where the
  levels come from the closed form instead, the same test applies after the
  fact: a returned $\alpha_k < \alpha$ is an extrapolation into unexplored
  spectrum and should be reported as a bound, not a value. The guard costs
  coverage rather than accuracy — levels genuinely below the probe become
  unreportable instead of wrong.
- *Reliability.* Later coefficients degrade as orthogonality is lost, worst
  on clustered spectra. At large $k$ a fresh direct measurement of $c_1$
  (§19) beats a drifting late coefficient.

Beyond the harvest, for $k$ small relative to the number of distinct nodes,
$\alpha_k \approx \alpha_\infty(\alpha_0/\alpha_\infty)^{1/(k+1)}$ with
$\alpha_\infty = \lim_k c_k^{1/k}$ the logarithmic capacity of the measure's
support. Substituting $c_1$ for $\alpha_\infty$ is the cheapest usable
approximation and errs by placing levels too high, so the policy believes a
level unavailable slightly before it is; the $(k+1)$-th root damps the
displacement.

### The ceiling levels

$\alpha_{\mathrm{wall}} = \lambda_Z \|x\| / \ell_d$ is the scale at which
refinement ceases to contract, with
$\lambda_Z = \lambda_{\min}(Z^\top A Z)$ on $Z = \ker B$ — the least
eigenvalue of the reduced Hessian, whose positive definiteness together with
full rank of $B$ is what makes $K$ nonsingular in the first place
[NW, Lem. 16.1; quoted in this setting at GG, p. 2077, where $Z$ spans the
null space of their $B^\top$, which is the present $\ker B$]. The ceiling's
accumulation scale therefore degenerates exactly as the system approaches
singularity. The same quantity bounds the parameter from *below*: for the
augmented block to be positive definite at all,
$\gamma$ must exceed a threshold inversely proportional to
$\lambda_{\min}(Z^\top K Z)$ and to the least nonzero singular value of the
constraint Jacobian [PS, Remark 2, p. 13, quoting Gould and Simoncini]; the
threshold is given directly as
$\gamma_{\min} = -\lambda_{\min}(A)/\lambda^*_{\min}(B^\top B)$ by
[RC, Thm. 1, p. 7]. It is $\le 0$ exactly when $A$ is positive semidefinite,
so §1's hypothesis is what makes this floor vacuous and leaves the parameter
constrained from below by cost alone — a consequence of that hypothesis, not
a property of the problem class, and the hypothesis is the one an
interior-point method loses as it proceeds [RC, Sec. 4.1, p. 7]. Where it
fails, $\gamma_{\min} > 0$ is a hard lower bound of a kind this note does not
produce, and it binds before any endpoint here does. That
is a floor on the parameter of a different kind from the one in §5 — imposed
by definiteness rather than by cost — and this note does not treat it. Since $\ell_d$ and
$\|x\|$ are already recorded, $\lambda_Z$ is the only unknown — and like
$c_1$ it depends on neither $\alpha$ nor $\varepsilon$, being a property of
the operator pair alone. It is therefore the ceiling's analogue of $c_1$:
measure it where the solve permits, carry it where it does not, and expect
it to move only as $A$ does.

**The primary route is an observed pass.** A refinement pass contracts the
dual residual by a factor that inverts to the wall:

    rho_d   = ||s_2|| / ||s_1||          # the FIRST ratio, not an average
    a_wall  = alpha / rho_d

Use $s_2/s_1$ and nothing further. The dual residual falls until it reaches
the rounding of the residual *evaluation*, after which the ratio tends to
one; a median or mean over later passes is dominated by that saturated tail
and is not a measurement of the contraction at all.

The reading inverts to the wall only where the dual mechanism is the one
being observed. The refinement contraction is the larger of two,

$$\varrho_{\mathrm{obs}} = \max\Big(\underbrace{\tfrac{1}{1+\alpha\sigma^2_{\min}}}_{\text{augmentation}},\ \ \underbrace{\alpha/\alpha_{\mathrm{wall}}}_{\text{dual}}\Big),$$

a shape falling as $1/\alpha$ on the first branch and rising in direct
proportion to $\alpha$ on the second, with a minimum where they cross. Only
the rising branch inverts. Equating the two gives the screen

    alpha^2 * sigma^2_min  >~  a_wall ,

with $\sigma^2_{\min}$ available from the harvest: a probe at or below the
crossing is reading the augmentation and its ratio carries no ceiling
information. The screen is necessary, and cheap; it is not sufficient, since
it cannot be evaluated when $\sigma^2_{\min}$ itself degenerates.

**The fallback is the factorization.** Where no pass ran, the factor already
in hand supplies the wall without one. On $\ker B$ the term $B^\top B$
vanishes, so $F$ reduces to $\beta A$ there and $\lambda_{\min}(F)$ is
$\lambda_Z/\alpha$ whenever the near-kernel of $F$ lies in $\ker B$ — the
situation $\lambda_Z$ describes. Hence

    lambda_Z  ~  alpha * lambda_min(F)
    a_wall    =  lambda_Z * ||x|| / l_d

and $\lambda_{\min}(F)$ needs no eigenvalue computation: the smallest pivot
bounds it, $\lambda_{\min}(F) \le \min_i L_{ii}^2$, each pivot being a
diagonal entry of a successive Schur complement arising in the factorization
of an SPD matrix, and so no smaller than
the whole matrix's least eigenvalue. The pivot is available wherever the
factor is; nothing is computed for it. Two riders: the bound's tightness
depends on the pivot ordering, only its direction being ordering-free; and
the identification $\lambda_Z \approx \alpha\lambda_{\min}(F)$ is an
approximation in its own right, so the total error is not one-signed.

**Do not average the two.** They are two derivations of one quantity rather
than two samples of it — they share $\ell_d$ and $\|x\|$, and their errors
move together. A mean lands between them, improving the weaker route and
degrading the better one. Where both are available, use the pass reading and
treat a large disagreement as evidence against it, since the screen above is
only necessary.

A third route, and not a marginal one: observing $j$ backsolves at $\alpha$
implies $\alpha \le \alpha_c^{(j)}$ by the definition of $j(\alpha)$, whence

$$\alpha_{\mathrm{wall}} \;\ge\; \big(\alpha^{\,j}/\alpha_c\big)^{1/(j-1)} .$$

It is a lower bound rather than an estimate, but it is not vacuous even at
$j = 2$, where it gives $\alpha^2/\alpha_c$: a pass fired, so
$\alpha > \alpha_c$, and hence $\alpha^2/\alpha_c > \alpha > \alpha_c$ —
strictly better than the fallback $\alpha_{\mathrm{wall}} \leftarrow
\alpha_c$ it replaces. §21's search uses it as its default for that reason,
so the bound does work rather than being recorded for completeness.

**When none is available, do not price $j \ge 2$.** Report the $j = 1$
answer: the least $k$ with $\alpha_k \le \alpha_c$, and the interval
$[\alpha_k, \alpha_c]$. This is conservative in both senses that matter.
The interval is contained in the true optimal region, since a larger $j$
would only raise the upper end and lower the required $k$; and the budget
$k+1$ is an over-estimate of the true $k'+j$. The estimate errs by
demanding more work over a narrower range, never by promising a range that
does not deliver.

This is the duality of the earlier part in operational form: a solve that
iterated sat below the ceiling and can price only the floor; a solve that
passed sat above the floor and can price only the ceiling. **Each prices
movement in the direction it already went**, and neither can recommend
reversing. The zero-work endpoints $\alpha_0$ and $\alpha_c$ are always
available and bound the region regardless; it is only the graded levels
that inherit the asymmetry.

### Exhaustion

A window exists at some small budget only while
$\alpha_\infty < \alpha_c^{(j)}$ for some affordable $j$. When it does not,
the criterion carries no information: a window becomes available only as $k$
approaches the number of distinct nodes — that is, only by solving the
correction system essentially exactly. This is a statement about the *levels*, not about the configuration, and
the distinction is not academic: the estimate's input converged, and so
exhibits a budget at which a window exists — its own, containing its own
$\alpha$. A configuration that a solve has solved cannot be infeasible.
When the criterion fires, therefore, what has been established is that the
level estimates cannot locate an interval, not that none exists; §21 returns
the witnessed one instead of a verdict.

The criterion retains its diagnostic value in that reading. It says the
tolerance, the preconditioner or the problem is where the difficulty lies
rather than $\alpha$ — but as a description of a regime the levels cannot
resolve, to be acted on by a caller, not as a claim the estimate is
entitled to make about feasibility.

## 19. Measuring $c_1$ directly

Only $n = 0$ leaves the floor family unknown, and that is exactly the case
in which the window is nonempty and no budget question arises: the window's
defining property is what makes it spectrally silent. Within the scope of a
single solve the levels are then not needed. Step 8 of §21 is the sole place
the question of buying them arises, and it declines — see the remark
following the algorithm, and the statement in §21 that $S$ is applied
nowhere. Nothing below is part of the estimate. A caller that reasons across
solves — a
tightening tolerance will raise $\alpha_0$ and may close the window — may
want them anyway. Evaluating quadratic forms in a matrix inverse by the
Lanczos/Stieltjes procedure, with the same remarks on inexact solves and
iterative refinement, is applied to this system at [GG, Sec. 3.2, p. 2086],
which credits it to Golub and Meurant. One application of
$S = BF^{-1}B^\top$ to $r_0$, a triangular-solve pair and two
multiplications by $B$, gives

    theta = (r0 . S r0) / (r0 . r0)          # the elasticity, = <mu>_w
    mu2   = (S r0 . S r0) / (r0 . r0)        # the second moment
    c1    = alpha * sqrt(mu2 - theta^2)

exactly, the two moments coming from the same application. $c_1$ is
$\alpha$ times the standard deviation of the measure — the $k = 1$ case of
§18, since $e_1$ is that same quantity — and the elasticity falls out of the
first moment at no extra cost. Being a standard deviation it is translation
invariant, so this is (21) at $k = 1$ exactly, with no asymptotic:
$\alpha\,\mathrm{sd}_w(\mu) = \mathrm{sd}(\tau)$. Routing instead through the
one-step minimal-residual reduction of (34) and inverting
(19) returns the same quantity divided by
$\langle\mu^2\rangle_w^{1/2}$; the discrepancy is precisely that equation's
$1 + o(1)$, and it vanishes as the measure concentrates near $\mu = 1$, which
is (A3). Prefer the form above, which carries no $o(1)$.

**The elasticity is free whenever any iteration ran.** $\theta$ is the first
moment $\langle\mu\rangle_w = r_0^\top S r_0 / r_0^\top r_0$, and the first
CG step length is $a_1 = (r_0^\top r_0)/(r_0^\top S r_0)$, so

    theta = 1 / a_1

exactly. No application of $S$ is needed for it. That leaves this section a
single purpose: the case $n = 0$, where no coefficient exists — and where,
the window being nonempty, the levels are not needed either. What may still
be wanted there is $\theta$ itself, as a check that the self-similarity the
floor formula assumes holds at the current $\alpha$; $\alpha_0$ rests on it,
and nothing else tests it.

**A $\theta$ bought at $n = 0$ must not be fed to the transport.** The
prohibition is worth stating because the purchase is available exactly where
the transport's one-sidedness fails. By (T1), $n = 0$ means
$\|r_0\| \le \varepsilon$, so the ratio $\|r_0\|/\varepsilon$ is at most $1$;
raising a quantity below $1$ to the power $1/\theta \ge 1$ *lowers* it, so
the corrected form returns a floor beneath the plain one and a window wider
than the tests support. §16's closure — that the plain form errs in one
direction only — holds because $\theta$ is unavailable on that branch, and
buying it removes the very fact that makes the closure exact. Use a measured
$\theta$ at $n = 0$ as a validity check on $\alpha_0$ and for nothing else;
it is the same refusal as §21's on $c_1$, one argument over, and for the
same reason: a quantity adequate for judging a reading is not thereby
adequate for moving it.

Since §16 now declines to apply $\theta$ at any $n$, this prohibition is the
$n = 0$ instance of a general rule rather than an exception to a general
licence. It is retained because the purchase is available only here, and
because the reason differs: elsewhere the correction fails by inflating the
distance, and here it would invert the sign of the bias.

Note the cancellation: inside the window the spectrum of $S$ clusters near
$1$, so $\theta$ and $\langle\mu^2\rangle_w$ are both near $1$ and their
difference is a subtraction of nearby quantities — the regime in which the
estimate is most likely to be wanted is also the one in which the
cancellation is worst. Form the residual vector $(S - \theta)r_0$ directly
and take its norm, rather than subtracting the scalars.

## 20. The shape of the problem

The two families are geometric, so each admits a real-valued count that
interpolates its staircase exactly. Inverting the level laws,

writing $\nu$ and $\vartheta$ for the real-valued counts (the latter
distinct from the elasticity $\theta$, with which it shares nothing but a
letterform):

$$\nu(\alpha)=\frac{\log(\alpha_0/\alpha)}{\log(\alpha/\alpha_\infty)},
\qquad
\vartheta(\alpha)=\frac{\log(\alpha_{\mathrm{wall}}/\alpha_c)}
                      {\log(\alpha_{\mathrm{wall}}/\alpha)},$$

with $n = \lceil\nu\rceil$ and $j = \lceil\vartheta\rceil$. Both are convex
on $(\log\alpha_\infty, \log\alpha_{\mathrm{wall}})$: writing $u$ for
$\log\alpha$, $\nu = L/(u - \log\alpha_\infty) - 1$ and
$\vartheta = M/(\log\alpha_{\mathrm{wall}} - u)$ have second derivatives
$2L/(u-\log\alpha_\infty)^3$ and
$2M/(\log\alpha_{\mathrm{wall}}-u)^3$, both positive there, so $C$ is
convex and its stationary point is its unique minimum. That is why the
closed form below exists; setting
$\nu = k$ recovers $\alpha_k = (\alpha_0 c_k)^{1/(k+1)}$ and setting
$\vartheta = j$ recovers $\alpha_c^{(j)}$. In the logarithm both are Möbius,
with poles at $\log\alpha_\infty$ and $\log\alpha_{\mathrm{wall}}$, and
$\nu + \vartheta$ is minimized in closed form. Writing
$L = \log(\alpha_0/\alpha_\infty)$ for the floor family's span,
$M = \log(\alpha_{\mathrm{wall}}/\alpha_c)$ for the ceiling's, and
$D = \log(\alpha_{\mathrm{wall}}/\alpha_\infty) = L + M - \log X$,

$$\log\alpha^\* =
\frac{\sqrt M\,\log\alpha_\infty + \sqrt L\,\log\alpha_{\mathrm{wall}}}
     {\sqrt L + \sqrt M},
\qquad
C^\* = \frac{(\sqrt L + \sqrt M)^2}{D} - 1. \tag{$\ast$}$$

Two things are worth reading off this, and one warning.

**The optimum is a weighted mean of the two poles, not of the two
endpoints.** $\alpha_0$ and $\alpha_c$ enter only through the spans, as
weights. The quantities that place the answer are $\alpha_\infty$ and
$\alpha_{\mathrm{wall}}$ — the two accumulation scales, which are also the
two quantities the estimate knows least well. That is a structural statement
about where the remaining uncertainty lives.

**$C^\*$ is a scalar difficulty.** It is computable before any $\alpha$ is
chosen and uses all four quantities, where the exhaustion ratio $X$ uses
two. $X$ answers whether a zero-work window exists; $C^\*$ answers what the
problem costs at its best.

**But $\alpha^\*$ is not the placement rule, and $C^\*$ is not the budget.**
Minimizing $\nu + \vartheta$ is not minimizing
$\lceil\nu\rceil + \lceil\vartheta\rceil$, and the reason is sharper than
rounding. The ceiling of a convex function is quasi-convex — for integer
$t$, $\{\lceil f\rceil \le t\} = \{f \le t\}$, an interval — so
$\lceil\nu\rceil$ and $\lceil\vartheta\rceil$ are each unimodal. But
quasi-convexity is not preserved under addition, and their sum is in general
neither convex nor unimodal: **the set of cost-minimizing $\alpha$ is a union
of intervals, not one interval.** The two staircases step at different
points, so between two equally-cheap plateaus the sum can rise where both
have risen and fall again beyond.

Note where the loss occurs. $\lceil \nu + \vartheta \rceil$ — the ceiling
taken last — is quasi-convex, being the ceiling of a convex function. It is
the *separate* rounding that costs unimodality, and it cannot be avoided: a
fractional iteration cannot be made up with a fractional pass. The smooth
minimum can therefore land not merely in a suboptimal region but in a gap
*between* optimal blocks, which is the failure the closed form is subject to
and the exact reason $\alpha^\*$ is a poorer placement than the search of
§21. What can be said exactly is two-sided: since
$\lceil\nu\rceil + \lceil\vartheta\rceil \ge \nu + \vartheta \ge C^\*$ at
every $\alpha$, and $\lceil\nu\rceil + \lceil\vartheta\rceil < \nu +
\vartheta + 2$,

$$\lceil C^\*\rceil \;\le\; \min_\alpha\big(\lceil\nu\rceil +
  \lceil\vartheta\rceil\big) \;<\; C^\* + 2 .$$

So $\lceil C^\*\rceil$ is a lower bound on the budget and never an
over-estimate of it, and the discrete search can exceed it by up to two
solve-pairs; the ceilings do not commute with the sum. The integrality is not incidental algebra to be
smoothed away: it is where the flatness lives, and the flatness is what the
estimate exploits. Use $(\ast)$ to understand the problem and to score it;
use §21 to answer it.

## 21. The estimate

Notation. $\alpha, \rho, \varepsilon$ are the parameter, shift and tolerance
of the completed solve; $n$ its Krylov iteration count and $j$ its backsolve
count; $r_0$ and $s_0$ the two recorded residuals; $(a_i, b_i)$ the CG
coefficients of the first $n$ iterations, from which $e_i = \alpha\sqrt{b_i}/a_i$
(§1); $s_i$ the dual residual recorded
after pass $i$ (only $s_1$ and, where a ceiling level is wanted, $s_2$ are
used); $\phi$ the ceiling transport exponent, supplied by the caller and
defaulting to $1$. Throughout, $S = BF^{-1}B^\top$ is the operator of §1, on
which the iteration runs and in whose variable every spectral quantity is
read. **The algorithm applies $S$ nowhere.** Every
quantity it uses is either a residual the solve already formed or a
coefficient the iteration already computed; the sole conditional extra work
is the bias correction of step 1, which is a correction rather than an
investment and is guarded by a free screen.

---

**Precondition.** The input solve **converged**: it reached the tolerance at
$\alpha$ in $n$ iterations and $j$ backsolves. Nothing below is defined for
an input that did not, and the precondition is what makes the return
contract of §15 exact rather than nearly so.

**Invariant.** The pair $([\alpha,\alpha],\ n+j)$ is admissible by that
hypothesis and is carried throughout as the incumbent: every branch returns
it unless it finds something cheaper, and no branch returns a budget above
$n+j$. In particular the algorithm has no infeasible outcome.

**Algorithm** ESTIMATE$(\alpha, \rho, \varepsilon; r_0, s_0, n, j, a, b, s, x)$

**returns** a triple $(\mathrm{lo}, \mathrm{hi}, \mathrm{cost})$ as specified in §15

**Order of computation.** The endpoint *values* depend on nothing but the
recorded residuals; their *qualities* depend on the treads, hence on the
wall and the levels. The steps are therefore grouped: values first, then the
scales the gates need, then the markings, then the comparison and the
search. Nothing below is a forward reference.

1. **[floor value]**
 **if** $\rho = 0$ **or** $\tfrac12\sqrt{\rho}\,\|x\| \ll \varepsilon$
 **then** *(the screen clears — free; no extra work)*
  **skip to** the assignment below
 **else** *(and only else)*
  $r_0 \leftarrow r_0 - \rho B F^{-1} x$
  *— **the one conditional solve-pair**: an $F^{-1}$ then a multiply by
  $B$, gated by the screen above and never executed when it clears —*
  **if** the reading lies on the bias floor **then return**
  $(\star, \alpha_c, n+j)$ *(floor invalid, ceiling intact; the budget is
  the witness, not a price on $(\star, \alpha_c]$)*.
 $\alpha_0 \leftarrow \alpha\,\|r_0\| / \varepsilon$;
 $\ \theta \leftarrow 1/a_1$ if $n \ge 1$, else unavailable *(free)*.

2. **[ceiling value]**
 $\hat\phi \leftarrow T_B/(T_A+T_B)$ with
 $T_A = \beta\max_i A_{ii}$, $T_B = \max_i (B^\top B)_{ii}$ *(free)*;
 **if** the reading and the endpoint lie on opposite sides of a shift
 engagement **then** $\alpha_c \leftarrow \star$ *(no value to transport;
 §16)* **and skip to** step 6, which will withhold.
 **otherwise** $\alpha_c \leftarrow \alpha\,(\varepsilon/\|s_0\|)^{1/\phi}$,
 with $\phi$ as supplied and $1$ by default.

3. **[wall]** — measured routes first, since they need only $\ell_d$,
 $\|x\|$ and the factor, none of which depends on $\alpha_c$:
 $\ \varrho_d \leftarrow \|s_2\|/\|s_1\|$ and
 $\alpha_{\mathrm{wall}} \leftarrow \alpha/\varrho_d$ if a pass ran and the
 screen of §18 is met; otherwise
 $\alpha_{\mathrm{wall}} \leftarrow \alpha\cdot\min_i L_{ii}^2\cdot\|x\|/\ell_d$
 from the factor, which is available at every $j$ including $j = 1$.
 If neither is available, mark $\alpha_{\mathrm{wall}}$ UNAVAILABLE; the
 pass-count default of step 6 will supply it.

4. **[levels]**
 $c_k \leftarrow \prod_{i\le k} \alpha\sqrt{b_i}/a_i$ for $k = 1,\dots,n$
 *(free; $n$ iterations give $n$ coefficients)*, and
 $c_k \leftarrow$ UNAVAILABLE beyond.
 **[direction guard]** clamp every level at the probe:
 $\alpha_k \leftarrow \max(\alpha_k,\ \alpha)$, marking any clamped level
 $\alpha_k^+$ — the true level lies at or below it (§18).
 $\ \Lambda \leftarrow \log(\alpha_0/\alpha_\infty)$ and
 $M \leftarrow \log(\alpha_{\mathrm{wall}}/\alpha_c)$, the two spans the
 gates below compare against; either is UNAVAILABLE if its pole is.

5. **[qualities]** Both gates compare a *magnitude* against a tread.

 - **Floor.** The bracket of §16 is on $\alpha_0$ itself, so the relevant
   tread is the first, $\Lambda/2$. Where
   $(1/\theta - 1)\log(\|r_0\|/\varepsilon) < \Lambda/2$, the value
   $\alpha(\|r_0\|/\varepsilon)^{1/\theta}$ may accompany $\alpha_0$ as the
   bracket's upper arm; otherwise the plain value stands alone, a bracket
   wider than a tread being uninformative. $\theta$ is reported in either
   case and **never applied as an exponent** (§16, §19).
 - **Ceiling.** Where
   $\big|1 - 1/\hat\phi\big|\cdot\big|\log(\varepsilon/\|s_0\|)\big|$
   exceeds $M/(j(j+1))$, mark $\alpha_c$ as $\alpha_c^-$ — a lower bound
   (§16). The magnitude is wanted here and the sign in §16's prose: since
   $\hat\phi \le 1$ the signed displacement is non-positive on the reachable
   branch, so a signed test would never fire. **If $M$ is UNAVAILABLE the
   test cannot be evaluated and the marking is applied**, stated rather
   than obtained from a degenerate tread.

6. **[compare]** Apply the rule of §15: act on the comparison only if the
 sign carried by $\alpha_c$ decides it, and otherwise withhold the
 conclusion and return the endpoints as the bounds they are.

 **if** $\alpha_0 \le \alpha_c$ **then return** $(\alpha_0, \alpha_c, 1)$,
 the ceiling carrying whatever quality step 5 gave it.
 *(Decided: with $\alpha_c$ exact, directly; with $\alpha_c^-$, because the
 true ceiling is higher still, so the interval returned is a subset of the
 true window and the budget holds. The solve stood inside the window iff
 $n=0$ and $j=1$.)*

 **else if** $\alpha_c$ is a lower bound **then return**
 $(\alpha_0, \alpha_c^-, n+j)$.
 *(Undecided, so the conclusion is withheld: an understated ceiling is
 precisely what produces a spurious emptiness. The floor stands, the
 ceiling is returned as the bound it is, the budget is the witness, and the
 remedy is a probe at larger $\alpha$ rather than a budget. The search of
 step 7 may not proceed, since it assumes emptiness.)*

 *Only the lower bound arises, by the switch–ceiling ordering of §16, so this
 is the whole of the rule. Were that ordering to fail — it rests on the
 flatness of $\ell_d$, not on a theorem — an $\alpha_c^+$ could occur, and
 §15's rule would then guard the opposite branch: the emptiness would be
 decided and the nonemptiness withheld.*

 — *from here the window is empty and known so: $\alpha_0 > \alpha_c$* —

7. **[the search]** Set two conservative defaults where the inputs are
 missing, and then minimize once.

 - **If $\alpha_{\mathrm{wall}}$ is unavailable**, set it to the bound the
   observed pass count already implies. The probe used $j$ backsolves, so
   $\alpha \le \alpha_c^{(j)}$ by the definition of $j(\alpha)$, whence
   $$\alpha_{\mathrm{wall}} \leftarrow \max\Big(\alpha_c,\
     \big(\alpha^{\,j}/\alpha_c\big)^{1/(j-1)}\Big) \quad (j \ge 2),
     \qquad \alpha_{\mathrm{wall}} \leftarrow \alpha_c \quad (j = 1).$$
   The floor $\alpha_c$ makes $\alpha_c^{(j')} = \alpha_c$ for every $j'$, so
   passes buy nothing and the minimization selects $j' = 1$; the bracket
   improves on that whenever a pass ran, and is the third route of §18 doing
   work rather than being noted for completeness.
 - **If $c_k$ is unavailable** — that is, $n = 0$ — set
   $\alpha_k \leftarrow +\infty$ for $k \ge 1$. Only $k = 0$ is then
   feasible, so the minimization selects the least $j'$ with
   $\alpha_c^{(j')} \ge \alpha_0$. This enforces structurally what §21's
   remark argues: with only $c_1$ in hand the extension
   $c_k \approx c_1^{\,k}$ places the levels too high, makes a floor climb
   appear cheaper than it is, and displaces an already-accurate ceiling
   answer. The rule is a default, not a prohibition to remember.

 With those in place, let $(k, j')$ minimize $k + j'$ subject to
 $\alpha_k \le \alpha_c^{(j')}$, **breaking ties by the largest $k$**, and
 **return** $(\alpha_k,\ \alpha_c^{(j')},\ k+j')$. If no pair qualifies,
 **return** the incumbent $([\alpha,\alpha],\ n+j)$.

 *(The tie-break selects the lowest of the optimal components, §15: since
 $\alpha_k$ falls with $k$ and $\alpha_c^{(j)}$ rises with $j$, trading a
 pass for an iteration at fixed total moves the whole interval down. It is
 the safe direction by §17 — an underspent $\alpha$ costs iterations the
 next solve reports, an overspent one risks a stall — and it costs nothing,
 the budget being equal by construction. It also selects the
 best-determined component: by the witness invariant $k \le n$, so the
 $c_k$ it uses was computed by this solve rather than extrapolated. Indeed
 the extrapolation beyond the harvest is unreachable: every level with
 $k \ge n$ clamps to $\alpha$, so $\alpha_\infty$ is never consulted by the
 search, and enters only §20's closed form.)*

 *(A clamped $\alpha_k$ is returned as the bound it is, $\alpha_k^+$: the
 true level lies at or below the probe, so the interval is a subset of the
 true one and the budget is attainable there, merely slack. Do not blank
 it — a $\star$ floor would forbid the caller from descending, when what is
 known is precisely that it may. Where the ceiling default was used and no
 pass ran, the upper endpoint is $\alpha_c$ itself and the answer is
 conservative twice over: a larger $j'$ would raise the upper end and lower
 the required $k$, so both the interval and the budget err inward.)*

---

**Cost.** The estimate is free. Both endpoints come from residuals the solve
already formed, the floor levels from coefficients the iteration already
computed, the ceiling levels from passes that already ran. The one
exception is conditional and is not an information purchase: step 1's bias
correction costs a solve-pair, is guarded by a free screen, and repairs a
floor that would otherwise be *wrong* rather than merely imprecise.

**Remark (why no level is bought).** When the floor family is unavailable
it is tempting to buy $c_1$ — one application of $S$ — and extend by
$c_k \approx c_1^{\,k}$. This is worse than not buying. The extension places
the levels too high, which is the conservative direction when *reporting* a
level, since one then believes a level unavailable slightly before it is;
but it is the *anti*-conservative direction when **choosing between
branches**, because a floor climb whose levels sit too high appears cheaper
than it is and displaces a ceiling answer that was already accurate. A
crude family is adequate for widening an interval and inadequate for a
decision. The same measurement remains the only source of $\theta$ (§19),
which is a validity check rather than a level, and that use is unaffected.

**Where the cases land.** With $(n,j)$ the counts of the completed solve:

    n = 0, j = 1    step 6 — window nonempty, nothing unknown
    phi_hat low     step 6 — nonempty reported as usual; apparent
                    emptiness withheld as (alpha_0, alpha_c^-, n+j)
    n >= 1, j = 1   step 7 — floor family free; ceiling defaults to alpha_c,
                    so the search selects j' = 1
    n = 0, j >= 2   step 7 — floor defaults to +infinity, so the search
                    selects k = 0 and climbs the ceiling family
    n >= 1, j >= 2  step 7 — both families known; the pass-count bracket
                    raises the ceiling default if no wall was measured
    stall           not an input to this algorithm; see §17

## 22. What to record

    triple       (lo, hi, cost) as returned, each endpoint carrying its
                 quality: exact, ^- (lower bound), ^+ (upper bound), or *;
                 an inversion marks a withheld comparison, never a verdict
    C*           the scalar difficulty of (*) in §20, when all four
                 quantities are available; a problem-level score that does
                 not depend on which alpha was used
    poles        alpha_inf and alpha_wall, with their provenance — these
                 place the optimum and are the least well known
    theta        the elasticity, = 1/a_1 on S (Sec. 1) whenever n >= 1.  It
                 is a meter, not an exponent: it reports A1's validity at the
                 probe and gates Sec. 16's bracket, and is not applied to the
                 floor transport (Sec. 16)
    phi_hat      T_B/(T_A+T_B) from the two diagonals; the ceiling's meter,
                 near 1 on the ramp and falling to 1/2 at the switch
    l_d          ||s0||/alpha^phi, the dual invariant at the exponent used;
                 two of these from different alpha determine phi for the next
                 estimate, and their agreement is the ramp check
    counts       Krylov iterations and refinement passes actually incurred
    rho          the shift applied by the factorization, 0 if none
    regime       ramp or bias-floor, if two readings are available
    coefficients (a_i, b_i), i = 1..n, for the levels and for later
                 comparison; a_1 alone gives theta = 1/a_1, and the levels
                 follow as c_k = prod e_i with e_i = alpha sqrt(b_i)/a_i
    Ldiag_min    the smallest pivot of the Cholesky factor, for lambda_Z
    ||x||        the assembled solution's norm, for a_wall and for the
                 regularized ceiling
    provenance   per endpoint: floor = harvested / extrapolated beyond the
                 harvest depth / clamped by the guard / unknown;
                 ceiling = observed pass / carried lambda_Z / single reading
    probe        the alpha the solve ran at — with it a clamped floor is
                 readable as "below alpha", which a bare * does not convey

## 23. Open items

The items divide by what a reader wants them for. The first group is what is
not known: questions live inside the note's frame, whose answers would make
it strictly better. The second is what is known and imperfect: things the
estimate does, each degrading in a stated direction. Nothing in the second
group invalidates a return; each item bounds what a return means, and that is
the group to read before trusting one.

### Not established

- **Whether $n_{\mathrm{eff}}$ can collapse.** The dominance closure of §17
  rests on assumption (D): that the backward error's participation ratio is
  large, so the single draw's dispersion stays small against the tread
  spacing. The rate $\operatorname{sd} \simeq c\,n_{\mathrm{eff}}^{-1/2}$ is a
  second-moment identity and is not at issue; what is at issue is whether
  $n_{\mathrm{eff}}$ itself can be small. (D) fails only if the perturbation
  concentrates on few directions,
  which operator ill-conditioning alone does not produce — it requires a
  spread in $F$ unmatched by the solution's. The open question is therefore
  narrower than it looks: not whether the dispersion obeys a law, but whether
  any operator class this note addresses makes $n_{\mathrm{eff}}$ collapse.
  Constructing such a case, or ruling one out, is the outstanding task; no
  attempt here has produced one.
- **How far the graded branch sits from the optimum.** The two branches
  differ in kind. The window branch reproduces the optimal interval exactly,
  even at a width of zero decades; the graded branch names an interval
  usually within one solve-pair of optimal and rarely the optimal one. The
  gap is in the level formula, not in the search. §20 bounds the *smooth*
  argument's gap, $\lceil C^\*\rceil \le \min(\lceil\nu\rceil +
  \lceil\vartheta\rceil) < C^\* + 2$, which says nothing about the interval
  steps 4–7 name. The witness of §15 caps the returned budget at $n+j$, hence
  the excess at $n + j - C^\*$ wherever $C^\*$ is available — but that is a
  bound inherited from the input, not one the levels earn, and it is the
  levels that would repay work. The regime where the question bites is deep
  emptiness with a badly conditioned $A$: both $X$ and $\kappa(A)$ are known
  before the estimate runs, so the case is detectable, but neither a bound on
  the excess nor what to report instead is established.
- **The ordering between the two wall routes.** §18 prefers the observed pass
  to the factorization, and no argument here ranks them: the pass route is
  screened only necessarily, and the pivot route compounds a bound whose
  tightness depends on the ordering with an identification
  $\lambda_Z \approx \alpha\lambda_{\min}(F)$ approximate in its own right,
  so neither error is one-signed. What is established is that either
  suffices, that the first ratio must be used rather than an average, and
  that the two must not be combined.
- **Whether the smooth form can replace the search.** $\nu$ and $\vartheta$
  interpolate the staircases exactly, are convex, and give a closed-form
  optimum, but the minimizer of $\nu + \vartheta$ is not the minimizer of
  $\lceil\nu\rceil + \lceil\vartheta\rceil$, which by §20 is in general not
  even unimodal. Any surrogate that is to replace the search must therefore
  reproduce a disconnected argmin, which no smooth convex function does;
  whether one respecting the flatness — penalising the ceilings rather than
  ignoring them — could do so in closed form is not known, and the obstacle
  is now identified rather than merely observed.
- **The optimal set may be disconnected, and the estimate returns one
  component.** §15 records this; what is not established is how often the
  gaps matter, how wide they are, or whether a caller holding one component
  can bound what lies in the others. The bound
  $\lceil C^\*\rceil \le \min(\lceil\nu\rceil+\lceil\vartheta\rceil)$
  applies to the global minimum over all components, so it at least says the
  returned budget is within two of every component's.

### Established, and imperfect

- **No level below the probe is determined.** A solve fixes no $\alpha_k$
  below the parameter it ran at, and the guard of §18 clamps the ones that
  would fall there. The case in which those levels would change the answer is
  exactly the case in which the solve performed no iterations — which is also
  the case in which the window is nonempty and no budget question arises, so
  the loss is smaller than it first appears. No single-solve remedy is known.
  The options are a $c_1$ calibrated across solves, a deliberate probe at
  lower $\alpha$, or accepting the conservatism; the estimate takes the
  third, and §21 declines to buy $c_1$ for the reason given there.
- **$\hat\phi$ meters the ceiling and does not correct it.** $\theta$ is
  measured from the solve; $\phi$ must still be supplied, a slope needing two
  readings. What $\hat\phi$ adds is a verdict on whether the default is right
  without measuring the slope, so the ceiling is no longer the only endpoint
  without a validity check. What remains unmetered is the size of the error
  once $\hat\phi$ falls: the meter sorts good readings from bad and repairs
  neither. The remedy is a probe at larger $\alpha$ rather than a supplied
  exponent, since a supplied $\phi$ fixes the arithmetic and not the fact
  that the sensor is reading its own noise floor.
- **The screen is necessary, not sufficient.** A probe failing
  $\alpha^2\sigma^2_{\min} \gtrsim \alpha_{\mathrm{wall}}$ is certainly
  reading the augmentation branch; one passing it may still be inaccurate.
  And $\sigma^2_{\min}$ degenerates when $A$ is singular, where the screen
  cannot be evaluated at all.
- **The bias correction uses the computed solution.** Step 1 forms
  $\rho B F^{-1} x$ with the solve's own $x$ in place of the exact one. This
  is sound while the solve is converging, which is not guaranteed in the
  regime where the shift engages.
- **The switch–ceiling ordering is derived; its converse is not.** That
  $\alpha_\star \le \alpha_c$ in any solvable problem — hence that an
  off-ramp reading always understates, and that the $\alpha_c^+$ marking is
  unreachable — follows from the upper bound (C3, (10b)), which is derived
  rather than assumed. What is not derived is the converse reading of §16's
  equivalence, that a ceiling below the switch means no ceiling at all: that
  needs (10b) attained and not merely bounding, which (C3) supplies only to
  within a factor of order one. Nothing in the algorithm depends on the
  converse.

## 24. Assumptions, organized by role

The premises fall into three groups, distinguished by what they do rather
than by which endpoint they belong to. They are all stated under one standing
hypothesis, which is not a member of any group because the groups presuppose
it.

| | statement | consequence | if it failed |
|---|---|---|---|
| **PSD** | $A$ is symmetric positive semidefinite, and $B$ has full row rank on the relevant subspace (§1) | $F$ is positive definite for every $\alpha > 0$ and admits Cholesky; the parameter is bounded below by cost alone; $\lambda_Z \ge 0$, so §10's wall and §18's fallback are meaningful | everything. An indefinite $A$ imposes $\gamma_{\min} = -\lambda_{\min}(A)/\lambda^*_{\min}(B^\top B) > 0$ [RC, Thm. 1, p. 7], a hard floor of a different kind from §5's and binding before it; no endpoint, family or cost below survives unaltered. The regime is not exotic — it is where an interior-point method spends part of its run [RC, Sec. 4.1, p. 7] and what [GS] is written to analyse |

### Group 1 — the convergence tests

These are read from the algorithm, not assumed about the problem. They are
what fix both constants to 1.

| | statement | consequence | if it changed |
|---|---|---|---|
| T1 | the Krylov stopping test is evaluated at entry, against an **absolute** tolerance | primal endpoint at $\|r_0\| = \varepsilon$ | a test on a multiple of $\varepsilon$ carries that multiple into (8); a *relative* test is not a multiple but a different family, $\alpha_k = (c_k/\varepsilon_{\mathrm{rel}})^{1/k}$, in which $\ell_p$ does not appear at all — see §1 |
| T2 | the refinement test admits no borrowing between components | dual endpoint at $\|s_0\| = \varepsilon$ | a combined test reduces the dual allowance by the other components' share |

### Group 2 — the scaling laws

One per endpoint; these and only these are consumed by the derivations. Each
has a local form, which is what is used, and a global form, which implies it
and is free of any reference to where the solve was probed.

| | statement | status |
|---|---|---|
| A1 | $\|r_0\| \propto 1/\alpha$ on the transport interval | **used** — gives the primal endpoint (8) |
| A2 | the same over the operating range | implies A1 |
| C1 | $\|s_0\| \propto \alpha$ on the transport interval | **used** — gives the dual endpoint (12) |
| C2 | $\|s_0\| \propto \alpha$ above the switch $\alpha_\star = \lambda_{\max}(A)/\sigma^2_{\max}$, (10a) | implies C1 when the transport interval lies above $\alpha_\star$; checkable from the operator alone |
| P | the correction of a refinement pass lies in $\ker B$, its size bounded by (26) | **used** — gives the pass contraction (27), hence the graded ceiling (28) and the wall | exact to the extent the pass's right-hand side has negligible primal component, which the primal condition guarantees; the bound is worst-case, so $\alpha_{\mathrm{wall}}$ errs conservative |

A1 fails as $\alpha$ falls and the tail populates, with error one-signed,
understating depth below $\alpha_0$; under a regularized factorization it
acquires a second failure at large $\alpha$, where the initial residual
meets a floor rather than continuing to fall. C1 fails below the switch $\alpha_\star$, where
$\|F\|$ ceases to grow with the parameter. Each therefore fails in the regime the
*other* endpoint governs — A as $\alpha$ falls, C as it rises — so neither
is asked to hold where it is least justified.

The two are mirror images: A asserts that a product $\alpha\|r_0\|$ is
$\alpha$-free, C that a quotient $\|s_0\|/\alpha$ is.

### Group 3 — the mechanism layer

Nothing here is consumed by a formula. This group explains why the scaling
laws hold, names their constants, and is what one would try to establish for
a class of operators.

| | statement | supports | status |
|---|---|---|---|
| A3 | resolved-mass: $\|P\delta\| \approx \|\delta\|$, $W$ negligible | A2, and names $\ell_p = \|\delta\|$ | quantified by the elasticity $\theta = \langle\mu\rangle_w \le 1$. That an alignment premise is needed at all is not peculiar to this route: the conditioning analysis requires one too, a minimum of $\kappa(A + wbb^\top)$ existing for $w>0$ iff $\|b^\top u_1\|/\|b^\top u_n\| < \kappa_2(A)$ — that is, iff the augmenting directions form a small angle with $A$'s smallest eigendirection [GGV, Prop. 3.3, p. 786]. Two analyses of the same parameter, by different criteria, each needing a statement about how $B$ meets the near-kernel of $A$ |
| C3 | $\|s_0\| \lesssim u\max(\lambda_{\max}(A), \alpha\sigma^2_{\max})\|x_0\|$, (10b) | **derived** from Prop. 2 and backward error; yields C2 and hence C1 | an upper bound, attained up to a factor of order one when the perturbation $\Delta F$ is not orthogonal to $x_0$ |
| S | two-scale spectrum: $A$-supported near-kernel of $B$ | explains the *cost* side, $\kappa(F) \propto \alpha$ (10); not used by the dual residual, which depends on $\|F\|$ alone | needed only where $\kappa(F)$ is; $h$ enters no endpoint. The split is exact in exact arithmetic: the spectrum of $A + \alpha B^\top B$ separates into intervals independent of $\alpha$, anchored on the spectrum of $Z^\top A Z$ — whose least eigenvalue is the $\lambda_Z$ of §18 — and intervals growing linearly in $\alpha$ on the complement of $\ker B$ [GS, Prop. 2.12, p. 1160]. That separation is why the conditioning worsens when the parameter is taken too large [GS, p. 1161], which is this note's ceiling seen in exact arithmetic rather than in rounding |
| D | no small subset of directions carries the backward error: the participation ratio $n_{\mathrm{eff}} = (\sum\eta_i^2)^2/\sum\eta_i^4$ is large, with $\eta_i \propto (\|L\|\|L^\top\|\|x\|)_i$. The condition is on $n_{\mathrm{eff}}$ absolutely, not on $n_{\mathrm{eff}}/n$ | the dominance closure of §17, through $\operatorname{sd}(\log\alpha_c) \simeq c\,n_{\mathrm{eff}}^{-1/2}$ with $c$ of order one — a second-moment identity, not a central limit, so the exponent survives the failure of Gaussianity while the constant does not. Correlation deflates the effective count, so the rate understates the dispersion and is used only as a lower bound on it | not consumed by any formula, and checkable free from the factor and the solution. Fails only if the perturbation concentrates on few directions, which operator ill-conditioning alone does not cause: it requires a spread in $F$ *unmatched* by the solution's. For the augmented block in an interior-point method the spread is structured rather than arbitrary: its large eigenvalues align with the range of the active Jacobian and its small ones with that Jacobian's null space [PS, Thm. 2, p. 17], and the resulting errors separate accordingly, so that quantities formed through the inverse retain full precision [PS, Prop. 3, p. 22] |
| R | the intermediate regime $k \ll N$ of §9 | the graded levels (22) and the accumulation form (24) | $R_k \propto \alpha^{-k}$ with $c_k = \|\pi_k\|$ is derived, not assumed (19), (21); what is assumed is only that $k$ is small relative to the number of distinct weighted eigenvalues, beyond which the family terminates rather than accumulating |

Note the asymmetry between the two mechanisms, and that it has narrowed.
A3 is direction-aware — a statement about where the gap's mass sits — and so
needs no companion condition; it remains an assumption, quantified but not
proved. C3 is now a derivation rather than an assumption, and needs no
companion either: because the recovery step reintroduces the $F$ that the
backsolve inverted, the dual residual is governed by $\|F\|$ alone, with no
condition number, no near-kernel scale, and no condition on which directions
the right-hand side excites. Of the two scaling laws, only the primal one
still rests on a hypothesis about the problem.

### The exact results

| | statement |
|---|---|
| Prop. 1 | $r_0 = -\frac{1}{\alpha}S(y_0-\lambda^*)$ — exact algebra, no assumption |
| Prop. 2 | the dual row of the assembled pair vanishes identically in exact arithmetic |

Proposition 2 gives the ceiling's mechanism an unusual status: $s_0$ is
rounding error in its entirety, so C3 is a statement about which rounding
process dominates rather than about any exact quantity.

### Why the local forms are the right ones

Both endpoints are *defined* by a residual meeting the tolerance, which
requires no assumption at all. The scaling laws enter only to carry a
reading from the parameter used to the endpoint sought, so their domain is a
bounded interval with computable endpoints, not the operating range. The
global forms are stronger, are what one would establish for a class of
operators, and are strictly more than the derivation needs — a distinction
that matters whenever the class in question has spectra that move over many
orders of magnitude.

## 25. Summary

For a saddle-point system solved by augmentation, Krylov correction and
refinement, the set of augmentation parameters at which one backsolve
suffices is the interval

$$W_0 = \Big[\frac{\ell_p}{\varepsilon},\ \frac{\varepsilon}{\ell_d}\Big]
= \Big[\frac{\alpha\|r_0\|}{\varepsilon},\
\frac{\alpha\varepsilon}{\|s_0\|}\Big],$$

whose endpoints are the two convergence tests read literally against two
residuals with opposite and exact dependence on $\alpha$: the initial
constraint residual, inversely proportional by the identity of §3 under
self-similarity; the dual residual, directly proportional by the
conditioning of §7 under attainability. The interval is nonempty iff
$\ell_p\ell_d \le \varepsilon^2$, and its width $\varepsilon^2/\ell_p\ell_d$
measures the precision available to the tolerance demanded. Admitting $k$
Krylov iterations extends the interval downward, to
$\alpha_k = (\ell_p c_k/\varepsilon)^{1/(k+1)}$, the levels crowding
together as $1/k^2$ about the logarithmic capacity $\alpha_\infty$ of the
weighted spectrum's support — not about unity, which the units forbid —
until they terminate, once $k$ reaches the number of distinct weighted
eigenvalues, in a window unbounded below. The constants $c_k$ are norms of
the orthogonal polynomials of the residual's spectral measure, of which the
iteration computes an exact $k$-node representation at the parameter it ran
at, and a systematically optimistic one below it — and none at all inside
the window, where no iteration runs. There the first level is available
either for one application of $S$, which returns $c_1$ as a standard
deviation of the spectrum, or from a constant calibrated over a family:
$c_1$ depends on neither $\alpha$ nor $\varepsilon$, and enters the
endpoints through a $(k+1)$-th root, so a calibration error is damped rather
than transmitted in full.

Admitting $j$ backsolves extends the interval upward by the mirror law,
$\alpha_c^{(j)} = \alpha_c^{1/j}\,\alpha_{\mathrm{wall}}^{1-1/j}$, the
levels crowding as $1/j^2$ about
$\alpha_{\mathrm{wall}} = \lambda_Z\|x\|/\ell_d$ — a scale at which the
passes cease to contract, and past which, unlike the primal family, no
amount of work suffices. The span the passes can buy is
$\log(\lambda_Z\|x\|/\varepsilon)$, the decades between the tolerance and
the near-kernel scale of $A$.

The two families are the two arms of the cost. Iterations and refinement
passes are paid in one currency, the triangular-solve pair, and the total
$n(\alpha) + j(\alpha)$ is a descending staircase plus an ascending one,
each with treads crowding as the inverse square of its index. Its minimum is
a broad plateau rather than a point, and the window $W_0$ — where both
staircases stand at their lowest tread at once — is a sufficient condition
for that minimum, not a necessary one: it may be empty while the cheapest
solve costs two or three solve-pairs.

## 26. Transfer along an interior-point trajectory

**This section is different in kind from the rest and is marked as such.**
Everything above concerns one solve and reports no measurements. What
follows concerns a *sequence* of solves — the KKT systems an interior-point
method generates as the barrier parameter falls — and it reports
measurements, because the questions it raises are not settled by derivation
alone. It is provisional. The derivable part is given first, so that what is
measured has something to be measured against.

The setting is the one an implementation faces. A bracket is computed at
barrier iteration $k$; at iteration $k+1$ the operator has changed, the
tolerance has changed, and the question is what of the previous answer
survives. In the runs examined the tolerance is tied to the barrier
parameter exactly, $\varepsilon = c\,\mu$ with $c$ constant within a
problem, so $\varepsilon$ and $\mu$ are collinear and no measurement below
can separate a law in one from a law in the other. That is a limitation of
the data and is flagged wherever it bites; it is not a limitation on using
the results to *predict*, since the two move together in the solver as well.

### 26.1 What transfers by derivation

Both endpoints depend on $\varepsilon$ through exactly one factor, and the
graded exponents damp it. Since $\ell_p$, $\ell_d$, $c_k$ and
$\alpha_{\mathrm{wall}}$ are all $\varepsilon$-free — the first two by
construction, $c_k$ because the Krylov coefficients do not depend on where
the iteration stops, $\alpha_{\mathrm{wall}}$ because it is spectral —

$$\mathrm{lo}' = \Big(\frac{\varepsilon}{\varepsilon'}\Big)^{1/(k+1)}\mathrm{lo},
\qquad
\mathrm{hi}' = \Big(\frac{\varepsilon'}{\varepsilon}\Big)^{1/j}\,\mathrm{hi},$$

exactly, for a fixed $(k,j)$ and a fixed operator. The window case is the
$k=0$, $j=1$ instance, where the exponents are $\mp 1$.

Three consequences. The exponents are **reciprocal** because the two tests
bound quantities that move oppositely in $\alpha$: tightening $\varepsilon$
demands more augmentation to satisfy the primal test and permits less before
the dual one fails, so the floor rises and the ceiling falls. The window's
width therefore loses *two* decades per decade of tightening —
$\log(\alpha_c/\alpha_0) = 2\log\varepsilon - \log(\ell_p\ell_d)$ — which is
the quadratic $\varepsilon$-dependence of the exhaustion ratio $X$ seen from
the trajectory, and the reason late iterations run out of window abruptly.
And a *graded* bracket is less sensitive to the tolerance than a zero-work
one, the exponents falling as $1/(k+1)$ and $1/j$: buying budget on either
side slows the closure.

The derivation holds within a branch. Across a branch change — the returned
$(k,j)$ shifting as $\alpha_0$ rises — it says nothing, and the change is
discrete.

### 26.2 The floor law holds

Measured against the observed levels, with the level index taken as the
base-solve Krylov count and grid-censored edges excluded, the floor exponent
reproduces $-1/(k+1)$:

| $k$ | measured | $-1/(k+1)$ |
|---|---|---|
| 0 | $-1.000$ | $-1.000$ |
| 1 | $-0.499$ | $-0.500$ |
| 2 | $-0.328$ | $-0.333$ |

Three decimal places at the first two levels. This is the clearest
confirmation of the level law (22) available, and it is worth saying what it
confirms: not merely that the floor moves with the tolerance, but that it
moves with the *exponent the graded family predicts*, which is a statement
about $c_k$ and the crowding rather than about $\alpha_0$ alone.

Two cautions. Censoring is not incidental — at $k = 3$ only about half the
pairs survive the filter and the fitted exponent collapses to zero, which is
residual censoring rather than a breakdown of the law. And the index must be
the Krylov count, not the logged solve-pair count, which exceeds it by one.

### 26.3 The dual level is not constant: $\ell_d \propto \mu^{1/2}$

The ceiling transfer assumes $\ell_d$ is a property of the operator and the
machine, hence stationary along the trajectory. It is not. Across the runs
examined,

$$\frac{d\log \ell_d}{d\log\mu} \;\approx\; +0.52
\qquad (\text{median over problems; interquartile } [0.43,\,0.68],
\text{ positive in } 40 \text{ of } 41),$$

so the dual level *falls* as the barrier tightens, roughly as $\sqrt{\mu}$.
The consequence is immediate: with $\varepsilon \propto \mu$,

$$\alpha_c = \frac{\varepsilon}{\ell_d} \;\propto\; \mu^{1-0.52}
= \mu^{0.48},$$

and the measured ceiling exponent at $j = 1$ — the one budget level free of
censoring — is $+0.417$. The ceiling does not transfer at the derived rate
because one of the quantities the derivation holds fixed is moving at half
the rate of the tolerance.

The drift is not explained by the obvious candidates. Correlations across
problems between the $\ell_d$ drift and the drift of $\|x\|$, of
$\alpha_\star$, and of $\alpha_{\mathrm{wall}}$ are $+0.07$, $+0.23$ and
$+0.19$ respectively — none of them. That the exponent is close to $1/2$ and
nearly the same on every problem argues for a common mechanism rather than a
problem-specific one, and finding it is the most interesting open question
this section raises. What can be said is that $\ell_d = c_u u \|F\|\|x\|$
has only two moving parts, and neither alone accounts for it.

For an implementation the practical form is

$$\mathrm{hi}' \approx \Big(\frac{\mu'}{\mu}\Big)^{1/2}\,\mathrm{hi},$$

an exponent near one half where the derivation gives one. This is a fitted
rule and is marked as such; it is stated because it is what the measurements
support, not because it is derived.

### 26.4 What is not established

The ceiling family beyond $j = 1$ is **unmeasured**, not measured and
refuted. The $j = 2$ top sits at the last point of the $\alpha$ sweep in
three-quarters of the cases examined and the $j = 3$ top almost always, so
the observed tops are the end of the grid rather than the end of the window.
Censoring compresses tops toward a constant and biases the fitted exponent
toward zero, which is consistent with the negative values obtained; no
conclusion about the graded ceiling's transfer should be drawn from them.

Two experiments would settle what this data cannot. Extending the $\alpha$
sweep several decades above its present top would uncensor the $j \ge 2$
windows — and the point at which the factorization fails is itself the wall,
worth recording rather than avoiding. And re-running a *fixed* iteration at
several forcing constants would hold the operator fixed while
$\varepsilon$ varies, separating the tolerance exponent from the barrier
drift that is presently collinear with it.

### 26.5 Two further observations

**Additivity holds.** §11 argues that the cost is $n + j$ rather than
$\sum_{\text{passes}}(1 + n_{\text{pass}})$, on the ground that a refinement
pass inherits an already-contracted $\ell_p$ and so needs no Krylov
iterations. In the runs examined the refinement Krylov count equals the pass
count in $99.5\%$ of some forty-seven thousand rows — that is, refinement
passes used *no* iterations, essentially always. The exceptions are worth
isolating, since the theory places them where $\lambda_Z$ is smallest.

**$\alpha_\star$ is a staleness detector, not only a probing floor.** It is
computable from two diagonals before any solve, and its position relative to
a carried bracket sorts three cases: below the bracket, the carried interval
is entirely on the ramp and usable; inside it, only the part above
$\alpha_\star$ can furnish a ceiling reading; above it, the operator has
moved enough that the carried bracket is stale. Since $\alpha_\star$ climbs
as $\lambda_{\max}(A)$ does while $\alpha_c$ falls with $\varepsilon$, the
two converge along the trajectory, and $\alpha_\star$ overtaking the ceiling
is the same event as the tolerance becoming unattainable (§16) — visible
several iterations in advance rather than at the point of failure.

## References

`[GG]` G. H. Golub and C. Greif, *On solving block-structured indefinite
linear systems*, SIAM J. Sci. Comput. **24** (2003), no. 6, 2076–2092.
DOI 10.1137/S1064827500375096.

`[GGV]` G. H. Golub, C. Greif, and J. M. Varah, *An algebraic analysis of a
block diagonal preconditioner for saddle point systems*, SIAM J. Matrix Anal.
Appl. **27** (2006), no. 3, 779–792. DOI 10.1137/04060679X. Their $B$ is the
transpose of the $B$ used here, their $W = \gamma I$ is the present $\alpha$,
and their $S(W)$ is this note's $S$ divided by $\alpha$.

`[GS]` N. I. M. Gould and V. Simoncini, *Spectral analysis of saddle point
matrices with indefinite leading blocks*, SIAM J. Matrix Anal. Appl. **31**
(2009), no. 3, 1152–1171. DOI 10.1137/080733413.

`[NW]` J. Nocedal and S. J. Wright, *Numerical Optimization*, Springer, New
York, 1999. **TODO** — Lem. 16.1 is cited as the source of the
nonsingularity condition, but is known to this note only through
[GG, p. 2077]; it has not been read in the original, and the numbering of
the 2006 second edition has not been checked against it.

`[PS]` F. Pacaud, S. Shin, A. Montoison, M. Schanen, and M. Anitescu,
*Condensed interior-point methods for scalable nonlinear programming on
GPUs*, arXiv:2405.14236v3 (2026). **TODO** — cited from the preprint; check
against the journal version for section and page numbering, and confirm the
statement on p. 20 that condition (a) implies $\|K\| = \Theta(\Xi)$, which
appears to be a typo for $\Theta(\Xi^{-1})$ given Thm. 2.

`[RC]` S. Regev, N.-Y. Chiang, E. Darve, C. G. Petra, M. A. Saunders,
K. Świrydowicz, and S. Peleš, *A hybrid direct-iterative method for solving
KKT linear systems*, arXiv:2110.03636v1 (2021). **TODO** — cited from the
preprint; check against the journal version for numbering, and note that
their $\tilde H$, $J$, $\gamma$ are the present $A$, $B$, $\alpha$, and that
their $S$ is $\gamma$ times this note's.

`[BO]` M. Benzi and M. A. Olshanskii, *An augmented Lagrangian-based approach
to the Oseen problem*, SIAM J. Sci. Comput. **28** (2006), no. 6, 2095–2113.
DOI 10.1137/050646421.
