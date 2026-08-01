# The augmentation window: a mathematical justification

This note concerns the solution of a single saddle-point (KKT) system by
augmentation and Krylov correction, inside a refinement loop, and the
choice of the augmentation parameter $\alpha$.

We identify the set of $\alpha$ for which the augmentation alone suffices —
for which the system is solved to tolerance by one backsolve, with no
Krylov iteration and no refinement pass — and show that this set is an
interval whose two endpoints are computable from quantities the solve
records, with no fitted constants. We then show that the graded family of
sets admitting at most $k$ iterations is a nested family sharing the same
upper endpoint, whose lower endpoints descend at a rate governed by the
spectrum.

The note is self-contained and purely mathematical: no measurements are
reported, and every constant is derived from a stated assumption. It
concerns one invocation of the solver. A computation requiring several
solves at a common $\alpha$ must intersect the intervals obtained for each;
that intersection is not treated here.

## 1. The solver

We solve

$$\begin{bmatrix} A & -B^\top \\ B & 0 \end{bmatrix}
\begin{bmatrix} x \\ y \end{bmatrix} =
\begin{bmatrix} f \\ g \end{bmatrix}$$

with $A$ symmetric positive semidefinite and $B$ of full row rank on the
relevant subspace. The system is indefinite and is not factored directly.
For a chosen $\alpha > 0$ (write $\beta = 1/\alpha$) form the **augmented
operator**

$$F = \beta A + B^\top B,$$

symmetric positive definite, admitting a Cholesky factorization.
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
   We write $n$ for the number of iterations performed.
3. **Recovery.** Assemble $x = x_0 + \delta x$, let $r = g - Bx$, and set
   $$y = y_0 + \alpha(\delta y + r).$$

**The stopping test.** The Krylov iteration is *entered* only if its
stopping criterion fails at the initial point:

> **(T1)** The iteration tests $\|g - Bx\| \le \varepsilon$ before its
> first step; if the test passes, $n = 0$.

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

## 2. The window

**Definition.** For $k \ge 0$ let

$$W_k = \{\alpha > 0 : n(\alpha) \le k \text{ and no refinement pass fires}\},$$

and call $W_0$ — the set of $\alpha$ at which one backsolve suffices — **the
window**.

This definition is absolute, not relative: it does not refer to a minimum
attainable cost, so it requires no prior knowledge of the cost curve.
Sections 3–7 determine $W_0$; section 8 treats $W_k$ for $k \ge 1$.

By (T1) and (T2), $\alpha \in W_0$ iff both entry conditions hold at the
initial backsolve:

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
splitting. Directions with $\alpha\sigma_i^2 \gg 1$ have $\mu_i \approx 1$
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

> **Assumption A (self-similarity).** *Over the range of $\alpha$ under
> consideration the gap's mass lies on resolved directions:*
> $\|P\delta\| \approx \|\delta\|$ *and* $W$ *is negligible.*

Under A,

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
residual can never fall faster than inversely with $\alpha$**, and A is the
statement $\theta \approx 1$, i.e. that the residual's weight sits on
resolved directions. As $\alpha$ decreases the boundary in (4) sweeps into
populated spectrum, $\theta$ decays, and the proportionality (6)
overstates the rate at which the residual falls.

## 5. The primal endpoint

By (T1) the primal condition of (2) is $\|r_0(\alpha)\| \le \varepsilon$.
Under A, $\|r_0\|$ is strictly decreasing in $\alpha$, so the condition
defines a half-line, whose endpoint is where the residual meets the
tolerance:

$$\boxed{\ \alpha \ge \alpha_p = \frac{\ell_p}{\varepsilon}
= \frac{\alpha\|r_0(\alpha)\|}{\varepsilon}\ } \qquad (8)$$

the second form expressing the endpoint through the recorded residual at
any $\alpha$ where A holds. Equivalently, in ratio form,

$$\frac{\alpha}{\alpha_p} = \frac{\varepsilon}{\|r_0(\alpha)\|}. \qquad (9)$$

**Remark (no free constant).** The number 1 implicit in (8) — the endpoint
is where the residual equals *one* tolerance — is (T1) read literally, not
a calibration: the iteration is entered exactly when the initial residual
fails the test it is given. Under a variant that tested a multiple of
$\varepsilon$, that multiple would appear in (8) unchanged in form.

**Remark (one-sidedness of the A-error).** By (7), $\theta \le 1$, so below
$\alpha_p$ the true residual is never smaller than (6) predicts. The
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

**Proposition 2 (exactness of the dual row).** *In exact arithmetic the
assembled pair of §1 — $x = x_0 + \delta x$, $y = y_0 + \alpha(\delta y +
r)$ with $r = g - Bx$ — satisfies the first block row identically,*
$$f - Ax + B^\top y = 0,$$
*for every $n \ge 0$, including $n = 0$.*

*Proof.* The defining equation of $x_0$ gives
$\beta(Ax_0 - f - B^\top y_0) = B^\top(g - Bx_0)$, i.e.
$f - Ax_0 + B^\top y_0 = -\alpha B^\top r_0$. The invariant
$F\delta x = B^\top\delta y$ gives
$\beta A \delta x = B^\top \delta y - B^\top B \delta x$, i.e.
$A\delta x = \alpha B^\top(\delta y - B \delta x)$. With
$r = r_0 - B\delta x$,

$$f - Ax + B^\top y
= -\alpha B^\top r_0 - \alpha B^\top(\delta y - B\delta x)
+ \alpha B^\top \delta y + \alpha B^\top r
= \alpha B^\top(r - r_0 + B\delta x) = 0. \qquad\square$$

At $n = 0$ the recovery $y = y_0 + \alpha r_0$ alone performs the
cancellation. Consequently the computed $s_0$ consists entirely of
accumulated rounding error: every exact-arithmetic contribution vanishes by
construction. This licenses treating it as a sensor of attainable accuracy.

## 7. Conditioning, attainability, and the dual endpoint

> **Assumption B (two-scale spectrum).** *$B$ has a nontrivial near-kernel
> on which $A$ acts with scale $h > 0$, while the top of $F$'s spectrum is
> pinned at $\sigma_{\max}^2$ independently of $\alpha$.*

On the near-kernel the only support is $\beta A$, so the corresponding
eigenvalues of $F$ scale as $h/\alpha$, and

$$\kappa(F) \approx \frac{\sigma_{\max}^2}{h/\alpha}
= \frac{\alpha \sigma_{\max}^2}{h} \propto \alpha. \qquad (10)$$

This is the mirror of (6): raising $\alpha$ improves the primal residual
and degrades the conditioning, each in direct proportion.

> **Assumption C (attainability carries the dual residual).** *Above the
> rounding level of the residual evaluation itself, the computed dual
> residual is carried by the attainability term of the factorization:*
> $$\|s_0\| \approx c_u u \Gamma \kappa(F) = \ell_d \alpha,$$
> *with $u$ the unit roundoff, $\Gamma$ a norm scale, $c_u$ of order one,
> and $\ell_d = c_u u \Gamma \sigma_{\max}^2 / h$ independent of $\alpha$.*

Because the computed quantity cannot fall below the evaluation's own
rounding level, what is available is

$$\|s_0\| = \max(\ell_d \alpha,\ \text{evaluation floor}), \qquad (11)$$

a proportional ramp above a constant floor.

By (T2) the dual condition of (2) is $\|s_0\| \le \varepsilon$, which by
(11) is a half-line bounded above:

$$\boxed{\ \alpha \le \alpha_d = \frac{\varepsilon}{\ell_d}
= \frac{\alpha\varepsilon}{\|s_0(\alpha)\|}\ } \qquad (12)$$

again with the second form expressing the endpoint through the recorded
residual. In ratio form,

$$\frac{\alpha_d}{\alpha} = \frac{\varepsilon}{\|s_0(\alpha)\|}. \qquad (13)$$

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
below the ramp, so (12) never overstates $\alpha_d$: when the evaluation
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

$$\frac{\alpha_d}{\alpha_p} = \frac{\varepsilon^2}{\ell_p \ell_d}. \qquad (16)$$

*Proof.* By §5 the primal condition holds exactly on $[\alpha_p, \infty)$
and by §7 the dual condition on $(0, \alpha_d]$; $W_0$ is their
intersection, which is $[\alpha_p, \alpha_d]$ when $\alpha_p \le \alpha_d$
and empty otherwise. Substituting (8) and (12) gives (15) and (16).
$\square$

Contiguity is thus derived, not assumed: it follows from the monotonicity
of the two residuals in $\alpha$, which in turn follows from (6) and (10).
Both quantities in (15) are properties of the solve — $\ell_p$ the size of
the multiplier gap the solve must close, $\ell_d$ the amplification the
factorization contributes per unit $\alpha$ — and both are measured by the
recorded residuals through (8) and (12).

**Precision exhaustion.** When (15) fails there is no $\alpha$ at which the
augmentation alone suffices: at every $\alpha$ either the backsolve leaves
too large a constraint residual or the factorization cannot deliver a dual
residual within tolerance. At equality the two endpoints coincide, and at
that point both residuals equal $\varepsilon$. Since $\ell_p$ and $\ell_d$
are properties of the problem and the arithmetic while $\varepsilon$ is
imposed, (15) states the precision available to the tolerance demanded: the
window closes from both ends as $\varepsilon$ tightens, and the quantity
$\ell_p \ell_d / \varepsilon^2$ is the natural measure of how close the
configuration is to exhaustion.

**Asymmetry of consequences.** The two endpoints fail differently. Below
$\alpha_p$ the solve is still correct: the iteration is entered and
converges, at a cost quantified in §9. Above $\alpha_d$ the refinement loop
fires, and by Assumption C its passes face the same attainability limit as
the base solve; when that limit exceeds $\varepsilon$ the passes cannot
reach the tolerance at all. Leaving the window downward buys iterations;
leaving it upward risks non-convergence. This asymmetry, not any difference
in the sharpness of (8) and (12), is what distinguishes the two endpoints in
use.

## 9. The graded family

For $k \ge 1$, $W_k$ relaxes only the primal condition: the dual condition
is unaffected by the iteration count, so all $W_k$ share the endpoint
$\alpha_d$. Let

$$R_k(\alpha) = \min_{\substack{p \in \Pi_k \\ p(0) = 1}}
\frac{\|p(S) r_0\|}{\|r_0\|} \qquad (17)$$

be the optimal $k$-step residual reduction of the Krylov iteration on $S$ —
the standard characterization of CRAIG/CG as a polynomial minimization over
the residual's spectral measure. Then $n(\alpha) \le k$ iff
$R_k(\alpha)\|r_0(\alpha)\| \le \varepsilon$, so

$$W_k = [\alpha_k, \alpha_d], \qquad
R_k(\alpha_k) \|r_0(\alpha_k)\| = \varepsilon, \qquad (18)$$

and $W_0 \subseteq W_1 \subseteq W_2 \subseteq \cdots$, a nested family
with a common upper endpoint.

**The rate of descent.** Under A the spectrum of $S$ concentrates near 1
with spread $O(1/\alpha)$: by (4), $1 - \mu_i = 1/(1+\alpha\sigma_i^2)$. A
polynomial of degree $k$ with $p(0) = 1$ can be chosen to annihilate the
first $k$ moments of that spread — e.g. $p(t) = (1-t)^k$ gives
$\|p(S)r_0\|/\|r_0\| \le \max_i (1-\mu_i)^k = O(\alpha^{-k})$ — so

$$R_k(\alpha) = \frac{c_k}{\alpha^{k}}\,(1 + o(1)), \qquad (19)$$

with $c_k$ determined by the spectral measure and independent of $\alpha$.
Combining (19) with (6) in (18),

$$\boxed{\ \alpha_k = \Big(\frac{\ell_p c_k}{\varepsilon}\Big)^{1/(k+1)}\ }
\qquad (20)$$

so that, writing $\Lambda = \ln(\ell_p/\varepsilon)$ and suppressing the
$c_k$ (which vary slowly with $k$ relative to the exponent),

$$\ln \alpha_k \approx \frac{\Lambda}{k+1},
\qquad
\ln \alpha_{k} - \ln \alpha_{k+1} \approx \frac{\Lambda}{(k+1)(k+2)}.
\qquad (21)$$

The levels of the graded family are therefore *not* equally spaced: each
successive iteration buys a smaller extension of the window than the last,
the increments falling as $1/k^2$. The first iteration is worth the most —
it halves $\ln\alpha_p$ — and the family accumulates towards $\alpha = 1$
in the limit of large $k$.

**Reading.** Equations (20)–(21) are where the problem's spectral content
enters. The window $W_0$ requires none: its endpoints are fixed by the two
tests and the two proportionalities. Every level below it is spectral,
through the constants $c_k$ of (19), which are properties of the residual's
spectral measure. Section 10 shows that the iteration itself determines
these constants.

**Validity.** The asymptotic (19) rests on A, i.e. on the spectrum
concentrating near 1, and on $k$ small relative to the number of distinct
clusters in the measure. It degrades for dense spectra, where no low-degree
polynomial is uniformly small on the support, and for $k$ large enough that
the polynomial must resolve individual eigenvalues; in the extreme case of
a spectrum with no clustering, $R_k$ falls only as fast as the classical
$\kappa(S)$-dependent bound permits, and the graded family compresses
accordingly.

## 10. The spectral content of the cost

The constants $c_k$ of (19) are not extra data: the iteration that produced
$r_0$ computes a finite representation of the measure that determines them.

**The quadrature.** Over $k$ iterations the Golub–Kahan recursion assembles
the Lanczos tridiagonal $T_k$ of $S$ with starting vector $r_0/\|r_0\|$.
Its eigenpairs give

$$(\theta_j,\ w_j), \qquad w_j = \|r_0\|^2 (e_1^\top v_j)^2,
\qquad j = 1,\dots,k, \qquad (22)$$

which are the nodes and weights of the $k$-point Gaussian quadrature of the
residual's spectral measure on $S$. Since (17) depends on $r_0$ only
through that measure, and since the Krylov space of dimension $k$ is
determined by its first $2k$ moments — exactly those the quadrature
integrates exactly — the pairs (22) reproduce the iteration's own residual
history for the first $k$ steps in exact arithmetic. The representation is
lossless at $\alpha$, not an approximation of it.

**Transport.** The dependence on $\alpha$ in (4) is a bijection between
$\mu_i$ and $\sigma_i^2$, so the nodes may be stripped of the parameter and
re-dressed at any $\alpha'$:

$$\hat\sigma_j^2 = \frac{\theta_j}{\alpha(1-\theta_j)},
\qquad
\hat\delta_j^2 = w_j\Big(\frac{\alpha}{\theta_j}\Big)^2, \qquad (23)$$

$$\hat\mu_j(\alpha') = \frac{\alpha'\hat\sigma_j^2}{1+\alpha'\hat\sigma_j^2},
\qquad
\hat w_j(\alpha') = \Big(\frac{\hat\mu_j}{\alpha'}\Big)^2 \hat\delta_j^2,
\qquad (24)$$

(24) being identity (3) applied per node. Evaluating (17) on the
transported pairs gives $R_k(\alpha')$, hence $c_k$ and, through (20), the
levels $\alpha_k$ — all from the single solve, with no assumption of
self-similarity spent on the transported directions, since (24) carries
each node's $\alpha$-dependence explicitly.

**What transport does and does not determine.** The representation (22) is
exact at $\alpha$ for the first $k$ steps. Transported, it is not: the map
(24) reweights *within* each node, but a node standing for a cluster of
distinct $\sigma^2$ transports as a single point, and the reweighting of
that cluster is only as accurate as the lumping. Moreover the measure is
known only where the Krylov space explored it — mass on directions the
starting vector did not excite, or that $k$ steps did not reach, is absent
from (22) altogether and cannot be recovered by any map. Both effects
concern spectrum at the low end, which is where descending $\alpha$ moves
weight; hence the transported estimate of $R_k(\alpha')$ for
$\alpha' \ll \alpha$ is systematically optimistic, and the levels (20)
computed from it are lower bounds on the true $\alpha_k$. This is a
statement about the representation, not about any particular
implementation: a $k$-node quadrature determines the measure's first $2k$
moments and nothing more.

## 11. Assumptions

| | statement | role | failure |
|---|---|---|---|
| T1 | the Krylov stopping test is evaluated at entry | fixes the primal endpoint at $\|r_0\| = \varepsilon$ | a test on a multiple of $\varepsilon$ carries that multiple into (8) |
| T2 | the refinement test admits no borrowing between components | fixes the dual endpoint at $\|s_0\| = \varepsilon$ | a combined test reduces the dual allowance by the other components' share |
| A | self-similarity: the gap's mass lies on resolved directions | gives $\|r_0\| \propto 1/\alpha$; quantified by $\theta = \langle\mu\rangle_w \le 1$ | degrades as $\alpha$ falls and the tail populates; error one-signed, understating depth below $\alpha_p$ |
| B | two-scale spectrum: $A$-supported near-kernel of $B$, $\alpha$-independent top | gives $\kappa(F) \propto \alpha$ | without such a near-kernel the upper endpoint has no mechanism |
| C | attainability carries the dual residual above the evaluation floor | gives $\|s_0\| \propto \alpha$, hence (12) | below the floor the computed value is a bound, not an equality; the estimate remains one-sided |
| — | Proposition 1 (identity), Proposition 2 (dual exactness) | exact algebra; no assumption | — |

Assumptions A and C are mirror images: A asserts that a product
$\alpha\|r_0\|$ is $\alpha$-free, C that a quotient $\|s_0\|/\alpha$ is.
Each is an assertion that one scaling dominates over the range considered,
and each fails in the regime the *other* endpoint governs — A as $\alpha$
falls, C as it rises — so that neither is asked to hold where it is least
justified.

## 12. Summary

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
Krylov iterations extends the interval downward only, to
$\alpha_k = (\ell_p c_k/\varepsilon)^{1/(k+1)}$, the levels crowding
together as $1/k^2$; the constants $c_k$ are properties of the residual's
spectral measure, of which the iteration computes an exact $k$-node
representation at the parameter it ran at, and a systematically optimistic
one below it.
