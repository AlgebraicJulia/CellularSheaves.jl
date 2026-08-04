# Pricing augmentation

## 0. Contract and notation

This note prices the augmentation parameter of a saddle-point solve: from
the records of one completed solve it estimates, with stated error, what
every other value of the parameter would have cost, and answers where that
cost is least. The estimate, not any single recommendation drawn from it,
is the object. Operationally, the parameter is *given* and the solve
at it is mandatory: this note prices what that solve pays, reads what
it records, and returns the estimate for the choice it could not
make — the next one. A solve's record is a set of readings at its
parameter, the **anchor**; the *designated solve* is the one the
caller ran, its anchor given and its readings free.

**The contract.** Every claim below is exactly one of four things:
*derived*, from the solver's definition and named hypotheses;
*assumed*, under a registered name **(A-·)**, collected with coverage and
breadth in §17; or *deferred*, under a registered name **TODO-·**,
collected with resolution criteria in §18; or **calibrated** — a
tuning constant admitted only with recorded provenance and
sensitivity, collected in §17's calibration block. The four
categories are closed: where the note would otherwise say that
something is observed to work, it names an assumption, a TODO, or a
calibration row. Empirical support beyond a row's provenance lives
in the companion technical report; the body cites verification
outcomes by tag, never by table.

Two local disciplines follow from past failures and are enforced at every
use, not only at first statement. A law with branches is cited together
with its active branch and the condition selecting it. A wall quantity is
cited together with its provenance: a *reading* (a measurement of the
quantity, admissible as a value) or a *bound* (a one-sided certificate,
admissible only as a constraint); the two are never interchangeable, and
the type system of §13 exists to make their confusion inexpressible.

**The system and the solver's letters.** We solve

$$\begin{bmatrix} A & -B^\top \\ B & 0 \end{bmatrix}
\begin{bmatrix} x \\ y \end{bmatrix} =
\begin{bmatrix} f \\ g \end{bmatrix}$$

with $A$ symmetric positive semidefinite of dimension $n$ and $B$ of full
row rank, $p \times n$. For $\alpha > 0$, $\beta = 1/\alpha$, the
augmented operator is $F = \beta A + B^\top B$, symmetric positive
definite; $\rho \ge 0$ — the shift enters through the factorization
alone and its perturbation is left to the refinement loop, priced in §15. The solver, its tests, and its halting
policy are fixed in §1. The bare letter $n$ names the dimension only
where the operator is described; as a count it is always subscripted, and
the two never share an expression.

**Model of computation.** Two conventions, stated once. *Arithmetic*
is the standard model of floating point (Higham's convention): machine
numbers with unit roundoff $u$, every elementary operation returning
$\mathrm{fl}(x \circ y) = (x \circ y)(1+\delta)$ with
$|\delta| \le u$. Every rounding statement in this note is quantified
adversarially over all $\delta$-patterns within that bound — which is
what makes the results independent of instruction order, fusion, and
code path: a particular machine's run is one realization, and the
note's inequalities bracket all of them. The application of $F^{-1}$
is an *oracle*: it returns the exact solve perturbed at a relative
error level $\eta$ — of order $u$ for a working-precision
factorization, larger under reduced precision or boundary-degraded
conditioning — and the factorization's interior stays behind that
boundary. *Cost* is oracle counting: programs act on machine numbers,
comparisons are exact and free, matrix–vector products are free by the
currency below, and the cost of a run is its number of oracle calls;
wall-clock time, memory traffic, and parallelism are outside the
model. The model supplies the *forms* of the laws — floors, levels,
poles, the $\eta$-proportional deposit — with constants pessimistic
and unpinned; the record supplies the realized constants, and no
returned quantity evaluates a model constant. Transfer across machines
is transfer of forms and brackets, each machine's $(u, \eta)$ and
record supplying its own numbers. Within the model the two arms'
platonic objects differ in type: the floor's exact-arithmetic
staircase is point-valued, floating point jittering around it; the
ceiling's are bracket-valued under the adversarial quantifier —
extrema over admissible histories, a realized machine one
deterministic sample inside — so at the wall the certified content is
natively a black/white pair, and every interior guess is a
measurement. Three grades recur for every counted quantity: the
measured value, the law's, and the smooth relaxation's; the signed
census bridges the first pair, corner displacement the second.

**Counts and cost.** The unit of cost is one application of $F^{-1}$: the
initial backsolve is one, each refinement pass is one, and each Krylov
iteration — preconditioned by $F^{-1}$ — is one. Writing $j$ for the
number of refinement passes, $n_0$ for the iteration count of the base
correction and $n_i$ for that of pass $i$,

$$\mathrm{cost} \;=\; (1 + j) \;+\; n_0 \;+\; \sum_{i=1}^{j} n_i ,$$

an integer, and the objective of the note. Its two-hyperbola relaxation
$C$ (§8) is real-valued and carries the guarantees.

**The estimate spends nothing.** Every returned quantity is computed
from the record and the already-held factors at zero additional
applications; the estimator holds no budget of its own. Extensions
that purchase information — a refinement pass minting a contraction
reading, any auxiliary spectral iteration — are a caller's options
and are outside this note's problem, which is exactly: the best
estimate at zero marginal cost. Where the record is silent the
output says so; silence priced honestly is the product, not a
deficiency a purchase would cure. One purchase is admitted, alone of
its kind: under regularization the recorded entry residual is one reading
in two unknowns — the decaying scale and the $\rho$-born bias — and no
manipulation of a single record separates them. One application of
$F^{-1}$ does. The exception is warranted by **identifiability**, not by
convenience: it is not that the free route is uninformative, but that
there is no free route. It fires on a zero-cost gate (§15) and returns a
bound where it does not.

**Logarithms and gaps.** All logarithms are natural; a *gap* is a log of
a ratio, measured in nats, and a *decade* is $\ln 10$ of them.
Throughout, $t = \log\alpha$ and $e = \log\varepsilon$.

**Parameters are tolerance-free.** Every constant of the cost model is a
property of the triple (operator, data, arithmetic) alone:
$\lambda_p = \log \ell_p$ and $\lambda_d = \log \ell_d$ (the two residual
scales, §4–§6), $t_\infty$ (the floor's accumulation point, §5), $T$ (the
realized wall, §7, always tagged reading-or-bound), the entry-gap
constants (§7), and $\varphi$ (the limiting accuracy, §7). The tolerance
enters every formula explicitly through $e$ and nowhere implicitly, and
that is the point of the discipline rather than a tidiness convention:
the parameters are measured at one $\varepsilon$ and the curves are
evaluable at any. Transfer in $\varepsilon$ is therefore available and
is used. What degrades with distance is the *box*, not the formula — the
harvest is read over the levels a solve explored, and carrying it far
past them is extrapolation, one-signed and governed by the same depth
discipline as transport in $\alpha$ (**TODO-G**).

**Fixed symbols.** $Z$ is a basis for $\ker B$;
$\lambda_Z = \lambda_{\min}(Z^\top A Z)$; $u$ is the unit roundoff;
$\tau$ is the stall threshold of test (T3); $c_k$ are the harvest
coefficients of §5; $r_0 = g - Bx_0$ and $s_0$ are the initial constraint
and dual residuals; $\varepsilon$ is the absolute tolerance of test (T1).
Assumption names in force are introduced where first used and indexed in
§17; the eleven TODO names are indexed in §18.
## 1. The solver

Multiplying the first block row by $\beta$ and adding $B^\top$ times the
second yields the equivalent single equation

$$F x \;=\; \beta f + B^\top(g + \beta y), \qquad (1)$$

with two unknowns, solved iteratively in $y$:

1. **Initial backsolve.** Given a starting multiplier $y_0$, solve
   $F x_0 = \beta f + B^\top(g + \beta y_0)$ and form the initial
   constraint residual $r_0 = g - B x_0$.
2. **Krylov correction.** A CRAIG iteration preconditioned by $F^{-1}$
   solves $B\,\delta x = r_0$ for the minimum-$F$-norm correction,
   maintaining the invariant $F\,\delta x = B^\top \delta y$ at every
   iterate and stopping when the constraint residual — **evaluated
   directly at the current iterate** — is within the invocation's
   tolerance $\varepsilon_j$.
3. **Recovery.** Assemble $x = x_0 + \delta x$, let $r = g - Bx$, and set
   $y = y_0 + \alpha(\delta y + r)$.

A **refinement pass** repeats the three steps with the current residuals
as data; $j$ counts these passes, and the base invocation is pass $0$.

**The interior discipline.** The base invocation and the exit run at
the caller's $\varepsilon$; an *interior* pass runs at
$\varepsilon_j = \max(\varepsilon,\ \theta\,\|r_0\|)$ with the named
constant $\theta = \tfrac1{10}$ — each pass reduces its own entry
residual tenfold and no further, the proportional discipline of the
inexact-multiplier literature, covered there for $\theta$ below a
threshold scaling as the inverse square root of the Schur
conditioning; (T3) is the net where that condition thins. When (T2)
first passes, one final invocation at full $\varepsilon$ enforces the
constraint, and **both residuals are verified directly**; a solve that
cannot meet them routes to (T3)'s halt. Exits are typed —
*reached-tolerance*, *floor* (the directly evaluated residual detaches
above the recursion's), or *declared error* — the type is recorded
(§16), and no other interior halting heuristic is admitted.

**The tests.** Four, and together with the threshold $\tau$ they are the
*halting policy*; every set and function in this note is conditional on
the policy as stated, a dependence made explicit where it acts (§7, §14).

> **(T1)** The Krylov iteration is entered only if
> $\|g - Bx\| \le \varepsilon_j$ fails at its initial point — the
> invocation's own tolerance, per the interior discipline; if the test
> passes, the invocation's count is $0$.

> **(T2)** A further pass is taken while the dual residual $s_j$ of the
> assembled iterate exceeds $\varepsilon$; the solve *converges* when it
> does not.

> **(T3)** The solve *stalls*, and halts, if a pass contracts the dual
> residual by less than the threshold: $s_{j+1} > \tau\, s_j$.

> **(T4)** *(assumption (A-term))* An entered Krylov invocation exits
> only by meeting its test or by declared error; silent returns are
> outside scope. Each invocation records a met-flag (§16), so (A-term)
> is checkable per solve.

(T4) is the load-bearing addition: it makes every recorded count a
*test-meeting* count, from which §7 derives that refinement iterations
are dual-born and that the pass-cost formulas price the right object.

**The tolerance is absolute**, and that is a restriction on the solves
addressed rather than a normalization: under a relative test the entry
check of (T1) can never pass, so no window exists and the accounting
below does not apply.

**Semidefiniteness of $A$** does quiet work: it makes $F$ positive
definite for every $\alpha > 0$, so the parameter is bounded below by
cost alone and by nothing else. An indefinite leading block imposes a
hard admissibility threshold [RC, Thm. 1] and is outside scope.

$\rho \ge 0$; the shift enters the factorization alone, and §15
prices its consequences: an additive dual scale, an entry bias, and a
ceiling arm that steepens toward an unmoved pole.
## 2. The cost

Fix the instance — operator, data, tolerance, policy — and let $\alpha$
vary. A solve that converges under the policy pays

$$\mathrm{cost}(\alpha) \;=\; (1 + j) \;+\; n_0 \;+\; \sum_{i=1}^{j} n_i
\qquad (2)$$

applications of $F^{-1}$; where the solve stalls, or where
$\varepsilon$ is unattainable (§7), set $\mathrm{cost}(\alpha) = +\infty$.
The convention is not decoration: cost so extended is the object every
later section estimates, bounds, and minimizes, and $+\infty$ is how the
estimate will say *infeasible* without a side channel (§13).

Three facts about (2) organize the note. It is a **staircase**: each
count is integer-valued and piecewise constant in $\alpha$, so the
minimizers form flat plateaus and no point inside a plateau is
distinguishable from its neighbors — the natural answer to "where is it
least" is an interval. It is **policy-conditional**: (T3)'s threshold
decides which low-$\alpha$ solves halt rather than pay, so the domain
and the values of (2) both move with $\tau$ (§7). And it has **two
arms**: $n_0$ grows without bound as $\alpha$ falls (§§4–5, the floor)
while the pass count $j$ grows as $\alpha$ approaches the wall from
below (§§6–7, the ceiling); everything expensive happens at the ends,
and the cheap middle is where the two arms are simultaneously small.

The **window**

$$W_0 = \{\alpha : \mathrm{cost}(\alpha) = 1\}$$

is the unit level set — one backsolve, entry test passed, no passes —
and §3 shows it is an interval with a computable nonemptiness criterion.
The old habit of treating $W_0$ as the subject survives here as a level
set of the actual subject.

**The objective.** From the records of one completed solve (§16),
produce an estimate of $\mathrm{cost}(\cdot)$ over the whole parameter
range, with stated error and stated domain of validity (§13), from which
the least-cost plateau and its budget are read as a query (§14). The
estimate is what a solve can know; the staircase is what it pays; the
distance between the two is priced, not ignored (§8, Appendix A).
## 3. The window

$\mathrm{cost}(\alpha) = 1$ requires exactly that both tests pass at the
base invocation: $\|r_0\| \le \varepsilon$ by (T1), so the correction is
never entered, and $s_0 \le \varepsilon$ by (T2), so no pass follows.
Sections 4–7 derive the two scale laws

$$\|r_0\| = \frac{\ell_p}{\alpha} \quad (§4\text{–}§5), \qquad
s_0 = \ell_d\,\alpha \quad (§6),$$

each exact in its stated regime and each monotone in $\alpha$, so each
condition holds on a half-line and their intersection is the interval

$$W_0 \;=\; \Big[\frac{\ell_p}{\varepsilon},\;
\frac{\varepsilon}{\ell_d}\Big], \qquad
W_0 \neq \emptyset \iff \varepsilon^2 \ge \ell_p\,\ell_d . \qquad (3)$$

In log coordinates the statement is cleaner and carries its own
sensitivity analysis: $t \in [\lambda_p - e,\; e - \lambda_d]$, an
interval of width

$$w(e) \;=\; 2\big(e - e_{\min}\big), \qquad
e_{\min} = \tfrac{1}{2}(\lambda_p + \lambda_d),$$

so the window loses two nats of width per nat of tightening, one from
each end, and closes at $\varepsilon_{\min} = \sqrt{\ell_p \ell_d}$ —
the geometric mean of the two residual scales, at the parameter value
$\alpha_g = \sqrt{\ell_p/\ell_d}$ where the closing ends meet.

Below $W_0$ the entry test fails and iterations are bought (§5);
above it the dual test fails and passes are bought (§7); the two
purchases are the floor and ceiling arms of the cost, and everything
after this section is an account of their prices. When (3) fails
outright the window is empty at this tolerance and the question is
purely one of pricing — the case the old treatment reached last and
this one treats as generic. Three one-line motions of the window, for
later use: tightening $\varepsilon$ moves both ends inward, one nat
each per nat; rescaling the data is indistinguishable from moving
$\varepsilon$; improving the start moves the left end alone (§4).
## 4. The initial residual

**Proposition 1.** *Let $(x^*, \lambda^*)$ solve the two-block system
exactly. Then, in exact arithmetic,*

$$r_0 = -\frac{1}{\alpha}\, S\, (y_0 - \lambda^*),
\qquad S = B F^{-1} B^\top. \qquad (4)$$

*Proof.* The exact solution satisfies (1) with $y = \lambda^*$:
$F x^* = \beta f + B^\top(g + \beta\lambda^*)$. Subtracting the defining
equation of $x_0$ gives
$x_0 - x^* = \beta F^{-1} B^\top (y_0 - \lambda^*)$; apply $B$ and use
$B x^* = g$. $\square$

$S$ is symmetric positive semidefinite with eigenvalues

$$\mu_i = \frac{\alpha\sigma_i^2}{1 + \alpha\sigma_i^2} \in (0,1),
\qquad \sigma_i^2 = \lambda_i\big(B A^{-1} B^\top\big), \qquad (5)$$

the identification being exact at every $\alpha$, not asymptotic: the
augmented Schur complement is $S_\bullet(I + \alpha S_\bullet)^{-1}$ in
terms of the unaugmented one [GGV, Prop. 2.7], and (5) is that statement
term by term. Where $A$ is singular, $A^{-1}$ is defined through the
perturbation $A + \epsilon B^\top B$ and the limit [GGV, §2.2]; the
$\sigma_i$ are the limit's. Directions with $\alpha\sigma_i^2 \gg 1$ are
**resolved** by the augmentation ($\mu_i \approx 1$); those with
$\alpha\sigma_i^2 \ll 1$ are the **tail**.

**Two objects, one letter avoided.** $\sigma$ with an index — and
$\sigma_{\min}$ — always denotes these *generalized* values, solving
$B^\top B v = \sigma^2 A v$; they carry the graded families and the
augmentation contraction. The ordinary largest singular value of $B$ is
written $\|B\|$ and never $\sigma_{\max}$; it arises from $\|F\|$ and so
appears in the conditioning, the switch, and the dual scale. The two
coincide only for $A = I$.

Identity (4) supplies, with no assumption, the prefactor $1/\alpha$.
Expanding in the eigenbasis of $S$, with $\delta = y_0 - \lambda^*$,

$$\|r_0(\alpha)\|^2 = \frac{\|P\delta\|^2}{\alpha^2} + W(\alpha),
\qquad (6)$$

$P$ the projector onto resolved directions and
$W = \sum_{\alpha\sigma_i^2 < 1} \sigma_i^4 \delta_i^2$ the tail's
contribution: resolved mass scales inversely, tail mass contributes
a step-increasing term — each entry is $\alpha$-free, but the index
set grows as $\alpha$ falls. The transport signs need no more than
the cap this puts on the decay: every direction decays at most as
$1/\alpha$, so the local exponent satisfies $q \le 1$ everywhere,
and slope-one transport errs optimistic upward — assuming decay at
least as fast as truth, the endpoint placed low, the census's sign —
and conservative downward, assuming growth at least as fast as
truth. The exponent itself is not monotone in $\alpha$: tail weight
grows faster than tail value, so $q$ can dip before recovering,
which is why the signs hang on the cap and on nothing finer.

**The scaling assumption, in three tiers**, each implying the ones below,
the weakest being what derivations consume:

> **(A1) — used.** Between the anchor $\alpha$ and the endpoint sought,
> $\|r_0\| \propto 1/\alpha$.

> **(A2).** The same over the whole operating range.

> **(A3).** The gap's mass lies on resolved directions:
> $\|P\delta\| \approx \|\delta\|$, $W$ negligible; then the constant is
> named, $\|r_0\| = \|\delta\|/\alpha$.

Only (A1) enters the floor endpoint; (A3) is the spectral account of
*why* it holds and supplies

$$\|r_0(\alpha)\| = \frac{\ell_p}{\alpha}, \qquad \ell_p = \|\delta\| .
\qquad (7)$$

Any mechanism producing inverse scaling serves (A1) equally; the tiers
are registered separately in §17 because their coverage differs. Read
as a statement about starts, (7) says the primal scale *is* the start
quality — a cold start is the special case $y_0 = 0$ — and every
appearance of $\ell_p$ below is a statement about the multiplier gap.
## 5. The primal arm

**The elasticity.** Differentiating the eigen-expansion (6) in
$\log\alpha$, using $d\mu_i/d\log\alpha = \mu_i(1-\mu_i)$,

$$-\frac{d \log\|r_0\|}{d \log\alpha} \;=\; \langle \mu \rangle_w
\;=:\; \vartheta \;\in\; (0, 1], \qquad (8)$$

the mean of the spectrum of $S$ against the residual weights
$w_i$ — the squared components of $r_0$ in the eigenbasis. Two facts are
carried by (8) alone, with no assumption: $\vartheta > 0$ strictly, so
$\|r_0\|$ is strictly decreasing and the entry condition
$\|r_0\| \le \varepsilon$ holds on an interval unbounded above — the
*structure* of the condition is Proposition 1's; and $\vartheta \le 1$, so
below the endpoint the true residual is never smaller than the scaling
law predicts. (A1) is needed only to say *where* the endpoint falls:

$$\alpha \ge \alpha_0 = \frac{\ell_p}{\varepsilon}
= \frac{\alpha\,\|r_0(\alpha)\|}{\varepsilon}, \qquad (9)$$

the second form reading the endpoint from the recorded residual at any
$\alpha$ where (A1) holds between there and $\alpha_0$. The implicit
constant $1$ is (T1) read literally, not a calibration. In finite
precision the computed $\|r_0\|$ bottoms out at the rounding of forming
$g - Bx_0$; monotonicity statements in §§4–5 concern the exact
quantity. By the one-sidedness above, (9) errs, when (A1) degrades below
the anchor, by *understating* how far below the endpoint a small $\alpha$
lies — the first member of a signed census kept at the end of this
section.

**The graded family.** For $k \ge 1$ let $W_k$ relax the entry condition
to $n_0 \le k$; the dual condition is untouched, so all $W_k$ share the
upper endpoint. Every iterate's constraint residual is a polynomial in
$S$ applied to $r_0$, so the achieved $k$-step reduction obeys

$$R_k(\alpha) \;\ge\; R_k^{\mathrm{opt}}(\alpha)
= \min_{p \in \Pi_k,\, p(0)=1} \frac{\|p(S)r_0\|}{\|r_0\|}, \qquad (10)$$

with inequality in general: CRAIG minimizes in its own norm, not the
tested one. Under (A3) the discrepancy is $1 + O(1/\alpha)$, but
$R_k^{\mathrm{opt}}$ under-prices the work, an optimism recorded in the
census. By (T1) at each step, $W_k = [\alpha_k, \alpha_c]$ with
$R_k(\alpha_k)\|r_0(\alpha_k)\| = \varepsilon$, a nested family. The interval form presupposes the reduction's
monotonicity in $\alpha$ — exact for $\|r_0\|$ by (8), asymptotic
for $R_k$, tagged **(A-mon)**, generic within (11)'s regime.

The reduction has an exact characterization and a usable asymptotic.
Exactly, $R_k^{\mathrm{opt}}(\alpha)^2 = (\sum_{j\le k}
\varphi_j(0)^2)^{-1}$, the Christoffel function at the origin of the
normalized residual measure, $\{\varphi_j\}$ its orthonormal
polynomials. Asymptotically, on resolved-weighted measures (A3),

$$R_k(\alpha) = \frac{c_k}{\alpha^k}\,(1 + o(1)), \qquad
c_k = \|\pi_k^\tau\|_{\hat w} = \prod_{j \le k} e_j, \qquad (11)$$

the $L^2$-norm of the degree-$k$ monic orthogonal polynomial of the
**limiting measure** — nodes $\tau_i = 1/\sigma_i^2$, weights
$\propto \delta_i^2$ — equivalently the running product of that
measure's Jacobi off-diagonals. What the recursion computes at finite
$\alpha$ is the affine node variable $\tau_i = \alpha(1 - \mu_i)$; the
identification with $1/\sigma_i^2$ is its resolved limit, immaterial
inside the window and not outside it. $c_k$ carries the units of
$\alpha^k$. The regime of (11) is a labeled triple — *domain*
(resolution $\alpha\sigma_i^2 \ge 1$, a partition moving with
$\alpha$, implied by the window at tight tolerance), *assumption*
((A3), an instance property), *depth* ($k \ll N$) — with per-axis
verification: metered, trusted, visible. Combining (11) with (7),

$$\alpha_k = \Big(\frac{\ell_p\, c_k}{\varepsilon}\Big)^{1/(k+1)} .
\qquad (12)$$

**Termination, and the floor's two parameters.** The limiting measure
has $N$ distinct weighted nodes; for $k \ge N$ the monic polynomial
vanishes, $c_k = 0$, and $W_k = (0, \alpha_c]$: the family terminates
rather than accumulates. The termination binds the curves: the
lower curve's returned floor contribution is capped at level $N$, a
clamp on returned values (the convex machinery object is unclamped),
since beyond termination the relaxation's pole is an envelope fiction
the exact staircase does not share. (The completed-family view of
these levels, with closures in the Golub–Meurant tradition, derives
both the cap and the anchored reading; the note ships the endpoints.) $N$ is structured: the preconditioned spectrum
carries eigenvalue $1$ with multiplicity $n - p +
\operatorname{nullity}(A)$ [GGV, Thm. 2.5], so a rank-deficient leading
block — the hard case for factorization — lowers $N$ and is the easy
case here. For $k \ll N$ the discrete measure is governed by its
continuum envelope and $c_k^{1/k} \to \alpha_\infty =
\operatorname{cap}(E)$, the logarithmic capacity of the support
($\tfrac14(\sigma_-^{-2} - \sigma_+^{-2})$ for an interval). In the
note's coordinates, with $t_\infty = \log\alpha_\infty$ and
$M_f = \lambda_p - e - t_\infty$,

$$\log\alpha_k = t_\infty + \frac{M_f}{k+1}, \qquad
n_0(t) \;\text{interpolated by}\;\; \nu_0(t) = \frac{M_f}{t - t_\infty}
- 1, \qquad (13)$$

the levels descending geometrically toward $t_\infty$, increments
falling as $1/k^2$, each iteration buying less than the last.
The interpolation's precision is exactly this: $\nu_0$ passes through
the corner $(t_k, k)$ iff $\log c_k = k\,t_\infty$, so
$\lceil \nu_0 \rceil$ recovers the staircase of (12) *identically*
iff the harvest is exactly geometric; in general each boundary is
displaced by $\delta t_k = (k+1)^{-1} \sum_{j \le k} (\log e_j -
t_\infty)$ — the deviation partial sums, $(k+1)$-damped — and the two
functions differ, by one level, only on intervals of that width around
the displaced corners. Sharp drops of $c_k$ at structure-resolving
degrees (clusters) are the large-displacement case. The displacement's
sign is that of the running deviation sum — zero for the Chebyshev-U
weight, one-sided for single-interval regular measures, oscillating
for gapped supports — hence excluded from the signed census (which
requires a universal sign) but free per instance from the recorded
harvest. Equally free is the ledger identity
$t_j = t_\infty + (M_f + D_j)/(j+1)$ with
$D_j = \log c_j - j\,t_\infty$: read at the horizon it anchors the
far field at the deepest witnessed corner,
$M_f + D_K = (K+1)(t_K - t_\infty)$, the reading §14 adopts.
Formula (12)'s boundary values are ordered within (11)'s regime;
where they disorder, the asymptotic has been left and the monotone
envelope is the reading — the disorder itself a recorded regime
signal. (An earlier version gated transported upper-curve floors on a
calibrated window width; the gate is retired — measurement showed it
thresholding a quantity that does not govern the failure. What does is
depth below the anchor: **TODO-G**.) This relates the
interpolant to the *law's* staircase; the realized count relates to
(12) through the signed census below, and §8's corridor prices the
remaining smooth-versus-integer gap when $C$ consumes $\nu_0$. Equation
(13) is the floor arm of §8, its two parameters $\varepsilon$-free
exactly as §0 declares, $\varepsilon$ entering through $M_f$'s explicit
$-e$. The budget when the window is empty is immediate:
$k + 1 \ge M_f / (\log\alpha_c - t_\infty)$, informative only when
$\alpha_\infty < \alpha_c$; otherwise no practical iteration count
delivers a window and the regime is exact solution of the correction
system, which finite precision withholds.

**The signed census.** The approximations above do not share a sign,
and §17's breadth column consumes this list. Optimistic — placing
$\alpha_k$ *below* truth, overstating the family's reach — and
compounding, since they share a sign: (A1) degrading below the anchor
($\vartheta \le 1$); the norm substitution in (10); transport of the
measure to parameters below the anchor; and the single-scale substitution
$c_k \approx c_1^k$, since $c_k^{1/k}$ ascends from the standard
deviation toward the capacity. Conservative: bounding $c_k$ by an
explicit polynomial; and the $o(1)$ of (11) near $t_\infty$, where the
resolved substitution fails — growing exactly where the optimistic group
grows, so no net sign is asserted near the accumulation point. With
$c_k$ from the recursion, what remains is the optimistic group, and a
computed $\alpha_k$ is properly read as a lower bound on the level.

**The harvest.** Over $n_0$ iterations the recursion assembles the
Lanczos tridiagonal of $S$ with starting vector $r_0/\|r_0\|$; its
eigenpairs are the nodes and weights of the $n_0$-point Gaussian
quadrature of the residual measure — a *lossless* representation of the
first $n_0$ steps at the anchor, and the source of the $c_k$
at no cost beyond recording (§16). The family's constants and the
quadrature are one object seen twice. Nodes strip and re-dress across
$\alpha$ by the bijection (5); the use and screening of that transport
belong to §11.

Exactly at the anchor the constants admit a second, identity-level
reading: the recorded residual ratios give
$c_k = \alpha^k\,\|r_k\|/\|r_0\|$, exact by the iteration's own
recurrence, and the level at $k = n_0$ then reproduces the anchor
identically — a self-consistency any harvest computation must pass.
Where the recursion's census carries slack, this form carries none;
both corners may consume it, and the endpoint's remaining error
enters level $k$ with sensitivity exactly $1/(k+1)$, so the
certified levels narrow toward the anchor by construction.

One caution attaches to depth, and only at the node level: in finite
precision a deep tridiagonal loses orthogonality precisely as Ritz
values converge, and clustered eigenvalues acquire ghost copies —
hazardous to the strip-and-re-dress transport of individual nodes,
and to nothing else here. Moments, quadrature values, $\vartheta$, and
the $c_k$ see only a nearby, slightly smeared measure, and the
residual-ratio form consults no spectral identification at all.

Two cheap exact quantities close the arm. The elasticity is free
whenever any iteration ran: $\vartheta = 1/a_1$, the reciprocal of the
first CG step length, exactly. And one application of $S$ — a
triangular-solve pair and two multiplications by $B$ — yields
$c_1 = \alpha\,\|(S - \vartheta)r_0\|/\|r_0\|$ exactly, the measure's
standard deviation with no asymptotic; form the residual vector rather
than subtracting scalars, since inside the window both moments sit near
$1$. At $n_0 = 0$ the solve is spectrally blind — the window's defining
property — and a $\vartheta$ bought there validates (9) and corrects
nothing: by (T1) the ratio $\|r_0\|/\varepsilon \le 1$ there, and the
correction would move the floor in the unsafe direction. A quantity
adequate for judging a reading is not thereby adequate for moving it.
## 6. The dual arm

For a pair assembled naively through $F$, the first-row residual is not a
finite-precision quantity: the defining equation forces it to equal
$\alpha B^\top(g - Bx)$ identically. The recovery step exists to cancel
exactly this coupling.

**Proposition 2 (the dual row).** *In exact arithmetic the assembled
pair — $x = x_0 + \delta x$, $y = y_0 + \alpha(\delta y + r)$,
$r = g - Bx$ — satisfies $f - Ax + B^\top y = 0$ for every $n_0 \ge 0$,
including $n_0 = 0$.*

*Proof.* The defining equation of $x_0$ gives
$f - Ax_0 + B^\top y_0 = -\alpha B^\top r_0$. The invariant
$F\delta x = B^\top\delta y$ gives
$A\delta x = \alpha B^\top(\delta y - B\delta x)$. With
$r = r_0 - B\delta x$, the two combine to
$f - Ax + B^\top y = \alpha B^\top(r - r_0 + B\delta x) = 0$. $\square$

(Under regularization the right side is $\alpha\rho\,x$ instead of $0$
— an exact term joining the rounding one, and the *signal* the
refinement loop exists to remove; §15.)

The computed $s_0$ therefore consists entirely of accumulated rounding:
every exact-arithmetic contribution vanishes by construction, at
$n_0 = 0$ through the recovery $y = y_0 + \alpha r_0$ alone. This
licenses treating $s_0$ as a *sensor* of attainable accuracy, and it is
the fact the whole ceiling side stands on. In exact arithmetic (T2)
therefore passes at the first test for *any* inner truncation — the
halting object is rounding-born in every detail — while the hand-over
quality $e_y$ keeps its own exact-arithmetic story; the two are
distinct objects and the note never trades one for the other.

**The mechanism, derived.** Let $\hat x_0 = x_0 + \Delta x_0$ be the
computed backsolve, exact for a perturbed operator with
$\|\Delta F\| \le c_u u \|F\|$. Expanding the dual row and killing the
exact part by Proposition 2, the introduced $F^{-1}$ is annihilated by
the $F$ the recovery reintroduces:

$$s_0 \;\approx\; \alpha\,\|\Delta F\, x_0\|
\;\lesssim\; c_u\, u\, \alpha\, \|F\|\, \|x_0\|
\;=\; c_u\, u \max\!\big(\lambda_{\max}(A),\ \alpha\|B\|^2\big)\|x_0\| .
\qquad (14)$$

No condition number appears — the cancellation removes it; no
near-kernel scale appears — the bottom of $F$'s spectrum is irrelevant;
no condition on the right-hand side appears — the perturbation acts on
whatever $x_0$ is. Equation (14) is a max-law and its branches are
named at every use: the **flat branch** below and the **ramp branch**
above the switch

$$\alpha_\star = \frac{\lambda_{\max}(A)}{\|B\|^2} . \qquad (15)$$

Below the switch $\beta A$ dominates $F$, the achievable accuracy does
not move with the parameter, and the dual residual is a flat signal
carrying no ceiling information; above it $B^\top B$ dominates and the
residual grows in proportion. A reading taken at or below $\alpha_\star$
extrapolates from the flat branch and the ceiling law fails accordingly;
§11 draws the exclusion zone.

**The scaling, in tiers**, parallel to §4's but with the mechanism
derived rather than assumed:

> **(C1) — used.** Between the anchor $\alpha$ and the endpoint sought,
> $s_0 \propto \alpha$.

> **(C2).** The same for every $\alpha > \alpha_\star$ — checkable, the
> switch being computable from two spectral quantities.

> **(C3).** The backward-error derivation (14), which yields (C2)
> outright and (C1) on any interval above the switch.

**The constant $c_u$** is the backward-error constant of the solve —
factorization plus triangular solves — a dimension-dependent quantity
this note does not pin. It is independent of $\alpha$, $\varepsilon$,
and the data, so (14)'s branches are $\alpha$-branches and never
$c_u$-branches; and *the estimate never uses it*: the ceiling reads
$s_0$ from the record rather than predicting it, so $c_u$ cancels from
every returned quantity and appears only in statements about the
mechanism. (Its one design use — placing a regularization ladder's first
rung — travels with **TODO-R**.)

**The ramp and the endpoint.** The computed residual cannot fall below
the evaluation's own rounding, so what is observable is

$$s_0 = \max\big(\ell_d\,\alpha,\ \text{evaluation floor}\big),
\qquad (16)$$

a two-branch law (ramp / floor; the floor's refinement-side analogue
$\varphi$ is §7's). On the ramp branch, (T2)'s condition
$s_0 \le \varepsilon$ is a half-line bounded above, its structure from
(14)'s monotonicity and its location from (C1):

$$\alpha \le \alpha_c = \frac{\varepsilon}{\ell_d}
= \frac{\alpha\,\varepsilon}{s_0(\alpha)}, \qquad (17)$$

the second form reading the endpoint from the record at any anchored
$\alpha$ on the ramp. $\lambda_d = \log\ell_d$ is $\varepsilon$-free as
§0 declares; $\varepsilon$ enters (17) explicitly. The conditioning
mirror is worth one line and no more: on a near-kernel of scale $h$,
$\kappa(F) \approx \alpha\|B\|^2/h$ grows in the same proportion the
primal residual shrinks — the cost of resolution, priced where it is
paid (§5, §7) rather than through the condition number, which (14) has
just shown the dual side never sees. The two-scale hypothesis behind
that mirror is registered **(S)** in §17's mechanism layer; no
derivation consumes it.
## 7. Refinement

Section 5 relaxed the entry condition by admitting iterations. The dual
condition relaxes too, by admitting passes, and the structure is the
mirror image with two qualifications the primal side does not need, both
finite-precision in origin: a pass must be *priced* before it can be
traded against an iteration, and its descent terminates at a floor
rather than continuing to the tolerance.

**What a pass is, and what it costs.**

> **(A-pass)** Each pass starts at $y_0 = 0$ and runs the full §1
> sequence — backsolve, entry test, correction if entered, recovery —
> on the recorded residual pair $(d, p) = (s_j, r_j)$.

> **(A-restart)** *(sibling style)* A pass may instead re-solve the
> original data with the current multiplier as start. Contraction and
> wall are identical — the kernel geometry is untouched — and the
> evaluation floor is (23) unchanged; what differs is the deposit
> branch of (29): under (A-restart) the recovery's deposit is
> proportional to the *current* gap, which contracts, where under
> (A-pass) the base's persists in the pass data. **The adjudication is closed: the deposit is a flow.** Measured
> pass-entry sequences follow (A-restart) — each pass's deposit forms
> at the *current* gap's scale, because the pass map is memoryless:
> residuals are recomputed from the state, so no stored error can
> persist except as the gap itself, which contracts (a lemma with an
> elementary three-step proof, leaning on the directly-evaluated stop
> of §1). The pricing below reads under (A-restart), with toll
> magnitudes keyed to the backsolve's error level. The remainder of this
> section reads under (A-restart).

Its cost is one application of $F^{-1}$ *plus* whatever iterations run
inside it; §0's accounting prices a pass at one application only if the
interior iterations vanish. That they do, in the regime addressed, is a
property of the entry test, in three steps.

*The skip lemma.* Proposition 1 concerns the operator pair, not the
right-hand side: for any data and any start, the pass's initial
constraint residual is $-\tfrac1\alpha S(y_0 - \lambda^*)$ with
$\sigma(S) \subset (0,1)$, so (T1) passes — the iteration does no work —
whenever

$$\|y_0 - \lambda^*\| \;\le\; \alpha\,\varepsilon . \qquad (18)$$

*The multiplier of a pass.* With $U, Z$ orthonormal bases of
$\operatorname{range}B^\top$ and $\ker B$, $A_{22} = Z^\top A Z$,
$\lambda_Z = \lambda_{\min}(A_{22})$, and $\omega = \|A_{12}A_{22}^{-1}\|$,
the null-space solution of the pass system gives, exactly,

$$\|\lambda^*_{\text{pass}}\| \le C_d\|d\| + C_p\|p\|, \quad
C_d = \frac{\sqrt{1+\omega^2}}{\sigma_{\min}}, \quad
C_p = \frac{\|A_{11} - A_{12}A_{22}^{-1}A_{21}\|}{\sigma_{\min}^2},
\qquad (19)$$

and semidefiniteness caps both by recorded spectral quantities:
$\omega \le \sqrt{\lambda_{\max}(A)/\lambda_Z}$ and
$C_p \le \lambda_{\max}(A)/\sigma_{\min}^2 = \alpha_\star \kappa^2(B)$.
*(Proof: $x = B^+p + Zv$; the $Z$-projection of the first row determines
$v$ through $A_{22}$; the multiplier's lead operator has norm
$1/\sigma_{\min}$ and the $d$-bracket norm $\sqrt{1+\omega^2}$ since
$Z^\top U = 0$.)*

*The pass is one application.* At entry to a firing pass,

> **(A-entry)** $\|p\| \le \varepsilon$ and $\|d\| \le \ell_d\,\alpha$.

The first clause is *derived* under (A-term): the prior invocation met
its test, so the primal component is within tolerance when refinement
re-tests; the clause remains registered only for implementations outside
(A-term). The second holds on the ramp branch while the loop progresses,
and (T3) disposes of the regime where it does not. Then (18) with
$y_0 = 0$ and (19) give: every firing pass is entered with (T1) already
satisfied provided

$$\alpha_c \;\ge\; C_d + C_p . \qquad (20)$$

> **(A-skip)** Condition (20) holds, evaluated through (19)'s caps from
> recorded quantities.

Under (A-pass), (A-entry), (A-skip), each pass costs exactly one
application — §0's currency, discharged rather than assumed. When the
window is nonempty, (20) is implied by the $\varepsilon$-free comparison
$\sqrt{\ell_p/\ell_d} \ge C_d + C_p$, which asks spectral quantities to
clear a $\sqrt{u}$-type scale — the same character of hypothesis as full
row rank. A structural margin keeps the check generous: passes fire only
above $\alpha_c$, while the skip threshold $\alpha\varepsilon$ grows
with exactly that $\alpha$. Outside (A-skip), the true cost of a solve
lies between the additive count and the envelope $j(1 + n_0)$.

*Uses:* (A-pass), (A-entry), (A-skip), (A-term), (T1), (T2), Prop. 1,
$\lambda_Z > 0$; the window discharge additionally uses (3). Failure
mode of (A-skip): pass cost exceeds one application; ceiling budgets
understate; optimum biased high.

**A pass and its contraction.** Above $\alpha_c$ the dual component
fails and a pass fires with data $(s_j, p)$, $\|p\| \le \varepsilon \ll
s_j$; to leading order the second row forces the correction into
$\ker B$, and projecting the first row there,

$$\|x_j\| \;\le\; \frac{s_j}{\lambda_Z} . \qquad (21)$$

The recovery makes the pass's exact dual residual vanish identically, so
the computed leftover is the rounding of the pass's own chain — (14)
read at the pass's solution — and composing with (21), the pass
contracts by at most

$$\varrho(\alpha) = \frac{\alpha}{\alpha_{\mathrm{wall}}}, \qquad
\alpha_{\mathrm{wall}} = \frac{\lambda_Z\,\|x\|}{\ell_d}, \qquad (22)$$

independent of any Krylov work inside the pass, since the correction
lies in $\ker B$. But contraction is not all a pass does: updating the
accumulated pair and evaluating the next residual costs rounding at the
accumulated scale, with no $\alpha$-amplification — no solve ever acts
on the accumulated pair. Two rounding sources feed the computed
residual, and one vector bounds both: the evaluation of
$f - Ax + B^\top y$ costs, per component, a modest multiple of
$u\,(|A||x| + |B^\top||y|)_i$ — the componentwise model of a
matrix–vector product and sum [H] — and the pass chain's leftover is
the backsolve's backward perturbation acting on the pass's solution,
$|\Delta F| \le c\,u\,|L||L^\top|$ componentwise, so its components
are bounded by the row scales $(|L||L^\top||x|)_i$, the $\eta_i$ of
§12, already recorded. Taking the vector norm of the common bound,

$$\varphi = c_\varphi\, u\,\big\|\,|A|\,|x| + |B^\top|\,|y|\,\big\|
\qquad (23)$$

for this **limiting accuracy** — a vector norm of componentwise
products, not a product of norms; the distinction is idle on uniformly
scaled operators and decisive on complementarity-scaled ones, where the
product form overestimates by the scaling — the recursion is affine,
not geometric:

$$\boxed{\; s_{j+1} \;\le\; \varrho(\alpha)\, s_j + \varphi \;}
\qquad (24)$$

— the standard limiting-accuracy form of iterative refinement, reached
from the note's own quantities. Near the wall the first term dominates;
far from it the second does; both regimes live in one inequality, and
every consumer below names which term it is reading.

*Uses:* (A-entry); $\lambda_Z > 0$; **(A-const)** — a single $c_u$
bounds base and pass chains alike, with $O(1)$ scatter; **(A-norms)** —
the accumulated norms in (23) are the converged pair's; **(A-cw)** —
the componentwise model's slack over the achievable floor is modest
and uniform within a decade across the pair's rows, metered by the
screens' own separation and by every verification outcome.

**Levels, cap, and the two walls.** Iterating (24), $j$ applications
leave $s \le \ell_d\alpha\,\varrho^{j-1} + \varphi/(1-\varrho)$; the
dual condition is satisfiable only if
$\varepsilon > \varphi/(1-\varrho)$, and within the policy's reach
($\varrho \le \tau$) the conservative requirement is
$\ell_d\alpha\,\varrho^{j-1} \le \varepsilon_{\mathrm{eff}}$,

$$\varepsilon_{\mathrm{eff}} = \varepsilon - \frac{\varphi}{1-\tau},
\qquad
\boxed{\;\alpha_c^{(j)}
= \big(\alpha_c^{\mathrm{eff}}\big)^{1/j}
\big(\alpha_{\mathrm{wall}}\big)^{1-1/j}\;},
\qquad
\alpha_c^{\mathrm{eff}} = \frac{\varepsilon_{\mathrm{eff}}}{\ell_d} .
\qquad (25)$$

For $\varepsilon \gg \varphi$ the superseded family is recovered
verbatim; the correction is $\varepsilon \to \varepsilon_{\mathrm{eff}}$
and, decisively, a termination the old family lacked: the descending
term meets the floor at

$$j_{\mathrm{floor}} = 1 +
\frac{\log(\ell_d\alpha/\varphi)}{\log(1/\varrho)}, \qquad (26)$$

$\varepsilon$-free outright, beyond which passes purchase nothing.

Where the record carries the measured floor of §16 — the compensated
pair's vector difference — the measurement **replaces**
$\varphi/(1-\tau)$ wholesale: the meter reads the achieved floor
$s_\infty$ itself, and the $(1-\tau)$ inflation modeled precisely
the approach that is now measured; stacking both double-counts. The
corner deductions become $\varepsilon - 3 f_m$ (pessimistic) and
$\varepsilon - \tfrac12 f_m$ (optimistic), the constants calibrated
and carried with their provenance in the companion report.

The wall statement splits, and the split is the honest one. Equation
(22) is a **bound** — worst case over alignment and realization. The
realized per-pass contraction $\hat\varrho \le \varrho$ defines the
**realized wall**

$$\tilde\alpha_{\mathrm{wall}}
= \sup\{\alpha : \hat\varrho(\alpha) < 1\}
\;\ge\; \alpha_{\mathrm{wall}}, \qquad (27)$$

a *reading-class* quantity: every consumer — (25), the span
$M = \log(\tilde\alpha_{\mathrm{wall}}/\alpha_c^{\mathrm{eff}})$, the
treads, the caps, the halt — takes $\tilde\alpha_{\mathrm{wall}}$
alone, and the note names no slack constant. The derived bound serves
twice: substituted into (25) it yields conservative (understated)
levels, and read against a measurement it is a validator — a reading
below the bound refutes a hypothesis, not the arithmetic. Measurement
belongs to §11. **TODO-U**: no certified *upper* bound on
$\tilde\alpha_{\mathrm{wall}}$ exists; its absence is why §9's
certificate ships in robust form.

*Uses:* the above, plus **(A-flat)** — the realized contraction is
proportional to $\alpha$ between a reading and the wall, as the bound
is. Consumed only by transported readings and the halt identity;
failure mode one-signed, overstating the wall from far readings.

**The wall's product form is a bracket.** The record supports two
tiers and the output declares both: the *certified floor* is the
witnessed frontier — the budget fact, the deepest converged
refinement — and the *reading-class top* is §11's inverted first pass
ratio. Equation (24) predicts the drift that contaminates it,
$s_{j+1}/s_j \approx \varrho + \varphi/s_j$, and a regression of ratios
on $1/s_j$ would separate the contraction (intercept) from the floor
(slope, an independent $\varphi$ meter); that deconvolution is the
route's destination and not its present method, for the reasons and under
the conditions of **TODO-C**. The raw maximum ratio is retired as a wall
reading — its floor inflation is orthogonal to fit quality and understates
the wall even from settled sequences. The white pole takes the certified floor;
the black pole takes the reading; the truth lies between, indexed by
the halting policy itself — the stall rule is part of the wall's
definition, and its placement carries implementation jitter of a few
tenths of a nat. Distinct from this stall wall, and unnamed until
now, is the *factorization wall*: the $\alpha$ beyond which
$\tfrac1\alpha A$ sinks below working precision against the
rank-deficient $B^\top B$ and the factorization itself fails — a
structural death the curves inherit as an outer bound on every pole.
**TODO-F**: (A-flat)'s breadth is the note's one named gap; families
exist on which no proportional regime is observed, and the pending
shape verdict resolves whether they lie outside the hypothesis class
or inside it.

**The stall test, classified.** By (24) the observed contraction is

$$\varrho_{\mathrm{obs},j} = \hat\varrho + \frac{\varphi}{s_j},
\qquad (28)$$

one formula, two regimes, and (T3) fires for two reasons that must not
be conflated. A **wall halt** — first term exceeds $\tau$ — occurs
under (A-flat) at $\alpha_{\mathrm{halt}} =
\tau\,\tilde\alpha_{\mathrm{wall}}$ and locates the realized wall with
an extrapolation of only $1/\tau$. A **floor halt** — second term
exceeds $\tau$, the residual having descended to
$s_j \sim \varphi/\tau$ — can occur at any $\alpha$ whose tolerance
sits at or below (25)'s cap and carries no wall information whatever.
The two are separated by comparing the halt's $s$ against $\varphi$;
§11 classifies before it reads, and a misclassified floor halt
understates the wall unboundedly. The policy truncates the family at
$j_{\max} = M/\log(1/\tau)$, subordinate to
$\min(j_{\max}, j_{\mathrm{floor}})$; under (A-flat) a different $\tau$
moves the halt and the truncation and nothing else — the realized wall
is a property of the arithmetic, not the policy.

**Iterations inside a pass.** Under (A-term) every entered invocation
meets its test, so the interior count of pass $j$ is determined by its
*entry gap* — the log of its initial constraint residual over
$\varepsilon$ — and by the rate; and the entry residual of a firing
pass is **dual-born**: the pass's data is $(s_j, p)$ with the primal
component already within tolerance, so what (T1) tests is the dual
channel's image under Proposition 1, of scale
$\tfrac1\alpha\,\chi\,C_d\,\|d\|$ with alignment $\chi \le 1$. The dual
data itself is two-branch — the ramp delivers $\ell_d\alpha$, and the
recovery *deposits* its own rounding, of scale
$c\,\eta\,\|B\|\,\|\lambda^* - y_{j-1}\|$ — $\eta$ the backsolve's
relative error level ($u$ for working-precision factors, larger under
reduced precision or boundary-degraded conditioning), $y_{j-1}$ the
*current* multiplier per the flow — which does not shrink with
$\alpha$ — so the entry gap is

$$a_j(\alpha) = \log\frac{\chi\,C_d\,
\max\big(\ell_d\,\alpha,\ c\,\eta\,\|B\|\,\|\lambda^*-y_{j-1}\|\big)}
{\alpha\,\varepsilon} \qquad (29)$$

with crossover between branches at

$$\alpha_x:\qquad \ell_d\,\alpha_x = c\,u\,\|B\|\,\|\lambda^*-y_0\| .
\qquad (30)$$

Above $\alpha_x$ the **ramp branch** governs, $\alpha$ cancels in (29),
and under (A-skip) the constant gap is safely negative: interior counts
vanish and the cost is all-or-nothing — the window corollary. Below
$\alpha_x$ the **deposit branch** governs and the gap grows as
$\alpha$ falls. What a pass *pays* against its gap is capped by §1's
interior discipline: every interior pass's count is a **constant by
construction** — $\log(1/\theta)/\bar\lambda$ plus the startup
transient $k_0$ — and the gap's size prices only the single exit
invocation at full $\varepsilon$; §8 prices both. The deposit's scale
is $\eta$: at working precision it sits at or below tolerance and the
tolls vanish; where the backsolve is genuinely inexact the tolls are
real, and the constant-toll discipline is what keeps their total
linear in the pass count rather than quadratic in the gap. The surcharge is **policy-coupled**: (T3) halts many
of the solves that would otherwise pay it, so its extent is conditional
on $\tau$ as stated — tightening the policy shrinks the region by
exclusion, not by economy. The deposit's scale is the same gap the
start controls: improving $y_0$ shrinks $\alpha_x$ and the surcharge
region in lockstep with the floor arm — one number, two nuisances.

*Uses:* (A-term); (A-entry) second clause; (A-skip) for the ramp
branch's sign; Prop. 1; (19). Branch of (29) named at every use;
failure mode of reading it on the wrong branch: below $\alpha_x$ the
all-or-nothing conclusion is false and the additive count silently
grows.

**Derived and defined.** (18)–(30) are derived under the tagged
hypotheses; $c_u, c_\varphi, c$ are bounded in form, unpinned in value,
and cancel from every returned quantity except $\varphi$'s role in
screens, where order of magnitude suffices and (A-cw) calibrates. The
realized quantities $\hat\varrho, \tilde\alpha_{\mathrm{wall}}$ are
defined here and measured in §11. The regularized case is §15's.
## 8. The relaxation

Assembling the arms gives the cost of (2) as three staircases:

$$\mathrm{cost}(\alpha) \;=\; 1 \;+\; j(\alpha) \;+\; n_0(\alpha)
\;+\; \Sigma(\alpha), \qquad (31)$$

with $n_0$ descending by §5's levels, $j$ ascending by §7's, and
$\Sigma$ the interior-iteration surcharge, identically zero on the ramp
branch of (29) under (A-skip) — the classical two-term account,
recovered as a theorem where it holds — and supported, where it exists,
on the deposit-branch region below $\alpha_x$, conditional on the
halting policy as §7 states. Both families crowd as $1/k^2$, so (31) is
a broad, flat-bottomed bowl: the cheapest $\alpha$ is a range, not a
point.

**The estimable part.** The staircases are what a solve pays; they are
not what a solve can know. Passing to the level laws' interpolants —
the floor inverts to (13); the ceiling law (25) inverts the same way,
with the same wall pole reweighted, whether or not the surcharge is
present — define, on $t \in (t_\infty, T)$ and $+\infty$ outside it,

$$C(t) \;=\; 1 \;+\; \Big(\frac{M_f}{t - t_\infty} - 1\Big)^{\!+}
\;+\; (1 + c_\theta)\,\Big(\frac{d_0 + t}{T - t}\Big)^{\!+}
\;+\; E(t)^{+}, \qquad (32)$$

with the parameters of §0's declaration: $M_f = \lambda_p - e -
t_\infty$ and $t_\infty$ from (13); $d_0 = \lambda_d - e_{\mathrm{eff}}$
from (17) and (25); and $T = \log\tilde\alpha_{\mathrm{wall}}$, a
*reading-class* pole — where only the bound of (22) is in hand, it
enters as a constraint, never as $T$ (§13 makes the confusion
untypable). Each parameter is $\varepsilon$-free with $\varepsilon$
entering affinely through $e$; the note stays at fixed $\varepsilon$.

The interior discipline of §1 makes every interior pass's toll a
constant by construction, $c_\theta = \log(1/\theta)/\bar\lambda +
k_0$ with $k_0$ the recorded startup transient (§16), folded into the
ceiling coefficient of (32). What remains of (29) is the single exit
invocation at full $\varepsilon$, whose gap declines with unit slope:

$$E(t) \;=\; \frac{(G_x - t)^{+}}{\bar\lambda}
\quad\text{on the deposit branch,}\qquad 0 \ \text{on the ramp branch},
\qquad (33)$$

with $G_x$ the log of the exit deposit over $\varepsilon$ —
$\eta$-parameterized as in §7, at or below zero for working-precision
factors, so that $E$ vanishes identically on clean machines.

$\bar\lambda$ a representative log-rate over the active region. The
operand is the *firing-interior* rate, recorded as its own field:
base-harvest rates run several times faster and under-price the toll,
by mode selection — late residuals concentrate on the slow modes, the
one phenomenon behind the ratio drift, the displacement slivers, and
this rate gap alike. **The
rate model is quarantined to the pass interior**: the base arm is
priced by the measured crowding law (13), whose shape requires no rate
model, and no rate expression touches $n_0$ anywhere in this note. The
quarantine's slack, together with (33)'s use of a single $\bar\lambda$,
is a constant $\kappa_r$ — **TODO-K**: bounded in form, unsized; its
criterion is a measured spread of realized rates over the deposit
region.

**Proposition 3 (convexity).** *On its finite domain every term of
(32) is convex — the floor and ceiling terms linear-fractional with
poles at the domain ends, the interior-toll coefficient a constant by
§1's discipline, and the exit toll piecewise-linear declining — so $C$
is convex.* The discipline is what abolished the former deficit: under
full per-pass enforcement the interior tolls formed a triangle,
quadratic in the starting gap, whose truncated corner carried a
concavity deficit; the constant toll removes the mechanism rather than
bounding it. The reading of $T$ absorbs the discipline's rate
perturbation ($\varrho \to \varrho + O(\theta)$): the wall is measured,
never modeled, as always. The enumeration on the exact layer (§9) is the backstop
for exactly this corner.

**The sandwich.** $C \le \mathrm{cost}$ pointwise — each staircase sits
within one rounding of its interpolant per application, and every term
of (32) is nonnegative — while in the other direction

$$C(t) \;\le\; \mathrm{cost}(t) \;<\; C(t) + J_p(t) + 2 + \kappa_r,
\qquad (34)$$

$J_p(t)$ the smooth pass count at $t$: the gap between what is paid and
what is knowable is at most one rounding per invocation plus the
rate-model slack $\kappa_r$, and it widens toward the wall exactly as
the pass count does.
Appendix A carries the accounting; nothing in the estimate depends on
its sharpness.

The decomposition $\mathrm{cost} = C + \text{rounding}$ is the useful
one: $C$ depends on a handful of recorded quantities and varies
smoothly and monotonically in each — the fact §13's bands are built
on — while the remainder is bounded and carries no structure one solve
can resolve. The minimization of §9 is conducted on $C$, with the
staircases retained as the exact layer wherever a tread boundary is in
doubt. **TODO-P**: the adoption of $C$ as the canonical objective is
provisional, its criterion the pre-registered comparison of §9's
guarantees against exact enumeration on deployment data; resolution
edits this paragraph and nothing else.

*Uses:* (31) from §0's currency + (A-skip) branch conditions; (32)'s
floor from (13) [(A-crowd): the crowding law's continuum form];
ceiling from (25) with $T$ reading-class [(A-flat) where transported];
(33) from (29)'s deposit branch [branch named] + the single-rate
quarantine [TODO-K]; Prop. 3 as stated; (34) from term-nonnegativity
and per-application rounding.
## 9. Where it is least

**The minimizer, in closed form.** Balancing (32)'s two poles,

$$t_\ast \;=\; \frac{T + r\,t_\infty}{1 + r}, \qquad
r = \sqrt{\frac{W_{\mathrm{eff}}}{M_f}}, \qquad
W = (1+c_\theta)(d_0 + T), \qquad (35)$$

the optimal $t$ dividing $[t_\infty, T]$ at the square root of the
ratio of the arm weights — with $W_{\mathrm{eff}}$ set by the exit
toll's branch. Write

$$h'(t) \;=\; -\frac{M_f}{(t - t_\infty)^2} \;+\; \frac{W}{(T-t)^2}$$

for the stationarity slope with the toll inactive. Since
$h''(t) = 2M_f/(t-t_\infty)^3 + 2W/(T-t)^3 > 0$ on the domain, $h'$ is
strictly increasing there, so the following three cases are **exhaustive
and the minimizer unique**:

*Toll inactive at the optimum* ($h'(G_x) < 0$): the classical two-pole
balance, $W_{\mathrm{eff}} = W$, exact outright.

*Toll active at the optimum* ($h'(G_x) > 1/\bar\lambda$): the
first-order condition is a two-pole balance with the self-consistent
weight $W_{\mathrm{eff}} = W - (T - t_\ast)^2/\bar\lambda$ — the
declining toll *lightens* the wall pole, the correction signed minus.
The identity is exact; it is solved by bisection, $h'(t) - 1/\bar\lambda$
being monotone, and *not* by resubstitution, which is a fixed point that
need not contract.

*The kink* (otherwise, i.e. $0 \le h'(G_x) \le 1/\bar\lambda$):
$t_\ast = G_x$ exactly — spend until the deposit is amortized, then stop.

Two earlier forms are withdrawn. The single-form folding of the toll into
$W_{\mathrm{eff}}$: a declining line does not reweight a pole in general,
and the toll-active case is a balance point the toll *creates*. And the
boundary $(T-G_x)^2 \le W\bar\lambda$, which omits $M_f$ and so cannot
sign a three-term competition; the condition above is the exact split.
The case boundaries are computable from the recorded quantities at zero
cost. Proposition 3 holds every term of (32) convex, so the sum is a
single bowl. Displacement under parameter error is second-order:
the bowl's curvature at $t_\ast$ is explicit from (32), and a
perturbation $\delta$ of a parameter first induces a cost gap
$\Delta \le |\partial C/\partial(\mathrm{param})|_{t_\ast}\,
\delta$ — the partial explicit from (32); for the wall,
$|\partial C/\partial T| = (1+c_\theta)(d_0+t_\ast)/(T-t_\ast)^2$
— and the minimizer then moves by at most $\sqrt{2\Delta/m}$ with
$m$ the curvature — pole errors additionally
enter log-compressed and damped by $1/(1+r)$, which is why wall
misestimates of even large factors displace the recommendation by
treads, not decades. Under warm starting $M_f$ falls and $r$ grows:
better starts move the recommendation down-parameter, buying distance
from the wall with the floor's savings.

**The answer is a plateau.** On the exact layer the minimizers of (31)
form flat treads; the recommendation is the tread containing $t_\ast$,
reported as the interval between the adjacent floor level below and
ceiling level above, together with its integer budget. Where two treads
tie — generic in the empty-window regime — **the tie breaks downward,
by arithmetic**: the surcharge and the rounding it rides on load the
high-$\alpha$ tread only, and the stall asymmetry of §12 points the
same way. The old preference for low $\alpha$ is here a consequence,
not a policy. Hump anatomy on the exact layer is rounding jitter of
depth at most one on the convex envelope; no structure deeper than one
application survives Proposition 3.

**Two guarantees, from two lines of counting.** Every term of (32) is
nonnegative, so the smooth pass count at the recommendation satisfies
$J_p(t_\ast) \le C(t_\ast) - 1 \le \min\mathrm{cost} - 1$; feeding this
into (34),

$$\mathrm{cost}(t_\ast) \;\le\; 2\,\min_\alpha \mathrm{cost}
\;+\; 1 + \kappa_r : \qquad (36)$$

**self-scaling** — the recommendation cannot lose more than the problem
was going to cost anyway, and the factor is tight on instances that pin
both arms. A posteriori, one verification solve at $t_\ast$ certifies
the realized loss:

$$\mathrm{cost}(t_\ast) - \min_\alpha \mathrm{cost}
\;\le\; \mathrm{cost}(t_\ast) - C(t_\ast) \;+\; e_C, \qquad (37)$$

where $\mathrm{cost}(t_\ast)$ is *recorded*, $C(t_\ast)$ is evaluated,
and $e_C$ is the estimation error of $C$'s parameters — **TODO-E**:
bounded in form through §13's band widths, unsized in value; its
criterion is the deployment distribution of band widths at the
minimizer. The certificate ships in this robust form and not the pure
one because no certified upper wall bound exists (**TODO-U**, §7): a
pure $C$-gap certificate is valid only when $T$ is a two-sided
reading, and the note does not assume one. The guarantee ladder is
therefore: (36) always; (37) after one solve; exact enumeration on the
recorded staircases when a tread boundary is genuinely in doubt.
Operationally these are promises to the next chooser: what following
the advice costs, if the operator it meets is the operator that
priced it — the conditional named in §14.

One warning the smooth layer owes the exact one, restated from the
original analysis: the sum of two ceilings of convex functions is
not quasi-convex — each staircase is unimodal, their sum need not
be — so the set of cost-minimizing parameters is in general a
**union of separated intervals**, and the smooth minimizer can land
in a gap between optimal blocks. The ceiling taken last is
quasi-convex; it is the *separate* rounding that costs unimodality,
unavoidably — a fractional iteration cannot be traded for a
fractional pass. Step 5's enumeration is the answer: the blocks
differ in height by at most the corridor, and the exact layer
visits them.

**Three demotions.** The window is a *sufficient* condition, not a
necessary one: $W_0$ is where both staircases stand at their minimum
simultaneously, and when (3) fails the cost is still finite and usually
small — one arm pays a step so the other need not. The endpoints
$\alpha_0, \alpha_c$ are demoted to poles: natural targets exactly when
simultaneously reachable, with (35) saying how far between them to
stand otherwise. And the conditioning of $F$, which an older treatment
took as the object of optimization, appears nowhere in (35): what is
balanced is work against work, and the condition number was only ever a
proxy for the floor arm's half of it.

*Uses:* (32) and Prop. 3 for (35) and the plateau statement; (29)'s
deposit branch + §12's stall asymmetry for the tie-break;
term-nonnegativity + (34) for (36); recorded-cost-vs-evaluated-$C$ +
band widths [TODO-E] for (37); TODO-U for the certificate's form.
Failure modes: the smooth minimizer is unique by Prop. 3's convexity,
but the exact staircase's plateau may still be non-unique beyond the
tie rule — the enumeration backstop governs; a
bound consumed as $T$ invalidates (35) and is untypable by
construction (§13).
## 10. What a solve determines

The estimate of §8 is generated by a short parameter vector: the two
scales $\ell_p, \ell_d$; the floor's pair $(t_\infty, M_f)$ through the
constants $c_k$; the wall $T$; the gap quantities of (29) together with
the toll constants $(\bar\lambda, k_0)$ and the exit gap $G_x$ that
(32) now consumes; the limiting accuracy $\varphi$; and $\lambda_Z$. They are not
equally available, and the availability structure is the section's
first half; the second half is the *type* each measurement carries.

**Determined at no cost.** $\ell_p = \alpha\|r_0\|$ and
$\ell_d = s_0/\alpha$ are read from residuals the algorithm forms
regardless — one as the correction's right-hand side before any work,
one as the refinement loop's own test after it — fixing the lowest
tread of each staircase exactly and without calibration. They differ in
kind: by Proposition 1, $\ell_p$ is an exact quantity; by Proposition
2, $\ell_d$ is a single realization of a rounding process, and
successive passes furnish further readings of *related* quantities —
each pass solving a smaller system — not repeated draws of one. §12
records what that costs. One caution travels with warm starts: a start
can shrink the gap while worsening its alignment with the resolved
directions; $\vartheta = 1/a_1$ is the meter, and size and alignment move
independently.

**Determined by the iteration, to its own depth.** The recursion's
coefficients yield $c_1, \dots, c_{n_0}$, one per iteration; the
informative ones are those *above* where the solve ran — the levels its
Krylov space explored. A single iteration fixes the first level; only
$n_0 = 0$ leaves the floor family unknown, and that is the spectrally
blind case of §5, in which the window is nonempty and no budget
question arises. Levels *below* the run are reachable only by
transport, which is one-signed (§5's census).

**Not determined by the iteration at all.** $\lambda_Z$ lives on
$\ker B$, exactly the subspace about which $B\,\delta x = r_0$ says
nothing; no Krylov work reveals it. But the wall need not be assembled
from it: one *observed pass* gives the ratio $s_{j+1}/s_j$, hence a
reading $T = \log(\alpha/\varrho_{\mathrm{obs}})$ — subject to the
screens of §11, since (28) shows the ratio reads the floor near it and
the augmentation's own contraction below the branch crossing, and
subject to (A-flat)'s one-signed high transport. The spectral formula
$\lambda_Z\|x\|/\ell_d$ is a *bound* on the **contraction** wall: it errs
conservative by worst-case alignment, by an uncontrolled factor, and it
enters as **validator** — never as $T$, and never as the clip, for which
§12 gives the reason. Under $\rho > 0$ the same quantity becomes readable
from the record: inverting (44),
$\lambda_Z \approx \alpha\rho\,(1-\varrho_d)/\varrho_d$ from one
recorded pass ratio, at no application. The reading is biased **high**,
so it too is a validator and not a clip source, and its usable range is
the neighbourhood of the knee — below it the ratio is floor-dominated,
above it $(1-\varrho_d)/\varrho_d$ is a difference of near-equal
numbers. So $\lambda_Z$ is not determined by the iteration at
$\rho = 0$ and is, one-sidedly, at $\rho > 0$.

**The duality.** Each family is observable exactly where it binds: the
floor is revealed by a solve that grinds, the ceiling by a solve that
passes, and each is silent about the other. At the cost minimum a solve
does neither — **the cheapest solve is the one that says least about
why it is cheapest** — and away from the minimum, in either direction,
the solve pays for information about the side it has moved toward.
This is why the estimate is banded rather than pointed, and why the
bands narrow with distance from the optimum in exactly the direction
the solve happened to sit.

**Every parameter is an interval with typed endpoints.** A measurement
is one of three kinds, and the kind travels with the value: a
**reading** — a measurement of the quantity itself, admissible as a
value, carrying its route and classification (a wall-type halt; a
screened ratio); a **fit** — carrying its residual and the width of its
support, disqualified when the support does not span the decay; a
**bound** — one-sided by construction, carrying its certified
direction, admissible only as a constraint. Multiple sources for one
parameter **compose by intersection**: each contributes its certified
side, nothing is averaged, and an empty intersection is not an error to
suppress but a fired consistency check — some tagged hypothesis is
false, and the estimate says which sources disagree. A one-signed
reading is a half-line anchored at its value; a two-sided reading is a
short interval; a bound is a half-line from its certificate; the
composed parameter is as narrow as its best sources and no narrower.

Since $C$ is monotone in every parameter — the ceiling term falls as
$T$ grows, the floor term rises with $M_f$, each entry one-signed — an
interval parameter vector pushes through (32) exactly: the two corner
evaluations bound $C$ from below and above with no slack introduced by
the composition. That fact is the whole of §13's construction, stated
here because it is a property of the *parameters*, not of the bands.

*Uses:* Props. 1–2 for the two scales' kinds; §5's harvest and census;
(28) and (A-flat) for the ratio route's tags; the worst-case alignment
of (21) for the spectral bound's direction; monotonicity of (32) in
each argument, by inspection. Failure mode of composition: an empty
intersection, which is reported, never resolved by averaging.
## 11. Measurement

The record of a solve (§16) contains the two residuals, the
per-invocation counts and met-flags, the harvest coefficients, the
per-pass residual pairs, and any halt with its final ratio. This
section turns the record into the typed intervals of §10, one parameter
at a time, with each route's screens and each screen's reason.

**The floor.** The harvest supplies $c_1, \dots, c_{n_0}$ exactly at
the anchor, hence the levels above the run; the direction
guard is absolute — **a level below the anchor is not reportable as a
value**: searches clamp at $\alpha' \ge \alpha$, and a closed-form
level below the anchor is reported as a bound, its optimistic direction
tagged (§5's census); depth attaches its caution at the node level
only, as §5 states. Where the
interpolant (13) is *fitted* to observed counts rather than read from
the harvest, the fit is meaningful only if its support spans the
decay — budgets differing substantially over decades — and the fit
residual is the meter: a residual that is an appreciable fraction of
the data's range **disqualifies** the interpolant in favour of the
enumerated levels — the fit consuming records the caller already
possesses, solves that were happening anyway, an interior-point
method's own iterations being the canonical case; the estimate
purchases none — which is a decline of §13's upper band on the floor
side, not a degraded fit to ship anyway.

**The exclusion zone.** On and below the switch (15) the flat branch
governs (14): the dual signal does not move with the parameter, the
ramp's constants are not asserted there, and ceiling readings from that
zone are excluded as a class — extrapolation from a flat signal, §6.
The zone's upper reach is problem-dependent and generous; a route whose
anchor sits at or below the switch reports itself closed.

**The wall, by three routes in order of extrapolation distance.**

*First, a classified halt.* A halt at $\alpha_{\mathrm{halt}}$ yields
$\tilde\alpha_{\mathrm{wall}} = \alpha_{\mathrm{halt}}/\tau$ with an
extrapolation of only $1/\tau$ — but only after classification, since
(28) gives halts two causes and only one is the wall. The comparison is
$s_{\mathrm{halt}}$ against the floor — the *measured* floor where the
compensated re-evaluation of §16 is recorded, $\varphi$ as fallback:
wall-type ($\gg$) is a reading; floor-type ($\sim$) is discarded,
carrying no wall information. The two
populations sit far apart, so the comparison is not delicate; using
floor-type halts as wall readings is the single largest corruption the
unclassified route admits. Under (A-term) a halt is a contract event,
which is what entitles its $\alpha$ to be read at all.

*Second, the first pass ratio.* A ratio $s_{k+1}/s_k$ inverts to the
wall. The route takes the **first** ratio and inverts it: by (28) the
observed ratio is the contraction plus the floor share $\varphi/s_k$, and
that share grows as $s_k$ shrinks, so the earliest reading is the least
contaminated. The *floor-share screen* remains — a reading whose share is
an appreciable fraction of the ratio measures the floor, not the
contraction — but no minimum sample and no composition rule is imposed;
both are **TODO-C**.

Every reading carries two contaminations of *opposite* sign, and which
dominates is set by the floor share. (A-flat)'s transport pushes the
inferred wall up — the unsafe side; the residual floor share pushes the
observed ratio up and hence the inferred wall *down*, the conservative
side. The estimate does not attempt to separate them; the clip exists for
the first and **TODO-C** for the second. Transport of a ratio across
$\alpha$ consumes (A-flat), with its one-signed high tendency; high is
the unsafe side, which the clip contains. **TODO-F** stands at this join.

*Third, an observed budget.* A **converged** solve with $j$ passes at
$\alpha$ implies $\alpha \le \alpha_c^{(j)}$, whence
$\tilde\alpha_{\mathrm{wall}} \ge
(\alpha^{j+1}/\alpha_c^{\mathrm{eff}})^{1/j}$ — a bound, not vacuous
even at one pass. A stalled solve's pass count implies nothing: its
passes were attempted, not sufficient, and the route consumes converged
records only.

The wall interval's **certified floor is the budget fact**: a solve
whose refinement contracted at $\alpha$ certifies the wall above
$\alpha$ — the one certificate among the routes, held by the thing
that happened. A classified halt, when present, overrides with its
two-sided reading.

**Feasibility is witnessed, not screened.** The input is a *completed*
solve at the caller's $\varepsilon$. A converged record is proof that the
tolerance is attainable on this instance — not estimated, witnessed — and
a halted record is classified already: a floor-type halt *is* the finding
that the arithmetic ran out. No separate screen is carried, and the
calibrated zone boundaries an earlier version used are retired; they
answered at the record's own $\varepsilon$ a question the record settles.
For a tolerance *tighter* than the record's, feasibility is genuinely
open, and the curves answer it directly: $\varepsilon$ enters explicitly
and infeasibility appears as the effective tolerance going non-positive
against $\varphi$.

**Composition is §10's intersection.** The routes are derivations of
one quantity sharing $\ell_d$ and $\|x\|$, not samples of it: each
contributes its certified side, a large disagreement is evidence
against the longer extrapolation, and an empty intersection is a fired
check naming its sources. Where no reading exists at all, the wall
interval remains a half-line and nothing further is done *here*: §13's
minimax on the wide box returns the conservative answer — in substance
the one-pass interval, demanding more work over a narrower range,
never promising a range that does not deliver — with no special case
in the pipeline.

**The cost below the crossover is policy-conditioned**, and so
therefore is its measurement: a pass that (T3) halts before its
interior work makes the solve fail there, while the same solve under a
laxer policy converges and pays the surcharge. Every level, budget,
and interval this section supplies is conditioned on the policy §1
states; a caller comparing across policies is comparing different
objective functions and should say so.

**Duality, operationally.** A solve that iterated sat below the
ceiling and can price only the floor; a solve that passed sat above
the floor and can price only the ceiling. **Each prices movement in
the direction it already went**, and neither can recommend reversing;
the zero-work endpoints bound the region regardless, and only the
graded levels inherit the asymmetry.

**Exhaustion has two causes**, separated. The levels may fail to
interleave at any affordable budget — the estimate cannot locate an
interval. Or the arithmetic forecloses the ceiling —
$\varepsilon_{\mathrm{eff}} \le 0$, or the cap reached at once — in
which case no policy and no budget recovers a window at this tolerance
and the difficulty is the tolerance against the limiting accuracy,
reportable as such. In either case what is established is that the
estimate cannot resolve, not that nothing exists: the input solve
converged, and a configuration a solve has solved cannot be
infeasible. The witnessed configuration of §13 is returned in place of
a verdict.

*Uses:* §5's harvest, guard, and census; (13)–(15), (25), (26), (28);
(A-term) for halt standing; (A-flat) with TODO-F for transport; the
worst-case-alignment direction of the spectral bound; §10's
intersection; the deposit branch of (29) with policy coupling; (A-cw)
for every $\varphi$-based screen's calibration.
Failure modes named per route above; the composed failure mode is the
fired intersection, which is output, not error.
## 12. The two ends differ in kind

The floor rests on an exact identity and an assumption with a
computable meter — a meter that also adjudicates warm starts, whose
gains in size can trade against alignment. The ceiling rests on residuals that are *entirely
rounding* — by the exactness of the dual row they have no
exact-arithmetic part at all — so its quantities are draws from a
distribution rather than readings of a constant. That difference in
kind persists; what this treatment adds is a measurement protocol on
the ceiling side in place of a single helpless reading, so the
dispersion that could not previously be calibrated becomes, in part, a
recorded spread.

$s_0$ is formed once, before any pass, and its scatter cannot be
calibrated from one draw; $\varphi$ sorts readings on the ramp from
readings off it without bounding the error of either, with the
componentwise scale (23) whose slack is modest and within-a-decade
uniform ((A-cw)). But at $j \ge 2$ — precisely where the ceiling
binds — the passes supply consecutive dual residuals, and under
(A-term) their ratios estimate the contraction directly. Each ratio is
screened by its floor share $\varphi/s_k$: readings with the floor a
small fraction of the residual lie on the geometric branch of (24) and
carry wall information; readings near the floor carry none and are
discarded rather than averaged in. A halt is classified before it is
used: $s_{\mathrm{halt}} \gg \varphi$ marks a wall halt whose
$\alpha/\tau$ *is* a wall reading; $s_{\mathrm{halt}} \sim \varphi$
marks exhaustion of the arithmetic and says nothing about the wall.
The realized wall §11 returns is taken from these routes in order of
extrapolation distance — a classified halt first, the inverted **first**
pass ratio second, the least floor-contaminated reading the record holds
— with the spectral formula typed as §10 types it, a bound with a
certified direction, admissible as a constraint and never as a value, and
the certified floor being §11's budget fact. The
screened readings' spread is then an honest per-solve proxy for the
ceiling-side dispersion: correlated draws still, normalized per pass,
with the residual one-signed tendency to read high charged to (A-flat)
rather than hidden in a constant. Less than a calibrated measurement;
considerably more than one draw.

What removes the remaining difficulty is not precision but
**dominance**. Cost is flat within a tread, so an error in a level is
free unless it crosses a tread boundary. Consecutive ceiling levels sit
$(1/j - 1/(j+1))M$ apart and consecutive floor levels
$(1/(k+1) - 1/(k+2))\Lambda$ apart, both spans now measured objects —
$\Lambda = M_f$ from the floor's poles, $M$ from the realized wall, a
wider and honester span than the spectral bound gives. The scatter is
irrelevant while $j \lesssim \sqrt{M/\mathrm{sd}}$ and
$k \lesssim \sqrt{\Lambda/\mathrm{sd}}$: quadratic in the budget, so
slack at small budget is large; both spans widen as the tolerance
tightens, in step; and $j$ stays small in the regime the ceiling
serves. The floor side runs the higher budgets and is the side that
can bind — deep in the empty-window branch, where adjacent levels are
in any case indistinguishable. The two known imprecisions coincide
rather than compound.

**Two inputs, not one.** The closure needs the dispersion small *and*
the spans non-degenerate, and these are different quantities: a
perfectly spread perturbation still fails the ceiling condition at
$j = 2$ if $M$ is a fraction of a decade. $M$ collapses exactly when
$\lambda_Z$ does — the near-singular regime — and the collapse is
*observed* rather than inferred: the realized wall is measured, so a
short span reports itself in what the estimate returns, the spectral
bound serving as explanation rather than messenger. The closure's
failure mode is a nearly singular reduced Hessian, not a misbehaviour
of the arithmetic.

*The dispersion, as commentary.* By Proposition 2 the dual residual is
the backsolve's backward perturbation acting on the solution, so its
squared norm is a sum over the perturbation's rows and its logarithm
concentrates at a rate governed by the participation ratio of the row
scales, $n_{\mathrm{eff}} = (\sum \eta_i^2)^2/\sum \eta_i^4 \in [1,n]$
with $\eta_i \propto (|L||L^\top||x|)_i$. Nothing returned depends on
this; it is recorded as an account of *why* the ceiling side scatters,
not as a computation, and it carries no registry row.

**Asymmetry of consequence.** Which direction an error takes is not
symmetric. Leaving the habitable range downward costs iterations,
which the next solve reports immediately. Leaving it upward fires
refinement against a level the arithmetic cannot reach, and under
(A-term) that is not a degraded return but a failed step by contract.
The asymmetry is why the recommendation is clipped. **The clip is not a
separate mechanism.** It is the wall interval's certified lower endpoint,
and it acts through the shape of the curve being minimized: $C''$ blows up
there, so (40) cannot recommend past it (§13). Only a source that *cannot
read high* may set that endpoint — the budget fact, and a classified halt.
The derived spectral bound may not: it cannot read above the *contraction*
wall, but the machine dies at the stall-policy wall, measured 0.8–2.4 nats
below the derived level on the test instances, so a clip there can
recommend into the dead band. The consequence is accepted rather than
patched: **a record with no halt and no pass has no certified endpoint
above its own $\alpha$**, and its clip is the anchor. A well-placed solve
buys the least room to move — the duality of §10, at its sharpest —
and why §9's ties break downward. The margin's direction is settled by
this asymmetry, not by the estimator's scatter.

*Uses:* (A-term); Prop. 2 and the dual row's exactness; (24) and (28)
for the ratio route and screens; (A-flat) for transported readings,
carrying the one-signed high tendency; the level laws (13), (25) for the treads; (A-cw) for $\varphi$'s
calibration.
## 13. The bands

Section 10 delivered the parameters as intervals with typed endpoints
and the monotonicity of (32) in each. This section assembles the
estimate those two facts determine. Nothing new is measured here; the
construction is bookkeeping made load-bearing.

**Two corner evaluations, exact.** Order the parameter intervals by the
direction each argument enters (32): for the *optimistic corner* take
$T$ at its upper endpoint and each cost-increasing parameter at its
lower; for the *pessimistic corner*, the reverse. Define

$$C'(t) = C\big(t;\ \theta_{\mathrm{opt}}\big), \qquad
C''(t) = C\big(t;\ \theta_{\mathrm{pess}}\big), \qquad (38)$$

each extended-real: $+\infty$ outside its own finite domain. By
monotonicity these are the exact pointwise infimum and supremum of $C$
over the parameter box — no relaxation, no sampling, two evaluations —
so wherever the tagged hypotheses hold,

$$C'(t) \;\le\; C(t) \;\le\; C''(t) \;\le\;
C'''(t) := C''(t) + J''_p(t) + 2 + \kappa_r, \qquad (39)$$

the last step stacking the sandwich (34) at the pessimistic corner: the
inner band $[C', C'']$ is *measurement ignorance* and shrinks with
information; the outer step to $C'''$ is *arithmetic* and shrinks with
nothing. The true cost lies in $[C', C''']$; the true $C$ in
$[C', C'']$; both statements are per-$t$ and testable.

**The extended-real convention does the domain's work.** Each curve's
finite domain is the open interval between its own poles, and the three
domains are nested: $C''$'s is the narrowest, since the pessimistic
wall is the nearest. The strip where $C'$ is finite under an infinite
$C''$ is precisely *the wall is somewhere in here* — the estimate can
neither promise nor exclude — and needs no annotation beyond the values
themselves. With a wall interval unbounded above ($T$ known only as a
constraint), $C'$'s ceiling term vanishes and $C''$'s domain ends at
the wall's certified lower endpoint: the bound-only regime is not a
second mode of the estimate but the same object with a wide box. **The
consequence deserves stating plainly.** A converged solve's certified
lower wall endpoint is the budget fact — the wall lies above the
anchor — so in the bound-only regime $C''$ has its pole *at the anchor*
and the promise does not extend above the parameter the caller already
ran. The upper curve prices the direction one came from. Only a wall
reading moves that pole, and §11's routes supply one exactly when the
solve was placed badly enough to make several passes; the duality of §10
is here at its sharpest.
Assertion trims further from the left — below the feasibility line
$\varepsilon \le \varphi/(1-\tau)$ nothing is finite, and on the flat
branch below the switch the ramp's parameters are not asserted (§6,
§11) — and the *anchor* $t_{\mathrm{a}}$, the parameter of the
designated solve, is part of the returned object: (A-flat) and its
kin are transport hypotheses, so the box's honesty degrades with
distance from $t_{\mathrm{a}}$, one-signedly where the tags say so.

No $-\infty$ arises, and the asymmetry is structural: cost is a count,
so nature supplies the lower edge $1$ before any measurement; total
ignorance is the proper pair $(1, +\infty)$, not an improper function.
Every scrap of floor information raises $C'$ from one; every scrap of
wall information lowers $C''$ from infinity; the band is the exact
account of what has been learned.

**The default query is minimax, and the width is part of the answer.**
The recommendation is

$$t_{\mathrm{rec}} = \arg\min_t C''(t), \qquad (40)$$

— which is also the minimax-regret query outright: the shipped
regret compares $C''(t)$ against $\min_s C'(s)$, a constant in $t$,
so $\arg\min_t [C''(t) - \min_s C'(s)] = \arg\min_t C''(t)$ and
(40) needs no separate justification — the safe query and the
regret-optimal query are the same choice.

the parameter with the best *guaranteed* cost. Three numbers accompany
it, none optional: the promise $C'''(t_{\mathrm{rec}})$, an upper bound
on what will be paid; the possibility $\min_t C'(t)$, a lower bound on
what the optimum could be; and their difference, which bounds the
regret —

$$\mathrm{cost}(t_{\mathrm{rec}}) - \min_\alpha \mathrm{cost}
\;\le\; \big[C''(t_{\mathrm{rec}}) - C'(t_{\mathrm{rec}})\big]
\;+\; J''_p(t_{\mathrm{rec}}) + 2 + \kappa_r, \qquad (41)$$

band width at the recommendation plus corridor. A wide band is not a
failure state; it is the estimate reporting low visibility, and the
width profile names its own remedy: the interval endpoint dominating
the width at $t_{\mathrm{rec}}$ identifies the *cheapest extension* —
a wall reading, a wider floor fit — and the recommendation carries that
identification. The pessimistic corner is also the safety: since $C''$
blows up at the wall interval's *lower* endpoint, (40) never
recommends past anything certified, and the clip of §12 is not a side
constraint but the shape of the curve being minimized. A bound
consumed as a pole is unwritable: bounds tighten interval endpoints,
and (38) has no other slot for them.

**The next solve closes the loop.** The next mandatory solve, read
against the transferred bands, records its cost and upgrades
(41) to the certificate (37) with the realized left side. Three
outcomes, all terminal: the solve converges within the promise — the
recommendation stands certified; it converges above the promise — some
tagged hypothesis failed, and the composed intervals say which sources
to distrust; it fails to converge — the promise is *withdrawn* and the
instance flagged, which is feasibility settled at a tolerance no record
had witnessed, by the one procedure that cannot be wrong about it. Where no finite $C''$ exists anywhere — a
fired fit check, an empty wall intersection — the estimate declines to
recommend, returns the finite $C'$ with the decline's named reason, and
offers the *witnessed configuration*: the designated solve that ran is
proof of its own feasibility and cost, and reporting it is the floor
below every failure mode.

**Black and white.** The exposition names the two outputs by shade:
**black** the certified-optimistic curve (nothing achievable lies below
it, *within the box's assertion range and to a tolerance owed by*
**TODO-G**), **white** the corridor-included promise (realized cost
never exceeds it, hypotheses named). Both are convex extended-real
functions — finite on an interval between their poles, $+\infty$
outside, convex there by Proposition 3, the lower curve's returned
floor additionally clamped at §5's termination level — ordered black
$\le$ white pointwise, with **the paid cost and the law's corner
set** inside the black–white bracket: the cost via the corridor, the
corners under their tagged hypotheses, the smooth $C$ the *instrument*
of these claims and not their referent. Black's infinity is the only verdict on
the problem: both bounds *identically* infinite — no finite window
anywhere on the axis — is a certificate of impossibility at the stated
tolerance, while pointwise double-infinity is mere abstention (a flank
the estimate declines to price, the entry clamp being the standing
example) and white's infinity alone is a shrug. On
every axis of decline the surrender order is white, then black.

**Derived guesses.** Between the certified curves a caller may want a
working number; the bracket supplies it *derivatively* — the midpoint,
or the companion report's calibrated quantile — with no independent
constitution, no claims, and one standing rule inherited from the
estimate: **no sentence containing a promise cites a guess**. The
hand-over carries the two curves; the certified recommendation carries
the promise; any interior guess travels tagged and uninsured — the
chooser's risk posture is priced, not chosen for them. On a zero-pass
record the certified window contains the global minimum outright and
the parameter is the chooser's free move, so *both flanks* are
decision-irrelevant there — the window is the report, complete in
itself; and a single purchased refinement pass would
mint the missing contraction reading — a caller's option, named as
the width's remedy and declined by §0's zero-cost discipline: the
note prices the silence rather than paying to end it.

*Uses:* §10's typed intervals and monotonicity; (32), (34);
Prop. 3 for both corners' convexity (each corner is a fixed parameter
vector, so $C', C''$, and $C'''$ are convex where finite); TODO-K in
(39)/(41) through $\kappa_r$; TODO-E as the deployment size of the
inner width; (A-cw) for $\varphi$'s calibration — verification retains the
feasibility duty at tolerances no record witnesses. Failure mode: hypotheses false *inside* their tags — transport
beyond (A-flat)'s class — moves $C$ outside the inner band; (39) is
conditional exactly as its constituents are, and the verification
solve is the detector of last resort.
## 14. The estimate

The algorithm, end to end. Its input is the record (§16) of one
*designated* solve and the policy in force; its output is the estimate — the typed parameter
box, the three curves it generates, and the answers to the queries
asked of them. Every step below is a citation; nothing is introduced
here.

**Step 1 — qualify.** Check the met-flags: a record violating (A-term)
is an implementation event, not data, and the algorithm declines it by
name. Classify any halt (§11). Feasibility at the record's own
$\varepsilon$ needs no screen: a converged record witnesses it, and a
halt is classified already (§11).

**Step 2 — measure.** Assemble the parameter box per §10–§11: the two
scales with their kinds; the floor pair read at the anchored corner — $\hat t_\infty$ from
the settled tail, $M_f = (K+1)(t_K - \hat t_\infty)$, the hyperbola
through the deepest witnessed corner, one free interval with derived
corners (sign flip at $t_K$) — with the guard and the residual meter; under $\rho > 0$ the entry reading decontaminated first — §15's one
purchase, its trigger the a-priori screen — before any floor claim;
the transported floor
admitted only within the box's assertion range (**TODO-G**); the wall
interval by
intersection of the routes' certified sides; the gap quantities with
their branch; $\varphi$ with its conservative-only tag; the anchor.

**Step 3 — assemble.** The two corner evaluations (38) and the corridor
give $C' \le C'' \le C'''$ per (39), extended-real, domains implicit in
the values. The toll quantities enter here with modes like every
parameter: $k_0$ read from the first invocation's iteration history at
no cost, $\bar\lambda$ a realized-rate reading carrying TODO-K's tag
(its corner spread is the measured rate spread), and the exit gap
$G_x$ read from the recorded entry gap — pessimistic corners take the
slow rate and the large gap, optimistic the reverse. If $C''$ is finite nowhere — a disqualified floor fit, an
empty wall intersection — skip to step 6's decline.

**Step 4 — query.** The query is tread-primary: the levels are
enumerable from the record, so the integer budget is evaluated on
the treads and minimized over them directly — the enumeration step
5 already performs, promoted from refinement to the query itself,
which removes the quasi-convexity gap mode §9 names (the smooth
minimizer landing between separated optimal intervals). The smooth
(40) supplies the seed and the analytics: $t_{\mathrm{rec}} =
\arg\min C''$, with the promise, the possibility, and the width, the
width's dominant contributor naming the cheapest extension. Other
queries evaluate the same object — the cost of holding a caller's
current $\alpha$ is $C''$ there against the possibility, and needs no
new machinery — but the default is what ships.

**Step 5 — refine on the exact layer.** The recommendation is the
*tread* containing $t_{\mathrm{rec}}$: the interval between the
adjacent enumerated levels, its integer budget attached, ties broken
downward by §9's arithmetic. Where the exit toll is
active — a positive recorded entry gap at the anchor — enumerate the
neighbouring treads outright, cheap insurance at the staircase's least
regular region; the smooth layer proposes, the exact layer disposes.

**Step 6 — report, in one of three shapes.** *A recommendation*: the
certified pair of §13 with the tread, the
budget, the promise, the possibility, the width with its extension,
the anchor, and the mode of every parameter — under the
policy, as stated. *A decline*: the named check that fired, the finite
$C'$ (never less informative than "at least one application"), and the
witnessed configuration — the designated solve is proof of its own
feasibility and cost, and reporting it is the floor below every
failure. *A withdrawal*: reached only through step 7.

**Step 7 — hand over, then read the next obligation.** The estimate —
bands, recommendation, widths, tags — is handed to the chooser of the
next parameter; no solve is performed here. The next *mandatory*
solve is then read against it, checkpointed in flight: after its
first application, $s_0$ against the transferred ramp; after its
first pass, the ratio against the transferred contraction — early
warning, delivered while an iteration remains to react in. Then
containment, in three outcomes: the solve lands within the
transferred bands — the hand-over is certified, by (37),
*conditionally on the smallness of the drift between iterates, which
is this step's named scope and not its theorem*; it lands above — a
tagged hypothesis failed or the drift did, the intersection and the
checkpoints apportioning suspicion, and the miss itself a reading; it
fails — the promise is *withdrawn*, the instance flagged
infeasible-at-tolerance, settled by the procedure that cannot be wrong
about it. In every outcome the new
solve's record refreshes the estimate, and the loop is the process
itself: obligation, record, estimate, hand-over, obligation.

**What is never returned.** No averaged routes; no value for any
unpinned constant; no level below the anchor as a value; no wall from an
unclassified halt or an unscreened ratio; no feasibility verdict
beyond the witnessed and withdrawn semantics; no recommendation past
the clip, which (40) enforces by shape rather than by exception; and
no point recommendation — the tread is the answer, because the cost is
flat across it and precision inside a tread is fiction.

Appendix B renders the solver itself as a listing, with the payment
sites and the recording arrays identified line by line.

**Cost of the estimate itself.** The estimate performs no solve.
Steps 1–6 are arithmetic on recorded quantities — zero applications
beyond the designated solve's own, plus any extensions purchased in
the currency; step 7 reads a solve that was mandatory anyway. The
estimate's marginal cost is therefore zero.

**Iteration zero, and the policy of silence.** The hand-over chain
needs a first link: with no record, the recordless default is generous
headroom above the switch (15) — the classical balance point, climbing
with the operator. While the default sits inside the window its solves
achieve cost $1$, the global minimum of $C$, so no estimate can improve
on the choice and the estimate's correct behavior is silence; the free
monitor is the pass count, zero passes being a witnessed in-window
certificate. The first nonzero passes are the corridor's arrival — the
descending edge has caught the climbing default — and the estimate
takes the steering: re-center through the transported window, price
the rest by the bracket's midpoint, derived and uninsured. The endgame closures — the tolerance
ladder against the measured floor, the switch crossing the corridor —
are re-scaling's territory, and the decline shapes say so honestly.

*Uses:* §10 (box, intersection, monotonicity), §11 (routes, screens,
guards), §13 (bands, query, decline, verification), §9 (tie-break,
enumeration backstop, certificate), (A-term) at step 1, (A-cw) at
steps 1 and 7, TODO-E through (37) at step 7. Failure mode of the
whole: hypotheses false inside their tags, detected at step 7 or not
at all — which is the honest boundary of a one-solve instrument.
## 15. Regularization

Production runs $\rho > 0$. The shift enters through the factorization and
nowhere else: with $F_\rho = \tfrac1\alpha A + B^\top B + \rho I
= \tfrac1\alpha A_\rho + B^\top B$ and $A_\rho = A + \alpha\rho I$, the
system solved is **always the original one**. $\rho$ degrades the oracle;
its perturbation is left to the refinement loop, which is what refinement
is for. Every statement of this note is the unregularized one with $A$
replaced by $A_\rho$, and every formula below reduces to its $\rho = 0$
form as $\rho \to 0$ — a check the section is built to pass.

Two consequences of the framing, stated once. **The tests read the
original residuals.** A convergence test that subtracted the $\rho$-born
term would be blind to exactly the error the loop exists to remove, and
would converge to the shifted solution while reporting success on the
original; Proposition 2's rider $\alpha\rho\,x$ is therefore *signal*, not
contamination, and Appendix B's raw `s[j]` is correct. Subtraction is
admissible at the **sensor** and nowhere else — the one rule, both arms:
*subtract at a fresh solve to read the floor, never along the chain.* And
**the caller supplies a $\rho$ at which the factorization succeeds**: a
precondition on the input, of the same kind as full row rank. The
escalating search implementations run to find such a $\rho$ costs nothing
in this note's currency and is outside its scope; but a $\rho$ arriving
from such a search is the smallest that survived, so any breadth claim
resting on deployed records is conditioned on that selection and must say so.

Parameter discipline extends accordingly: the measured quantities are
**tolerance- *and* shift-free**, with both dials entering the formulas
explicitly.

**The dual arm: $\rho$ enters where machine epsilon does.** By
Proposition 2 the dual residual carries the exact term $\alpha\rho x$
alongside the rounding term of (14), so

$$\ell_d = \big(\rho + \rho_\star\big)\|x\|, \qquad
\rho_\star = c_u u \|B\|^2. \qquad (42)$$

Below $\rho_\star$ the shift is invisible; above it, $\alpha_c$ falls as
$1/\rho$. The threshold is **derived, not assumed** — it is a comparison
of two computable quantities, named as a branch at every use — and the
two scales are separately recorded (§16). They are *not* separable on
scalars: the converged original dual residual equals $\alpha\rho\|x\|$ to
four digits, so $\rho_\star\|x\|$ sits far below the noise in
$\ell_d - \rho\|x\|$. The separation is the **vector** subtraction of A6,
$\|w - \alpha\rho x\|/\alpha$, and only that. The
ceiling law (25) is *unchanged in form*: it reads the recorded $s_0$ as
before and returns a smaller number. The scaling row **(C1)** of §6 — not
the listing's payment line — is *strengthened*,
since $\alpha\rho\|x\|$ is exactly linear where the rounding term is flat
below $\alpha_\star$, so the switch (15) moves down to
$u\lambda_{\max}(A)/\rho$ and proportionality holds over a wider range.

Because $\alpha\rho x$ is an exact vector identity rather than a bound, it
holds in whatever norm (T2) employs, and the **zero-pass** dual endpoint

$$\alpha_c^{\rho} = \frac{\varepsilon}{\rho\,\|x\|_{\mathrm{rec}}}
\qquad (43)$$

is computable before $\alpha$ is chosen — the note's only endpoint that
need not be measured. It omits the rounding term and so overstates the
endpoint; it is sound for $\rho \gg \rho_\star$, and it is an endpoint,
not a wall.

**The wall's pole is unmoved; its arm is not.** The shift enters the
contraction twice. It raises the dual scale by (42), and it raises the
operator bounding the correction, since $Z^\top A_\rho Z = Z^\top A Z +
\alpha\rho I$ and hence $\lambda_{\min} = \lambda_Z + \alpha\rho$:

$$\varrho_d(\alpha) \;=\; \frac{\alpha\,(\rho + \rho_\star)}
{\lambda_Z + \alpha\rho}. \qquad (44)$$

Numerator and denominator are both $\Theta(\alpha)$; neither wins. So

$$\frac{d\varrho_d}{d\alpha}
= \frac{(\rho+\rho_\star)\,\lambda_Z}{(\lambda_Z + \alpha\rho)^2} > 0,
\qquad
\lim_{\alpha\to\infty}\varrho_d = 1 + \frac{\rho_\star}{\rho} \;>\; 1 :$$

the contraction **degrades monotonically** with $\alpha$, as it does
without the shift. What $\rho$ changes is that the degradation *saturates*
rather than growing linearly, and the saturation point is the **knee** at
$\alpha = \lambda_Z/\rho$, where $\varrho_d \approx \tfrac12$. The old
analysis's usable condition $\alpha\rho \lesssim \lambda_Z$ — a condition
on the *product*, so that the large $\alpha$ the primal side favours is
what makes a given $\rho$ expensive — is exactly this, and stands.

Setting $\varrho_d = 1$ gives $\alpha\rho_\star = \lambda_Z$: **the $\rho$
terms cancel exactly and the pole is unmoved from its $\rho = 0$ position**
$\lambda_Z/\rho_\star$. That is a statement about (44) and not about where
the arm becomes unusable — $\varrho_d$ has saturated near one decades
earlier. The pole is where it was; the arm approaching it is not.

**The arm, and it is one law.** Section 7 stopped at the contraction and
did not convert it into passes. It needs nothing new. Refinement must
carry the dual residual from $s_0 = \ell_d\alpha$ down to $\varepsilon$, a
log-range of $\log(s_0/\varepsilon) = \lambda_d + t - e = d_0 + t$ — the
ceiling numerator already in hand — so

$$\boxed{\;j(t) \;=\; \frac{d_0 + t}{\log\big(1/\varrho_d(t)\big)}\;}
\qquad (45)$$

At $\rho = 0$, $\varrho_d = \alpha\rho_\star/\lambda_Z =
\alpha/\alpha_{\mathrm{wall}}$, so $\log(1/\varrho_d) = T - t$ and (45) is
$(d_0+t)/(T-t)$ verbatim: the hyperbola of §8 is this law's unregularized
reading. Near the knee, where $\log(1/\varrho_d) \approx
\lambda_Z/(\lambda_Z + \alpha\rho)$, the same law reads

$$j(t) \;\approx\; (d_0 + t)\Big(1 + \frac{\rho}{\lambda_Z}e^{t}\Big),$$

exponential in $t$. **No constant is introduced**; the two shapes are one
formula evaluated in two regimes, and the linearisation is a reading of
(45) near the knee, not a substitute for it. And $F_\rho \succeq \rho I$,
so the *factorization* wall of §7 is gone — the shift removes that death
and steepens this arm.

**The primal arm: a bias, and the floor law it induces.** The identity of
§4 acquires a term that does not vanish with $\alpha$,

$$r_0 = -\frac{1}{\alpha}S_\rho\,\delta + \rho\,BF^{-1}x^\ast,
\qquad
b_\infty = \big\|B\,(B^\top B + \rho I)^{-1} x^\ast\big\|, \qquad (46)$$

so $\|r_0\|$ does not fall to zero but onto the floor $\rho b_\infty$, and
the entry test is unsatisfiable at any $\alpha$ once $\rho b_\infty \ge
\varepsilon$. **This empties the window; it does not fail the solve.** The
correction solves $B\,\delta x = r_0$ for whatever $r_0$ is, and it
removes the bias along with everything else; the cost stays finite and the
graded family survives. (A1) thus fails at *large* $\alpha$ as well as
small, and the window's nonemptiness condition (3) degrades on both
factors, becoming $\ell_p\ell_d \le \varepsilon(\varepsilon - \rho
b_\infty)$.

The graded family takes the substitution and nothing more. The level
condition $R_k(\alpha_k)\|r_0(\alpha_k)\| = \varepsilon$ reads

$$c_k\,\frac{\ell_p + \rho b_\infty\,\alpha}{\alpha^{k+1}} = \varepsilon,$$

which is the $\rho = 0$ condition with $\ell_p \to \ell_p + \rho
b_\infty\alpha$ — and that is exactly the recorded quantity
$\alpha\|r_0\|$. Hence (13) holds verbatim with a $t$-dependent span,

$$\nu_0(t) = \frac{M_f(t)}{t - t_\infty} - 1, \qquad
M_f(t) = \log\big(\ell_p + \rho b_\infty e^{t}\big) - e - t_\infty.
\qquad (47)$$

The pole is unmoved: $t_\infty$ is the logarithmic capacity of the
measure's *support*, and the bias changes which vector starts the Krylov
space — the weights — not the nodes. The term is convex on the domain, so
Proposition 3 and the corner construction of §13 are untouched. Where
$\rho b_\infty > \varepsilon$ the interpolant tends to a positive limit
decaying as $(\log\rho b_\infty - e)/(t-t_\infty)$, so $\lceil\nu_0\rceil
= 1$ at every large $\alpha$: the entry test never passes and one
iteration always suffices, which is what the machine does.

**What does not transport.** $S_\rho$'s spectrum caps at
$\sigma^2/(\sigma^2+\rho)$, so directions with $\sigma^2 \ll \rho$ are
never resolved, the regime triple's domain axis truncates there, and the
harvest's limiting measure is $\rho$-capped. The $c_k$ are therefore
**$\rho$-conditioned readings**: valid at the anchor's $\rho$, not
transportable across it. This is a transport restriction, not a new
parameter, and it is the binding constraint on (46) rather than the law
itself.

**The one purchase.** $\ell_p$ and $b_\infty$ enter the record only
through their sum: $\alpha\|r_0\| = \ell_p + \alpha\rho b_\infty$ is one
reading in two unknowns, and no manipulation of a single record separates
them. The zero-cost route is the a-priori bound $\rho b_\infty \le
\tfrac12\sqrt{\rho}\,\|x\|$, admissible as a bound and carried as one. The
reading-class route is one application of $F^{-1}$ forming $\rho BF^{-1}x$
at the **accumulated** pair — the raw backsolve is the wrong operand — and
this is the note's single admitted purchase, warranted by identifiability
rather than by convenience. It fires on the gate

$$\gamma \;=\; \frac{\tfrac12\sqrt{\rho}\,\|x\|}{\|r_0\|}, \qquad
\text{purchase iff } \gamma \ge \gamma_\star, \qquad (48)$$

one division and no application; below $\gamma_\star$ the bound suffices and
the split does not move $\ell_p$ appreciably. The threshold is a
calibration row (§17). Like every reading, $b_\infty$ is an $\alpha \to
\infty$ limit and a reading taken below the switch (15) is biased — the
same exclusion the dual arm already carries, now on both.

**In the cost**, the ascending arm's *slope* is proportional to $\rho$ by
(45) — not merely its position — while (46) raises the descending arm and
neither is lowered. The habitable interval closes from both ends, and the
minimum cost is monotone increasing in $\rho$. There is no interior optimum, and the cost model therefore has no
opinion on $\rho$ beyond *as small as the factorization tolerates* — which
is a decision made outside it, on grounds it does not price.

*A recorded instance.* On a power-flow sequence [RC, Fig. 4] a
Ruiz-scaled problem required $\delta_1 \approx 1$ to make the augmented
block definite and was abandoned there; in the terms above that is both
failures at once — a $\rho$ of order one collapses (43), and
$\rho b_\infty$ overruns $\varepsilon$ so no window exists at any
$\alpha$.

*A shift on the other block is not treated.* Implementations also
regularize the Schur complement itself, solving with $S + \delta_2 I$
when $B$ loses row rank. That shift moves every node of the measure on
which §5 reads its harvest, so it is not the $\rho$ of this section
under another name, and nothing here transports across it.

*Uses:* Prop. 2's rider ($\alpha\rho x$ exact); (13), (25); the purchase
discipline of §0; §17's calibration rows. Failure modes: a reading of
$b_\infty$ below the switch (biased, screened); the harvest transported
across $\rho$ (excluded); the bound consumed as a value (untypable, §13).

## 16. What to record

The estimate consumes nothing outside this list; an implementation
recording it makes every route of §11 available and every check of
§14 executable. Quantities per solve unless marked per pass or per
invocation.

    anchor       the alpha the solve ran at, and epsilon, and the policy
                 (tau) in force; every returned object is conditioned on
                 the last and anchored at the first
    residuals    ||r_0|| before any correction; the dual s_j after each
                 pass, as the two-block residual's first component —
                 per pass, since the ratios and the floor shares of §11
                 are formed from consecutive pairs
    counts       per invocation: the Krylov iterations run, and the
                 MET-FLAG — whether the invocation exited by meeting its
                 test; the flag is (A-term)'s per-solve check, and a
                 record without it cannot be qualified (§14, step 1)
    entry gaps   per pass, theta-truncated under §1's interior
                 discipline: the initial constraint residual the pass's
                 (T1) tested, for the branch of (29) and the surcharge
    halt         if (T3) fired: the final ratio, the final s, and the
                 alpha — the classification of §11 needs the middle one
    coefficients (a_i, b_i), i = 1..n_0: the harvest; a_1 alone gives
                 the elasticity vartheta = 1/a_1, and the levels follow as
                 c_k = prod e_i with e_i = alpha sqrt(b_i)/a_i
    norms        ||x|| and ||y|| of the assembled pair, for phi (23),
                 the spectral bound, and the amplifier caps (19)
                 feeding validator and clip
    vartheta     recorded when n_0 >= 1; a meter on (A1) at the anchor,
                 never an exponent applied to transport (§5); read at the
                 anchor.  NOT theta, which is §1's interior factor
    s_comp       a compensated (per-row fsum) re-evaluation of the
                 final dual residual, zero applications; the floor
                 meter is the VECTOR difference ||w - w_comp|| (a
                 difference of norms cancels at large residuals),
                 read at the earliest large-residual pair and invalid
                 below the resolution u*s; s_comp reads the pair's
                 true residual, and the pair is the measured floor
                 that classification, the corners (§7), and (A-cw)
                 consume
    k0           the inner startup transient — iterations to the first
                 decade from a fresh start — consumed by (32)'s
                 interior-toll constant
    demand_j     the realized per-pass reduction log(||r0_j||/eps_j),
                 consumed by the interior discipline's verification
    lam_fire     the firing-interior log-rate, operand of c_theta
    rho, purchase   the regularization and whether the one purchase fired
    lp_rho       the decontaminated entry scale, when purchased
    ld_sub       the SUBTRACTED dual scale ||w - alpha*rho*x||/alpha, the
                 floor scale rho_star*||x||.  Distinct from ell_d, which
                 is the RAW s_0/alpha that (T2) tests and (43) consumes;
                 the two differ by rho*||x|| and are not separable on
                 scalars (§15), so both are recorded
    ratio_1      the first pass ratio s_1/s_0, the wall route's reading
                 (§11) and, under rho > 0, the lambda_Z validator of §10
    exit_type    the invocation's typed exit (reached-tolerance, floor,
                 declared error), consumed by classification and the
                 decline shapes

Each recorded wall or level quantity carries its provenance as §10
types it — reading, fit, or bound, with direction — and the floor's
provenance distinguishes harvested from extrapolated from clamped by
the guard; with the anchor recorded, a clamped floor is readable as
*below the anchor*, which a bare flag does not convey. Nothing here
requires an extra solve or an application of $S$; the list is what the
algorithm already computes, kept instead of discarded.
## 17. Assumptions

This catalog is *generated*: the consumer column is assembled by scan
from the uses-lists of §§1–15 (`gen_catalog.py`, part of the note's
build), so a derivation cannot consume a hypothesis without acquiring a
row, and a consumed row cannot outlive its last consumer; rows kept
without a consumer — (A2), (S) — are marked *unconsumed*,
documentation of the road not taken. Hand edits to the
consumer column are prohibited; the remaining columns are maintained
here. **Coverage** classifies each row: *discharged* — proved from
others in a stated regime; *checked* — carries a per-solve meter the
algorithm evaluates; *trusted* — consumed without a per-solve check,
the honesty residue. **Breadth** states for which problems the
assumption is believed to hold, as a characterization or as a named
gap, never as silence.

| Name | Statement (short) | Coverage | Breadth | Failure mode | Consumers |
|---|---|---|---|---|---|
| (A-term) | An entered invocation exits by meeting its test or by declared error | **checked** — met-flag per invocation | property of the implementation, not the problem; violations are defects | counts stop meaning test-meeting work; passes lose their dual-born standing; §7's chain collapses | §1, §7, §11, §12, §14, §16 |
| (A1) | $\|r_0\| \propto 1/\alpha$ between anchor and endpoint | **checked** — $\vartheta = 1/a_1$ meters it at the anchor | generic under (A3)'s mechanism; degrades below the anchor, one-signedly | floor endpoint placed low; census sign 1 | §4, §5, §16 |
| (A2) | (A1) globally | trusted (unconsumed by derivations) | as (A1); stated for cleanliness | — | §4 |
| (A3) | gap mass on resolved directions | trusted; explanatory tier | operator-and-data property; fails when the multiplier gap excites the tail | (A1) inherits the failure; constants unnamed | §4, §5 |
| (A-mon) | the $k$-step reduction non-increasing in $\alpha$ (up-set form of $W_k$) | trusted | exact for $\|r_0\|$ by (8); asymptotic for $R_k$ within (11)'s triple | tread intervals mis-shaped; enumeration repairs locally | §5 |
| (A-dec) | the decontamination operand is the accumulated pair, whose distance to $x^*$ is $O(\rho)$ | **checked** — the raw-backsolve operand measured failing at twice the floor | any implementation following §15's purchase discipline | floor under-corrected; entry endpoint misplaced | §15 |
| (C1) | $s_0 \propto \alpha$ between anchor and endpoint | **discharged** above the switch by (C3) | all problems, above (15) | ceiling endpoint from a flat signal; §11 excludes the zone | §6 |
| (C2) | (C1) for all $\alpha > \alpha_\star$ | **discharged** by (C3) | as above | — | §6 |
| (C3) | backward-error mechanism (14) | **discharged** — derivation | backward-stable solves | — | §6 |
| (S) | two-scale spectrum on a near-kernel | trusted; mechanism layer, consumed by no derivation | the conditioning picture's hypothesis; retained for exposition | none — nothing returned depends on it | §6 |
| (A-pass) | passes start at $y_0 = 0$ and run the full sequence | **checked** — implementation inspection, once | property of the implementation | multiplier bound (19) misapplies; skip lemma misfires | §7 |
| (A-entry) | at a firing pass, $\|p\| \le \varepsilon$ and $\|d\| \le \ell_d\alpha$ | first clause **discharged** under (A-term); second **checked** on the ramp, (T3) disposing of the rest | all in-policy solves above the switch | pass entry residual mispriced; (29) wrong branch | §7 |
| (A-skip) | $\alpha_c \ge C_d + C_p$ through (19)'s caps | **checked** — one inequality per solve; discharged outright when the window is nonempty | spectral quantities clearing a $\sqrt{u}$ scale; same character as full row rank | pass cost exceeds one application; budgets understate; optimum biased high | §7, §8 |
| (A-const) | one $c_u$ bounds base and pass chains | trusted — $O(1)$ scatter absorbed | uniform-arithmetic implementations | recursion constants off by $O(1)$; folded into $\kappa_r$ | §7 |
| (A-norms) | (23)'s norms are the converged pair's | trusted | corrections small against the pair — in-policy convergent solves | $\varphi$ mis-scaled; screens and classification degrade | §7 |
| (A-cw) | componentwise slack of (23) modest, within-a-decade uniform | **checked** — the measured floor of §16 (compensated re-evaluation, $\varphi$-independent) and every verification outcome | floating-point implementations with standard elementwise arithmetic | $\varphi$ mis-scaled; classification and the corner deductions degrade | §7, §11, §12, §13, §14 |
| (A-restart) | sibling pass style: re-solve original data at the current multiplier | **checked** — implementation inspection, once | implementations exposing a start argument | surcharge mispriced conservatively if (33) applied unhedged | §7 |
| (A-flat) | realized contraction $\propto \alpha$ between reading and wall | trusted — **the note's one named gap** (TODO-F) | families with a proportional regime; families exist without one, and whether they lie inside or outside the class is the pending verdict | one-signed: transported readings overstate the wall — the unsafe side; the clip exists for this | §7, §8, §10, §11, §12, §13 |
| (A-crowd) | the crowding law's continuum form (13) | **checked** — the fit residual meters it; disqualification is a decline | $k \ll N$, interval-like support; fails on clustered spectra, which the meter catches | floor arm misplaced; §11 declines the upper band | §8 |

Three reading aids. The **trusted residue** is seven rows —
(A-const), (A-norms), (A-flat), (A-mon), and the unconsumed
(S)/(A2)/(A3) tier — and of these only (A-flat) has a named breadth
gap; the other trusted rows are
uniform-arithmetic or smallness statements whose failure is bounded and
folded into constants. Every *checked* row's meter is in the record of
§16, so checking costs no extra solve. And the census of §5 is this
table's signed shadow: where several rows fail together, their
directions compound as recorded there, and no cancellation is assumed
anywhere in the note.

**Calibration block.** The tuning constants, each with provenance
and sensitivity; a constant absent from this block is not in the
note.

| constant | value | provenance | sensitivity | consumers |
|---|---|---|---|---|
| interior factor $\theta$ | $1/10$ | §1's design: one decade of closure per interior | $c_\theta$ linear in $\log(1/\theta)$ | §1, §8 |
| stall threshold $\tau$ | policy | the wall is policy-indexed; placement jitter a few tenths of a nat under implementation detail | one-signed: stricter policies lower the operative wall | §1, §7, §12 |
| corner deductions | $\varepsilon - 3f_m$, $\varepsilon - \tfrac12 f_m$ | companion report | calibrated against the measured floor | §7 |
| purchase gate $\gamma_\star$ | $1/10$ | 192 configurations, four operator families, six $\rho$: at $1/10$ no case whose $\ell_p$ error exceeds one nat is missed, and no case whose error is below $0.1$ nat fires | looser thresholds add waste only; the free bound alone reports an empty window on 44% of configurations where the true bias is negligible | §15 |
| derived-clip gap | 0.8–2.4 nats | the derived bound exceeds the operative (policy) wall by this margin on the test instances | consumed as a labeled caveat, never subtracted | §12 |
## 18. Open items

Eleven debts, registered by name where incurred, each with its
resolution criterion and its edit locus — the text a resolution
changes, and nothing else.

**TODO-K** ($\kappa_r$, §8, §13). The rate quarantine and the
single-$\bar\lambda$ surcharge carry a slack constant, bounded in form,
unsized. *Criterion:* the measured spread of realized interior rates
over the deposit region — of which the explored-levels component is
computable now, the harvest's coefficients bounding realized rates
over the levels the solve visited; the unexplored remainder is the
production half. *Edit:* a number in (34), (39), (41).

**TODO-E** ($e_C$, §9, §13, §14). The certificate's estimation term is
bounded through band widths, unsized in deployment. *Criterion:* the
distribution of inner-band widths at the minimizer over deployed
records. *Edit:* a number in (37); possibly a tightened default in
step 7.

**TODO-F** ((A-flat)'s breadth, §7, §11). The note's one named gap:
families exist on which no proportional contraction regime is
observed. *Criterion:* the contraction-shape reading across the
deployed family — flat transported walls resolve it inward, cliff-like
shapes resolve those families outward, as scope. *Edit:* (A-flat)'s
breadth cell in §17, and possibly the ratio route's standing in §11.

**TODO-U** (certified upper wall bound, §7, §9). Only the lower bound
(22) is certified; a certified upper bound would upgrade the
certificate (37) from robust to unconditional. *Criterion:* a
derivation. *Edit:* the certificate paragraph of §9; the wall row of
§17.

**TODO-P** (the objective's primacy, §8). $C$'s adoption as canonical
is provisional against exact enumeration. *Criterion:* the
pre-registered comparison of §9's guarantees against enumeration on
deployed records. *Edit:* one paragraph of §8, as stated there.

**TODO-R** (regularization residue; the main debt is paid in §15).
Remaining: the grade-(iii) harvest theory on the $\rho$-capped limiting
measure — saturation's effect on the $c_k$, which is what bounds (46)'s
transport rather than the law itself; and the skip analysis under
$A_\rho$. *Criterion:* each derived with the $\alpha\rho x$ term
carried. *Edit:* §15's transport paragraph; §7's skip lemma.

**TODO-C** (ratio composition, §7, §11). The route takes the first pass
ratio and inverts it; §7's reading-class top names the deconvolution as
this debt's destination, not as the present method. Where several screened ratios exist, the drift
$s_{j+1}/s_j \approx \varrho + \varphi/s_j$ of (28) admits a
deconvolution — regress the ratios on $1/s_j$; the intercept is the
contraction and the slope an independent floor meter — which returns
$\varphi$ as a by-product and does not merely screen the contamination
but identifies it. Adopting it requires settling the minimum sample and
the admissibility of one- and two-ratio records, neither of which is
done; a two-parameter fit on correlated points, worst-conditioned exactly
near the wall where the residuals cease to spread, is not obviously
better than the first reading. Four points is the smallest sample at
which such a fit can be *checked* rather than merely computed — two
determine the line exactly — which is where that figure comes from and
why it is a property of the deferred method rather than a calibrated
constant of the shipped one. *Criterion:* the minimum sample and the
admissibility of one- and two-ratio records, together. *Edit:* §11's
second route and §7's bracket.

**TODO-G** (containment of the lower curve, §5, §11, §13). The
certified-optimistic curve is not, in general, a lower bound, and two
distinct mechanisms are responsible.

*Below the anchor.* The signed census records transport of the measure as
optimistic and consumes that sign universally, which is what licenses the
lower curve to take a transported floor ungated. The sign is not
universal: the census's optimism dominates near the anchor, and the
crowding law's continuum form breaking against the *terminating* family
dominates at depth, erring the other way. Where the second dominates the
licence lapses. The remedy is an **assertion range** on the box — a limit
on how far a transport hypothesis is carried, shipped with the parameters
and outside which the estimate declines. It must be a statement about the
parameters, not a patch on either curve: a substituted value would be
inexpressible in the box, would break convexity, and would give a
consumer evaluating the formula a different function than the estimator's.

*Above the anchor.* No transport occurs in that direction, so no
assertion range reaches it. The mechanism is the wall: read from a ratio
whose floor share is appreciable, it comes out low, which puts a pole
where the exact cost is flat and takes the lower curve with it. The
direction of that error is determined by the floor share, hence by
$\varphi$; the estimate does not carry $\varphi$ into the wall route, and
so cannot distinguish a reading a little high from one badly low.

A third consequence, accepted rather than deferred: with the derived
spectral bound excluded as a clip source (§12), a record with neither a
halt nor a pass has no certified wall endpoint above its own $\alpha$,
and its clip is the anchor. That is a loss of function, not a defect —
the estimate declines to recommend up-parameter from a solve that
witnessed nothing up-parameter — and it is stated in §12 rather than
patched.

*Criterion:* a derived or metered depth threshold for the first; for the
second, either a route that is one-signed by construction or a typed
admission of $\varphi$. *Edit:* §5's census sign becomes local; §11's
floor route gains the assertion range; §13's black-curve claim is
qualified there already.

**TODO-T** (the surcharge term, §7, §8, §9, §13, §17). The interior
surcharge and its apparatus — the two-branch entry law, the crossover,
the toll constant, the firing-interior rate, the startup transient, the
exit-toll term of (32), and the slack constant $\kappa_r$ — price a
quantity that has not been observed to be nonzero. A backward-stable
factorization commits its error where the recovery step cancels it, which
is Proposition 2's own content. *Criterion:* whether any regime in scope
makes the term nonzero; the untested candidate is genuinely
reduced-precision factorization, as against a badly conditioned working
precision. *Edit if it resolves toward deletion:* (31) loses $\Sigma$;
(32) loses $E$ and $c_\theta$; (35) collapses to the two-pole balance;
§17 loses the rate rows and the interior-factor row; **TODO-K and
$\kappa_r$ disappear, so (39) and (41) lose their unsized constant** —
the largest simplification available to the note. *Standing note
meanwhile:* the surcharge is priced but has not been observed, and every
quantity it introduces is one the note would otherwise not carry.

**Extended precision, and what it would buy.** The limiting accuracy
$\varphi$ of (23) is a componentwise *model* of the evaluation's
rounding. It can instead be measured. Evaluating $f - Ax + B^\top y$ with
the multiply-accumulate over all $n + p + 1$ terms carried in compensated
(double-double) form, and differencing against the ordinary evaluation as
a *vector*, returns the rounding that evaluation committed — no
applications of $F^{-1}$, and the same construction applies to the
constraint residual $g - Bx$ for the primal-side floor. Three things are
worth recording about it. Its cost is a dependency chain, not a flop
count: compensated accumulation cannot be reassociated, so the operation
does not reduce to a library matrix–vector product, and the constant is
accordingly large — of order a hundred plain evaluations at moderate
dimension, which is nothing against a factorization and not nothing
against a matvec. Recovering the individual products' errors as well as
the accumulation's changes the measured value by under two percent even
under fourteen orders of cancellation, so the accumulation term dominates
and the partial form captures essentially all of it. And what is measured
is the *evaluation's* rounding; its use as the refinement floor is
(A-cw)'s bridge, not an identity. The estimate as specified does not
consume it, and none of §14's steps depend on it; it is recorded because
the quantity it would supply is the one standing between an observed
ratio and a wall reading (**TODO-G**), and because a caller who can
afford it can check (A-cw) directly rather than trusting it.

**Beyond the registry**, two mathematical opens that change no
returned quantity. Proposition 3's convexity is now unconditional for
the smooth layer; the residual open is the exact staircase's
block structure under the toll terms, and the enumeration backstop
covers it, so the open is aesthetic. And on operators with an isolated near-null reduced
direction, the realized wall is conjecturally a bulk-spectrum object —
early passes blind to a direction of generic weight $1/\sqrt{n}$ — in
which case (27) is loose there by exactly the isolation; the bound's
validity is unaffected, and the conjecture is recorded because it
would explain the looseness rather than merely bound it.

**TODO-W** — the warrant program: for each empirically-warranted
reading, a target theorem (circumstances + bound); rows: the
wall-intercept margin under (A-flat), the bracket-midpoint quality,
the retired pivot (re-entry only by proven warrant).

**TODO-RS** (stall rescue; §1's appendix boundary). A stalled solve may
be retried with the backsolves themselves refined —
factor-preconditioned CG on $F$ — which revives wall-adjacent stalls at
bounded cost and *moves the wall*; it is outside the ideal machine,
flagged in the record where used, and unpriced here. *Criterion:* the
rescued-solve cost model and the halt-record retyping, at production
scale. *Edit locus:* §1 (appendix note), §11 (halt semantics), §16
(rescue flag).
## 19. Summary

For a saddle-point system solved by augmentation, Krylov correction,
and refinement under the stated policy, the cost of a solve is the
count (2) of applications of $F^{-1}$, an extended-real staircase in
the parameter. Its two arms are priced by §§4–7: the floor by the
exact identity (4), the elasticity (8), and the crowding law (13) with
parameters $(t_\infty, M_f)$; the ceiling by the exactness of the dual
row, the backward-error ramp (14), and the affine recursion (24),
whose level family (25) is corrected at the effective tolerance,
capped at (26), and bounded by a wall that splits honestly into a
derived lower bound (22) and a measured realized wall (27). Interior
iterations obey the two-branch entry law (29): absent above the
crossover — the classical two-term account, recovered as a theorem —
and a priced, policy-conditioned surcharge below it.

What a solve can know is the relaxation (32): two hyperbolas and a
surcharge in tolerance-free parameters, convex unconditionally by
Proposition 3, under-counting pointwise, within the sandwich (34) of what
is paid. Its minimizer is closed-form (35); the recommendation is a
tread, ties breaking downward by arithmetic; and two lines of counting
give the guarantees — the self-scaling bound (36), and the one-solve
certificate (37) in robust form.

What the algorithm returns is the estimate itself: parameters as
intervals with typed endpoints, composed by intersection; the exact
corner envelopes (38) with the corridor (39); the minimax
recommendation (40) carrying its promise, its possibility, and its
width with the width's own remedy; a verification solve whose three
outcomes — certified, hypothesis-failed with sources named, withdrawn
— close every path; and, beneath every failure mode, the witnessed
configuration. The window $W_0$ of (3) survives as the unit level set
of the object the note is actually about, and the old question —
where is the window — survives as the first query one asks of the
answer to the better one: what does every parameter cost, and how
certainly.

Every claim above is derived under the registered rows of §17 — of
which the trusted residue is small and only (A-flat) carries a named
breadth gap — or is deferred under the registered debts of §18. Nothing else is
assumed, and nothing observed is reported here.

---

## Appendix A. The corridor and the regret accounting

**A.1 The sandwich (34).** Each staircase sits within one rounding of
its interpolant per application — the count is the ceiling of the
real-valued level crossing — contributing $1$ each through $n_0$ and
$j$; the surcharge's triangle adds one more through its own ceiling;
and the quarantine and deficit contribute $\kappa_r$. Hence
$C \le \mathrm{cost} < C + J_p + 2 + \kappa_r$, with the corollaries:
hump depth above the convex envelope is bounded by the local budget
$J_p + 2 + \kappa_r$ — tight far from the wall, loose near it, where
the enumeration is mandatory; and on the ramp branch the sandwich
collapses to the additive $+2$, hump depth one.

**A.2 A-priori regret.** With $t_C = \arg\min C$ over the feasible
domain, $\Delta(t) = J_p(t) + 2 + \kappa_r$, and $t_\ast$ the true
minimizer:
$\mathrm{cost}(t_C) \le C(t_C) + \Delta(t_C) \le C(t_\ast) +
\Delta(t_C) \le \mathrm{cost}(t_\ast) + \Delta(t_C)$, so the regret is
at most the corridor width *at the recommendation* — which lands in
the valley where $J_p$ is small; near-wall corridor bloat never
enters.

**A.3 Displacement.** Where $C$ is $m$-strongly convex on the bracket,
$|t_\ast - t_C| \le \sqrt{2\Delta(t_C)/m}$, with the curvature in
closed form, $C'' = 2M_f/(t - t_\infty)^3 + 2W_{\mathrm{eff}}/(T-t)^3$,
in nats. The bound is global now that Proposition 3 is unconditional;
the exact layer's enumeration still governs the final tread.

**A.4 The certificate, pure and robust.** The pure form
$R \le \mathrm{cost}(t_C) - C(t_C)$ uses only the lower arm
($\min\mathrm{cost} \ge \min C = C(t_C)$) and one verification solve —
but it presumes the lower arm trustworthy at $t_C$, which without a
certified upper wall bound (TODO-U) it need not be. The robust form
adds the components' estimation error: $R \le \mathrm{cost}(t_C) -
C(t_C) + e_C$, with $e_C$ bounded through the inner band width
(TODO-E). The note ships the robust form (37).

**A.5 Self-scaling (36).** Term-nonnegativity gives
$J_p(t_C) \le C(t_C) - 1 \le \min\mathrm{cost} - 1$; feeding this into
A.2, $\mathrm{cost}(t_C) \le 2\min\mathrm{cost} + 1 + \kappa_r$. One
inequality, two readings: additive and sharp on cheap problems, factor
two and safe on expensive ones; the factor is attained. Per instance,
balancing the hyperbolas gives $J_p(t_C) = (M/S)(1 + 1/r)$ with $M$
the ceiling span, $S = T - t_\infty$, $r$ as in (35) — a certificate
computable from the recorded parameters at recommendation time.

*Uses throughout:* term-nonnegativity of (32); the sandwich's
hypothesis stack (§8); strong convexity on the bracket for A.3;
TODO-U for A.4's form; TODO-K and TODO-E where their constants
appear. Estimation-robust versions of A.2 and A.5 carry $+2e_C$.

---

## Appendix B. The algorithm as a listing

The solver of §1 and the loop of §7, rendered as one listing. Payments
— applications of $F^{-1}$ — occur at lines A1, A4, B2, B3, and the
optional P1; everything else is free, and the declared arrays are
exactly the record of §16. The $\rho$-guarded lines are
*algorithmically* complete and included here because the control flow
is settled even where the pricing is not: the refinement loop is
identical for $\rho = 0$ and $\rho > 0$ by the following, whose cost
analysis is §15's.

**Proposition (P-tel).** *Let $(dx, dy)$ be the machinery's exact
output on the residual pair $(d_j, p_j)$ of an accumulated pair, and
$(x_{j+1}, y_{j+1}) = (x_j + dx,\, y_j + dy)$. Then
$f - Ax_{j+1} + B^\top y_{j+1} = \alpha\rho\,dx$.* Proof: Proposition
2 applied to the pass's own invocation gives
$d_j - A\,dx + B^\top dy = \alpha\rho\,dx$; substitute
$d_j = f - Ax_j + B^\top y_j$ and collect. $\square$ Each pass's
leftover is proportional to its *own* increment, so the raw residual
is the convergence object in both regimes and (T2) needs no guard;
the subtraction at A6 is the sensor's tool only, valid at a fresh
solve.

```
ALGORITHM AW.  (Augmented solve with refinement; shifted
     factorization handled in-line.)

     Given A (n×n, PSD), B (p×n), data (f, g), α > 0, ε > 0,
     stall ratio τ, pass cap jcap, iteration cap nmax, and a
     shift ρ ≥ 0 chosen at factorization time.  β := 1/α;

          F := β·A + BᵀB + ρ·I,

     factored once, outside the currency.  Payments are the applications of
     F⁻¹: lines A1, A4, B2, B3, and the one optional line P1.
     All else is free.

     declarations
        real array  a[1:nmax], b[1:nmax];      comment harvest;
        real array  s[0:jcap];        comment dual residuals,
                                       one per pass, s[0] base;
        integer array  niter[0:jcap];  comment F⁻¹ counts per
                                       invocation, niter[0]=n_0;
        real  ℓ_p, ℓ_d, φ, ϑ;          comment scales, limiting
                                       accuracy, elasticity;
        boolean array  met[0:jcap];    comment (A-term) flags;
        real array  gap[1:jcap];       comment pass entry gaps;
        real  halt_s, halt_ratio, halt_α;  comment the halt triple;
        real  xnorm, ynorm;            comment accumulated norms;
        real array  x, dx, r, w, v[1:n];  y, dy, p[1:p];
        integer  j, n_0, k;


procedure SOLVE (f, g, y_0, ρ);
  begin
     comment ── base invocation ──────────────────────────────;

A1:  r := g;  if y_0 ≠ nil then r := r + β·y_0 fi;
     x := F⁻¹(β·f + Bᵀr);

A2:  r := g − B·x;
     ℓ_p := α·‖r‖;
     comment  under ρ > 0 this ℓ_p is the SUM ℓ_p + α·ρ·b_∞ of
     comment  §15: one reading in two unknowns.  The split, if
     comment  bought at all, is bought at C1 — AFTER A4, since
     comment  the operand is the accumulated pair and this x is
     comment  still the raw backsolve;

A3:  if ‖r‖ ≤ ε then
        n_0 := 0;  met[0] := true;
        comment  an unentered invocation satisfies (A-term)
        comment  vacuously;
        go to A5 fi;

A4:  n_0 := 0;
     while ¬(‖g − B·(x+dx)‖ ≤ ε) do
        begin  craigstep(dx, dy);  n_0 := n_0 + 1;
               store (a[n_0], b[n_0])
        end;
     x := x + dx;   niter[0] := n_0;
     met[0] := (‖g − B·x‖ ≤ ε);   ϑ := 1/a[1];
     if ρ > 0 then
        comment  the ρ-bias split, at the ACCUMULATED pair.
        comment  gate γ := ½·√ρ·‖x‖ / ‖r‖, zero applications.
        comment  NOT g, which is the constraint right-hand side:
        if γ < γ⋆
           then widen the ℓ_p interval by ½·√ρ·‖x‖
                comment  bound-typed, the default;
        else
P1:          r := (g − B·x) − ρ·B·F⁻¹(x);
             ℓ_p := α·‖r‖;   lp_rho := ℓ_p
             comment  the one purchase: ONE payment, reading-
             comment  typed, operand the accumulated pair —
             comment  the raw backsolve leaves twice the floor
        fi
     fi;
     comment  CRAIG unchanged by ρ: it solves B·dx = r for
     comment  whatever r is; the shift lives inside F⁻¹ only;

A5:  r := g − B·x;
     y := y_0 + α·(dy + r);

A6:  w := f − A·x + Bᵀy;
     s[0] := ‖w‖;
     ℓ_d := s[0] / α;              comment RAW scale: what (T2)
                                   comment tests, and what the ceiling
                                   comment endpoint and (43) consume;
     if ρ > 0
        then ℓ_d_sub := ‖w − α·ρ·x‖ / α
             comment  the sensor subtracts the exact part
             comment  α·ρ·x — an axpy, x in hand, free.  This is
             comment  the FLOOR scale ρ⋆·‖x‖, NOT ℓ_d: the two
             comment  differ by exactly ρ·‖x‖ and are not
             comment  separable on scalars (§15).  Record both.
             comment  valid at a fresh solve only (P-tel):
             comment  do NOT subtract along the chain;
        else ℓ_d_sub := ℓ_d
     fi;
     xnorm := ‖x‖;  ynorm := ‖y‖;
     comment  a compensated re-evaluation of w here — zero
     comment  payments — supplies §16's measured floor, as a
     comment  vector difference;

     comment ── refinement loop ──────────────────────────────;
     comment  identical for ρ = 0 and ρ > 0: the raw residual
     comment  is the convergence object in both — each pass
     comment  leaves α·ρ·(its own increment), which contracts
     comment  (P-tel), so (T2) needs no guard;

A7:  j := 0;
     while s[j] > ε do
        begin
           if j ≥ min(jcap, j_floor) then go to HALT fi;
B1:        w := f − A·x + Bᵀy;   p := g − B·x;
B2:        dx := F⁻¹(β·w + Bᵀp);
B3:        gap[j+1] := ‖p − B·dx‖;
           if gap[j+1] > ε then
              while test not met do craigstep;
                             niter[j+1] := niter[j+1] + 1 od
           fi;
           met[j+1] := true;  comment or declared ERROR (A-term);
B4:        dy := α·(dy_c + (p − B·dx));
           comment  dy_c is B3's CRAIG dual, zero when no
           comment  craigstep ran.  §1's recovery is
           comment  y_0 + α·(δy + r); dropping δy is correct
           comment  only inside (A-skip), where B3 does no work;
           x := x + dx;  y := y + dy;  j := j + 1;
           s[j] := ‖f − A·x + Bᵀy‖;
B5:        if s[j] > τ·s[j−1] then go to HALT fi
        end;
     cost := 1 + j + Σ_k niter[k];   return CONVERGED;

HALT: halt_s := s[j];  halt_ratio := s[j]/s[j−1];  halt_α := α;
     classify halt_s against the measured floor (φ as fallback);
     comment  wall-type (s[j] ≫ φ) vs floor-type (s[j] ∼ φ);
     comment  if ρ > 0 and α·ρ ∼ λ_Z: knee regime, contraction
     comment  saturates below one — expect smeared wall-type
     comment  halts; classification unchanged, its reading is;
     return STALLED
  end SOLVE;


     Remarks.
     (i)   ρ appears in the factorization and in one demand-priced
           payment (P1), which is a METER and fires after the
           accumulated pair exists — not on the raw backsolve of A1.
           The tests read the raw residual; see §15.
     (ii)  With ρ > 0 the zero-pass ceiling endpoint is
           computable before any α is chosen:
                α_c(ρ) = ε / (ρ·‖x‖).
     (iii) The declarations are the record of §16 entire:
           (a, b) the harvest; s[·] and gap[·] the ratio and
           surcharge routes' food; niter[·] and met[·] the
           counts with their (A-term) flags; the halt triple;
           (ℓ_p, ℓ_d, ϑ, xnorm, ynorm) the scales,
           meter, and bound routes.  The procedure arguments
           (α, ε, τ, ρ) supply §16's anchor row; nothing is
           recorded twice.  The estimate reads the record and
           pays nothing; verification pays its own bill.
```

Two unconditional payment sites (A1, B2), one conditional (the
craigsteps of A4 and B3, wherever the entry gap is positive), one
optional and demand-priced (P1). Recovery, tests, measurements,
screens, bands, and the recommendation itself ride free on vectors
the payments already produced.

## Appendix C. Unpacked readings

The body is written compressed. That is a choice, and it has a cost: a
sentence carrying three claims in apposition is hard to *check*, and
checkability is the whole point of §0's contract. This appendix pays the
cost back for the sentences that carry the most weight. Each entry quotes
the compressed form and then says the same thing at length. Nothing here is
new; where an unpacked reading and the body disagree, the body governs and
the disagreement is a defect to report.

**C.1 — §0, the two arms.** *"everything expensive happens at the ends, and
the cheap middle is where the two arms are simultaneously small."*

Turn the parameter down and the initial guess is bad, so the iterative
correction has to work: you pay Krylov iterations. Turn it up and the
arithmetic can no longer deliver the accuracy in one shot, so you pay
refinement passes. Both costs blow up at their own end, and neither is
large in between. The cheapest parameter is therefore an interval in the
middle, not a point, and finding it means balancing two growing costs
rather than minimizing one.

**C.2 — §0, the adversarial quantifier.** *"Every rounding statement in this
note is quantified adversarially over all $\delta$-patterns within that
bound — which is what makes the results independent of instruction order,
fusion, and code path."*

Floating-point results depend on the order operations happen in, which the
compiler and hardware may change. Rather than model that, every statement
here holds for *any* pattern of rounding errors within the standard bound.
So a real machine's run is one case among those covered, and the results
transfer between machines without re-deriving them.

**C.3 — §0, black and white.** *"at the wall the certified content is
natively a black/white pair, and every interior guess is a measurement."*

Near the upper end the estimate cannot produce a single number, only a
range: a lower curve nothing can beat, and an upper curve nothing will
exceed. Any number picked from inside that range is not a deduction; it is
someone's guess, and it must be labelled as one.

**C.4 — §2, the staircase.** *"the minimizers form flat plateaus and no
point inside a plateau is distinguishable from its neighbors."*

Cost is a count of solves, so it is an integer and it changes in steps. Over
a whole range of parameter values the count is identical. Reporting a single
best parameter would therefore be false precision: the honest answer is the
range over which the count is the same.

**C.5 — §4, the resolved and the tail.** *"Directions with
$\alpha\sigma_i^2 \gg 1$ are **resolved** by the augmentation; those with
$\alpha\sigma_i^2 \ll 1$ are the **tail**."*

The augmentation fixes some directions of the problem and leaves others
untouched, and which is which depends on the parameter. Raise it and more
directions come under control. The ones still outside are what the initial
residual is made of.

**C.6 — §5, why the levels crowd.** *"the levels descending geometrically
toward $t_\infty$, increments falling as $1/k^2$, each iteration buying less
than the last."*

Each extra Krylov iteration lets you get away with a smaller parameter, but
by a shrinking amount. The parameter values at which the count drops by one
bunch up as the count grows, approaching a limit they never pass. So a few
iterations buy a lot of room and many iterations buy little.

**C.7 — §5, the harvest.** *"The family's constants and the quadrature are
one object seen twice."*

The numbers needed to predict how much each future iteration would buy are
already produced by the iteration that ran — they are its own coefficients,
read differently. Nothing extra is computed; something already computed is
kept.

**C.8 — §6, why the dual residual is rounding.** *"every exact-arithmetic
contribution vanishes by construction."*

The recovery step is built so that in exact arithmetic the first block row
of the system is satisfied exactly. Whatever is left over when you measure
it is therefore entirely the arithmetic's own error. That makes the leftover
a *sensor*: it reads how accurately this machine can solve this problem, and
it is the only direct reading of that quantity available.

**C.9 — §7, the two-sided wall.** *"a *reading* (a measurement of the
quantity, admissible as a value) or a *bound* (a one-sided certificate,
admissible only as a constraint)."*

Some numbers are measurements of the thing you want; others only say the
thing is above or below some level. Using a one-sided fact as if it were a
measurement is how you get an answer that looks certified and isn't. The
type system exists so that the confusion cannot be written down.

**C.10 — §7, the limiting accuracy.** *"a vector norm of componentwise
products, not a product of norms."*

Two ways to bound the same rounding error. One multiplies the sizes of the
whole matrix and the whole vector; the other multiplies entry by entry and
then measures. On problems whose entries span many orders of magnitude —
which is the interior-point case — the first is far too pessimistic, and the
difference matters because this quantity decides feasibility.

**C.11 — §8, what a relaxation is for.** *"The staircases are what a solve
pays; they are not what a solve can know."*

The real cost is a step function, and no single solve can locate every step.
What a solve can pin down is a handful of smooth quantities. So the estimate
works with a smooth curve, and separately accounts for the distance between
that curve and the steps it approximates.

**C.12 — §9, self-scaling.** *"the recommendation cannot lose more than the
problem was going to cost anyway."*

The guarantee is a factor, not a fixed number of solves. On a cheap problem
you lose at most a little; on an expensive one you may lose more in absolute
terms but never more than roughly what you were already paying. A bound that
scales with the problem is worth more than a constant that is tight on some
instances and vacuous on others.

**C.13 — §10, the duality.** *"the cheapest solve is the one that says least
about why it is cheapest."*

A solve that grinds through iterations shows you the lower end; a solve that
pays refinement passes shows you the upper end; a solve sitting at the
optimum does neither, so it reveals almost nothing. The better your
parameter, the worse your information about how good it was. This is why the
answer is a band and not a number, and why the band is narrow in the
direction you happened to have moved.

**C.14 — §10, composition by intersection.** *"each contributes its
certified side, nothing is averaged, and an empty intersection is not an
error to suppress but a fired consistency check."*

When two sources bound the same quantity, take the tighter bound from each
side rather than splitting the difference — averaging a certified fact with
an uncertain one throws away the certification. If the sides cross, some
stated hypothesis is false, and that is a finding to report, not a number to
fix up.

**C.15 — §11, why a closed route is not zero.** *"a closure the route must
*report*, since the rising branch can fail to resolve within the anchored
range, and a closed route is not a route that reads zero."*

If a measurement cannot be taken, the estimate must say so rather than
supplying a default. A route that silently returns nothing looks the same as
a route that measured nothing, and only one of those is information.

**C.16 — §12, the asymmetry.** *"Leaving the habitable range downward costs
iterations, which the next solve reports immediately. Leaving it upward
fires refinement against a level the arithmetic cannot reach."*

Guessing the parameter too low is a bounded mistake: the solve is slower and
you find out at once. Guessing too high is not: the solve fails to reach the
tolerance at all. Since the two failures cost differently, the estimate
should be pushed toward the cheap one, and that is why ties break downward
and why the recommendation is clipped from above.

**C.17 — §13, extended-real domains.** *"total ignorance is the proper pair
$(1, +\infty)$, not an improper function."*

Cost is a count of solves, so it is at least one before anything is
measured. Knowing nothing is therefore a perfectly well-formed answer —
between one and infinity — and every scrap of information narrows it from
one side or the other. The estimate never has to invent a value to have
something to say.

**C.18 — §13, the promise stops at the probe.** *"the pessimistic wall
endpoint *is* the anchor and the promise does not extend above the parameter
the caller already ran."*

A solve that converged proves the upper limit lies somewhere above the
parameter it used, and nothing more. So the guaranteed-cost curve runs out
exactly where the caller already is, and the estimate can price moving down
but not moving up. Only a solve that made refinement passes, or one that
stalled, says anything about how much further up there is to go.

**C.19 — §15, the shift degrades the oracle.** *"the system solved is
**always the original one**. $\rho$ degrades the oracle; its perturbation is
left to the refinement loop."*

The shift is a crutch to make the factorization succeed. It does not change
the problem being solved — it makes the solver's inner tool slightly wrong,
and refinement's job is to clean that up. So the convergence tests must
measure the original problem's residual. A test that subtracted the shift's
contribution would be blind to exactly the error the loop exists to remove,
and would report success while converging to the wrong answer.

**C.20 — §15, the pole and the arm.** *"The pole is where it was; the arm
approaching it is not."*

Without the shift, refinement's effectiveness degrades steadily as the
parameter rises and eventually stops working at a definite point. With the
shift, that point is unchanged — but the degradation *saturates* on the way
there, so the number of passes needed grows explosively well before it. The
limit is in the same place; what has changed is how expensive it gets on the
approach, and that is what closes the usable range.

**C.21 — §15, one reading, two unknowns.** *"$\ell_p$ and $b_\infty$ enter
the record only through their sum."*

Under the shift, the residual you measure at the start is two things added
together: a part that shrinks as the parameter grows and a part that does
not. The cost model needs the first on its own. One measured number cannot
be split into two, at any parameter, ever. So either you accept a range —
using a worst-case bound on the stuck part — or you spend one extra solve to
measure it. There is no third option, and that is why the one purchase this
note admits is admitted.

---
