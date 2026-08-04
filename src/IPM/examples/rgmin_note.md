# The regularization floor: `rgmin` as a scaled quantity

Status: candidate. Nothing here is decided.

This note concerns one setting — the lower bound of the solver's $\rho$-shift
ladder — and recommends replacing its fixed value with a problem-scaled one.
The recommendation is independent of the augmentation controller; it is
argued from the size of the object the shift is correcting.

## 1. What the shift is fixing

When the Cholesky factorization of $F = \tfrac1\alpha A + B^\top B$ fails,
the solver retries with $F + \rho I$, starting at $\rho = \texttt{rgmin}$ and
multiplying by 8 until it succeeds or exceeds `rgmax`. Present values:
$\texttt{rgmin} = 10^{-9}$, $\texttt{rgmax} = 10^{-6}$.

Measured at the point where the ladder engages, on captured solve states:

    capture                lambda_min(F)     lambda_max(F)   lambda_min/lambda_max
    07 ipm  alpha 1e10.9    +2.2e-15           3.6e+01            6.2e-17
    07 ipm  alpha 1e12.0    -2.6e-15           3.6e+01           -7.2e-17
    07 hsd  alpha 1e10.0    +3.5e-15           3.6e+01            9.6e-17
    e15 ipm alpha 1e11.9    -1.4e-15           5.8e+00           -2.4e-16

The matrix is not meaningfully indefinite. Its smallest eigenvalue sits at
the unit-roundoff scale relative to its largest — in two of the four cases
still positive, with the factorization failing only because accumulated
rounding tipped a computed pivot to the wrong side of zero. Assembling $F$
two different ways changes $\lambda_{\min}$ in the last digit and agrees to
$10^{-25}$ absolute, which is the signature of a decision being made at the
roundoff boundary rather than on the matrix's actual definiteness.

The job the shift is being asked to do is therefore of size $u\|F\|$, and
$\|F\| \approx \sigma_{\max}^2(B)$ once $B^\top B$ dominates. Define

$$\rho_\star = u\,\sigma_{\max}^2(B).$$

## 2. How far the present value overshoots

$\sigma_{\max}^2(B)$ is a property of the constraint matrix alone and does
not change over a solve, so $\rho_\star$ is a per-problem constant computable
once, at construction, alongside the norms the augmentation default already
uses. Across the fine-grid fleet:

    median rho_star                    5.2e-16
    range                              3.1e-22  ..  5.8e-14
    median rgmin / rho_star            1.3e+06

So the ladder's first rung is a median of six orders of magnitude above the
size of the perturbation required, and the overshoot varies by eight orders
across problems — which is what a fixed absolute constant does to a quantity
that scales with the problem.

Direct check, factoring the captured operators at a range of shifts:

    capture                smallest shift for which Cholesky succeeds
    07 ipm  alpha 1e10.9     none required
    07 ipm  alpha 1e12.0     1e-14        (rho_star = 8.0e-15)
    07 hsd  alpha 1e10.0     none required
    e15 ipm alpha 1e11.9     1e-16        (rho_star = 1.3e-15)

In every case a shift at or below $\rho_\star$ sufficed. (Caveat: this is a
dense factorization of the captured matrices; the solver's chordal
factorization performs the same operations in a different grouping, and the
two disagree exactly at the boundary — which is the point of §1, not an
objection to it.)

## 3. Why the overshoot is not free

The shift is not a neutral repair. With $F_\rho = \tfrac1\alpha A + B^\top B
+ \rho I$, the assembled solution satisfies

$$f - Ax + B^\top y = \alpha\rho\, x$$

*exactly, at every iteration count* — verified to seven digits on every
capture. So the dual residual, which is pure rounding error when $\rho = 0$,
acquires a deterministic term proportional to $\rho$ that no Krylov work and
no refinement pass removes. Correspondingly the initial constraint residual
acquires a constant bias $\rho\,b_\infty$ that does not scale with $\alpha$.

Both effects are proportional to $\rho$, so the overshoot is transmitted in
full. Measured on the captures at the alphas where the ladder engaged:

    rho          bias as a fraction of ||r0||
    rho_star           0.1% – 1.8%
    1e-12              5% – 87%
    1e-9              100% – 103%

At $\rho = 10^{-9}$ the recorded primal residual is *entirely* regularization
bias, carrying no information about the quantity it is supposed to measure.
At $\rho_\star$ the bias is negligible.

The dual side shows the same thing as a discontinuity. Because
$\alpha\rho\|x\|$ at $\rho = 10^{-9}$ is roughly $10^6$ times the residual's
natural level, the moment the ladder engages the dual residual jumps by that
factor — measured jumps of 2 to 7 decades in a single grid step — rather
than continuing its smooth climb. At $\rho_\star$ the added term is
comparable to what is already there and the climb is uninterrupted.

## 4. Recommendation

    rgmin = u * sigma2max(B)        # computed once, at construction

with the ladder otherwise unchanged. Two adjustments are worth considering
alongside it.

**Increase the multiplier.** The ×8 step suits a ladder spanning
$10^{-9}$ to $10^{-6}$ — three decades, three rungs. Starting at $\rho_\star$
spans eight or nine decades, and ×8 would need six additional
factorizations in the worst case. A multiplier of $10^3$ covers the same
range in three rungs, and overshooting the minimal shift by $10^3$ still
leaves it $10^3$ below the present starting point.

**Reconsider `rgmax`.** A factorization that genuinely requires $10^{-6}$ is
not being repaired; it is being replaced by a different problem whose primal
and dual residuals are dominated by the substitution. Whether to continue
there or to fail is a design judgement this note does not settle, but the
current ladder makes the two indistinguishable, since it starts high enough
that every engagement looks alike.

## 5. Expected effect, and what is untested

**Expected.** On the fine-grid fleet the shift sets the window's upper edge
on 13% of iterations. On those, the ceiling becomes readable from below
rather than being set by an invisible discontinuity, and the primal
residual's high-$\alpha$ saturation moves out of the operating range. No
change on the other 87%.

**Not expected.** Any iteration-count improvement. The affected windows are
already a median 7.8 decades wide, and cost is flat across a window, so
widening them buys nothing directly. The value is that two estimates become
correct where they were previously wrong by 1 to 5 decades.

**Untested.** Whether $\rho_\star$ suffices for the solver's own
factorization, on problems beyond the two captured; whether the ladder
climbs more often at the lower starting rung, and how much that costs; and
whether the wider windows that result are as safe to sit in as the present
narrower ones. A sweep at the new setting answers all three.
