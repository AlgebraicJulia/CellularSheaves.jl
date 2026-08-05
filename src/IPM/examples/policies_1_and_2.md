# Choosing the augmentation parameter: policies 1 and 2

*Self-contained. Describes the solve, exactly what must be recorded, and two
policies for choosing the parameter of the next solve. Policy 1 needs two
residual norms. Policy 2 additionally needs the correction iteration's
residual history.*

---

## 1. The solve

Each system is

$$\begin{bmatrix} A & -B^\top \\ B & 0\end{bmatrix}
\begin{bmatrix} x \\ y\end{bmatrix} =
\begin{bmatrix} f \\ g\end{bmatrix}$$

with $A$ symmetric positive definite, $B$ full row rank $p\times n$. For a
chosen $\alpha>0$, write $\beta = 1/\alpha$ and

$$F \;=\; \beta A + B^\top B .$$

$F$ is factored once. **Cost is the number of applications of $F^{-1}$**;
matrix–vector products and vector work are free.

```
SOLVE(f, g, alpha, eps, tau):

  --- base invocation ---
A1  x  = F^-1 ( beta*f + B'*g )                        <- 1 application
A2  r  = g - B*x
A3  if ||r|| > eps:                            (T1) entry test
       CG on S = B F^-1 B'  for  B*dx = r,
       stopping when ||g - B*(x+dx)|| <= eps   <- 1 application per iteration
       x = x + dx
A4  r  = g - B*x
    y  = alpha*(dy + r)                        recovery
A5  s  = || f - A*x + B'*y ||

  --- refinement ---
A6  while s > eps:                             (T2) dual test
       w  = f - A*x + B'*y ;  p = g - B*x      residuals from current state
       dx = F^-1 ( beta*w + B'*p )             <- 1 application
       gap = || p - B*dx ||
       if gap > eps_j:  CG until met           eps_j = max(eps, gap/10)
       dy = alpha*( dy_c + (p - B*dx) )
       x += dx ;  y += dy
       s_new = || f - A*x + B'*y ||
       if s_new > tau * s        -> HALT, STALLED      (T3)
       s = s_new
    on first pass with s <= eps, run one final pass at full eps and
    verify BOTH residuals directly
```

Fixed choices assumed throughout: cold start $y_0=0$; no shift; absolute
tolerance; $\tau = 1/2$; refinement solves the **residual** system from a zero
multiplier (not a warm-started re-solve of the original data — that injects
error proportional to $\alpha$ and raises the achievable accuracy floor by
orders of magnitude).

---

## 2. What to record

Everything below is already computed by the solve. Nothing extra is
evaluated and no additional application is spent.

### Required for policy 1

| symbol | what | where |
|---|---|---|
| `alpha` | the parameter used | input |
| `eps` | the tolerance | input |
| `r0_norm` | $\|g - Bx\|$ **before any correction** | line A2 |
| `s0` | $\|f - Ax + B^\top y\|$ **after the base recovery**, before any refinement | line A5 |
| `n0` | correction iterations in the base invocation | A3 |
| `j` | refinement passes | A6 |
| `exit` | converged / stalled / cap | — |

Derived, both free:

```
ell_p = alpha * r0_norm          the primal scale
ell_d = s0 / alpha               the dual scale
```

**Note `s0` must be the dual residual after the base invocation only.** If
refinement runs, later values are contracted and must not be used here.

### Additionally required for policy 2

| symbol | what | where |
|---|---|---|
| `res_hist[0..n0]` | the norm $\|g - B(x+\delta x)\|$ at **every** correction iteration, including iteration 0 (which is `r0_norm`) | A3 |

This is the directly evaluated constraint residual the stopping test already
forms at every iteration. Store the sequence instead of discarding it.

**Do not use the CG coefficients** $(\alpha_k,\beta_k)$ for this. They give the
same quantities in exact arithmetic but drift badly at depth in floating
point: on a 42-iteration solve the coefficient route failed the
self-consistency check of §3.2 by 1.28 nats, while the residual history
passed it to 0.006.

If `n0 = 0` the history is empty and policy 2 falls back to policy 1's
behaviour for that solve (see §5).

---

## 3. The two boundaries

Both are computed from the record at zero cost, and both are *sharp* — they
are the parameters at which a test changes state, not fitted quantities.

### 3.1 The dual boundary

$$\alpha_c \;=\; \frac{\varepsilon}{\ell_d}$$

Above $\alpha_c$ the dual test (T2) fails and refinement passes are paid;
at or below it, no passes.

*Basis.* The base dual residual is entirely rounding-born — in exact
arithmetic the recovery step makes it vanish identically — and it grows in
proportion to $\alpha$. So $s_0 = \ell_d\,\alpha$ with $\ell_d$ a property of
the operator and arithmetic, and (T2) holds iff $\alpha \le \varepsilon/\ell_d$.

*Validity.* The proportionality requires $\alpha$ above the switch
$\alpha_\star = \lambda_{\max}(A)/\|B\|^2$; below it $s_0$ is flat in
$\alpha$ and $\alpha_c$ is not meaningful. **Check $\alpha \ge \alpha_\star$
before using it.**

### 3.2 The entry boundary

$$\alpha_0 \;=\; \frac{\ell_p}{\varepsilon}$$

Below $\alpha_0$ the entry test (T1) fails and correction iterations are
paid; at or above it, none.

*Basis.* $\|r_0\| = \ell_p/\alpha$ exactly, from the identity
$r_0 = -\alpha^{-1}S(y_0-\lambda^*)$. Verified to three significant figures
over 2.5 decades.

*Safe downward.* The residual's elasticity in $\log\alpha$ lies in $(0,1]$
with no assumption, so below the anchor the true residual is *smaller* than
$\ell_p/\alpha$ predicts — the entry test keeps passing further down than
claimed. So $\alpha_0$ is a conservative estimate of where iterations begin.

*Measured sharpness.* At $0.9\,\alpha_0$ the machine needed 1 iteration; at
$1.2\,\alpha_0$ it needed 0. Two instances, exact to grid resolution.

### 3.3 The window

$$W_0 = [\alpha_0,\ \alpha_c], \qquad
\text{nonempty} \iff \varepsilon^2 \ge \ell_p\ell_d .$$

Inside $W_0$ the solve costs exactly **1 application**: one backsolve, entry
test passed, no refinement. Its width in logs is

$$\text{width} \;=\; \log\alpha_c - \log\alpha_0
\;=\; 2\log\varepsilon - \log(\ell_p\ell_d),$$

which **narrows by two nats for every nat the tolerance tightens**, one from
each end. Width is measurable every solve with no transport and is the
natural early-warning signal: it shrinks as the sequence proceeds and
reaches zero when the window closes.

---

## 4. Policy 1 — window, then hold

The baseline. Uses `ell_p` and `ell_d` only.

```
POLICY_1(record):
    ell_p = alpha * r0_norm
    ell_d = s0 / alpha
    a0    = ell_p / eps
    ac    = eps  / ell_d

    if a0 <= ac:                       # window open
        return sqrt(a0 * ac)           # geometric mean
    else:                              # window closed
        return alpha                   # hold
```

**Rationale for the geometric mean.** Inside $W_0$ every parameter costs 1,
so there is nothing to optimise; the only question is robustness. The two
boundaries close toward each other at one nat per nat, so the midpoint in
log coordinates is the last point standing and is the maximin choice.

**Behaviour.** Correct and cheap while the window is open. Once it closes,
policy 1 stops adapting; the cost then rises on its own as the operator
drifts. This is the thing to beat.

---

## 5. Policy 2 — stay below the dual boundary

Refinement is the expensive and dangerous side: a pass costs one application
plus any interior iterations, and far enough above $\alpha_c$ the solve
**stalls** and the whole solve is lost. Correction iterations cost one each
and the overshoot is reported immediately by the next solve. Policy 2
therefore keeps the parameter on the correction side and treats $\alpha_c$
as a ceiling, not a target.

### 5.1 Levels

When the correction iteration runs, the recorded residual history gives, for
$k = 1,\dots,n_0$,

```
c_k     = alpha^k * res_hist[k] / res_hist[0]
alpha_k = ( ell_p * c_k / eps ) ^ (1/(k+1))
```

$\alpha_k$ is the parameter at which **$k$ correction iterations suffice**.
The levels descend, $\alpha_1 > \alpha_2 > \dots > \alpha_{n_0}$, and
$\alpha_{n_0}$ reproduces the parameter actually solved at — a
self-consistency check every implementation should assert (agreement to
$10^{-2}$ in logs; failure indicates the coefficient route was used instead
of the residual history).

The **level spacing** $\log\alpha_{k} - \log\alpha_{k+1}$ is the distance in
$\log\alpha$ that buys exactly one more iteration. It is the natural unit of
margin: a fixed factor such as a decade buys an unknown number of iterations
because the levels crowd, tightening roughly as $1/k^2$.

Levels are known only for $k \le n_0$. **Nothing below $\alpha_{n_0}$ is
priced**, and it cannot be recovered from this record.

### 5.2 The policy

```
POLICY_2(record, state):
    ell_p = alpha * r0_norm ;  ell_d = s0 / alpha
    a0 = ell_p / eps ;  ac = eps / ell_d

    if a0 <= ac:                        # window open
        return a0                       # sit at the LOWER edge

    # window closed
    if n0 >= 2:
        gap = log(alpha_{n0-1}) - log(alpha_{n0})    # one level spacing
        state.spacing = gap                          # carry it forward
    gap = state.spacing  (default: one third of a decade if never set)

    return exp( log(ac) - gap )         # one level below the ceiling
```

**Why the lower edge while the window is open.** Every point in $W_0$ costs
1, so the choice is purely about which failure to guard against. Sitting at
$\alpha_0$ puts the entire window width between the parameter and $\alpha_c$
— the maximum available margin against an upward jerk — while an
undershoot below $\alpha_0$ costs one iteration and is reported at once.

**Why one level below the ceiling once it closes.** $\alpha_c$ itself has
zero margin: any drift puts the next solve into refinement. One level
spacing is the cheapest nonzero margin that can be stated in units of cost.

**Bootstrapping.** While the window is open the solve pays no correction
iterations, so there is no history and no spacing. This is harmless: the
policy does not need a spacing until the window closes, and once it does,
the parameter sits below $\alpha_0$, so the correction iteration runs and the
history appears. The first closed-window solve uses the default and every
subsequent one uses a measured spacing.

**On overshoot.** If refinement fires anyway (the boundary moved), the next
solve's record contains $\alpha_c$ *witnessed* rather than predicted, and the
policy simply reapplies. No special case is needed. The pass ratios that
overshoot produces are the only route to pricing the refinement side; policy
2 does not use them, but they should be recorded for later work.

---

## 6. What to measure when evaluating

Per solve: applications, `n0`, `j`, exit type, `alpha`, and the derived
`a0`, `ac`, window width.

Per run, for each policy:

- **cumulative applications** over the endgame sequence — the headline
- **stalls** — policy 2 should have none; any stall is a policy failure, not
  a solver failure
- **overshoots** — solves where refinement fired unintentionally; the rate
  indicates whether one level of margin is enough
- **window width against iteration index** — the drift signal, and the thing
  that says when each policy stops being applicable

A useful control is an oracle that solves at the true cheapest parameter,
found by sweeping. The gap between policy 2 and that oracle is what a third
policy would have to recover.

---

## 7. Assumptions and limits

- $\alpha_c$ requires $\alpha \ge \alpha_\star = \lambda_{\max}(A)/\|B\|^2$;
  below the switch the dual residual does not scale with $\alpha$ and the
  boundary is meaningless. Assert this.
- Everything is conditional on $\tau$. The stall fires when the *observed*
  ratio crosses $\tau$, so the point of no return sits where the contraction
  reaches $\tau$, not 1. Measured, moving $\tau$ from $0.5$ to $0.1$ moved
  the failure point by more than half a decade on a fixed problem.
- Nothing below $\alpha_{n_0}$ is priced, ever, from one record.
- The refinement side is not priced by either policy. Pricing it requires a
  contraction ratio, which requires a pass to have fired.
- Both policies assume the operator drifts smoothly between solves. A large
  jerk can put policy 2's parameter above $\alpha_c$; the consequence is
  refinement passes, not a stall, provided the jerk is smaller than the
  distance to the stall point.
