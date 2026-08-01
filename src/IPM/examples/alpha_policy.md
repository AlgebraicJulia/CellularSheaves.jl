# The augmentation policy

Companion to *The augmentation window: a mathematical justification*, which
derives everything asserted here. This note states what is measured after
each solve, what is recorded, and how placement is decided.

Status: candidate. Nothing here is decided.

Supersedes the earlier controller note, which targeted the min-cost+1
window and carried calibrated constants. Both are gone: the target is now
the zero-work window, and every quantity below is either recorded by the
solve or computed from it without fitting.

## 1. Shape of the policy

Measurement is per role and mechanical. Combination is policy and
revisable. The two are kept apart:

    after each solve:
        for each role r:  measure the role's window, record it
    before each solve:
        policy reads the recorded history, chooses alpha

Nothing is aggregated at measurement time. Three windows go into the
history; the policy decides what to do with them. This preserves
attribution (which role binds which edge), keeps provisional readings
local to the role that produced them, and lets per-role drift be seen —
the roles are warm-started differently and do not move together.

## 2. Per-role measurement

Each role's solve records $\alpha$, $\varepsilon$, its initial constraint
residual $\|r_0\|$, its dual residual $\|s_0\|$, its CRAIG iteration count,
and its refinement pass count. From these:

    alpha_0 = alpha * ||r0|| / eps        # zero-work floor
    alpha_c = alpha * eps  / ||s0||       # ceiling

Both are free: each residual is already computed by the test that consumes
it. Both constants are 1 — the two convergence tests read literally, not
calibrated quantities.

**If $\alpha_0 \le \alpha_c$** the role's zero-work window is nonempty.
Record the interval and stop; there is nothing further to compute.

**If $\alpha_0 > \alpha_c$** the window is empty and the solve necessarily
did work. Which kind identifies where the step sat:

- *CRAIG iterations ran* — the step was below the ceiling. The Krylov
  recursion has built a spectral representation; harvest it (§3) to obtain
  the graded levels $\alpha_1, \alpha_2, \dots$ and record the deepest
  interval that is nonempty.
- *Refinement passes ran instead* — the step was above the ceiling, and no
  Krylov space was built. Record $\alpha_c$ (which is sharp here, since
  $\|s_0\| > \varepsilon$ puts the sensor well onto its ramp) and mark the
  floor **unknown**. The state clears in one step: descend below
  $\alpha_c$, and the next solve will iterate and produce a harvest.

## 3. The graded levels

Admitting $k$ Krylov iterations relaxes the floor only; the ceiling is
common to all levels. The floors are

$$\alpha_k = \big(\alpha_0 \cdot c_k\big)^{1/(k+1)},$$

where $c_k$ is the norm of the $k$-th orthogonal polynomial of the limiting
spectral measure — a property of the operator pair and the gap direction,
independent of $\alpha$ and of $\varepsilon$.

**Obtaining $c_k$.** The Golub–Kahan recursion computes the Jacobi
coefficients $b_j$ from which $c_k = \big(\prod_{j \le k} b_j\big)^{1/2}$,
as a running product. A solve that ran $k$ iterations therefore determines
$c_1, \dots, c_k$ at no cost beyond recording them. Three points:

- **Log the coefficients, not the Ritz pairs.** Nodes and weights are
  downstream of the coefficients (the recursion builds the bidiagonal, then
  diagonalizes); recovering $c_k$ from them is a round trip that loses a
  little and costs more storage. The coefficients also arrive
  incrementally, so the level sequence can be truncated as soon as
  $\alpha_k \le \alpha_c$.
- **Depth is the limit, not method.** At $k = 1$ the recursion, a direct
  two-moment computation, and a harvest all agree — one step has no
  orthogonality to lose. Deeper coefficients degrade as orthogonality is
  lost, worst when the spectrum is clustered.
- **The levels you need are the shallow ones you did not run.** A step that
  ran $k$ iterations already witnesses that $k$ suffices; what is wanted is
  whether $k-1$ or $k-2$ would have. Those come from the later coefficients
  in the sequence, which are the less reliable ones. When $k$ is large,
  a fresh direct measurement of $c_1$ (and possibly $c_2$) may be preferable
  to reading a drifting $c_{k-1}$.

**When no level works.** If $\alpha_c \le c_1$ the equation
$\alpha_k \le \alpha_c$ has no positive solution. Since $c_1$ approximates
the accumulation point of the level sequence, this says the ceiling has
descended below the point the levels converge towards: no small iteration
budget opens a window. Record the condition; it is a signal to change the
tolerance, the preconditioner, or the problem — not $\alpha$.

## 4. Measuring $c_1$ directly

When a role's window is nonempty, no iterations run and no harvest exists —
the window's defining property is what makes it spectrally silent. If the
levels below are wanted anyway, one application of
$S = BF^{-1}B^\top$ to $r_0$ — one triangular-solve pair against the
factorization in hand, two multiplications by $B$ — gives

    m1 = (r0 . S r0) / (r0 . r0)
    m2 = (S r0 . S r0) / (r0 . r0)
    R1 = sqrt(1 - m1^2 / m2)
    c1 = alpha * R1

exactly. The same computation returns $\theta = m_1$, the elasticity, which
reports whether the self-similarity assumption — and therefore the floor
formula itself — still holds at the current $\alpha$.

Note the cancellation: inside the window the spectrum of $S$ clusters near
1, so $m_1$ and $m_2$ are both near 1 and $R_1^2$ is a difference of nearby
quantities. Compute the residual vector directly rather than subtracting
the scalars.

## 5. What the history records

Per solve, per role:

    interval        [alpha_0, alpha_c], or [alpha_k, alpha_c], or
                    [unknown, alpha_c]
    level           k, the iteration budget the interval assumes
    binding         which edge was tight, if placement was clamped
    theta           if measured
    exhausted       if alpha_c <= c1
    counts          CRAIG iterations and refinement passes actually incurred

## 6. Placement

The policy reads the history and chooses $\alpha$ for the next solve. The
default: intersect the recorded intervals, place at the geometric midpoint,
clamp a margin below the lowest ceiling.

The margin is one-sided by design. Leaving the window downward costs Krylov
iterations, which the counts report immediately. Leaving it upward fires
refinement, which faces the same attainability limit as the base solve and
may not converge at all. The asymmetry is qualitative, and the clamp is
where it is accounted for.

Intersection is a default, not a requirement — with the three windows in
hand the policy can do otherwise, for instance accepting an iteration in a
cheap role to buy conditioning margin for an expensive one. That option
does not exist once the windows are collapsed into one, which is the
argument for keeping them separate.

Degenerate readings: if a role's gauges are unavailable, do not estimate —
reduce $\alpha$ by a fixed factor and re-read. If a role's floor is
unknown, the intersection is a lower bound on the true floor; placement is
provisional and the next step's counts will correct it.

## 7. Measured behaviour

Replay across the fine-grid sweep (43 problem/solver pairs), comparing
three policies at matched starting points: this policy; the earlier
controller targeting min-cost+1 with calibrated constants; and a fixed
$\alpha$ held at its initial value.

    total CRAIG iterations      new 1198    old 1224    fixed 1258
    steps placed outside the
      operative window          new   15    old   33

The margin over the earlier controller is small on problems whose windows
stay wide — any reasonable policy sits inside them — and concentrates where
windows narrow or close. On the collapse cases the difference is
structural: the earlier controller aims at a window defined relative to a
cost minimum it cannot observe, and drifts high as that window's upper edge
outruns the zero-work one.

Endpoint accuracy on the same sweep, against directly observed transitions:

    alpha_0     median error -0.04 dec,  MAE 0.07,  96-98% within one lattice step
    alpha_c     median error +0.05 dec,  MAE 0.46 (ipm) / 0.26 (hsd)
    alpha_k     medians within 0.06 dec through k = 3; spread grows with k

The floor is measurement-limited — the observing lattice is 0.1 decades.
The ceiling's median is equally good but its spread is larger; this is not
the sensor's peg (the error is flat in probe distance) and not the
combination rule (`min` over roles beats the alternatives decisively). It
is the gap between the continuous dual crossing and the discrete event that
the observation records.

## 8. Open items

- The ceiling's spread. Assumption C's own claim measures well — the ramp's
  slope is 0.996 and $\|s_0\|/\alpha$ varies by 0.13 decades within a
  solve — so the residual scatter lies between the crossing and the
  observed cliff, not in the sensor.
- What happens above the ceiling has never been measured: whether
  refinement passes fire and converge cheaply, stall, or diverge. The
  policy's asymmetry argument assumes the second or third. A sweep through
  and past the ceiling, recording passes and their residuals, would settle
  it — and would also test whether the ceiling should be defined by "no
  pass fires" or by "passes converge".
- Per-pass residuals are not recorded, only the first. The ratio
  $\|s_1\|/\|s_0\|$ would measure directly whether refinement shares the
  base solve's attainability floor.
- The margin below the ceiling has not been re-derived under the
  constant-free endpoints.
- Whether the per-role `max`/`min` remains correct when roles carry
  different tolerances has not been examined.
