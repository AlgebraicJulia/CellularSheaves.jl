# The augmentation KKT solver

*A formal description of the linear-system solver at the heart of the interior-point
method: the augmentation backend, its inner Krylov correction, the outer refinement
loop, the inexactness knob θ, and the tolerance scaling. No code — enough to
reimplement. The ambient IPM is sketched first so the interface makes sense.*

---

## 1. Ambient problem and interior-point loop

The solver targets the conic quadratic program

$$\min_p \ \tfrac12 p^\top Q p + c^\top p \quad\text{s.t.}\quad Bp = g,\ \ p \in K,$$

with $Q \succeq 0$ symmetric, $B \in \mathbb{R}^{m\times n}$ of full row rank, and $K$ a
product of convex cones (nonnegative, second-order, semidefinite, exponential, power).
The dual carries a free multiplier $y \in \mathbb{R}^m$ and a conic slack $d \in K^\*$.

**Optimality (KKT).** A triple $(p, y, d)$ is optimal iff

- primal feasibility: $Bp = g$, $p \in K$;
- dual feasibility: $r_d := d - c - Qp + B^\top y = 0$, $d \in K^\*$;
- complementarity: $p^\top d = 0$.

The central-path parameter is $\mu := p^\top d / \nu$, with $\nu$ the total cone degree.

**One IPM iteration.** Linearize about the current $(p, y, d)$. The cone is replaced by
its Nesterov–Todd scaling — a symmetric positive-definite matrix $H$ (barrier Hessian at
the scaling point; a Tunçel matrix for nonsymmetric cones). With $A := H + Q$ (SPD), the
Newton direction solves a **saddle system**

$$\begin{bmatrix} A & -B^\top \\ B & 0 \end{bmatrix}\begin{bmatrix}\Delta p\\ \Delta y\end{bmatrix} = \begin{bmatrix} f \\ g \end{bmatrix}. \tag{KKT}$$

Mehrotra's scheme solves (KKT) **twice per iteration through the same $A,B$** — an affine
*predictor* and a centered *corrector*, differing only in the right-hand side $(f,g)$ —
then takes a fraction-to-the-boundary step $p \leftarrow p + \tau\,\Delta p$, etc.

**Stopping.** Convergence is tested on residuals scaled by the data norms
$n_g := \lVert g\rVert$ and $n_c := \lVert c\rVert$:

$$\mathrm{pres} := \frac{\lVert r_p\rVert}{1+n_g} < \texttt{feas\_tol}, \qquad
\mathrm{dres} := \frac{\lVert r_d\rVert}{1+n_c} < \texttt{feas\_tol},$$

together with a gap condition — either $\mu < \texttt{gap\_tol}$, or the duality gap
$\mathrm{pobj}-\mathrm{dobj} < \texttt{gap\_tol}\,(1+|\mathrm{pobj}|+|\mathrm{dobj}|)$, where
$\mathrm{pobj} = \tfrac12 p^\top Qp - c^\top p$ and $\mathrm{dobj} = g^\top y - \tfrac12 p^\top Qp$.

Everything below is the black box that solves one instance of (KKT) for a given $(f,g)$ to a requested accuracy.

---

## 2. Augmentation

Fix a penalty $\alpha > 0$, set $\beta := 1/\alpha$, and form the **augmented (1,1)
operator**

$$F := \beta A + B^\top B \quad(+\,\rho I \text{ if a diagonal shift is needed to factor}).$$

$F$ is SPD; **factor it once** (sparse chordal $LDL^\top$). The entire cost of the solve
is counted in **applications of $F^{-1}$**; matrix–vector products and vector arithmetic
are free.

The idea: eliminating the constraint by an augmented Lagrangian replaces the indefinite
saddle system with repeated SPD solves against $F$.

---

## 3. Base solve — backsolve plus dual correction

Given a multiplier estimate $y_0$ (a warm start, or $0$):

**Backsolve.**
$$x \leftarrow F^{-1}\!\big(\beta f + B^\top(g + \beta y_0)\big).$$
This $x$ is exact only if $y_0$ is; in general it leaves a primal residual
$r := g - Bx \neq 0$.

**Dual correction.** Remove $r$ by the minimum-norm correction

$$\min_{\delta x}\ \lVert \delta x\rVert_{FF^\top} \quad\text{s.t.}\quad B\,\delta x = r,$$

which is exactly CG on the Schur complement $S := B F^{-1} B^\top$: solve
$S\,\delta y = r$, then $\delta x = F^{-1}B^\top \delta y$. This is run as an inner Krylov
iteration (CRAIG / Golub–Kahan bidiagonalization of $BF^{-\top}$); **each inner iteration
costs one application of $F^{-1}$**. Accumulate $x \leftarrow x + \delta x$.

**Recovery.**
$$y \leftarrow y_0 + \alpha\,(\delta y + r).$$
$\delta y$ is the dual correction; the leftover primal residual $r$ folds into the same
recovery, so no extra matvec is spent to form $y$.

The inner solve stops when the directly evaluated constraint residual
$\lVert g - Bx\rVert_2 \le \texttt{atol}$ (see §6). It writes that residual back in place,
so the caller reads off $y$ without recomputing $Bx$.

```
BASESOLVE(f, g, α, y0, atol):
    β ← 1/α
    x ← F⁻¹(β·f + Bᵀ(g + β·y0))          # 1 application
    r ← g − B·x
    (δx, δy) ← CRAIG(B, F, rhs=r, atol)   # k applications; solves B·δx = r
    x ← x + δx
    y ← y0 + α·(δy + r)
    return x, y, r        # r is the exit residual g − B·x
```

---

## 4. Iterative refinement — the outer loop

The base solve is inexact (finite $\alpha$; inexact inner solve). Refine the **full 2-row
residual** of (KKT). With the accumulated pair $(x, y)$:

```
REFINE(x, y, f, g):
    repeat:
        s_d ← f − A·x + Bᵀ·y                  # dual residual
        s_p ← g − B·x                          # primal residual
        if max( ‖s_d‖/(1+n_c), ‖s_p‖/(1+n_g) ) ≤ tol:  return CONVERGED
        if that ratio rose since the last pass:         return STALLED
        (dx, dy) ← BASESOLVE(f=s_d, g=s_p, α, y0=0, atol_interior)   # §6 for atol
        x ← x + dx ;  y ← y + dy
```

Two invariants make this correct:

- **Residual system from a zero multiplier.** The correction is a base solve on the
  residual pair $(s_d, s_p)$ with $y_0 = 0$ — *not* a warm-started re-solve of the
  original $(f,g)$. A warm re-solve injects an error proportional to $\alpha$ and raises
  the achievable accuracy floor by orders of magnitude.
- **Each pass's leftover is proportional to its own increment**, so the raw residual is
  the convergence object in every regime; the loop needs no special guard.

---

## 5. The nesting, and the cost

Two loops, nested:

- **Outer** — iterative refinement over the true 2-row KKT residual (§4). Each pass is one
  base solve.
- **Inner** — within each base solve, the CRAIG iterations of the dual correction (§3).
  Each iteration is one application of $F^{-1}$.

For a solve that runs $j$ refinement passes with $n_i$ inner iterations in pass $i$
(and $n_0$ in the base), the total cost is

$$\text{cost} \;=\; \underbrace{1}_{\text{base backsolve}} + \; n_0 \;+\; \sum_{i=1}^{j} \big(1 + n_i\big) \quad\text{applications of } F^{-1}.$$

---

## 6. Inexactness (θ) and tolerance scaling

**Target.** The outer loop drives to

$$\texttt{tol} := \max(\texttt{force\_tol},\ \texttt{floor\_tol}),$$
$$\texttt{force\_tol} := \min\!\Big(\texttt{forcing\_frac}\cdot \tfrac{\mu}{\mu_1},\ \texttt{forcing\_ceil}\Big), \qquad \texttt{floor\_tol} := 100\,u\,(1 + \max(\mathrm{pres},\mathrm{dres})),$$

where $\mu_1$ is the initial centrality and $u$ the unit roundoff. `force_tol` is the
inexact-Newton **forcing sequence** — loose early, shrinking with $\mu$ so the linear
solve is never solved tighter than the current Newton step warrants; `floor_tol` is the
attainable accuracy floor.

**Scaling the inner tolerance.** The outer test uses residuals scaled by the data norms,
but CRAIG measures the *unscaled* 2-norm $\lVert g - Bx\rVert_2$ — and its exit residual
*is* the next pass's $s_p$. So the conversion is done once, at the call site:

$$\texttt{atol} := \texttt{tol}\cdot(1 + n_g).$$

Only $n_g$ enters: CRAIG owns the primal row; the dual row is CRAIG-independent. This
makes the inner solve stop exactly at the outer loop's scaled target — no over-solving to
a bare 2-norm.

**θ-inexactness.** An interior refinement pass need not solve its correction all the way
to `atol` — only enough for the outer loop to keep contracting. So an interior pass caps
its inner tolerance at

$$\texttt{atol} \leftarrow \max\big(\texttt{atol},\ \theta \cdot \lVert r\rVert_{\text{entry}}\big), \qquad \theta = 0.1,$$

i.e. cut the pass's entry residual by roughly $10\times$ and let the outer loop finish the
job. The **base solve and the final exit pass run exact** ($\theta = 0$). This is inexact
Uzawa / inexact Newton: intermediate corrections are cheap; the outer loop supplies the
accuracy.

---

## 7. Exit status

Every solve returns a typed exit and its cost record (base and per-pass inner counts, pass
count):

- **REACHED_FORCE** — converged to `force_tol`;
- **REACHED_FLOOR** — converged only to `floor_tol` (the accuracy floor bounded it);
- **STALLED** — the scaled residual ratio rose above a stall threshold between passes
  (the penalty $\alpha$ is too aggressive for this system);
- **ITMAX** — pass budget exhausted.

---

## 8. Homogeneous (HSD) extension

The homogeneous self-dual embedding appends a scalar $(\tau,\kappa)$ row, giving a 3-row
system

$$\begin{bmatrix} A & -B^\top & -c \\ B & 0 & -g \\ c^\top & g^\top & \delta \end{bmatrix}\begin{bmatrix} x \\ y \\ \sigma \end{bmatrix} = \begin{bmatrix} f \\ h \\ e \end{bmatrix}.$$

It is solved by **bordering**: the 2-row augmentation of §2–6 is the leading block, and the
extra column/row is absorbed by a Woodbury (capacitance) update — one additional backsolve
to price the border — reusing the *same* factorization $F$ and the *same* refinement loop.
Now three solves share the factorization per iteration (predictor, corrector, and the
Woodbury/border solve). The primal, dual, and $\tau$ residuals are scaled by $1+n_g$,
$1+n_c$, and $1+n_g+n_c$ respectively.
