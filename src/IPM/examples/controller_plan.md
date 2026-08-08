# Plan: the α controller, in stages

*Stage 1 is a complete controller reading the predictor only. It also logs
what Stage 2 needs, so Stage 2 is a decision made on data rather than a
guess.*

---

## Background: the two boundaries

Both come from one solve's record, at zero cost.

```
ell_p = alpha * ||g - B*x||     # r0, measured BEFORE any CRAIG
ell_d = s0 / alpha              # s0 AFTER the base recovery, BEFORE refinement

alpha_0 = ell_p / ptol          # below this, CRAIG iterations
alpha_c = ytol  / ell_d         # above this, refinement passes
```

Inside `[α₀, α_c]` the solve costs **1 application**. That is the object the
controller sits in.

**Why this is now the whole problem.** Measured on the tested trajectories,
HSD needs no refinement at the optimal α. So the cost model is `1 + n₀`,
`α_c` is a hard ceiling rather than the start of a trade, and there is no
ceiling law to build.

Three properties that shape the design:

- **`α_c` transports well.** `s₀` is entirely rounding-born and grows in
  proportion to α, so `ℓ_d` is α-invariant — *provided* α is above the switch
  `α_⋆ = λmax(A)/‖B‖²`. Below it `s₀` is flat and `ℓ_d` is meaningless.
- **`α₀` transports badly.** `ℓ_p = ‖S_α(y₀ − λ*)‖` and `S` depends on α.
  Measured: reading `α₀` from 1.5 decades above overestimated it by **14×**;
  from just below, by 1.5×. It is a local quantity.
- **`α₀` is warm-start dependent, `α_c` is not.** `ℓ_p` *is* the multiplier
  gap, so a warm start shrinks it directly — measured 3–6 decades between
  predictor and warm-started corrector on the same operator. `ℓ_d` is
  factorization backward error and does not care.

---

## Stage 1 — predictor-only controller

### The policy

```
if the predictor refined (j > 0):
    alpha <- alpha_c / RETREAT            # RETREAT ~ 3 (half a decade)
else:
    i     <- smallest i <= n0 with alpha_i <= alpha_c
    alpha <- sqrt(alpha_i * alpha_c)      # i = 0 gives the cost-1 window
```

with

```
alpha_k = ( ell_p * c_k / ptol )^(1/(k+1)) ,   log c_k = k log(alpha) + log(||r_k||/||r_0||)
```

from the correction iteration's **residual history** — not the CG
coefficients, which drift badly at depth (measured: 1.28 nats of
self-consistency error at 42 iterations, versus 0.006 for the residual form).

### Why the retreat targets `α_c`, not "current α minus a factor"

`α_c` is measured on the solve that just refined, and it is where the dual
test starts passing. How far you overshot is not known; where you need to be
is. Go there, then apply the factor for margin.

The factor's size is not critical — below `α_c` there is no information
anyway, so it is drift margin, not an optimisation. A level spacing would be
the principled unit, but levels require `n₀ ≥ 2` and a solve that refined has
`n₀ = 0`. So a fixed factor is the honest choice here.

### Guards

- **Assert `α ≥ α_⋆`** before using `α_c`. Below the switch the formula is
  not valid.
- **Rate-limit the climb.** `α₀` is only locally valid, so a target read a
  decade away is wrong. Cap the upward step at ~half a decade and re-measure;
  the downward step needs no cap, since retreating is cheap and self-reporting.
- **`s₀` must be the base-solve dual residual.** If refinement ran, later
  values are contracted and `ℓ_d` comes out wrong.
- **Never step above the last α that converged.** Down is recoverable — you
  pay CRAIG iterations and find out next solve. Up risks the wall.

### Baseline

Compare against the current default (hold at `1e7·α_⋆`) on the same
trajectories, best-α oracle where a sweep is affordable. Report **two
columns, not one**: total applications and lost solves. There is no
defensible exchange rate between them.

---

## Stage 1b — instrument for Stage 2 (same commit)

While Stage 1 runs, log per solve, for **all** roles — predictor, corrector,
and the Woodbury column:

| field | why |
|---|---|
| `ell_p`, `ell_d`, `alpha_0`, `alpha_c` | the per-role windows |
| `n0`, `j`, exit status | cost and which test failed |
| `alpha`, `alpha_star` | position, and the validity assert |
| `ptol`, `ytol` | the boundaries move with these |

**The questions Stage 2 has to answer, and this log answers:**

1. **Does the predictor bind?** Earlier measurement said the shared-α optimum
   equals the split-α optimum — the corrector's window is wide enough
   (warm-started) to contain the predictor's optimum. If the log confirms
   `α₀_corr ≪ α₀_pred` throughout, aggregation over those two is a no-op and
   Stage 2 is only about the column.
2. **What is the column's window now?** Pre-change-2 it separated from the
   Newton windows by ~3 decades over a run — its right-hand side is fixed
   data while theirs shrink — and at e08 iteration 10 the three-way
   intersection was down to a factor of 1.1. **Post-change-2 the column is
   solved to a target derived from `Δτ`, not to `ptol`/`ytol`, so its
   `wfres`/`wdres` may no longer mean the same thing.** Establish what they
   mean before aggregating them.
3. **How often does the retreat fire?** If refinement is rare at the
   controller's α, the retreat branch is nearly dead code and its factor does
   not matter.

---

## Stage 2 — aggregation

Decided by the Stage 1b log, not in advance. The candidate shapes:

- **Predictor only** (i.e. Stage 1 unchanged) if (1) confirms the corrector
  never binds and the column is out of scope post-change-2.
- **`max(ℓ_p)`, `max(ℓ_d)` over roles** — the intersection of the per-role
  windows, which is what `augwindow1` already does. Correct if every role
  must be inside its own window.
- **Predictor + column only**, if the corrector never binds but the column
  does.

The measured fact that makes this cheap: one α must serve all roles, since a
second α means a second factorization, which costs far more than anything
being saved.

---

## Stage 3 — open, not planned

- Anticipation from the secant (how the boundaries moved between solves)
  rather than reaction. Not attempted; the boundaries move ~1–2 decades per
  system, monotonically with a single turning point, so it is plausible.
- The cold start: system 0 has no record.
