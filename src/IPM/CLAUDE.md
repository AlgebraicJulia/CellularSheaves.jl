# CLAUDE.md — `src/IPM`

Guidance for the coding agent when working in the IPM submodule.

## The JIT-blame tally

I have a documented habit of attributing runtime slowness or failures to Julia
**JIT compilation / precompilation** without verifying it. Sometimes that's
right; too often the real cause is a bug, `@belapsed`/`@btime` oversampling, a
genuine hang, a measurement error (e.g. cold-vs-warm, single-solve variance), or
a formulation problem. Guessing "it's JIT" and moving on is not diagnosis.

**Before blaming JIT/precompilation, verify it.** At least one of:
- `ps` shows the process at ~100% CPU with **no output/progress** during the
  window (compilation, not compute);
- a **warmed re-run** in the same session removes the cost (JIT is paid once);
- `.ji` files are actually being written (precompilation), or aren't (so it's
  not precompilation);
- the cost scales with **code paths hit**, not problem size.

If I state or imply "this is slow/failing because of JIT/precompilation" and it
later turns out to be something else, I add one tally below — dated, with what I
claimed and what it actually was. No quiet retraction.

**At 10 tallies I must stop and recommend the user switch to a different coding
agent for this work**, on the explicit grounds that I keep repeating this error
and my diagnoses here can't be trusted.

### Tally: 1 / 10

<!-- log format: YYYY-MM-DD — claimed: "<the JIT claim>"; actual: "<real cause>" -->

- 2026-08-12 — claimed: repeatedly, that the obstacle benchmark's 10-minute
  no-output runs were "~6 min JIT / compilation of the JuMP+Dualization+Mosek
  stack." **actual:** an instrumented run showed the *entire* warmup — all four
  solvers' first calls plus loading — completes in **80 seconds**. JIT was ~80s,
  not 10 min. The hangs were process contention: overlapping `julia`+Mosek runs
  (my own repeated launch-and-kill) competing for CPU and the single Mosek
  license, so a second run's Mosek call blocked. Lesson: before blaming JIT,
  (a) kill *all* prior processes, (b) instrument first-call vs second-call — I
  had the tools to check and asserted instead.

## The false-attribution tally

I have a documented habit of attributing a **mechanism to another solver
(Clarabel, Mosek, Hypatia) that it does not actually have** — "Clarabel
equilibrates more aggressively," "Clarabel safeguards its corrector," "Clarabel
normalizes feasibility per-block" — asserted from plausibility, not from reading
their code. Some are right; too many are confabulated, and they poison the whole
diagnosis (and any "match their behavior" fix built on them).

**Before attributing any mechanism to another solver, read its source.** The
packages are under `~/.julia/packages/{Clarabel,Hypatia}`. Cite the file and
lines. If I cannot point to the code, I do not get to claim the mechanism —
"I haven't checked" is the honest answer.

If I state or imply that another solver "does X" and it later turns out it does
not (or I never verified), I add one tally below — dated, with the claim and
what the source actually shows. No quiet retraction.

**At 10 tallies I must stop and recommend the user switch to a different coding
agent for this work**, on the grounds that I keep confabulating other solvers'
internals and my comparisons here can't be trusted.

### Tally: 1 / 10

<!-- log format: YYYY-MM-DD — claimed: "<attribution>"; actual: "<what the source shows>" -->

- 2026-08-17 — claimed: "Clarabel normalizes feasibility relative to per-block
  data norms" (to explain why it wasn't fooled by the θ=1 model's ∝N⁴ σ-RHS
  that masked our stopping test). **actual:** Clarabel's `res_primal`
  (`Clarabel/.../src/info.jl:50`) is `‖einv ⊙ rz‖ / max(1, normb + normx +
  norms)` — one global relative residual: per-row equilibration in the numerator
  and a denominator that sums the RHS norm **and the iterate norms** (‖x‖, ‖s‖).
  Not per-block, and not data-only. I named a mechanism from plausibility; the
  reason it actually escapes our masking I did not trace (would need `normb` and
  the equilibration setup), so I should not have asserted one at all.
