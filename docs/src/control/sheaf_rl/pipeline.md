# The pipeline: problem and proposed solution

A one-page map of what we are solving, how the architecture follows from it, and where each piece
stands. The [architecture page](index.md) has the details; the demonstrations are the evidence.

## The problem

Everything in this guide is a relaxation of a single optimal control problem. For a team of agents
whose coordination requirements are encoded by a cellular sheaf ``\mathcal F``:

```math
\min_{u(\cdot)} \int_0^T \Bigl[\; \underbrace{\tfrac12\, z^\top L_{\mathcal F}\, z}_{\text{coordination}}
\;+\; \underbrace{\tfrac12\, \lVert u \rVert_R^2}_{\text{effort}} \;\Bigr]\, \mathrm{d}t
\qquad \text{s.t.}\quad
\begin{aligned}
&\dot x_i = f_i(x_i, t) + g_i(x_i)\, u_i, \quad f \text{ unknown}\\
&\lVert u_i \rVert \le \bar u \quad \text{(actuator cone)}\\
&h_{ij}(x) \ge 0 \quad \text{(safety)}
\end{aligned}
```

The coordination cost is the **sheaf Dirichlet energy** — total disagreement across every edge,
measured through the restriction maps. This is the point of the formalism: coordination cost *is*
Dirichlet energy, so the sheaf appears in the objective rather than being imposed by hand. Targets
are boundary data; agents are the free cells.

As stated it is intractable — ``f`` is unknown, safety is non-convex in ``x``, and solving it over a
horizon every step is the model-predictive control we are trying to avoid.

## The proposal: four relaxations

| | Relaxation | Assumes | Gives |
|---|---|---|---|
| **R1** | Coordination | no constraints, ``f`` known, ``g_i`` full row rank | the **harmonic extension** ``q^\star`` — and the optimal feedback is gradient flow on the Dirichlet energy |
| **R2** | Execution | ``f`` is **unknown** | a learned **residual** on that base law, searching near the known-``f`` solution |
| **R3** | Actuation | keep the actuator cone | ``\lVert u\rVert \le \bar u`` enforced without an online solver |
| **R4** | Safety | horizon-wide safety → **pointwise** CBF condition | a small convex program each step — no horizon, no MPC |

```math
u_i \;=\; \underbrace{\text{actuation}}_{\text{R3}}\;\circ\;\underbrace{\text{safety}}_{\text{R4}}\;\circ\;
\Bigl(\;\underbrace{u_{\text{base}}(x_i, b^\star_i)}_{\text{R1}} \;+\; \alpha\,\underbrace{\Delta_\theta(o_i)}_{\text{R2}}\;\Bigr)
```

**The claim this buys us.** All topological reasoning lives in R1, so the learned policy sees only
its own state and its own reference — never a restriction map, never the graph. The *same*
decentralized policy therefore runs on any sheaf, at any team size.
[Demonstration 7](generalization.md) is the evidence: rotations and subspace projections change
``b^\star``, and nothing in the policy changes at all.

## Where each layer stands

| | What we have | Next capability |
|---|---|---|
| **R1** Coordination | `harmonic_extension` on a chordal sparse factorization | **distributed solve** — the multifrontal solve is already scheduled by the clique tree (bottom-up reduce, top-down broadcast), so disjoint subtrees are independent. The per-node-addressed serial rewrite is **done and verified**; threading and benchmarking are next |
| **R2** Execution | residual TD3 + behaviour cloning; 8 demonstrations; cluster training | scale the structure-blindness claim across many sheaves and team sizes |
| **R3** Actuation | closed-form projection onto the actuator ball | **constraint-aware RL** — train against a conic oracle rather than projecting after the fact |
| **R4** Safety | [decentralized CBF filter](cbf.md) on the sheaf IPM; holds the safe distance exactly | **high-order CBF** for the underactuated quadrotor |
| Visualization | Isaac Sim pipeline (3-D rollouts + 2-D companions) | hardware |

## The open question

R3 and R4 are currently solved **separately** — the barrier program does not know about the actuator
cone, and the projection is applied afterwards. That is only sound while the cap is slack. When the
actuator limit binds, the projection can move the command off the barrier constraint the filter just
solved for.

The resolution is to solve **one conic program** — minimize ``\lVert u - u_{\text{nom}}\rVert^2``
subject to the barrier constraints *and* the actuator cone. This is exactly the problem class the
IPM already handles, and it is the point at which the conic machinery becomes necessary rather than
convenient.

Doing so surfaces the real question: **safety under bounded actuation can be infeasible.** With a
hard cap there exist states from which no admissible input keeps the system in the safe set.
Characterizing that viable set — and deciding what the controller does when it is empty — is where a
theoretical contribution is available to us rather than borrowed.
