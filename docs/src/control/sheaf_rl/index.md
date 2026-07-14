# Sheaf-Coordinated Learned Control

Multi-agent target tracking can be posed as coordination on a **cellular sheaf**: agents and
targets are vertices of a graph, an edge carries a *restriction map* relating the states it
joins, and the harmonic extension of the sheaf Laplacian assigns each agent a coordination
reference. This guide presents a **layered, learned** controller for that problem. The sheaf
provides the coordination reference *analytically*; a learned policy provides the per-agent
control, taking only an agent's own state and its sheaf reference — with **no linear-time-invariant
``(A,B)`` model encoding**, and no model-predictive control (no prediction horizon, no online
optimization). The inner controller is learned *feedback*.

## The two-layer architecture

| Layer | Role | Learned? | Uses ``(A,B)``? |
|-------|------|----------|-----------------|
| **Coordination** | sheaf harmonic extension → per-agent reference ``b^\star_i`` | no — topology / linear algebra | no |
| **Execution** | learned residual on a stabilizing base controller, ``u_i = u_{\text{base}} + \alpha\,\Delta_\theta(o_i)`` | only the residual ``\Delta_\theta`` | no (policy); the base may use a fixed gain |

### Layer 1 — coordination (analytic)

Let ``\mathcal F`` be a cellular sheaf over the graph of agents and targets, with restriction
maps ``\mathcal F_{i\trianglelefteq ij}`` and sheaf Laplacian ``L_{\mathcal F}``. Partition
``L_{\mathcal F}`` into the agent block ``H`` (the agent interaction map ``\mathcal H = L_q + D``)
and the agent–target coupling ``B``, and pin the
current target positions ``p`` as boundary data. **The solution of Layer 1 is the harmonic
extension** — the cochain that assigns each agent its reference. (Targets are *uncontrolled*:
they evolve by their own autonomous motion ``\dot p_k = f_k(p_k,t)``, with no control input;
in the demonstrations below this is instantiated by prescribed trajectories whose current
positions are pinned as the boundary data.)

```math
q^\star = H^{-1} B\, p, \qquad b^\star_i = q^\star_i,
```

equivalently the unique ``q`` with ``(L_{\mathcal F}\,x)\big|_{\text{agents}} = 0`` for the
targets held fixed. It is recomputed each step as the targets move; nothing in it is learned.

### Layer 2 — execution (learned residual)

```math
u_i = \operatorname{clip}_{\bar u}\!\Bigl( u_{\text{base}}(x_i, b^\star_i) \;+\; \alpha\,\Delta_\theta(o_i) \Bigr),
\qquad o_i = \bigl[\,x_i;\; b^\star_i - p_i;\; \text{mem}_i\,\bigr].
```

- ``\pi_\theta`` — **the policy**, ``u_i = \pi_\theta(x_i, b^\star_i)``: state and sheaf reference
  in, control out (no explicit ``(A,B)`` prediction in the learned map).
- ``u_{\text{base}}`` — an analytic stabilizing controller. When the control-effectiveness
  ``g_i`` is full row rank, it is the harmonic feedback law of the source construction,
  ``u_{\text{base}} = -k\,g_i^{+}\,\eta_i``, where ``\eta = \mathcal H q - B p =
  (L_{\mathcal F} z)\big|_{\text{agents}}`` is the sheaf **disagreement** and ``g_i^{+}`` the
  pseudoinverse; for the fully-actuated integrators ``g_i = I``, so this reduces to
  ``-k\,(L_{\mathcal F} z)_i``. The underactuated quadrotor has ``g`` *not* full row rank, so
  ``g^{+}`` cannot right-invert it and we use an LQR about hover instead.
- ``\Delta_\theta`` — the learned **residual**, a single weight-shared ``\tanh`` network;
  ``\alpha`` its scale. *Only this is trained;* at initialization ``\Delta_\theta\approx 0``, so
  the policy is exactly the base controller.
- ``\operatorname{clip}_{\bar u}(\cdot)`` — actuator **saturation** to ``\lVert u\rVert\le\bar u``.
- ``o_i`` — the observation; ``\text{mem}_i`` an optional history window ``\{(x,u)_{t-H:t-1}\}``.

The residual is trained by TD3 under an anchor,
``\max_\theta \mathbb E[Q_\phi(o,\Delta_\theta(o))] - \beta\,\mathbb E\lVert\Delta_\theta(o)\rVert^2``,
with reward ``r_i = -\lVert p_i - q^\star_i\rVert``. The anchor keeps the policy at the base
controller unless the critic clearly rewards departing from it. Per control step, one
coordination solve produces every ``b^\star_i``; the same network is then evaluated once per
agent on that agent's own observation (a single batched forward pass), so the inner controller is
**decentralized** and the *same* policy runs at any team size.

---

The demonstrations in this guide increase in difficulty. Each states its **dynamics**, its **sheaf
structure** (vertex stalks, edge families, restriction maps), and its **coordination
requirements** before the result.

Throughout, the reported **tracking error** is the mean of ``\lVert x_i - q^\star_i\rVert`` over
the trajectory and the fleet — each agent's distance to its *sheaf reference* ``q^\star_i`` (the
harmonic extension of Layer 1), **not** to a raw target. For an unpinned agent ``q^\star_i`` is
defined purely by consensus, so a small error means the *whole coordinated configuration* is being
held, not just individual target-following. Lower is better. The "analytic" curve is the sheaf
feedback law ``u_{\text{base}}`` alone (no learned residual); "sheaf+RL" is the full policy.

## In this guide

**Start here:** [The pipeline: problem, architecture, and status](pipeline.md) states the optimal
control problem the whole stack is a relaxation of, walks each layer, and gives an honest status —
what is implemented, what is verified, and what is still open. The demonstrations below are the
evidence for it.

1. [Replicating the analytic controller (simple → complex sheaf)](replicating.md) — the learned
   policy first reproduces the analytic sheaf controller, on progressively richer sheaves
   (Demonstrations 1 and 4).
2. [Where learning wins: drift, underactuation](rl_wins.md) — unknown time-varying drift, the
   underactuated quadrotor, and the multi-agent structure under drift, plus data collection and
   training (Demonstrations 2, 3, 5).
3. [Generalization and genuinely heterogeneous sheaves](generalization.md) — zero-shot
   generalization to unseen drifts and a genuinely heterogeneous sheaf (Demonstrations 6, 7).
4. [Hard actuator constraints, without MPC](constrained.md) — second-order-cone actuator limits
   as a closed-form projection, and what the learning adds overall (Demonstration 8).
5. [A hard safety layer: control barrier functions](cbf.md) — a decentralized CBF safety filter
   that guarantees inter-agent collision avoidance, composed on top of either the analytic or
   learned execution layer.
