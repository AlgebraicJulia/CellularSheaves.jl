# 018 — Assume–guarantee certificates over relative-mode summaries

**Status:** parked. We are doing cooperative multi-agent planning first
(optimistic relative-mode bounds, studied with the path-lifting paper).
This issue records the non-cooperative extension so it is not lost. Revisit if
the agents must be certified separately, or if some of them are not controlled
by the planner.

## Mathematical Background

In cooperative planning one planner controls every agent. Agent ``i``'s coarse
problem keeps its own state ``x_i`` and a small *mode* summary
``m = m(x_{-i})`` of the others: for example, the sector, distance band and
closing/opening bit of its nearest neighbour. The **optimistic** quotient along
``\mu_i(x) = (x_i, m)`` assumes that the others move favourably, so its value is
an admissible lower bound
([`TropicalAbstraction`](../src/features/tropical_abstraction.md)).

Without a single planner, agent ``i`` must be certified against *any* behaviour of
the others allowed by an **assumption**. Write the assumption as a relation
``A_i \subseteq M \times M`` on modes (``m \to m'`` is allowed). The certified
value is the min–max fixed point

```math
V(x_i, m) = \min\Big(b(x_i),\ \min_{x_i \to x_i'} \Big[ W_{x_i x_i'} + \max_{(m, m') \in A_i} V(x_i', m') \Big]\Big),
```

where a pair ``(x_i', m')`` that violates the safety constraint (a collision mode)
has value ``\infty``. Where ``V < \infty`` it is a control Lyapunov function for
agent ``i`` that holds for every behaviour of the others consistent with ``A_i``.
This is the ∀∃ ("must") side of the may/must abstraction, with the other agents
as the universally quantified player.

This is the "``M_i \parallel A_i \models G_i``" step of assume–guarantee
reasoning (Misra & Chandy 1981; Jones 1983; Abadi & Lamport 1993, 1995;
Henzinger, Qadeer & Rajamani 1998; Benveniste et al. 2018). The guarantee
``G_i`` is agent ``i``'s effect on its neighbours' modes under the certified
policy.

The operator is monotone but not min-plus *linear*: the ``\max`` over mode
successors is a superior function. Knuth's generalization of Dijkstra's
algorithm (Knuth 1977) computes it with one priority queue. A state is final
once all of its mode successors under some own-move are final. This is
retrograde analysis for a game in which the agent minimizes and the
environment maximizes.

**Composition** (out of scope here) puts the contracts on the interaction graph
as a sheaf of lattices. Vertex stalks hold the contracts ``(A_i, G_i)``, edge
stalks hold sets of relative-mode transitions, and compatibility is a lax-section
condition, with negotiation as a Tarski fixed point (Ghrist & Riess 2022).
Circular safety contracts compose by induction on time. Circular *liveness*
does not.

## Codebase Orientation

| File | Why it matters |
|---|---|
| `src/network_sheaves/TropicalAbstraction.jl` | `robust_quotient` (the ∀∃ quotient when the universal quantifier ranges over the states of a block), `is_control_lyapunov`, `greedy_policy`, `is_lyapunov` (the certificate checks to reuse), `_astar` and the heap helpers (the priority-queue code to extend). |
| `src/network_sheaves/GraphHomomorphisms.jl` | `GraphHomomorphism` for the mode map ``\mu_i``. |
| `scripts/tropical_multiagent_study.jl` of the path-lifting paper | the cooperative experiment, with the relative-mode maps this issue would reuse. |
| `docs/src/features/tropical_abstraction.md` | the lax/oplax framework; add a section if implemented. |
| `test/network_sheaves/TropicalAbstraction.jl` | tests go here (one test file per module). |

## Requested Implementation

One function in `TropicalAbstraction.jl`, exported:

```julia
"""
    robust_mode_value(W, assumption, collision, b) -> (V, policy)

The certified cost-to-go of one agent whose state graph is `W` (``n × n`` min-plus
costs) while the rest of the world is summarized by a mode in `1:K` that evolves
nondeterministically within `assumption` (a ``K × K`` Boolean matrix,
`assumption[m, m′]` when the others may move the mode from `m` to `m′`):

    V(x, m) = min(b[x], min_{x → x′} W[x, x′] + max_{m′ : assumption[m, m′]} V(x′, m′)),

with `V(x, m) = ∞` where `collision[x, m]`. Returns `V` as an ``n × K`` matrix and
the feedback `policy[x, m]`, the agent's next state (0 where `V` is infinite). Where
`V` is finite it is a control Lyapunov function for every behaviour of the others
allowed by `assumption`; it is at least the optimistic value of the same mode graph
(the ∃∃ quotient).
"""
function robust_mode_value(W::SparseMatrixCSC, assumption::AbstractMatrix{Bool},
        collision::AbstractMatrix{Bool}, b::AbstractVector) end
```

Algorithm sketch (Knuth 1977):

1. Validate with `@argcheck`:
   - `W` is square with nonnegative costs;
   - `size(assumption) == (K, K)`;
   - `size(collision) == (n, K)` and `length(b) == n`;
   - every row of `assumption` is nonempty (the world always has a next mode).
2. Initialize `V[x, m] = b[x]` for non-colliding states, and `∞` for colliding
   ones and elsewhere. Push the finite values.
3. For each pair (own move ``x → x'``, mode ``m``), keep a counter of the mode
   successors ``m'`` of ``m`` that are not yet final, and the largest final value
   among them.
4. Pop the smallest ``(x', m')``. For every predecessor ``x → x'`` and every ``m``
   with `assumption[m, m′]`, decrement the counter and update the running max.
   When a counter reaches zero, the candidate
   ``W_{xx'} + \max`` is final for that own-move. Relax ``V(x, m)`` with it,
   recording `policy[x, m] = x′`.
5. Never finalize colliding states.

The cost is ``O((\mathrm{nnz}(W) \cdot \mathrm{nnz}(A)) \log(nK))``.

## Tests to Write

```julia
# A line of 5 cells, goal at 5; one mode, always allowed: the robust value is the distance.
W = spdiagm(1 => ones(4), -1 => ones(4))
V, pol = robust_mode_value(W, trues(1, 1), falses(5, 1), [Inf, Inf, Inf, Inf, 0.0])
@test V[:, 1] == [4, 3, 2, 1, 0]
@test pol[1, 1] == 2

# Two modes: mode 2 blocks cell 3 (collision there). If the world can switch to mode 2 at
# any time, cell 2 can never safely move to 3, so the goal is unreachable from the left.
coll = falses(5, 2); coll[3, 2] = true
A = trues(2, 2)
V, _ = robust_mode_value(W, A, coll, [Inf, Inf, Inf, Inf, 0.0])
@test all(isinf, V[1:2, :])
@test V[4, 1] == 1 && V[5, 2] == 0

# If mode 2 can only be left (assumption: 2 → 1 only, 1 → 1 only), the agent waits in mode 2
# (needs a stay edge) and then proceeds: the robust value is finite.
Ws = W + sparse(1:5, 1:5, fill(1.0, 5))      # stay costs 1
A2 = Bool[1 0; 1 0]
V2, pol2 = robust_mode_value(Ws, A2, coll, [Inf, Inf, Inf, Inf, 0.0])
@test V2[2, 2] == V2[2, 1] == 3                # every successor of mode 2 is mode 1

# Lyapunov: the closed loop lowers V by the move's cost for every allowed mode successor.
for x in 1:4, m in 1:2
    (isfinite(V2[x, m]) && pol2[x, m] != 0) || continue
    y = pol2[x, m]
    @test all(V2[y, m′] + Ws[x, y] <= V2[x, m] for m′ in 1:2 if A2[m, m′])
end

# Robust ≥ optimistic: on random instances, V is at least the cost-to-go of the ∃∃ mode graph
# (own move and mode successor chosen together), built with cost_to_go on the product graph.
```

## Verification Checklist

- [ ] `robust_mode_value` is exported from `TropicalAbstraction` and listed in its module docstring.
- [ ] The docstring cites Knuth (1977) and the may/must abstraction, and explains the min–max
      fixed point for non-specialists.
- [ ] Tests pass, including robust ≥ optimistic on random instances and the Lyapunov decrease
      check for every allowed successor.
- [ ] On the multi-agent study lattice, certify one agent against a "nearest neighbour closes
      by at most one band per step" assumption, and report the certified fraction.
- [ ] The feature guide gains a short "Non-cooperative agents" section.

## Out of Scope

- Contract negotiation, compatibility checks, and the sheaf of contracts on the interaction
  graph (the Tarski fixed point).
- Circular liveness and any small-gain-style composition argument.
- Synthesizing guarantees ``G_i`` automatically from a certified policy.
- Continuous-time or stochastic environments.
