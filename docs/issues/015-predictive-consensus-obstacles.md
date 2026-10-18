# 015 — Predictive consensus with obstacles (Stage C)

Stages A (unbounded LQ) and B (Euclidean thrust bound) live in
`src/ControlSheaves/PredictiveConsensus.jl`. This issue is the next seam: one
function that adds obstacle avoidance to each agent's local trajectory problem.

## Mathematical Background

Each agent is a planar double integrator, `x = (q, v) ∈ ℝ⁴`, coupled to its
neighbours only through the pinned coordination sheaf: the team minimizes

    Σₜ [ (w/2)‖q(t) − q⋆‖²_H + (c/2)‖v(t)‖² ] + Σₜ (r/2)‖u(t)‖²

subject to per-agent dynamics and `‖uᵢ(t)‖ ≤ ū`. Agent `i` re-plans against
its neighbours' *predicted trajectories* `q̂ⱼ(·)`, which enter only the linear
term `−w Hᵢⱼ q̂ⱼ(t)` of its QP (block coordinate descent = zero-overlap
multiplicative Schwarz on trajectory space).

A circular obstacle `‖q − o‖ ≥ ρ` is nonconvex. The standard convexification
(sequential convex programming; Schulman et al., *Motion planning with
sequential convex optimization and convex collision checking*, IJRR 2014;
Augugliaro, Schoellig & D'Andrea, *Generation of collision-free trajectories
for a quadrocopter fleet: a sequential convex programming approach*, IROS
2012) linearizes the constraint about the current plan `q̄(t)`:

    nₜᵀ (q(t) − o) ≥ ρ,     nₜ = (q̄(t) − o)/‖q̄(t) − o‖,

a half-plane, i.e. one `PositiveCone` slack per step. The half-plane lies
inside the feasible set, so every SCP iterate is collision-free (an inner
approximation), and a trust region `‖q(t) − q̄(t)‖ ≤ δ` (one SOC per step)
keeps the linearization honest. Inter-agent collision avoidance uses the same
cut with `o` replaced by the neighbour's predicted position `q̂ⱼ(t)`, which is
exactly the information already exchanged. For nonconvex problems the sweeps
converge to a Nash point of the agents' games, not the team optimum.

## Codebase Orientation

| File | Why |
|---|---|
| `src/ControlSheaves/PredictiveConsensus.jl` | `_trajectory_qp` builds the local conic QP; `TrajectoryLayout` indexes stalks and rows; `_replan!` is the local update |
| `src/ControlSheaves/CoordinationBenchmarks.jl` | `CoordinationScenario` supplies `H`, `Bmat`, neighbour lists |
| Mumblebee.jl `src/IPM` | `PositiveCone`, `SecondOrderCone`, `reinit!` (f, g, warm start) |
| `test/ControlSheaves/PredictiveConsensus.jl` | existing reference tests (Riccati, centralized QP) to extend |

## Requested Implementation

```julia
"""
    CircularObstacle(center, radius)

A disc the agents' positions must avoid.
"""
struct CircularObstacle
    center::Vector{Float64}
    radius::Float64
end

"""
    ConsensusLQ(...; obstacles = CircularObstacle[], separation = 0.0, trust_radius = Inf)

`separation > 0` adds pairwise collision cuts between neighbours.
"""
```

Algorithm sketch (inside `_replan!`, no new public solver):
1. Linearize each obstacle (and neighbour, if `separation > 0`) about the
   agent's previous plan; append one `PositiveCone` slack stalk per active cut
   (`nᵀq − s = ρ + nᵀo`).
2. Because the cut normals change the constraint matrix, the local
   `IPMSolver` must be rebuilt when the active set changes; keep the solver
   when only the right-hand side moves (fixed set of cuts per obstacle/step).
3. Optional trust region: SOC stalk `(δ, q(t) − q̄(t))`.

## Tests to Write

```julia
obs = CircularObstacle([2.0, 1.5], 0.5)
lq = ConsensusLQ(scenario, plant, 25; targets, obstacles = [obs])
plan = solve(ConsensusProblem(lq, x0), PredictedTrajectorySweeps())
@test all(norm(plan.states[1:2, t, i] - obs.center) >= obs.radius - 1e-6 for t in 1:26, i in 1:9)
@test plan.converged
# no obstacle in the way ⇒ same as Stage A
far = CircularObstacle([100.0, 100.0], 1.0)
@test solve(ConsensusProblem(ConsensusLQ(scenario, plant, 25; targets, obstacles = [far]), x0),
    PredictedTrajectorySweeps()).cost ≈ solve(prob, JointRiccati()).cost rtol = 1e-6
# collision cuts keep neighbours apart
@test minimum(norm(plan.states[1:2, t, i] - plan.states[1:2, t, j])
    for t in 1:26, (i, j) in neighbour_pairs) >= separation - 1e-6
```

## Verification Checklist

- [ ] Every iterate is collision-free (inner approximation), not only the last.
- [ ] Rebuild of the local solver happens only when the cut set changes.
- [ ] Receding-horizon rollout around an obstacle reaches the formation.
- [ ] Docstrings cite the SCP references above.

## Out of Scope

- Nonlinear dynamics (unicycle, quadrotor): Stage D, via linearization of the
  dynamics about the predicted trajectory (same SCP loop).
- HJB value-function methods on the 4-D state grid (`32⁴` points per step);
  revisit only for small fleets or as a terminal cost.
- Polygonal obstacles (union of half-planes needs integer choices).
