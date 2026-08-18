# 006 — Closed-loop MPC for multi-agent tracking sheaf problems

## Mathematical Background

Current examples solve a single open-loop quadratic program on a fixed horizon. We want a receding-horizon controller.

Let the aggregate agent state at MPC step `t` be

```math
x_t = (x_t^{(1)}, \dots, x_t^{(n_a)}).
```

At each step, solve a finite-horizon sheaf-constrained quadratic program over a prediction window of length `W`:

```math
\min_{z_t} \; J_t(z_t)
```

subject to

```math
z_t \in \mathcal{T}_t(x_t, \hat{y}_{t:t+W}),
```

where `\mathcal{T}_t` is the feasible affine trajectory space induced by the sheaf constraints (agent dynamics, consensus, tracking edges, and pinned target vertices on the local window), and `\hat{y}` is the target forecast over the local window.

Use the same quadratic objective family already used in open-loop mode:

```math
J_t(z) = \tfrac12 z^\top H_t z + f_t^\top z,
```

with control-penalty structure assembled by existing `QuadraticCosts` utilities.

After solving, apply only the first control block:

```math
u_t = \Pi_u z_t^\star,
```

advance the true system one step,

```math
x_{t+1}^{(i)} = A_i x_t^{(i)} + B_i u_t^{(i)},
```

and repeat.

This is the standard receding-horizon MPC loop (Rawlings–Mayne–Diehl), specialized to the sheaf-feasible trajectory set. In sheaf terms, this turns one global open-loop solve into repeated local solves with updated boundary data, while preserving the same algebraic structure used in existing trajectory and harmonic-extension workflows (Hansen & Ghrist style sheaf-constrained dynamics).

## Codebase Orientation

| File | Why it matters |
| --- | --- |
| `src/ControlSheaves/MultiAgentTracking.jl` | Current `run_scenario` open-loop pipeline and data extraction. MPC API should live here. |
| `src/ControlSheaves/QuadraticCosts.jl` | Reuse `build_control_cost_matrix` and `solve_quadratic_on_basis` for each MPC solve. |
| `src/ControlSheaves/TrackingProblems.jl` | `TrackingProblem` structure and indexing helpers needed for state/control slicing. |
| `src/ControlSheaves/ScenarioResults.jl` | Closed-loop rollout should return this existing result type (or a closely related companion type). |
| `docs/literate/control/multi_quadrotor_target_tracking.jl` | Natural first documentation consumer for MPC closed-loop example. |
| `test/network_sheaves/MultiAgentTracking.jl` | Existing tests for `run_scenario`; use style and fixtures for MPC tests. |
| `test/network_sheaves/QuadraticCosts.jl` | Existing quadratic-cost behavior to preserve under repeated solves. |

## Requested Implementation

Add one public closed-loop routine plus a small helper for target-window extraction.

### Exact API

```julia
"""
    run_mpc_scenario(label, prob, x0_agents, target_trajs, times;
                     window,
                     y_col=1,
                     z_col=2,
                     cost=1.0,
                     tracking_weight=nothing,
                     include_target_dynamics=nothing)
        -> ScenarioResult

Run a receding-horizon MPC controller for a `TrackingProblem`.

At each simulation step `t`, solve a local sheaf-constrained quadratic problem on
`[t, min(t+window, k)]` with the current agent state pinned at the local initial
vertices and the provided target trajectory samples pinned on the same window.
Apply only the first control from the local optimum, advance agent dynamics, and
repeat until the end of `times`.

Returns a `ScenarioResult` containing the closed-loop agent trajectories and the
provided target trajectories.
"""
function run_mpc_scenario(
    label::String,
    prob::TrackingProblem,
    x0_agents::Vector{Vector{Float64}},
    target_trajs::Vector{Vector{Vector{Float64}}},
    times::AbstractVector{<:Real};
    window::Int,
    y_col::Int=1,
    z_col::Int=2,
    cost::Union{Number,Function}=1.0,
    tracking_weight=nothing,
    include_target_dynamics=nothing,
)

"""
    window_targets(target_trajs, t0, t1)

Return target trajectories restricted to integer timesteps `t0:t1`.
"""
window_targets(target_trajs, t0::Int, t1::Int)
```

### Implementation sketch

1. Validate inputs with `@argcheck`.
- `window >= 1`
- `length(x0_agents) == prob.n_agents`
- `length(target_trajs) == prob.n_targets`
- trajectory lengths match `length(times)`

2. Keep a mutable closed-loop state vector `x_now[i]` for each agent, initialized from `x0_agents`.

3. For each MPC step `t = 0:(k-1)`:
- Set `t_end = min(t + window, k)`.
- Build a local problem on horizon `t:t_end` using the same edge maps/weights as `prob`.
- Pin local initial agent states to `x_now`.
- Pin local target vertices with `target_trajs[:, t:t_end]`.
- Solve local open-loop via existing flow:
  - sheaf build
  - `harmonic_extension`
  - quadratic reduction with `build_control_cost_matrix` + `solve_quadratic_on_basis`
- Extract first control `u_t^{(i)}` from each agent local trajectory.
- Advance real state with `A_i`, `B_i` from `prob`.

4. Accumulate closed-loop agent trajectories as matrices `((k+1) × nx_i)` and return `ScenarioResult(label, times, ...)` with supplied targets.

5. Residual reporting:
- Define residual as RMS of per-step local sheaf residuals or final-step residual from the last MPC solve.
- Document clearly which convention is used.

6. Keep everything sparse and dependency-free (no new optimization packages).

## Tests to Write

Create `test/network_sheaves/MultiAgentTrackingMPC.jl` and include from `test/runtests.jl`.

Use a small two-agent, two-target setup with short horizon (e.g. `k=8`) to keep tests fast.

Required tests:

1. `run_mpc_scenario` returns `ScenarioResult` with expected lengths and dimensions.
2. Closed-loop initial states equal `x0_agents` exactly.
3. With stable simple dynamics and moderate control cost, all states remain finite.
4. `window=1` runs and acts as one-step greedy MPC.
5. `window=k` runs and produces a trajectory consistent with one-shot open-loop solve (within tolerance in a deterministic fixture).
6. Bad `window` (`0`), mismatched agent counts, and mismatched target lengths throw `ArgumentError`.
7. A regression test where mirrored-consensus map `(R_y, -R_y)` is active confirms lateral symmetry tendency (`y1 + y2` small near constrained steps).

## Verification Checklist

- [ ] New API is declared/implemented in `MultiAgentTracking` and exported consistently.
- [ ] Implementation reuses existing quadratic/sheaf solve path instead of introducing a second optimizer stack.
- [ ] No new package dependencies are added.
- [ ] Input validation uses `@argcheck` at public boundaries.
- [ ] MPC behavior is documented in a literate control example.
- [ ] Tests in `MultiAgentTrackingMPC.jl` pass locally and in CI.

## Out of Scope

- Nonlinear dynamics and nonlinear MPC.
- Inequality constraints (state/control bounds, obstacle avoidance).
- Robust/stochastic MPC, tube MPC, or chance constraints.
- Distributed/decentralized MPC across agents.
- New QP solver integrations beyond current linear-algebra stack.
- Real-time code generation and warm-start benchmarking.

## References

- J. B. Rawlings, D. Q. Mayne, and M. M. Diehl, *Model Predictive Control: Theory, Computation, and Design*, 2nd ed., Nob Hill Publishing.
- J. Hansen and R. Ghrist, "Toward a Spectral Theory of Cellular Sheaves," *Journal of Applied and Computational Topology* (2019). DOI: 10.1007/s41468-019-00038-7, arXiv:1808.01513.
- T. Hanks, H. Riess, S. Cohen, T. Gross, M. Hale, and J. Fairbanks, "Distributed Multi-agent Coordination over Cellular Sheaves," arXiv:2504.02049 (2025).
