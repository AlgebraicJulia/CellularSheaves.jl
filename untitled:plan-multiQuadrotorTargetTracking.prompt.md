# Multi-Quadrotor Target Tracking — Literate Example Plan

## Goal

Create a new literate example `docs/literate/control/multi_quadrotor_target_tracking.jl`
demonstrating two planar quadrotors each tracking a distinct 2-D target using:
- A `ControlledTrajectorySheaf` with a two-vertex base sheaf
- LQR optimal control solved via the trajectory affine space
- The sparse LDL⊤ solver (`ChordalLDLt`) for the reduced Hessian
- Visualization of the LQR state/control trajectory
- Visualization of the nullspace basis as a trajectory family

## Key Design Decision (pending fix)

The restriction maps on the agent–target coupling edges must **project onto the y-coordinate only**.
The quadrotor state is `[y, z, φ, ẏ, ż, φ̇]` (6D) but the target "state" is just a 2-D position `(y, z)`.
Coupling the full 6D stalk to a full 6D terminal boundary condition is infeasible unless the
terminal velocity/roll states are exactly zero — which they are (hover), but z may differ.

**Resolution:** Keep the boundary conditions standard (both agents start and end at hover with
desired lateral displacements). The issue was setting `xk1[2] = 0.5` (a nonzero z target) while
the *feasibility check* tests reachability under the ZOH dynamics. For the planar quadrotor the
coupling term `−g φ` in `ẍy` means a z-displacement requires net vertical thrust, which is
feasible only if the trajectory has nonzero z-dynamics throughout. The `ControlledTrajectorySheaf`
encodes feasibility via the sheaf coboundary, so **any** endpoint pair that is reachable under
the k-step ZOH dynamics is fine — the problem was a mismatch between the horizon length `k=16`
and the chosen endpoint z-offset.

Two options:
1. **Simplest fix:** Use only lateral (y) displacement targets; set z to 0 for both agents.
   This makes the example directly analogous to the single-quadrotor example.
2. **More interesting:** Build a separate tracking sheaf where the edge between agent vertex
   and target vertex uses a projection restriction map `Pᵧ : ℝ⁶ → ℝ¹` (picks out y).
   This is the "correct" sheaf-theoretic framing for target tracking.

## Current State

- `multi_quadrotor_target_tracking.jl` created but fails with:
  `ArgumentError: Endpoint states are not mutually reachable` because `xk1[2] = 0.5`
  (z-target for Q1) was set but z-reachability fails with the chosen horizon.
- `docs/make.jl` already updated to include the new page.

## Fix Required

Edit `multi_quadrotor_target_tracking.jl`:

```julia
# Change endpoint conditions to y-only displacements (z stays 0):
x1  = zeros(n)
xk1 = zeros(n)
xk1[1] = 0.5   # Q1 moves 0.5 m laterally
xk1[7] = -0.5  # Q2 moves -0.5 m laterally (opposing direction for visual interest)
```

This matches the single-quadrotor example pattern exactly and is guaranteed feasible.

## Visualization Plan

1. **Heatmap of cost Hessian H** — shows block-diagonal band structure.
2. **Time-series plot** — position/roll for both quadrotors + 4 control inputs (2 panels).
3. **Physical (y, z) plane** — spatial trajectory for each quadrotor with time colormap.
4. **Spy plot** of nullspace basis matrix.
5. **Nullspace trajectory family** — lateral position `y₁(t)` and `y₂(t)` for the first
   `min(r, 6)` basis directions overlaid on the optimal trajectory.

## Files to Edit

| File | Change |
|------|--------|
| `docs/literate/control/multi_quadrotor_target_tracking.jl` | Fix endpoint z-offset; verify runs end-to-end |
| `docs/make.jl` | Already updated ✓ |

## Verification Checklist

- [ ] `julia --project=docs docs/literate/control/multi_quadrotor_target_tracking.jl` runs without error
- [ ] All endpoint and dynamics assertions pass (no warnings)
- [ ] Optimality residual `< 1e-10`
- [ ] Plots are generated without error (no `nothing` returned)
- [ ] `julia --project=docs docs/make.jl --no-literate` builds without cross-reference warnings for the new page

## Out of Scope

- Adding a separate agent–target tracking sheaf with projection restriction maps (that is a
  separate feature and should be its own issue/example).
- Any changes to the `ControlledTrajectorySheaf` struct or `TrajectorySheaves` module.
- Animation/GIF output (already done in plotter.jl on the multi-agent-target-tracking branch).
- Multi-agent consensus sheaf (separate from this controlled-trajectory example).
