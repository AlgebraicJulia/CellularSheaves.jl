# Trajectory and Control

This page highlights the externally facing APIs for trajectory sheaves and control-oriented examples.

## When to use this page

Use this guide when you want to:

- Encode a discrete-time dynamical system as a sheaf.
- Solve boundary-value trajectory problems via sheaf methods.
- Work with controlled trajectory models and platoon-style examples.

## Trajectory Sheaf Construction

- Main APIs: `TrajectorySheaf`, `build_trajectory_sheaf`, `colocation_trajectory`.
- See: [Trajectory Sheaves](../api/trajectory_sheaves.md).

## Controlled Dynamics

- Main APIs: `ControlledTrajectorySheaf`, `controlled_vehicle_platoon_with_herding`.
- See: [Trajectory Sheaves](../api/trajectory_sheaves.md) and [Herding Platoon Utilities](../api/herding_platoon.md).

```julia
# Sketch: discrete-time trajectory sheaf from x_{t+1} = A * x_t
S = TrajectorySheaf(A, horizon)

# Colocation solve for endpoint constraints
traj = colocation_trajectory(S, x0, xT)
```

For end-to-end worked examples, start in:

- `generated/trajectory_sheaf.md`
- `generated/control/controlled_vehicle_platoon.md`
- `generated/control/controlled_mass_spring_damper_chain.md`
