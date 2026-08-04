# Usage and API

The API preserves the architecture's isolation boundary: solve the harmonic
problem first, then filter its reference. `TikhonovFilter` does not own ``H``
and does not know whether ``q^\star`` came from a dense solve, cached clique
tree, distributed message passing, or differentiable conic layer.

## Basic use

```julia
using CellularSheaves.ControlSheaves.Tikhonov

q0 = tikhonov_equilibrium(H, rhs0)
q1 = tikhonov_equilibrium(H, rhs1)
filter = TikhonovFilter(q0; epsilon=0.2)

tikhonov_step!(filter, q0, q1, dt)
delivered_reference = filter.x
```

For a callable reference, the convenience method evaluates both interval
endpoints:

```julia
tikhonov_step!(filter, qstar_at, t, dt)
```

## Exact first-order-hold update

Between consecutive solved references, approximate

```math
q^\star(t_k+s)=q_k+\frac{s}{\Delta t}(q_{k+1}-q_k).
```

The normalized linear ODE can then be integrated exactly. With
``z=\Delta t/\epsilon`` and ``\rho=e^{-z}``, the implementation uses

```math
x_{k+1}=\rho x_k+(1-\rho)q_k+
\left[1-\frac{1-\rho}{z}\right](q_{k+1}-q_k).
```

This is not an RK4 approximation. It has three practical advantages:

- no ODE stability restriction couples ``\Delta t`` to ``\epsilon``;
- each interval needs only its two harmonic endpoint references;
- as ``\epsilon\to0^+``, ``x_{k+1}\to q_{k+1}`` without numerical blowup.

The small-``z`` coefficient is evaluated by series expansion and `expm1` to
avoid cancellation. ``\epsilon=0`` is rejected; use the direct harmonic
reference when the algebraic limit is desired exactly.

## Feedforward

For known or estimated harmonic-reference velocity:

```julia
qdot0 = tikhonov_reference_rate(H, rhs_rate0)
qdot1 = tikhonov_reference_rate(H, rhs_rate1)

u0 = tikhonov_feedforward_reference(q0, qdot0, epsilon)
u1 = tikhonov_feedforward_reference(q1, qdot1, epsilon)
tikhonov_step!(filter, u0, u1, dt)
```

With a cached tree factor, solve `rhs_rate` through the same workspace used for
position rather than refactorizing or forming ``H^{-1}``.

## Diagnostics

For stationary-target experiments,

```julia
error = filter.x - qstar
Vdot = tikhonov_dissipation(error, filter.epsilon)
```

returns ``-norm(error)^2 / epsilon``. For moving targets this value is only the
homogeneous part of the derivative; include the
``-e^\top\dot q^\star`` input term when evaluating the full Lyapunov balance.

## Tuning checklist

1. Bound or measure ``\lVert\dot q^\star\rVert`` over representative target motion.
2. Choose ``\epsilon`` from the acceptable uncompensated lag tube.
3. Use feedforward when reference velocity is sufficiently accurate.
4. Choose the sample period for target and actuator bandwidth, not RK4 stability.
5. Rebuild the harmonic provider if topology, restrictions, or pinned cells change.

## API reference

Full docstrings are in the [Tikhonov](@ref) API reference.
