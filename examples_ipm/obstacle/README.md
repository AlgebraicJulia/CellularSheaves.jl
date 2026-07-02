# Obstacle / American Option Benchmarks

Obstacle problems and American options — where **Q is the problem**, not regularization.

## The Setup

The objective is the PDE energy (dense-banded per-patch stiffness), cells are overlapping
subdomain patches (Schwarz domain decomposition), and the orthant constraint carries the
option's **early-exercise premium** w = u - ψ ≥ 0.

Key characteristics:
- **Dense Q per cell**: weighted stiffness + lumped mass (PD natively, ε = 0)
- **Macroscopic active set**: ~60% of variables on the cone boundary at optimum
- **No relaxation gap**: pure-Q, pure-glue flagship

## Files

| File | Problem | Cones |
|------|---------|-------|
| `obstacle_option.jl` | Perpetual/American put, game options | PositiveCone, CofreeCone |

## Running

```bash
# Default: OSQP (open-source)
julia --project=examples_ipm examples_ipm/obstacle/obstacle_option.jl

# With Mosek
julia --project=examples_ipm examples_ipm/obstacle/obstacle_option.jl --mosek
```

## Modes

- **`:dirichlet`** — standard obstacle (perpetual American put)
- **`:game`** — two-sided obstacle ψ ≤ u ≤ ψ + δ (Israeli/game options)
- **`:minimal_surface`** — SOC leaves for min ∫√(1+u'²), equals Dirichlet in 1D

## Time Stepping

The finite-horizon American put is a chain of 200 obstacle QPs with **fixed Q** —
only c changes. This is the cleanest factorization-reuse showcase: same Q, same B,
same cones, new linear term, two hundred times.

## References

Merton (1973), Cryer SIAM J. Control (1971), Jaillet-Lamberton-Lapeyre (1990),
Kifer Finance Stoch. (2000), Lions (1988), Badea SINUM (1991, 2006).
