# Network Utility Maximization Benchmarks

Kelly's NUM framework as a sheaf on the source-link bipartite graph — where
**the dual is deployed engineering**: dual decomposition is literally distributed
TCP congestion control.

## The Fairness-Cone Ladder

| Mode | α | Utility | Cones | Protocol |
|------|---|---------|-------|----------|
| `:throughput` | 0 | w·x | LP | — (starves long paths) |
| `:proportional` | 1 | w·log x | ExponentialCone | TCP Vegas |
| `:delay` | 2 | −w/x | SecondOrderCone | TCP Reno |
| `:maxmin` | ∞ | shared floor | LP (Cofree + slack) | — |

## Files

| File | Cones | Description |
|------|-------|-------------|
| `num_sheaf.jl` | LP/Exp/SOC | Four-mode builder, price reader, demos |

## Running

```bash
# Default: Clarabel (open-source)
julia --project=examples_ipm examples_ipm/num/num_sheaf.jl

# With Mosek
julia --project=examples_ipm examples_ipm/num/num_sheaf.jl --mosek
```

## Kelly's Identity

At the proportional-fair optimum: w_s/x_s = Σ_{l∈path} λ_l

The capacity duals λ_l are **link congestion prices** — the signal that TCP uses
for distributed rate control.

## References

Kelly-Maulloo-Tan (1998), Low-Lapsley (1999), Mo-Walrand (2000),
Chiang-Low-Calderbank-Doyle "Layering as Optimization Decomposition" (2007).
