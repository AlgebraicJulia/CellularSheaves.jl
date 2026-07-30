# Optimal Transport & Quantum Benchmarks

Wasserstein transport, quantum marginals, and related SDP problems. These problems feature dense per-cell Gram matrices and/or SDP constraints.

## Problems

### LP + Exponential Cone

**map_lp.jl** — Discrete optimal transport maps with entropy regularization.

```bash
julia --project=. examples_ipm/map/map_lp.jl
```

**wasserstein_graph.jl** — Graph Wasserstein distance with Laplacian Grams. Features dense Q from the graph Laplacian.

```bash
julia --project=. examples_ipm/map/wasserstein_graph.jl
```

### SDP Problems

**map_sdp.jl** — Quantum channel capacity and unital channel problems.

```bash
julia --project=. examples_ipm/map/map_sdp.jl
```

**quantum_marginal.jl** — Quantum marginal consistency (N-representability). Tests whether marginal density matrices are compatible.

```bash
julia --project=. examples_ipm/map/quantum_marginal.jl
```

**wasserstein_gauss.jl** — Gaussian (Bures-Wasserstein) optimal transport. Joint covariance blocks as SDP cells, principal-block compression as restriction maps.

```bash
julia --project=. examples_ipm/map/wasserstein_gauss.jl
```

## Mathematical Background

### Wasserstein Distance

The 2-Wasserstein distance between measures μ and ν is:
```
W₂(μ,ν)² = min_{π} ∫ c(x,y) dπ(x,y)
```
where π ranges over couplings with marginals μ, ν.

### Quantum Marginals

Given a multipartite quantum state ρ_ABC, the marginal consistency problem asks whether given ρ_AB, ρ_BC, ρ_AC are compatible with some global ρ_ABC.

### Gaussian Transport (Bures-Wasserstein)

For Gaussian measures N(m₁,Σ₁) and N(m₂,Σ₂):
```
W₂² = |m₁-m₂|² + Tr(Σ₁ + Σ₂ - 2(Σ₁^½ Σ₂ Σ₁^½)^½)
```
The coupling is characterized by a joint covariance matrix.

## Performance Summary

| Problem | IPM vs Mosek | Notes |
|---------|--------------|-------|
| Discrete transport LP | IPM wins 2-3x | Sparse structure |
| Graph Wasserstein (Higham weights) | IPM wins 5-50x | Dense Laplacian Gram |
| Quantum marginals | Mosek wins 2-5x | Large SDP cells |
| Gaussian transport | Mixed | Depends on d and graph |

**Key finding**: Dense per-cell Q (e.g., Laplacian Grams, Higham weights) shifts the balance toward IPM, even with SDP cells.
