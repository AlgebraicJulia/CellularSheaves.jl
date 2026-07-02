# IPM Benchmark Suite

Benchmarks comparing the CellularSheaves IPM solver against commercial (Mosek) and open-source (Clarabel, OSQP) solvers.

## Quick Start

```bash
# Run with open-source solvers (default)
julia --project=. examples_ipm/spline/nonnegative_spline_lp.jl

# Run with Mosek
julia --project=. examples_ipm/spline/nonnegative_spline_lp.jl --mosek

# Customize runs
julia --project=. examples_ipm/opf/opf_sheaf.jl --nruns=10 --nwarmup=3
```

## Command-Line Options

| Flag | Description |
|------|-------------|
| `--open` | Use open-source solvers (default) |
| `--mosek` | Use Mosek (commercial) |
| `--nruns=N` | Number of benchmark runs (default: 5) |
| `--nwarmup=N` | Number of warmup runs (default: 2) |

All benchmarks run both primal and dual forms (via Dualization.jl).

## Solver Selection

| Problem Type | Open-Source | Commercial |
|--------------|-------------|------------|
| LP/QP only | OSQP | Mosek |
| SOC/SDP/Exp cones | Clarabel | Mosek |

## Benchmark Groups

### 1. Spline Fitting (`spline/`)

Shape-constrained spline estimation problems. Tests LP, SOC, SDP, and exponential cones.

| File | Cones | Description |
|------|-------|-------------|
| `nonnegative_spline_lp.jl` | LP | Nonnegativity via control-point constraints |
| `nonnegative_spline_exact.jl` | SDP+SOC | Exact nonnegativity via SOS |
| `tjunction_spline.jl` | LP | T-junction corridor fitting |
| `tjunction_sdp.jl` | SDP | T-junction with CVX constraints |
| `mle_spline.jl` | Exp | MLE density estimation |
| `tensor_spline.jl` | Exp | Tensor-product splines |

See [spline/README.md](spline/README.md) for details.

### 2. Optimal Transport (`map/`)

Wasserstein transport, quantum marginals, and related SDP problems.

| File | Cones | Description |
|------|-------|-------------|
| `map_lp.jl` | LP+Exp | Discrete transport maps |
| `map_sdp.jl` | SDP | Channel capacity problems |
| `quantum_marginal.jl` | SDP | Quantum marginal consistency |
| `wasserstein_graph.jl` | LP+Exp | Graph Wasserstein with Laplacian Grams |
| `wasserstein_gauss.jl` | SDP | Gaussian (Bures-Wasserstein) transport |

See [map/README.md](map/README.md) for details.

### 3. Optimal Power Flow (`opf/`)

AC-OPF relaxations as a ladder of covers. Tests SOC, SDP, and mixed cones with native quadratic objectives.

| File | Cones | Description |
|------|-------|-------------|
| `opf_sheaf.jl` | SOC+SDP | Full OPF benchmark suite |

See [opf/README.md](opf/README.md) for details.

## What We Measure

Each benchmark reports:
- **IPM(ms)**: CellularSheaves IPM solve time
- **Primal**: Reference solver (primal form) solve time
- **Dual**: Reference solver (dualized via Dualization.jl) solve time
- **P/IPM**: Speedup ratio (primal / IPM)
- **D/IPM**: Speedup ratio (dual / IPM)

## Key Findings

| Structure | IPM vs Commercial |
|-----------|-------------------|
| LP/QP with dense Q | IPM wins 2-10x |
| SOC with dense Q | IPM wins 2-5x |
| Large SDP cells | Mosek wins 2-5x |
| Mixed cones + dense Q | IPM competitive |

The sweet spot for our solver is **problems with dense per-cell Gram matrices** (e.g., state estimation, weighted completion).
