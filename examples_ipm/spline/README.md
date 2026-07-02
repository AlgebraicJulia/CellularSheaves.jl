# Spline Fitting Benchmarks

Shape-constrained Bézier spline estimation problems. The domain is partitioned into pieces; each piece has polynomial coefficients living on a cell, and pieces are joined by C^k continuity constraints (edges in the sheaf).

## Problems

### LP Problems (OSQP compatible)

**nonnegative_spline_lp.jl** — Nonnegativity via control-point constraints (sufficient but not exact). The shape cone sits directly on the segment stalk.

```bash
julia --project=. examples_ipm/spline/nonnegative_spline_lp.jl
julia --project=. examples_ipm/spline/nonnegative_spline_lp.jl --mosek
```

**tjunction_spline.jl** — T-junction corridor fitting with positivity constraints.

```bash
julia --project=. examples_ipm/spline/tjunction_spline.jl
```

### SDP Problems (Clarabel compatible)

**nonnegative_spline_exact.jl** — Exact nonnegativity via SOS (Gram matrices). The SDP cone enforces polynomial nonnegativity exactly.

```bash
julia --project=. examples_ipm/spline/nonnegative_spline_exact.jl
```

**tjunction_sdp.jl** — T-junction with CVX (convex-over-variable) constraints via SDP.

```bash
julia --project=. examples_ipm/spline/tjunction_sdp.jl
```

### Exponential Cone Problems

**mle_spline.jl** — Maximum likelihood density estimation. Uses exponential cones for the log-likelihood.

```bash
julia --project=. examples_ipm/spline/mle_spline.jl
```

**tensor_spline.jl** — Tensor-product splines with entropy regularization.

```bash
julia --project=. examples_ipm/spline/tensor_spline.jl
```

## Performance Summary

| Problem | IPM vs Mosek | Notes |
|---------|--------------|-------|
| LP (nonneg, monotone, convex) | IPM wins 1.5-4x | Block-diagonal Q |
| SDP (exact nonneg) | Mosek wins 2-5x | Large PSD cells |
| Exp (MLE) | Mosek wins 2-10x | Exponential cone |
| T-junction LP | IPM wins 2-3x | Many small cells |
| T-junction SDP | Mosek wins 3-5x | Large PSD cells |

The IPM wins on LP problems with many small cells and block-diagonal objectives.
