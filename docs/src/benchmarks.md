# Benchmarks

The benchmark suite is defined once and then executed in three different environments:

- **Local desktop** for quick smoke runs
- **GitHub Actions** for the `small` tier on pushes and pull requests
- **SLURM** for the `large` tier and full cluster-scale runs

The suite is sharded by operation family and deployment tier:

- `assembly-small`, `solver-small`, `extension-small`
- `assembly-large`, `solver-large`, `extension-large`

Only the `small` tier is intended for GitHub Actions. The `large` tier is reserved for SLURM or explicit local runs.

## Local commands

Instantiate the benchmark environment once:

```julia
julia --project=bench -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
```

Run all small shards locally:

```julia
julia --project=bench bench/run_benchmarks.jl
```

Run one shard explicitly:

```julia
BENCHMARK_PROFILE=small BENCHMARK_SHARD=assembly-small julia --project=bench bench/run_benchmarks.jl
```

Render a merged report from prior artifacts:

```julia
BENCHMARK_INPUT_DIR=bench/results BENCHMARK_OUTPUT_DIR=bench/results julia --project=bench bench/render_report.jl
```

Compare the current checkout against `main` with `PkgBenchmark.jl`:

```julia
BENCHMARK_PROFILE=small BENCHMARK_BASELINE_REF=main julia --project=bench bench/compare_benchmarks.jl
```

## GitHub Actions

The dedicated `benchmarks.yml` workflow runs the `small` shards as a matrix and then aggregates the uploaded artifacts into a merged HTML and Markdown report.

## SLURM

The `bench/slurm_benchmarks.sh` launcher submits one `large` shard per array task and then runs a dependent aggregation job.

## Generated report

The latest checked-in benchmark report placeholder lives here:

- [Generated benchmark report](generated/benchmark_report.md)
- [Generated benchmark report](benchmark_report.md)
