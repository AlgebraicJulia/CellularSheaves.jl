# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

**Run all tests:**
```julia
julia --project=. -e 'using Pkg; Pkg.test()'
```

**Run a single test file:**
```julia
julia --project=. -e 'include("test/network_sheaves/SheafLaplacian.jl")'
```

**Setup benchmark environment (first time only):**
```julia
julia --project=bench -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
```

**Run benchmarks:**
```
julia --project=bench bench/run_benchmarks.jl
```

**Run one benchmark shard:**
```
BENCHMARK_PROFILE=small BENCHMARK_SHARD=assembly-small julia --project=bench bench/run_benchmarks.jl
```

**Run the full suite (all shards, all sizes):**
```
BENCHMARK_PROFILE=full julia --project=bench bench/run_benchmarks.jl
```

**Compare against a baseline ref:**
```
BENCHMARK_PROFILE=small BENCHMARK_BASELINE_REF=main julia --project=bench bench/compare_benchmarks.jl
```

**Render benchmark reports:**
```
BENCHMARK_INPUT_DIR=bench/results BENCHMARK_OUTPUT_DIR=bench/results julia --project=bench bench/render_report.jl
```

**Run large tiers on SLURM:**
```
bench/slurm_benchmarks.sh submit
```

**Install docs dependencies (first time only):**
```julia
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
```

**Build docs:**
```julia
julia --project=docs docs/make.jl
```

**Build docs without regenerating literate examples:**
```julia
julia --project=docs docs/make.jl --no-literate
```

**Regenerate asynch experiment data:**
```
julia --project=docs docs/scripts/acc26_experiments.jl
```

**Regenerate asynch doc figures from data:**
```
julia --project=docs docs/scripts/generate_figures.jl
```

## Architecture

`CellularSheaves` re-exports `NetworkSheaves`, which aggregates all submodules:

```
src/network_sheaves/
  SheafInterface.jl      # Abstract declarations — new operations go here first
  EuclideanSheaves.jl    # EuclideanSheaf type, coboundary map, Laplacian, solvers
  Morphisms.jl           # SheafMorphism, ComplexMorphism, compose, is_morphism
  Pushforwards.jl        # pushforward_sheaf, fiber operations
  Pushouts.jl            # Sheaf pushouts over graph spans
  GraphHomomorphisms.jl  # GraphHomomorphism type and helpers
  TrajectorySheaf.jl     # Trajectory/controlled trajectory sheaves for LTI systems
  ADT.jl + Parser.jl     # @cellular_sheaf DSL macro
  BlockSparseArrays.jl   # blocksparse helper used by coboundary_map
  PotentialSheaves.jl    # Potential sheaf constructions
```

The primary concrete type is `EuclideanSheaf{T}`, a struct with `vertex_stalks::Vector{Int}`, `edge_stalks::Dict{UnorderedPair{Int},Int}`, `underlying_graph::Graph`, and `restriction_maps::Dict{Pair{Int},Matrix{T}}`. The restriction map from vertex `v1` to edge `(v1, v2)` is keyed `v1 => v2`.

**Adding a new operation:** declare it in `SheafInterface.jl`, implement on `EuclideanSheaf` in the appropriate module, add the export to `NetworkSheaves.jl`.

**Test layout:** one test file per source module in `test/network_sheaves/`, included from `test/runtests.jl`. Do not add `CellularSheaves` to `test/Project.toml` — this breaks CI.

## Benchmark architecture

The benchmark suite lives entirely in `bench/src/`:

```
bench/src/
  BenchmarkSuite.jl    # build_suite(): defines all @benchmarkable leaves
  BenchmarkShards.jl   # SHARD_ORDER, SHARD_MANIFEST, shard/profile selection logic
  BenchmarkReports.jl  # report rendering; run/compare/render entry points
  CellularSheavesBenchmarks.jl  # module aggregator
```

Key env vars read by `run_benchmarks_from_env!()`:

| Variable | Default | Meaning |
|---|---|---|
| `BENCHMARK_PROFILE` | `small` | Size tier: `small`, `large`, or `full` |
| `BENCHMARK_SHARD` | `all` | Shard name or `all` to run the whole profile |
| `BENCHMARK_RESULT_DIR` | `bench/results` | Output directory for artifacts |
| `BENCHMARK_SECONDS` | `5` | Min seconds per leaf |
| `BENCHMARK_SAMPLES` | `25` | Max samples per leaf |

**Adding a new benchmark and shard:**
1. Add `@benchmarkable` leaves to `build_suite()` in `BenchmarkSuite.jl` under a descriptive group key.
2. In `BenchmarkShards.jl`, add a branch to `benchmark_ids()` (or specify IDs inline) and add entries to `SHARD_MANIFEST`.
3. Add the new shard names to `SHARD_ORDER` (small before large for the same category).
4. In `.github/workflows/hpc.yml`, add the new names to the `SHARDS` arrays for the relevant profiles.
5. In `bench/slurm_benchmarks.sh`, add them to `ALL_SHARDS`.



Vertex `v`'s DOF slice in a flat cochain vector `x`:
```julia
offsets = [0; cumsum(vertex_stalks(s))]
x[offsets[v]+1 : offsets[v+1]]
```
Return cochains as `BlockArray(x, vertex_stalks(s))`.

## Solver conventions

Use `CliqueTrees.Multifrontal.ChordalLDLt` (not dense `ldlt`) for all new sparse symmetric PSD linear systems. The factorization convention is `X = P' * L * D * L' * P`:

```julia
using CliqueTrees.Multifrontal
M = ldlt!(ChordalLDLt(A), RowMaximum())
# solve M*v = b:
c = M.P' \ b; z = M.L \ c; w = M.D \ z; y = M.L' \ w; v = M.P \ y
```

`nullspace_ldlt` in `EuclideanSheaves.jl` is the authoritative reference for this pattern. Avoid `Matrix(A)` or `Array(A)` conversions in hot paths — keep matrices sparse throughout.

## Dependency philosophy

No new top-level dependencies without checking `Project.toml` first. Implement graph-combinatorial operations directly on `Graphs.SimpleGraph` (union-find, etc.) rather than pulling in Catlab.jl. Do not add `LDLFactorizations.jl` — `ChordalLDLt` is already present.

## Code conventions

- Input validation at public boundaries via `@argcheck` (ArgCheck.jl, already a dependency).
- No Java-style getter/setter names (`get_foo` / `set_foo`); prefer idiomatic Julia names.
- Do not put line comments inside function definitions in literate examples (`docs/literate/`) — this breaks the Literate.jl pipeline. Place explanatory comments above the function or as end-of-line comments.

## Documentation

- Keep `makedocs(modules=...)` in `docs/make.jl` aligned with all documented modules.
- Add a separate `@autodocs` block per module in `docs/src/api/*.md`.
- When adding a new module, update both `docs/make.jl` and the relevant `docs/src/api/` file in the same change.
- Check for missing cross-references after building (they appear as warnings, not errors).

## Feature requests

Store feature request drafts in `docs/issues/` (numbered markdown files). Each issue must include: **Mathematical Background**, **Codebase Orientation** (table of files + why), **Requested Implementation** (exact signature + docstring + algorithm sketch), **Tests to Write** (concrete `@test` statements), **Verification Checklist**, and **Out of Scope**. Split issues at natural seams — one self-contained function per issue.

## Mathematical context

The codebase targets researchers in sheaf theory, category theory, and applied linear algebra (Hansen, Ghrist, Curry, Riess). Docstrings for public-facing functions should explain the math for users who are not specialists. When a request is mathematically ambiguous, ask for clarification before implementing. Cite relevant papers and use their nomenclature.
