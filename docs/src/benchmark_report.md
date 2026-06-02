# Benchmark report

This placeholder page is kept in the repository so the docs can link to a benchmark report without running the benchmark suite during the docs build.

Refresh it by generating benchmark artifacts first and then copying the rendered Markdown here:

```julia
julia --project=bench bench/run_benchmarks.jl
BENCHMARK_INPUT_DIR=bench/results BENCHMARK_OUTPUT_DIR=bench/results BENCHMARK_DOCS_MARKDOWN=docs/src/benchmark_report.md julia --project=bench bench/render_report.jl
```
