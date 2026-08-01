#!/usr/bin/env bash

set -euo pipefail

ALL_SHARDS=("assembly-small" "solver-small" "extension-small" "assembly-large" "solver-large" "extension-large")
MODE="${1:-submit}"
RESULT_DIR="${BENCHMARK_RESULT_DIR:-bench/results/slurm}"
PROFILE="${BENCHMARK_PROFILE:-full}"
EXPECTED_SHARDS="$(IFS=,; echo "${ALL_SHARDS[*]}")"

run_shard() {
  local shard="$1"
  export BENCHMARK_PROFILE="${PROFILE}"
  export BENCHMARK_SHARD="${shard}"
  export BENCHMARK_RESULT_DIR="${RESULT_DIR}"
  export BENCHMARK_RENDER_REPORT="false"
  julia --project=bench bench/run_benchmarks.jl
}

case "${MODE}" in
  submit)
    count="${#ALL_SHARDS[@]}"
    shard_job_id="$(
      sbatch --parsable \
        --array="1-${count}" \
        --export=ALL,BENCHMARK_PROFILE="${PROFILE}",BENCHMARK_RESULT_DIR="${RESULT_DIR}" \
        "$0" run-array
    )"
    sbatch \
      --dependency="afterok:${shard_job_id}" \
      --export=ALL,BENCHMARK_PROFILE="${PROFILE}",BENCHMARK_RESULT_DIR="${RESULT_DIR}",BENCHMARK_REPORT_MODE=summary \
      "$0" aggregate
    ;;
  run-array)
    : "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required for run-array mode}"
    shard="${ALL_SHARDS[$((SLURM_ARRAY_TASK_ID - 1))]}"
    run_shard "${shard}"
    ;;
  aggregate)
    export BENCHMARK_INPUT_DIR="${RESULT_DIR}"
    export BENCHMARK_OUTPUT_DIR="${RESULT_DIR}"
    export BENCHMARK_EXPECTED_SHARDS="${BENCHMARK_EXPECTED_SHARDS:-${EXPECTED_SHARDS}}"
    julia --project=bench bench/render_report.jl
    ;;
  *)
    echo "Usage: $0 [submit|run-array|aggregate]" >&2
    exit 1
    ;;
esac
