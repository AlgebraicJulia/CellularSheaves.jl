include("def.jl")

run()

# =============================================================================
# Sample run: 2026-07-22 (--quick)  [commit 53a237e9]
# -----------------------------------------------------------------------------
# Did not complete: MLE gate (test_mlem_vs_ipm, N=512) aborts with
# AssertionError "MLE solve: NUMERICAL_FAILURE" before the timing table.
# (Poisson-TV exp-cone MLE, known-fragile — not a regression from this refactor.)
# =============================================================================
