include("def.jl")

run()

# =============================================================================
# Sample run: 2026-07-24 (--quick --mosek). Mosek is benchmarked DUAL-only —
# its best form on the dense-Q tracking objective (dual ~1.4-1.6x over primal).
# Ratios below are vs HSD (the affine IPM stalls at N >= 240 — see def.jl banner).
# -----------------------------------------------------------------------------
# N=240   dof=3603    IPM —        HSD   28.3ms  Cla   33.1ms (1.17x)  MskD   53.7ms (1.90x)
# N=480   dof=7203    IPM —        HSD   59.0ms  Cla   73.6ms (1.25x)  MskD   88.0ms (1.49x)
# N=960   dof=14403   IPM —        HSD  116.7ms  Cla  157.5ms (1.35x)  MskD  454.5ms (3.89x)
# Slopes: HSD DOF^1.02  Clarabel DOF^1.13  Mosek(dual) DOF^1.54
#
# Dense-Q verdict: HSD < Clarabel < Mosek-dual at every N, and the gap GROWS
# with N (Mosek's slope is worst — 1.9x at N=240 to 3.9x at N=960). The
# quadratic tracking Q folds cheaply into the block Hessian for the sheaf
# solver but is a reformulation tax for the conic solvers. This INVERTS the
# Q ≡ 0 ranking, where Mosek led. The affine IPM stalls past N=240 (a
# near-degenerate-face limit cycle; mechanism open, candidates refuted by
# experiment — see def.jl banner). HSD is the solver at scale here.
# =============================================================================
