include("def.jl")

run()

# =============================================================================
# Sample run: 2026-07-25 (--quick --mosek). Mosek is benchmarked PRIMAL —
# its best form on the dense-Q tracking objective (primal ~2x faster than the
# dualized problem, and far better-behaved at scale; the earlier dual-only
# policy was an un-remeasured guess and is retired).
# Ratios below are vs HSD (the affine IPM stalls at N >= 240 — see def.jl banner).
# -----------------------------------------------------------------------------
# N=240   dof=3603    IPM —        HSD   27.9ms  Cla   32.8ms (1.18x)  Msk   47.5ms (1.70x)
# N=480   dof=7203    IPM —        HSD   57.5ms  Cla   73.3ms (1.27x)  Msk  100.1ms (1.74x)
# N=960   dof=14403   IPM —        HSD  114.4ms  Cla  157.8ms (1.38x)  Msk  197.1ms (1.72x)
# Slopes: HSD DOF^1.02  Clarabel DOF^1.13  Mosek DOF^1.03
#
# Dense-Q verdict: HSD < Clarabel < Mosek at every N, gap steady (~1.7x for
# Mosek). The quadratic tracking Q folds cheaply into the block Hessian for the
# sheaf solver but is a reformulation tax for the conic solvers. This INVERTS the
# Q ≡ 0 ranking, where Mosek led. (Dualizing Mosek was tried and rejected — it
# was competitive at N=240 but blew up to 454ms by N=960, slope 1.54.) The affine
# IPM stalls past N=240 (a near-degenerate-face limit cycle; mechanism open,
# candidates refuted by experiment — see def.jl banner). HSD is the solver at scale here.
# =============================================================================
