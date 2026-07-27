include("def.jl")

run()

# =============================================================================
# Sample run: 2026-07-25 (--quick --mosek)  [Mosek benchmarked DUAL — ~3.7x faster than primal]
# -----------------------------------------------------------------------------
# IPM stalls at N=512 and N=2048 (converges at N=1024) — instance-dependent exp-cone
# centrality stall: one exp cone is ejected off the central path on the first step and the
# affine method has no centrality-aware step control to recover (H2, see e07_diagnostic_spec.md).
# HSD solves all sizes. No slope (IPM has one finite point). Mosek is dualized here
# (Dualization.jl): dual is ~3.7x faster than primal and edges past HSD at N=2048.
# N=512   dof=8514   n1=7477  blk=13    IPM —         HSD  188.7ms       Cla  276.0ms       Msk  118.0ms
# N=1024  dof=17028  n1=14957 blk=13    IPM  512.2ms  HSD  411.9ms (0.80x)  Cla  891.3ms (1.74x)  Msk  319.0ms (0.62x)
# N=2048  dof=34068  n1=29925 blk=13    IPM —         HSD  852.4ms       Cla 1478.7ms       Msk  684.2ms
# slopes: HSD DOF^1.09  Clarabel DOF^1.21  Mosek DOF^1.27
# =============================================================================
