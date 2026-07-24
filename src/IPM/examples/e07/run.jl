include("def.jl")

run()

# =============================================================================
# Sample run: 2026-07-24 (--quick --mosek)  [rhs/rebuild — anchored raug default]
# -----------------------------------------------------------------------------
# IPM stalls at N=512 and N=2048 (converges at N=1024) — instance-dependent exp-cone
# centrality stall: one exp cone is ejected off the central path on the first step and the
# affine method has no centrality-aware step control to recover (H2, see e07_diagnostic_spec.md).
# HSD solves all sizes. No slope (IPM has one finite point).
# N=512   dof=8514   n1=7477  blk=13    IPM —         HSD  190.4ms       Cla  275.6ms       Msk  503.5ms
# N=1024  dof=17028  n1=14957 blk=13    IPM  517.0ms  HSD  414.4ms (0.80x)  Cla  886.8ms (1.72x)  Msk 1185.9ms (2.29x)
# N=2048  dof=34068  n1=29925 blk=13    IPM —         HSD  869.2ms       Cla 1498.5ms       Msk 3649.0ms
# =============================================================================
