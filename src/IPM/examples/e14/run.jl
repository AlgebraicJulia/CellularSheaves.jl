include("def.jl")

run()

# =============================================================================
# Sample run: 2026-07-25 (--quick --mosek). IPM = affine primal-dual, HSD =
# homogeneous self-dual (both the sheaf solver). Mosek benchmarked PRIMAL
# (measured faster than dual at every point on both dials).
# -----------------------------------------------------------------------------
# K-dial (grid size K, p=3, T=6):
# K=8    dof=2307   n1=864    IPM 11.7ms  HSD  9.8ms  Cla 13.1ms  Msk  34.0ms
# K=12   dof=4995   n1=1872   IPM 26.0ms  HSD 22.7ms  Cla 32.6ms  Msk  67.2ms
# K=16   dof=8707   n1=3264   IPM 45.9ms  HSD 46.3ms  Cla 64.6ms  Msk 103.8ms
# Slopes: IPM DOF^1.03  HSD DOF^1.17  Clarabel DOF^1.20  Mosek DOF^0.84
#
# T-dial (block fatness T at fixed K=12 — the fat-block kernel experiment):
# T=6    dof=4995   n1=1872   IPM  25.5ms  HSD  22.7ms  Cla  32.6ms  Msk   67.5ms
# T=12   dof=8739   n1=3744   IPM  51.3ms  HSD  42.4ms  Cla  56.0ms  Msk  106.6ms
# T=24   dof=16227  n1=7488   IPM 112.8ms  HSD 106.9ms  Cla 162.2ms  Msk  134.9ms
# T=48   dof=31203  n1=14976  IPM 418.6ms  HSD 334.5ms  Cla 442.3ms  Msk 1124.0ms
# T-slopes: IPM T^1.32  HSD T^1.30  Clarabel T^1.28  Mosek T^1.25
#
# Dense-Q verdict: both sheaf solvers beat Clarabel at every K/T (HSD fastest),
# and Mosek is slowest. The dense GP precisions W_e fold into the block Hessian
# for free while taxing the conic solvers, and the harmonic stalk is a dense
# COLUMN of B — eliminated last, staying sparse.
#   T-dial (registered prediction — margins widen vs both baselines as blocks
# fatten): CONFIRMED at the fat end (T=48: HSD leads Mosek 3.36x, Clarabel 1.32x)
# but NON-MONOTONE — both baselines close in at T=24 (Mosek to just 1.26x) before
# the margin reopens. That mid-range flat spot makes the printed T-slopes
# misleading (they under-report Mosek's true growth: the T=24→48 segment alone is
# HSD x3.1 vs Mosek x8.3 for a 2x T step). Read the endpoint RATIOS, not the
# slope. An earlier 3-point extrapolation (T<=24) wrongly predicted a Mosek
# crossover; T=48 refuted it — hence the sweep runs to 48.
#
# Oracle reference (tools/ocean_gyre_oracle.py, executed 2026-07-24, K=8):
#   dim H1 = 2; rank[d0|d1'|eta] = 144 = m (coker d0 alone: 64)
#   gyre LS err 1.309e-1 (floor 0.125); island-blind ablation 5.06e-2 →
#   2.72e-1 with 62.1% residual on top-10% harmonic edges; screening margin
#   14.1×; support recovery 6/6 exact; debiased ≡ mask-oracle (0.0e0)
# =============================================================================
