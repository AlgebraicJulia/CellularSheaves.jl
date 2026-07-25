include("def.jl")

run()

# =============================================================================
# Sample run: 2026-07-25 (--quick --mosek). IPM = affine primal-dual, HSD =
# homogeneous self-dual (both the sheaf solver). Mosek benchmarked DUAL-ONLY.
# -----------------------------------------------------------------------------
# K=8    dof=2307   n1=864    IPM 11.5ms  HSD  9.5ms  Cla 13.1ms  Msk 30.8ms
# K=12   dof=4995   n1=1872   IPM 25.9ms  HSD 22.2ms  Cla 33.1ms  Msk 64.1ms
# K=16   dof=8707   n1=3264   IPM 45.8ms  HSD 46.2ms  Cla 64.7ms  Msk 93.5ms
# Slopes: IPM DOF^1.04  HSD DOF^1.18  Clarabel DOF^1.20  Mosek DOF^0.84
#
# Dense-Q verdict: both sheaf solvers beat Clarabel at every K (HSD fastest at
# small K), and Mosek-dual is slowest. The dense GP precisions W_e fold into the
# block Hessian for free while taxing the conic solvers, and the harmonic stalk
# is a dense COLUMN of B — eliminated last, staying sparse. (Non-quick reaches
# K=24, where IPM ~132ms leads Clarabel ~231ms, Mosek-dual ~206ms.)
#
# Oracle reference (tools/ocean_gyre_oracle.py, executed 2026-07-24, K=8):
#   dim H1 = 2; rank[d0|d1'|eta] = 144 = m (coker d0 alone: 64)
#   gyre LS err 1.309e-1 (floor 0.125); island-blind ablation 5.06e-2 →
#   2.72e-1 with 62.1% residual on top-10% harmonic edges; screening margin
#   14.1×; support recovery 6/6 exact; debiased ≡ mask-oracle (0.0e0)
# =============================================================================
