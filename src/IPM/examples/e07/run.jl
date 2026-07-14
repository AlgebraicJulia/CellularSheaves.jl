include("def.jl")

run()

# =============================================================================
# Sample run (julia --project e8.jl --quick --mosek):
#
#   Gate tests (N = 512):
#     [PASS] counts: med 3.0, max 41.0, 106 zero bins (linear-only), 406 exp cones
#     [PASS] compressed MLE: IPM == Clarabel (rel 2.4e-9); nodal MLEM sits 0.012 below (compression cost)
#     [PASS] objective identity: F(u) = -5560.2241 (rel 5.7e-8)
#     [PASS] model value: MSE Poisson-TV 0.0083 < 0.75 × best weighted-Gauss-TV 0.0178
#     [PASS] positivity diagnostic (compressed space): min u = -0.08 without the cone, MSE ratio 1.0x; the dipole blow-up is nodal-only (oracle [5])
#     [PASS] IPM vs Clarabel (same conic program): ‖Δu‖∞ = 0.0002
#
#   N=512   dof=8514   n1=7477  blk=13    IPM  323.6ms  Cla  297.7ms (0.92x)  Msk  523.9ms (1.62x)
#   N=1024  dof=17028  n1=14957 blk=13    IPM  570.2ms  Cla  954.4ms (1.67x)  Msk 1273.0ms (2.23x)
#   N=2048  dof=34068  n1=29925 blk=13    IPM 1183.5ms  Cla 1602.2ms (1.35x)  Msk 4188.6ms (3.54x)
#
#   Fitted slopes: IPM DOF^0.94, Clarabel DOF^1.21, Mosek DOF^1.50
# =============================================================================
