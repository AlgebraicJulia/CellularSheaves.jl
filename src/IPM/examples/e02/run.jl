include("def.jl")

run()

# =============================================================================
# Sample run: 2026-07-14 (--quick --mosek, M3 Max)
# -----------------------------------------------------------------------------
#   Gate tests:
#     [PASS] ED oracle N=6: E₀=-2.8028 (true ring energy; fixed 2026-07-14)
#     [PASS] Lower bound: E_sdp=-3.0 ≤ E_ed=-2.8028 (gap=0.197)
#
#   N=8    dof=1088   IPM  27ms  Cla  145ms (5.37x)  Msk  56ms (2.09x)
#   N=10   dof=1360   IPM  32ms  Cla  187ms (5.84x)  Msk  61ms (1.91x)
#   N=12   dof=1632   IPM  40ms  Cla  246ms (6.22x)  Msk  80ms (2.02x)
#
#   Slopes: IPM DOF^0.94, Clarabel DOF^1.30, Mosek DOF^0.85
# =============================================================================
