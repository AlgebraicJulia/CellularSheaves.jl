include("def.jl")

run()

# =============================================================================
# Sample run (--quick --mosek, M3 Max, raug = 1e4, METIS ordering):
#
#   Gate tests (M = 8, Lh = 1024, Lf = 64, ε_LS = 0.7314):
#     [PASS] dense by nature: Θ 1.0, G_m ≥ 1.0, L_m' ≥ 1.0, joint G 1.0
#     [PASS] polyphase exclusion (D = 4 control): off-class Gram EXACTLY 0
#     [PASS] MINT threshold: ε_LS 0.009 → 5.2e-14 at Lf* = 40
#     [PASS] whitening: joint 0.0, split 2.0e-16
#     [PASS] TRS: rel Δv 4.7e-9, ν = 1209.65
#     [PASS] certificate: ball use 1.0 (exact form)
#     [PASS] conservatism: 1.664 → 1.126 over εx = 0.05 → 0.0125
#     [PASS] noise shaping: move 12.2%, excess ×1.51
#
#   M=4    dof=521    IPM  11ms  Cla   42ms (3.62x)  Msk   40ms (3.51x)
#   M=8    dof=1041   IPM  52ms  Cla  286ms (5.48x)  Msk  116ms (2.22x)
#   M=16   dof=2081   IPM 239ms  Cla 1987ms (8.30x)  Msk  428ms (1.79x)
#
#   Slopes: IPM DOF^2.19, Clarabel DOF^2.79, Mosek DOF^1.71
#
#   READING: METIS ordering drops the IPM slope from 2.30 → 2.19 (vs AMF) —
#   better elimination order for the star-over-simplex topology. IPM wins vs
#   both baselines at all M. The M² coupling mass is still visible in the
#   slope; Mosek's 1.71 is better because its KKT doesn't pay the triangular
#   slice hyperedges. The --budget form (O(M) coupling) should recover the
#   slope at the price of its measured conservatism (~1.2x at εx = 0.0125).
# =============================================================================
