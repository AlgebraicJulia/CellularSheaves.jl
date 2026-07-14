include("def.jl")

run()

# =============================================================================
# Sample run (raug = 1e4, --quick --mosek, 2026-07-13):
#
#   Gate tests (n = 12, m = 4, T = 100):
#   [PASS] dense by nature: density@1e-8 A 0.986 W 1.0 V 1.0 M 1.0
#   [PASS] analytic references: RTS recursion vs information form 1.5e-14
#   [PASS] objective identity: F(x) = 2230.0244 (rel 2.7e-9)
#   [PASS] deformation to RTS: ‖x(δ) − RTS‖∞ 0.055 → 0.0036 over δ = 4 → 16 (rate δ^-1.97)
#   [PASS] robustness: clean-region RMSE ×1 robust 0.2642 RTS 0.2567; ×32 robust 0.3071 RTS 0.9513
#   [PASS] support recovery: top-24 residuals are exactly the corrupted steps (margin ×3.0; influence cap 61.0×)
#   [PASS] IPM vs Clarabel (same conic program): ‖Δx‖∞ = 0.00022
#
#   T=64    dof=1926   n1=1081  blk=6   IPM   16.7ms  Cla   14.0ms (0.84x)  Msk   41.9ms (2.51x)
#   T=128   dof=3846   n1=2169  blk=6   IPM   34.8ms  Cla   28.2ms (0.81x)  Msk   77.1ms (2.21x)
#   T=256   dof=7686   n1=4345  blk=6   IPM   78.2ms  Cla   62.2ms (0.80x)  Msk  166.6ms (2.13x)
#
#   Fitted log-log slopes (time vs DOF):
#     IPM: DOF^1.12  Clarabel: DOF^1.08  Mosek: DOF^1.00
# =============================================================================
