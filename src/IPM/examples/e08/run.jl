include("def.jl")

run()

# =============================================================================
# Sample run (julia --project e9.jl --quick --mosek):
#
#   Gate tests (n = 8, T = 64):
#     [PASS] analytic (Almgren power-law ODE, E = 0.00448774): sup|x_dt − x(t)| = 0.00037 at N = 64
#     [PASS] objective identity: F(x) = 233.9868 (rel 5.0e-8)
#     [PASS] IPM vs Clarabel (same conic program): ‖Δx‖∞ = 5.3e-5
#
#   T=64    dof=2040   n1=1024  blk=3     IPM  108.3ms  Cla   30.7ms (0.28x)  Msk   43.1ms (0.40x)
#   T=128   dof=4088   n1=2048  blk=3     IPM  150.2ms  Cla   83.7ms (0.56x)  Msk  123.2ms (0.82x)
#   T=256   dof=8184   n1=4096  blk=3     IPM  498.0ms  Cla  283.9ms (0.57x)  Msk  630.4ms (1.27x)
#
#   Fitted slopes: IPM DOF^1.10, Clarabel DOF^1.60, Mosek DOF^1.93
# =============================================================================
