include("def.jl")

run()

# =============================================================================
# Sample run: 2026-07-24 (--quick --mosek)  [rhs/rebuild — anchored raug default]
# -----------------------------------------------------------------------------
# T=64    dof=2040   n1=1024  blk=3     IPM   56.2ms  HSD   41.9ms (0.75x)  Cla   30.4ms (0.54x)  Msk   42.3ms (0.75x)
# T=128   dof=4088   n1=2048  blk=3     IPM  107.1ms  HSD   81.7ms (0.76x)  Cla   81.8ms (0.76x)  Msk  121.2ms (1.13x)
# T=256   dof=8184   n1=4096  blk=3     IPM  225.0ms  HSD  163.5ms (0.73x)  Cla  275.6ms (1.23x)  Msk  604.7ms (2.69x)
# IPM: DOF^1.00  HSD: DOF^0.98  Clarabel: DOF^1.59  Mosek: DOF^1.92
# (ALMOST_OPTIMAL at 1e-8 — accuracy-ceiling cell; times finite and on the saved line.)
# =============================================================================
