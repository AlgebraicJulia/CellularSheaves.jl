include("def.jl")

run()

# =============================================================================
# Sample run: 2026-07-24 (--quick)  [rhs/rebuild — anchored raug default]
# -----------------------------------------------------------------------------
# T=64    dof=2040   n1=1024  blk=3     IPM   56.0ms  HSD   41.7ms (0.74x)  Cla   30.2ms (0.54x)  Msk —        (—)
# T=128   dof=4088   n1=2048  blk=3     IPM  106.4ms  HSD   82.0ms (0.77x)  Cla   81.8ms (0.77x)  Msk —        (—)
# T=256   dof=8184   n1=4096  blk=3     IPM  225.6ms  HSD  163.8ms (0.73x)  Cla  275.4ms (1.22x)  Msk —        (—)
# IPM: DOF^1.00  HSD: DOF^0.98  Clarabel: DOF^1.59  Mosek: n/a
# (ALMOST_OPTIMAL at 1e-8 — accuracy-ceiling cell; times finite and on the saved line.)
# =============================================================================
