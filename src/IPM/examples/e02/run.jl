include("def.jl")

run()

# =============================================================================
# Sample run: 2026-07-25 (--quick --mosek)  [Mosek benchmarked DUAL — ~1.15x faster than primal]
# -----------------------------------------------------------------------------
# N=8    dof=1088   n1=289   blk=136   IPM   18.4ms  HSD   25.7ms (1.40x)  Cla  146.3ms (7.94x)  Msk   43.6ms (2.37x)
# N=10   dof=1360   n1=361   blk=136   IPM   23.5ms  HSD   33.4ms (1.42x)  Cla  194.8ms (8.28x)  Msk   54.5ms (2.32x)
# N=12   dof=1632   n1=433   blk=136   IPM   28.6ms  HSD   40.5ms (1.42x)  Cla  245.1ms (8.56x)  Msk   66.3ms (2.32x)
# IPM: DOF^1.09  HSD: DOF^1.12  Clarabel: DOF^1.27  Mosek: DOF^1.03
# =============================================================================
