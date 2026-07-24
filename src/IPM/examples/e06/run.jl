include("def.jl")

run()

# =============================================================================
# Sample run: 2026-07-24 (--quick --mosek)  [rhs/rebuild — anchored raug default]
# -----------------------------------------------------------------------------
# P=32    dof=864    n1=125   blk=27    IPM    3.1ms  HSD    3.6ms (1.16x)  Cla    3.0ms (0.95x)  Msk    7.0ms (2.26x)
# P=64    dof=1728   n1=253   blk=27    IPM    6.9ms  HSD    7.4ms (1.07x)  Cla    6.2ms (0.89x)  Msk   20.3ms (2.94x)
# P=128   dof=3456   n1=509   blk=27    IPM   12.2ms  HSD   14.1ms (1.16x)  Cla   11.8ms (0.97x)  Msk   41.1ms (3.37x)
# IPM: DOF^0.99  HSD: DOF^0.99  Clarabel: DOF^1.00  Mosek: DOF^1.27
# =============================================================================
