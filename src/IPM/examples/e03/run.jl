include("def.jl")

run()

# =============================================================================
# Sample run: 2026-07-24 (--quick --mosek)  [rhs/rebuild — anchored raug default]
# -----------------------------------------------------------------------------
# M=2    dof=100    n1=60    blk=25    IPM    0.5ms  HSD    0.7ms (1.24x)  Cla    1.5ms (2.74x)  Msk    1.7ms (3.08x)
# M=4    dof=400    n1=360   blk=25    IPM    3.7ms  HSD    4.5ms (1.21x)  Cla   11.3ms (3.03x)  Msk   14.7ms (3.94x)
# M=8    dof=1600   n1=1680  blk=25    IPM   25.1ms  HSD   28.0ms (1.12x)  Cla   91.7ms (3.66x)  Msk   69.8ms (2.78x)
# IPM: DOF^1.38  HSD: DOF^1.35  Clarabel: DOF^1.49  Mosek: DOF^1.35
# =============================================================================
