include("def.jl")

run()

# =============================================================================
# Sample run: 2026-07-24 (--quick --mosek)  [rhs/rebuild — anchored raug default]
# -----------------------------------------------------------------------------
# K=4   dof=1312   n1=672    blk=36    IPM   38.9ms  HSD   46.5ms (1.19x)  Cla   92.1ms (2.37x)  Msk   50.9ms (1.31x)
# K=6   dof=2952   n1=1512   blk=36    IPM  120.0ms  HSD  145.7ms (1.21x)  Cla  339.2ms (2.83x)  Msk  142.7ms (1.19x)
# K=8   dof=5248   n1=2688   blk=36    IPM  260.2ms  HSD  316.4ms (1.22x)  Cla  700.6ms (2.69x)  Msk  288.5ms (1.11x)
# IPM: DOF^1.37  HSD: DOF^1.39  Clarabel: DOF^1.47  Mosek: DOF^1.25  (K-dial)
# d=4   dof=1312   n1=672    blk=36    IPM   38.7ms  HSD   46.0ms (1.19x)  Cla   93.2ms (2.41x)  Msk   50.9ms (1.32x)
# d=6   dof=2832   n1=1488   blk=78    IPM  223.6ms  HSD  267.4ms (1.20x)  Cla  544.0ms (2.43x)  Msk  245.3ms (1.10x)
# d=8   dof=4928   n1=2624   blk=136   IPM  672.4ms  HSD  812.1ms (1.21x)  Cla 2333.3ms (3.47x)  Msk  674.9ms (1.00x)
# =============================================================================
