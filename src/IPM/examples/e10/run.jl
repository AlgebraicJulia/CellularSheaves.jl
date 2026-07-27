include("def.jl")

run()

# =============================================================================
# Sample run: 2026-07-25 (--quick --mosek)  [Mosek benchmarked DUAL — ~1.4x faster than primal]
# -----------------------------------------------------------------------------
# Mosek is dualized (Dualization.jl): dual ~1.4x faster than primal and, on the
# fatness dial, flips Mosek from a tie with the affine IPM (primal was 1.00x at d=8)
# to a win (0.73x). The sheaf solver still leads the K-dial at every K.
# K=4   dof=1312   n1=672    blk=36    IPM   38.4ms  HSD   46.0ms (1.20x)  Cla   94.1ms (2.45x)  Msk   57.2ms (1.49x)
# K=6   dof=2952   n1=1512   blk=36    IPM  118.1ms  HSD  142.5ms (1.21x)  Cla  339.5ms (2.87x)  Msk  146.8ms (1.24x)
# K=8   dof=5248   n1=2688   blk=36    IPM  260.1ms  HSD  311.5ms (1.20x)  Cla  763.9ms (2.94x)  Msk  280.7ms (1.08x)
# IPM: DOF^1.38  HSD: DOF^1.38  Clarabel: DOF^1.52  Mosek: DOF^1.15  (K-dial)
# d=4   dof=1312   n1=672    blk=36    IPM   38.2ms  HSD   45.9ms (1.20x)  Cla   93.8ms (2.46x)  Msk   57.2ms (1.50x)
# d=6   dof=2832   n1=1488   blk=78    IPM  221.3ms  HSD  262.9ms (1.19x)  Cla  600.0ms (2.71x)  Msk  177.0ms (0.80x)
# d=8   dof=4928   n1=2624   blk=136   IPM  669.5ms  HSD  796.7ms (1.19x)  Cla 2544.6ms (3.80x)  Msk  488.7ms (0.73x)
# =============================================================================
