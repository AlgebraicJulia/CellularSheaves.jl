include("def.jl")

run()

# =============================================================================
# Sample run: 2026-07-24 (--quick --mosek)  [rhs/rebuild — anchored raug default]
# -----------------------------------------------------------------------------
# M=4    dof=521    n1=261    blk=64    IPM    9.1ms  HSD   10.8ms (1.18x)  Cla   41.6ms (4.56x)  Msk   38.9ms (4.26x)
# M=8    dof=1041   n1=521    blk=64    IPM   38.9ms  HSD   48.7ms (1.25x)  Cla  286.6ms (7.37x)  Msk  115.1ms (2.96x)
# M=16   dof=2081   n1=1041   blk=64    IPM  218.6ms  HSD  244.7ms (1.12x)  Cla 1979.5ms (9.05x)  Msk  419.1ms (1.92x)
# IPM: DOF^2.29  HSD: DOF^2.26  Clarabel: DOF^2.79  Mosek: DOF^1.72  (elim = METIS)
# =============================================================================
