include("def.jl")

run()

# =============================================================================
# Sample run: 2026-07-24 (--quick)  [rhs/rebuild — anchored raug default]
# -----------------------------------------------------------------------------
# K=4   dof=1312   n1=672    blk=36    IPM   38.3ms  HSD   45.6ms (1.19x)  Cla   93.6ms (2.44x)  Msk —        (—)
# K=6   dof=2952   n1=1512   blk=36    IPM  117.5ms  HSD  143.8ms (1.22x)  Cla  335.2ms (2.85x)  Msk —        (—)
# K=8   dof=5248   n1=2688   blk=36    IPM  256.6ms  HSD  311.2ms (1.21x)  Cla  674.5ms (2.63x)  Msk —        (—)
# IPM: DOF^1.37  HSD: DOF^1.39  Clarabel: DOF^1.43  Mosek: n/a  (K-dial)
# d=4   dof=1312   n1=672    blk=36    IPM   38.6ms  HSD   46.2ms (1.20x)  Cla   94.1ms (2.44x)  Msk —        (—)
# d=6   dof=2832   n1=1488   blk=78    IPM  221.4ms  HSD  264.9ms (1.20x)  Cla  595.2ms (2.69x)  Msk —        (—)
# d=8   dof=4928   n1=2624   blk=136   IPM  665.1ms  HSD  790.1ms (1.19x)  Cla 2621.6ms (3.94x)  Msk —        (—)
# =============================================================================
