include("def.jl")

run()

# =============================================================================
# Sample runs (MacBook Pro M-series, Clarabel + Mosek):
# -----------------------------------------------------------------------------
#
# EXACT FORM (default):
#   julia --project e11/run.jl --quick  (2026-07-16)
#
#   M=4    dof=521    blk=64   IPM   9.8ms  HSD  13.4ms (1.37x)  Cla   42.1ms (4.30x)  Msk —
#   M=8    dof=1041   blk=64   IPM  48.0ms  HSD  46.8ms (0.97x)  Cla  291.5ms (6.07x)  Msk —
#   M=16   dof=2081   blk=64   IPM 232.5ms  HSD 298.8ms (1.29x)  Cla 1981.7ms (8.52x)  Msk —
#
#   julia --project e11/run.jl --mosek --quick  (2026-07-14)
#
#   M=4    dof=521    blk=64   IPM  11.3ms  HSD  13.8ms (1.22x)  Cla   41.2ms (3.65x)  Msk   39.5ms (3.49x)
#   M=8    dof=1041   blk=64   IPM  49.3ms  HSD  53.4ms (1.08x)  Cla  284.2ms (5.77x)  Msk  117.3ms (2.38x)
#   M=16   dof=2081   blk=64   IPM 244.9ms  HSD 269.9ms (1.10x)  Cla 1952.5ms (7.97x)  Msk  437.1ms (1.78x)
#
#   Slopes: IPM DOF^2.22, HSD DOF^2.15, Clarabel DOF^2.79, Mosek DOF^1.74
#
# BUDGET FORM (--budget, star hub):
#   julia --project e11/run.jl --mosek --budget
#
#   M=4    dof=520    blk=65   IPM  10.2ms  HSD   8.7ms (0.84x)  Cla   11.6ms (1.13x)  Msk   17.1ms (1.66x)
#   M=8    dof=1040   blk=65   IPM  31.8ms  HSD  25.4ms (0.80x)  Cla   24.1ms (0.76x)  Msk   33.0ms (1.04x)
#   M=16   dof=2080   blk=65   IPM 142.0ms  HSD 106.5ms (0.75x)  Cla   52.3ms (0.37x)  Msk   65.8ms (0.46x)
#   M=32   dof=4160   blk=65   IPM 585.3ms  HSD 407.4ms (0.70x)  Cla   98.1ms (0.17x)  Msk  125.6ms (0.21x)
#
#   Slopes: IPM DOF^1.97, HSD DOF^1.87, Clarabel DOF^1.04, Mosek DOF^0.96
#   NOTE: IPM LOSES 6x to Clarabel at M=32 — the star hub bottleneck.
#
# TREE FORM (--tree, binary hypot combiners):
#   julia --project e11/run.jl --mosek --tree
#
#   M=4    dof=529    blk=64   IPM   4.6ms  HSD   4.7ms (1.01x)  Cla    6.6ms (1.43x)  Msk   12.2ms (2.64x)
#   M=8    dof=1061   blk=64   IPM   7.4ms  HSD   9.1ms (1.24x)  Cla   15.9ms (2.16x)  Msk   25.1ms (3.40x)
#   M=16   dof=2125   blk=64   IPM  14.9ms  HSD  22.3ms (1.50x)  Cla   32.3ms (2.17x)  Msk   54.7ms (3.67x)
#   M=32   dof=4253   blk=64   IPM  36.5ms  HSD  40.5ms (1.11x)  Cla   64.8ms (1.77x)  Msk   93.3ms (2.56x)
#
#   Slopes: IPM DOF^0.99, HSD DOF^1.06, Clarabel DOF^1.09, Mosek DOF^0.99
#   NOTE: Tree fix restores the win — 1.77x vs Clarabel at M=32 (was 0.17x).
#
# REGIME RULE: Flagship (exact, zero conservatism) for M <= 16; --tree
# (conservatism 1.04) beyond. The budget-star is dominated: same
# conservatism as tree, worse topology.
# =============================================================================
