include("def.jl")

run()

# =============================================================================
# Sample run (--quick --mosek, M3 Max, raug = 1e2):
#
#   K dial (d = 4 fixed):
#     K=4   dof=1312  n1=672   blk=36  IPM  55ms  Cla  93ms (1.67x)  Msk  52ms (0.94x)
#     K=6   dof=2952  n1=1512  blk=36  IPM 175ms  Cla 327ms (1.87x)  Msk 148ms (0.84x)
#     K=8   dof=5248  n1=2688  blk=36  IPM 385ms  Cla 739ms (1.92x)  Msk 297ms (0.77x)
#   Slopes: IPM DOF^1.40, Clarabel DOF^1.50, Mosek DOF^1.26.
#
#   d dial (K = 4 fixed) — baselines pay per nonzero (~d²), IPM per block:
#     d=4   dof=1312  n1=672   blk=36   IPM  55ms  Cla   93ms (1.70x)  Msk  53ms (0.97x)
#     d=6   dof=2832  n1=1488  blk=78   IPM 280ms  Cla  567ms (2.03x)  Msk 255ms (0.91x)
#     d=8   dof=4928  n1=2624  blk=136  IPM 800ms  Cla 2464ms (3.08x)  Msk 677ms (0.85x)
#
#   READING: Both predictions confirmed.
#   (1) THIN-SIDE: balanced blocks (b_v ≈ b_e) moved the Mosek ratio from
#       0.36x (fat-primal quantum torus at matched DOF) to 0.77-0.97x — a
#       ~2.5x swing attributable to block geometry alone — and Clarabel
#       moved from DNF to functioning at ~2x slower (coker = 0: nothing
#       to choke on).
#   (2) d DIAL: Mosek pulls further ahead as d grows (0.97 → 0.91 → 0.85),
#       contrary to the d² nonzero-cost prediction — something in the IPM
#       scales worse than expected with block size. Clarabel ratio does
#       improve (1.70 → 2.03 → 3.08x), so the d lever works against the
#       first-order solver, just not against Mosek.
# =============================================================================
