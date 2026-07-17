# =============================================================================
# e13/run.jl — Distributed sound-zone control (see def.jl for the full
#              story, oracle figures, and first-run checklist)
# =============================================================================

include("def.jl")

run()

# =============================================================================
# Sample runs (MacBook Pro M-series, Clarabel + Mosek):
# -----------------------------------------------------------------------------
#
# FLAGSHIP (default):
#   julia --project e13/run.jl --quick  (2026-07-16)
#
#   Z=4    dof=3040   blk=65   IPM   69ms   HSD   68ms (0.98x)  Cla  204ms (2.96x)  Msk —
#   Z=8    dof=7004   blk=65   IPM  201ms   HSD  198ms (0.99x)  Cla  539ms (2.69x)  Msk —
#   Z=16   dof=14932  blk=65   IPM  470ms   HSD  511ms (1.09x)  Cla 1190ms (2.53x)  Msk —
#
#   julia --project e13/run.jl --mosek  (2026-07-14)
#
#   Z=4    dof=3040   blk=65   IPM   75ms   HSD   79ms (1.05x)  Cla  200ms (2.67x)  Msk  167ms (2.22x)
#   Z=8    dof=7004   blk=65   IPM  214ms   HSD  226ms (1.06x)  Cla  541ms (2.53x)  Msk  494ms (2.31x)
#   Z=16   dof=14932  blk=65   IPM  513ms   HSD  567ms (1.11x)  Cla 1154ms (2.25x)  Msk 1453ms (2.83x)
#   Z=32   dof=30788  blk=65   IPM 1126ms   HSD 1222ms (1.09x)  Cla 2409ms (2.14x)  Msk 6331ms (5.62x)
#
#   Slopes: IPM DOF^1.17, HSD DOF^1.19, Clarabel DOF^1.07, Mosek DOF^1.55
#
# JOINT ABLATION (--joint, monolithic per-ball SOCs dim ~194):
#   julia --project e13/run.jl --mosek --joint
#
#   Z=4    dof=2988   blk=65   IPM   90ms   HSD  103ms (1.15x)  Cla  161ms (1.80x)  Msk  141ms (1.58x)
#   Z=8    dof=6880   blk=65   IPM  266ms   HSD  309ms (1.16x)  Cla  436ms (1.64x)  Msk  336ms (1.26x)
#   Z=16   dof=14664  blk=65   IPM  680ms   HSD  768ms (1.13x)  Cla  958ms (1.41x)  Msk  987ms (1.45x)
#   Z=32   dof=30232  blk=65   IPM 1498ms   HSD 2005ms (1.34x)  Cla 1999ms (1.33x)  Msk 4281ms (2.86x)
#
#   Slopes: IPM DOF^1.22, HSD DOF^1.28, Clarabel DOF^1.08, Mosek DOF^1.46
#   NOTE: IPM ratio degrades 2.14x -> 1.33x (balanced-blocks lesson at the constant).
#
# NO-BUDGET CONTROL (--nobudget, problem decomposes into Z programs):
#   julia --project e13/run.jl --mosek --nobudget
#
#   Z=4    dof=2376   blk=65   IPM   56ms   HSD   60ms (1.07x)  Cla  153ms (2.73x)  Msk  142ms (2.52x)
#   Z=8    dof=5544   blk=65   IPM  152ms   HSD  173ms (1.14x)  Cla  432ms (2.85x)  Msk  399ms (2.63x)
#   Z=16   dof=11880  blk=65   IPM  364ms   HSD  425ms (1.17x)  Cla  966ms (2.65x)  Msk 1010ms (2.78x)
#   Z=32   dof=24552  blk=65   IPM  842ms   HSD  986ms (1.17x)  Cla 2029ms (2.41x)  Msk 2952ms (3.50x)
#
#   Slopes: IPM DOF^1.16, HSD DOF^1.20, Clarabel DOF^1.10, Mosek DOF^1.29
#   NOTE: Only Mosek shows decomposition-awareness (1.55 -> 1.29); Clarabel does not.
#
# BLOCK-SIZE DIAL (--ldial, Lf 16 -> 48 at Z = 8):
#   julia --project e13/run.jl --mosek --ldial
#
#   Lf=16  dof=7004   blk=65   IPM  204ms   HSD  218ms (1.07x)  Cla  527ms (2.59x)  Msk  484ms (2.37x)
#   Lf=32  dof=13788  blk=129  IPM  951ms   HSD 1048ms (1.10x)  Cla 3158ms (3.32x)  Msk 1515ms (1.59x)
#   Lf=48  dof=20572  blk=193  IPM 2626ms   HSD 3191ms (1.22x)  Cla10561ms (4.02x)  Msk 3843ms (1.46x)
#
#   Slopes: IPM DOF^2.36, HSD DOF^2.47, Clarabel DOF^2.77, Mosek DOF^1.90
#   NOTE: Ratios NOT stable — Clarabel steepest, Mosek shallowest (count-vs-size rule).
#
# PRE-REGISTRATION SCORECARD (predictions written before first run):
#   * Z-slope ~1.1-1.3 vs E11's 2.29:            CONFIRMED (1.17)
#   * win vs both at every Z:                    CONFIRMED (constant-factor
#     vs Clarabel — its slope is 0.10 shallower, crossover ~1e7+ DOF;
#     exponent win vs Mosek)
#   * --joint degrades the ratio:                CONFIRMED (2.14x -> 1.33x,
#     at the constant, slope ~unchanged)
#   * --nobudget: baselines gain via presolve:   PARTIAL (Mosek 1.55 ->
#     1.29 yes; Clarabel 1.10 no; IPM indifferent at 1.16 as predicted)
#   * --ldial ratios roughly stable:             REFUTED — they diverge
#     (Cla 2.77, Msk 1.90): the count-vs-size regime rule, see examples.md
#
# Oracle figures: see def.jl header.
# =============================================================================
