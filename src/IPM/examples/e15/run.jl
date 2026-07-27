include("def.jl")

run()

# =============================================================================
# Sample run: 2026-07-25 (--quick --mosek), after the Cofree-dual fix (see note).
# All 7 gates green. IPM = affine primal-dual, HSD = homogeneous self-dual (both
# the sheaf solver). Mosek shown primal/dual; primal is its better form here.
# -----------------------------------------------------------------------------
# V-sweep (matched density — the mission config; N = 60) — THE HEADLINE:
# V=4    dof=4128   n1=3168    IPM  29.7ms  HSD  30.3ms  Cla  43.0ms  Msk  73.4ms/103.3ms
# V=8    dof=8256   n1=6336    IPM  65.2ms  HSD  60.7ms  Cla  88.5ms  Msk 279.4ms/245.8ms
# V=16   dof=16512  n1=12672   IPM 186.3ms  HSD 196.2ms  Cla 435.9ms  Msk 368.6ms/548.1ms
# V-slopes:  IPM V^1.32  HSD V^1.35  Clarabel V^1.67  Msk(p) V^1.16  Msk(d) V^1.20
#
# N-sweep (V = 4; stall watch):
# N=60   dof=4128   n1=3168    IPM  29.2ms  HSD  29.7ms  Cla  43.1ms  Msk  73.1ms/103.3ms
# N=120  dof=8208   n1=6288    IPM  63.1ms  HSD  58.8ms  Cla 122.8ms  Msk 177.5ms/279.6ms
# N=240  dof=16368  n1=12528   IPM 142.4ms  HSD 136.0ms  Cla 208.7ms  Msk 289.0ms/433.1ms
# N-slopes:  IPM N^1.14  HSD N^1.10  Clarabel N^1.14  Msk(p) N^0.99  Msk(d) N^1.03
#
# Verdict: both sheaf solvers beat the conic baselines at every V and N — at V=16
# HSD is 2.2x faster than Clarabel and 1.9x than Mosek, and the gap still widens
# with V (HSD^1.35 < Clarabel^1.67). On the N-dial HSD leads and the affine
# long-horizon stall did NOT recur on CW chains (clean to N=240).
#   CAVEAT ON THE SLOPE: an earlier (pre-fix) run reported HSD V^0.97 — that
# apparent sublinearity was partly an artifact of the Cofree-dual bug below: at
# larger V, HSD terminated early (μ<0 → NEAR_OPTIMAL) and so looked faster/flatter
# than honest convergence. With the bug fixed HSD converges fully and the true
# V-slope is ~1.35. The moat is a constant-factor + shallower-slope win, not
# sublinear scaling.
#
# RESOLVED (2026-07-25): at κ = 1e-4 HSD had been leaving the cone (μ<0) and
# mislabeling NEAR_OPTIMAL. Root cause: the Cofree (free) dual is definitionally 0
# (dual cone {0}) but was drifting — the free-row dual residual (~1e-6, the Woodbury
# floor amplified by the collapsing border capacitance 1/S) times the large free
# primals polluted μ = (pᵀd+τκ)/(ν+1), whose numerator summed free cones though ν
# excludes them (degree 0). Fix: pin the Cofree dual direction to 0 at the recovery
# sites in step! (hsd.jl/ipm.jl). HSD now reports OPTIMAL and lands within ~5e-4 of
# the oracle at κ=1e-4. The frontier gate still uses IPM for the tightest reference.
#
# Frozen references (tools/cw_oracle.py, fully executed): obj(base) = 1.30292563;
# κ = 0 per-vehicle fuels [0.175486, 0.001035, 0.148128, 0.013406] m/s; hover rate
# ω²√(9·80² + 40²) = 3.106832e-4 (oracle match 1.6e-4); frontier fuels 0.3885 /
# 0.8456 / 1.8628 m/s at κ = 1e-4/1e-3/1e-2; STM analytic vs expm 5.2e-14.
# =============================================================================
