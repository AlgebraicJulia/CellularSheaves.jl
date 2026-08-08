# =============================================================================
# soft_landing_oracle.py — self-contained validation for E14 (powered descent)
#
# Minimum-fuel planetary soft landing, 3-DoF, fixed mass, by lossless
# convexification (LCvx): the convex (SOC) relaxation of
#
#     min  ∫ ||T(t)||/m dt                     (fuel; thrust acceleration u = T/m)
#     s.t. r̈ = u + g,       r(0), ṙ(0) given,  r(tf) = ṙ(tf) = 0,
#          ρ1 <= ||u(t)|| <= ρ2                (nonconvex lower bound),
#          z(t) >= tan(γ_gs) ||(x(t), y(t))||  (glideslope cone),
#
# via the slack σ:  ||u|| <= σ,  ρ1 <= σ <= ρ2  (Açıkmeşe & Ploen, JGCD 2007).
# Discretized by exact ZOH on a uniform grid (double integrator: A, B exact).
# Fixed mass is the published simplification of Luo–Echigo–Açıkmeşe
# (arXiv:2410.09748, §6.1, after Sostaric & Rea 2005); it keeps the program
# exactly conic (the varying-mass log-mass treatment makes the thrust bounds
# only approximately SOC). ALL solvers receive this identical convex
# relaxation; LCvx validity (||u|| >= ρ1 nodewise) is a MEASURED claim, not an
# assumption — in discrete time it can fail at up to n_x - 1 = 5 grid points
# (Luo–Echigo–Açıkmeşe 2024, Thm 20; typically 0–1 in practice).
#
# The example's reason to exist: t_f dials the program across a PHYSICALLY
# MEANINGFUL feasibility boundary (minimum landing time). Below t_f^min the
# program is primal infeasible — the HSD solver must return a Farkas
# certificate, and this oracle certifies where the boundary is by three
# independent routes (bisection, min-deviation root, 1-D closed form).
#
# Checks:
#   [1] ANALYTIC (1-D vertical subfamily, ρ1 = 0):
#       (a) closed-form optimal value: fuel* = g·tf + |v0z| exactly (L1
#           impulse identity: ∫||u|| >= ||∫u|| with equality iff u_z >= 0);
#       (b) the discrete value is DISCRETIZATION-EXACT (ZOH integrates the
#           double integrator and the fuel integral exactly): identical value
#           at N = 15 and N = 240;
#       (c) closed-form minimum landing time (suicide burn: free-fall coast,
#           then max-thrust brake; one quadratic) vs conic bisection.
#   [2] LCvx VALIDITY CENSUS (base 3-D instance): count nodes with
#       ||u_t|| < ρ1; assert <= n_x - 1 = 5 across an N-sweep; verify the
#       instance is in the NORMAL case (some σ_t > ρ1), i.e. not long-horizon.
#   [3] FEASIBILITY BOUNDARY (3-D): t_f^min by Clarabel status bisection;
#       independently as the root of the min-deviation problem
#       d(tf) = min ||(r_N, v_N)||  s.t. all other constraints (always
#       feasible; d = 0 iff the landing problem is feasible — the
#       minimum-landing-error construction of Blackmore–Açıkmeşe–Scharf,
#       JGCD 2010); d monotone decreasing below the boundary.
#   [4] LIMITS: (a) tf -> tf_min+ : thrust saturation fraction -> 1
#       (minimum time = maximum thrust); (b) ρ1 = 0, generous tf:
#       bang-off-bang — σ concentrates on {0, ρ2}; (c) fuel(tf) is convex
#       with an interior minimum (the Açıkmeşe–Ploen line-search structure).
#   [5] GLIDESLOPE LIVE: the cone is active at the base instance (margin ~ 0
#       at some node) and load-bearing (removing it strictly lowers fuel and
#       moves the trajectory).
#   [6] COSTATE (discrete PMP): dynamics-equality duals satisfy the adjoint
#       recursion η_{t-1} = A'η_t on glideslope-inactive interior nodes
#       (Luo–Echigo–Açıkmeşe 2024, Thm 10); thrust anti-parallel to B'η on
#       burn arcs.
#   [7] CROSS-SOLVER: Clarabel vs SCS on states and value.
#   [8] REFINEMENT (base 3-D): fuel(N) converges under N-refinement at fixed
#       tf; measured rate (control ZOH parameterization error).
#
# Usage: python3 soft_landing_oracle.py [--quick]
# =============================================================================

import sys
import numpy as np
import cvxpy as cp
import warnings
warnings.filterwarnings("ignore", category=UserWarning)

QUICK = "--quick" in sys.argv
np.set_printoptions(precision=4, suppress=True)

CLARABEL_TIGHT = dict(solver=cp.CLARABEL, tol_gap_abs=1e-11, tol_gap_rel=1e-11,
                      tol_feas=1e-11, max_iter=5000)

FEASIBLE = ("optimal", "optimal_inaccurate")
INFEASIBLE = ("infeasible", "infeasible_inaccurate")


# -----------------------------------------------------------------------------
# Instance
# -----------------------------------------------------------------------------
# Physical constants after Luo–Echigo–Açıkmeşe (2024) §6.1 / Sostaric–Rea
# (2005): m = 1000 kg, T_max = 7500 N, thrust range [0.3, 0.8] T_max, lunar
# gravity 1.62 m/s². Thrust ACCELERATION bounds ρ1 = 2.25, ρ2 = 6.0 m/s².
# Initial state adds a lateral divert (so the glideslope is live); target is
# the origin at rest.

class Instance:
    def __init__(self, r0, v0, rho1=2.25, rho2=6.0, gmag=1.62, gs_deg=25.0):
        self.r0 = np.asarray(r0, float)
        self.v0 = np.asarray(v0, float)
        self.rho1, self.rho2, self.gmag = rho1, rho2, gmag
        self.gvec = np.array([0.0, 0.0, -gmag])
        self.tan_gs = np.tan(np.deg2rad(gs_deg))


# gs = 55 deg: the unconstrained trajectory's approach angle dips to 51.6 deg,
# so the 55-deg cone is genuinely live mid-trajectory (measured in check [5]).
BASE = Instance(r0=[600.0, -400.0, 2000.0], v0=[30.0, -10.0, -60.0], gs_deg=55.0)
TF_BASE = 60.0
N_BASE = 60

# 1-D vertical subfamily for the analytic gates: no lateral state, ρ1 = 0.
VERT = Instance(r0=[0.0, 0.0, 2000.0], v0=[0.0, 0.0, -60.0], rho1=0.0)


def zoh(tf, N):
    """Exact ZOH discretization of the double integrator."""
    dt = tf / N
    A = np.block([[np.eye(3), dt * np.eye(3)],
                  [np.zeros((3, 3)), np.eye(3)]])
    B = np.vstack([dt**2 / 2 * np.eye(3), dt * np.eye(3)])
    return dt, A, B


# -----------------------------------------------------------------------------
# Conic solves
# -----------------------------------------------------------------------------

def build(inst, tf, N, mode="fuel", glideslope=True):
    """mode: 'fuel' (landing problem) or 'dev' (min terminal deviation)."""
    dt, A, B = zoh(tf, N)
    X = cp.Variable((N + 1, 6))
    U = cp.Variable((N, 3))
    s = cp.Variable(N)
    cons = [X[0] == np.hstack([inst.r0, inst.v0])]
    for t in range(N):
        cons.append(X[t + 1] == A @ X[t] + B @ (U[t] + inst.gvec))
    cons += [cp.norm(U, 2, axis=1) <= s, s <= inst.rho2]
    if inst.rho1 > 0:
        cons.append(s >= inst.rho1)
    else:
        cons.append(s >= 0)
    if glideslope:
        cons.append(cp.norm(X[:, :2], 2, axis=1) * inst.tan_gs <= X[:, 2])
    if mode == "fuel":
        cons.append(X[N] == np.zeros(6))
        obj = cp.Minimize(cp.sum(s) * dt)
    else:  # 'dev': terminal state free, minimize distance to target
        obj = cp.Minimize(cp.norm(X[N]))
    return cp.Problem(obj, cons), X, U, s


def solve(prob, tight=True):
    try:
        if tight:
            prob.solve(**CLARABEL_TIGHT)
        else:
            prob.solve(solver=cp.CLARABEL)
    except (cp.SolverError, Exception):
        pass
    return prob.status


def solve_fuel(inst, tf, N, glideslope=True, tight=True):
    prob, X, U, s = build(inst, tf, N, "fuel", glideslope)
    st = solve(prob, tight)
    return st, prob, X, U, s


# -----------------------------------------------------------------------------
# Boundary machinery
# -----------------------------------------------------------------------------

def bisect_tfmin_status(inst, N, lo, hi, tol=1e-3):
    """Smallest feasible tf by Clarabel status bisection. lo infeasible, hi feasible."""
    st_lo, *_ = solve_fuel(inst, lo, N, tight=False)
    st_hi, *_ = solve_fuel(inst, hi, N, tight=False)
    assert st_lo in INFEASIBLE and st_hi in FEASIBLE, (st_lo, st_hi)
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        st, *_ = solve_fuel(inst, mid, N, tight=False)
        if st in FEASIBLE:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def deviation(inst, tf, N):
    """Min-deviation value d(tf) >= 0; d = 0 iff the landing problem feasible."""
    prob, X, U, s = build(inst, tf, N, "dev")
    st = solve(prob, tight=False)
    assert st in FEASIBLE, st
    return prob.value


def bisect_tfmin_dev(inst, N, lo, hi, dev_tol=1e-5, tol=1e-3):
    """Smallest tf with d(tf) <= dev_tol, by bisection (d continuous, monotone)."""
    assert deviation(inst, lo, N) > dev_tol and deviation(inst, hi, N) <= dev_tol
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if deviation(inst, mid, N) <= dev_tol:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


# -----------------------------------------------------------------------------
# 1-D closed forms (VERT subfamily, ρ1 = 0)
# -----------------------------------------------------------------------------

def fuel_1d_closed_form(inst, tf):
    """L1 impulse identity: ∫||u|| >= ||∫u dt|| = ||v_f - v_0 - g tf|| with
    equality iff u_z >= 0 throughout (attainable when a monotone braking
    profile exists). For descent to rest: fuel* = g tf + |v0z|."""
    return inst.gmag * tf + abs(inst.v0[2])


def tfmin_1d_closed_form(inst):
    """Minimum landing time, 1-D, ρ1 = 0: POWERED DIVE, not coast-then-brake.
    With ||u|| <= ρ2 free in direction, the min-time profile is bang-bang with
    a thrust reversal: phase 1 thrusts DOWN (a_z = -(ρ2 + g), diving faster
    than free fall), phase 2 max-brakes up (a_z = +(ρ2 - g)). PMP: the
    switching function is linear in the costate — one switch. This is the
    standard min-time double integrator with asymmetric authority
    a_z ∈ [-(ρ2+g), ρ2-g] (gravity biasing).
      v: 0 = v0 - a1 t1 + a2 t2,  a1 = ρ2 + g, a2 = ρ2 - g
      z: braking from v1 = v0 - a1 t1 needs altitude v1²/(2 a2):
         h + v0 t1 - a1 t1²/2 = (v0 - a1 t1)²/(2 a2)     (quadratic in t1)
    (The naive coast-then-brake ansatz overestimates tf_min — measured 35.50
    vs the true 32.22 on the VERT instance; the conic solver found the dive.)"""
    h, v0 = inst.r0[2], inst.v0[2]
    a1, a2 = inst.rho2 + inst.gmag, inst.rho2 - inst.gmag
    A2 = -a1 / 2 - a1**2 / (2 * a2)
    A1 = v0 * (1 + a1 / a2)
    A0 = h - v0**2 / (2 * a2)
    disc = A1**2 - 4 * A2 * A0
    assert disc >= 0
    roots = [(-A1 + np.sqrt(disc)) / (2 * A2), (-A1 - np.sqrt(disc)) / (2 * A2)]
    t1 = min(r for r in roots if r >= 0)
    t2 = (a1 * t1 - v0) / a2
    assert t2 >= 0
    return t1 + t2

# -----------------------------------------------------------------------------
# Checks
# -----------------------------------------------------------------------------

def check1_analytic_1d():
    print("[1] ANALYTIC (1-D vertical, rho1 = 0)")
    tf = 60.0
    ref = fuel_1d_closed_form(VERT, tf)
    st, prob, X, U, s = solve_fuel(VERT, tf, 60)
    assert st in FEASIBLE, st
    err = abs(prob.value - ref) / ref
    print(f"    (a) closed-form fuel* = g tf + |v0z| = {ref:.6f};"
          f" conic = {prob.value:.6f}; rel err = {err:.2e}")
    assert err < 1e-7
    vals = []
    for N in (15, 240):
        stN, probN, *_ = solve_fuel(VERT, tf, N)
        assert stN in FEASIBLE
        vals.append(probN.value)
    dis = abs(vals[0] - vals[1]) / ref
    print(f"    (b) discretization-exact: |fuel(N=15) - fuel(N=240)| / fuel* = {dis:.2e}")
    assert dis < 1e-7
    ref_tf = tfmin_1d_closed_form(VERT)
    N = 40 if QUICK else 80
    tfmin = bisect_tfmin_status(VERT, N, lo=0.5 * ref_tf, hi=1.5 * ref_tf,
                                tol=1e-3)
    print(f"    (c) suicide-burn tf_min closed form = {ref_tf:.4f} s;"
          f" conic bisection = {tfmin:.4f} s; |Δ| = {abs(tfmin - ref_tf):.2e} s")
    assert abs(tfmin - ref_tf) < 0.1
    return ref_tf


def check2_lcvx_census():
    print("[2] LCvx VALIDITY CENSUS (base 3-D instance, n_x - 1 = 5)")
    grid = (20, 40, 60) if QUICK else (20, 40, 60, 120)
    worst = 0
    for N in grid:
        st, prob, X, U, s = solve_fuel(BASE, TF_BASE, N)
        assert st in FEASIBLE, st
        unorm = np.linalg.norm(U.value, axis=1)
        bad = int(np.sum(unorm < BASE.rho1 * (1 - 1e-6)))
        worst = max(worst, bad)
        normal = np.max(s.value) > BASE.rho1 + 1e-3
        gap = np.max(s.value - unorm)
        print(f"    N = {N:4d}: invalid nodes = {bad}  (max sigma - ||u|| = {gap:.2e};"
              f" normal case: {normal})")
        assert normal, "long-horizon case — choose smaller tf"
    print(f"    worst census = {worst} <= 5 (Luo-Echigo-Acikmese 2024, Thm 20)")
    assert worst <= 5
    return worst


def check3_boundary():
    print("[3] FEASIBILITY BOUNDARY (base 3-D instance)")
    N = 40
    tf_status = bisect_tfmin_status(BASE, N, lo=20.0, hi=45.0, tol=1e-3)
    tf_dev = bisect_tfmin_dev(BASE, N, lo=20.0, hi=45.0, dev_tol=1e-5, tol=1e-3)
    print(f"    status-flip bisection: tf_min = {tf_status:.4f} s")
    print(f"    min-deviation root:    tf_min = {tf_dev:.4f} s"
          f"   |Δ| = {abs(tf_status - tf_dev):.2e} s")
    assert abs(tf_status - tf_dev) < 5e-3
    tfs = np.array([0.6, 0.7, 0.8, 0.9, 0.96]) * tf_status
    ds = [deviation(BASE, tf, N) for tf in tfs]
    mono = all(ds[i] > ds[i + 1] for i in range(len(ds) - 1))
    print("    d(tf) below boundary: " +
          ", ".join(f"{tf:.1f}s: {d:.3f}" for tf, d in zip(tfs, ds)))
    print(f"    monotone decreasing: {mono};  d(1.1 tf_min) = "
          f"{deviation(BASE, 1.1 * tf_status, N):.2e}")
    assert mono
    assert deviation(BASE, 1.1 * tf_status, N) < 1e-6
    return tf_status


def check4_limits(tfmin3d):
    print("[4] LIMITS")
    N = 40
    fracs = []
    for mult in (1.002, 1.02, 1.1, 1.5):
        st, prob, X, U, s = solve_fuel(BASE, mult * tfmin3d, N)
        assert st in FEASIBLE, st
        frac = np.mean(s.value >= BASE.rho2 * (1 - 1e-4))
        fracs.append(frac)
        print(f"    (a) tf = {mult:.3f} tf_min: saturated fraction = {frac:.3f}")
    assert fracs[0] > 0.9 and all(fracs[i] >= fracs[i + 1] - 1e-9
                                  for i in range(len(fracs) - 1))
    relax = Instance(r0=BASE.r0, v0=BASE.v0, rho1=0.0)
    st, prob, X, U, s = solve_fuel(relax, TF_BASE, 60)
    assert st in FEASIBLE
    sv = s.value
    onoff = np.mean((sv < 1e-3 * BASE.rho2) | (sv > BASE.rho2 * (1 - 1e-3)))
    print(f"    (b) rho1 = 0, tf = {TF_BASE:.0f}: bang-off-bang fraction = {onoff:.3f}")
    assert onoff > 0.8
    tf_grid = np.array([40, 45, 50, 55, 60, 70, 85, 100.0])
    fuels = []
    for tf in tf_grid:
        stf, probf, *_ = solve_fuel(BASE, tf, N)
        assert stf in FEASIBLE
        fuels.append(probf.value)
    fuels = np.array(fuels)
    kmin = int(np.argmin(fuels))
    d2 = np.diff(fuels, 2)
    print(f"    (c) fuel(tf) interior minimum at tf ~ {tf_grid[kmin]:.0f}"
          f" (fuel = {fuels[kmin]:.2f}); convex on grid: {np.all(d2 > -1e-6)}")
    assert 0 < kmin < len(tf_grid) - 1 and np.all(d2 > -1e-6)


def check5_glideslope():
    print("[5] GLIDESLOPE LIVE (base instance)")
    st, prob, X, U, s = solve_fuel(BASE, TF_BASE, N_BASE)
    assert st in FEASIBLE
    # exclude the final node: the target IS the cone apex, margin 0 trivially
    margin = (X.value[:-1, 2]
              - BASE.tan_gs * np.linalg.norm(X.value[:-1, :2], axis=1))
    st0, prob0, X0, U0, s0 = solve_fuel(BASE, TF_BASE, N_BASE, glideslope=False)
    assert st0 in FEASIBLE
    dtraj = np.max(np.abs(X.value - X0.value))
    dfuel = prob.value - prob0.value
    print(f"    min cone margin = {np.min(margin):.2e} m (active);"
          f" without cone: fuel Δ = {dfuel:+.4f}, ||Δtraj||_inf = {dtraj:.2f} m")
    assert np.min(margin) < 1e-4
    assert dfuel > 1e-3 and dtraj > 1.0


def check6_costate():
    print("[6] COSTATE (discrete PMP; Luo-Echigo-Acikmese 2024, Thm 10)")
    st, prob, X, U, s = solve_fuel(BASE, TF_BASE, N_BASE)
    assert st in FEASIBLE
    dt, A, B = zoh(TF_BASE, N_BASE)
    margin = X.value[:, 2] - BASE.tan_gs * np.linalg.norm(X.value[:, :2], axis=1)
    inactive = margin > 1e-3   # nodes where the glideslope is slack
    etas = [prob.constraints[1 + t].dual_value for t in range(N_BASE)]
    res, den, align = [], [], []
    for t in range(1, N_BASE):
        if inactive[t]:
            res.append(np.linalg.norm(etas[t - 1] - A.T @ etas[t]))
            den.append(np.linalg.norm(etas[t - 1]))
        # cvxpy's equality-dual sign makes u parallel to +B'eta_cvxpy; the PMP
        # statement u ∥ -B'eta holds for eta = -eta_cvxpy (recursion unaffected).
        bu = B.T @ etas[t]
        u = U.value[t]
        if np.linalg.norm(u) > 0.5 * BASE.rho1 and np.linalg.norm(bu) > 1e-9:
            align.append(np.dot(u, bu) / (np.linalg.norm(u) * np.linalg.norm(bu)))
    rel = max(res) / max(den)
    print(f"    adjoint recursion eta_(t-1) = A' eta_t on {len(res)}"
          f" cone-inactive nodes: rel residual = {rel:.2e}")
    print(f"    thrust anti-parallel to B'eta on {len(align)} burn nodes:"
          f" min cos = {min(align):.6f}")
    assert rel < 1e-5
    assert min(align) > 1 - 1e-4


def check7_cross_solver():
    print("[7] CROSS-SOLVER (Clarabel vs SCS)")
    st, prob, X, U, s = solve_fuel(BASE, TF_BASE, N_BASE)
    assert st in FEASIBLE
    xa, va = X.value.copy(), prob.value
    prob2, X2, U2, s2 = build(BASE, TF_BASE, N_BASE)
    prob2.solve(solver=cp.SCS, eps=1e-9, max_iters=200_000)
    assert prob2.status in FEASIBLE, prob2.status
    dx = np.max(np.abs(xa - X2.value))
    dv = abs(va - prob2.value) / va
    print(f"    ||ΔX||_inf = {dx:.2e} m; rel Δfuel = {dv:.2e}")
    assert dx < 1e-3 and dv < 1e-6


def check8_refinement():
    print("[8] REFINEMENT (base 3-D instance, fixed tf)")
    Ns = (20, 40, 80) if QUICK else (20, 40, 80, 160)
    Nref = 320
    stR, probR, *_ = solve_fuel(BASE, TF_BASE, Nref)
    assert stR in FEASIBLE
    ref = probR.value
    errs = []
    for N in Ns:
        stN, probN, *_ = solve_fuel(BASE, TF_BASE, N)
        assert stN in FEASIBLE
        errs.append(abs(probN.value - ref))
        print(f"    N = {N:4d}: fuel = {probN.value:.6f}  |err| = {errs[-1]:.2e}")
    rates = [np.log2(errs[i] / errs[i + 1]) for i in range(len(errs) - 1)]
    print(f"    measured order (log2 err ratios): "
          + ", ".join(f"{r:.2f}" for r in rates)
          + f"   (reference N = {Nref}, fuel = {ref:.6f})")
    assert errs[-1] < errs[0]


def main():
    print("=" * 77)
    print("E14 soft-landing oracle" + ("  [--quick]" if QUICK else ""))
    dt, A, B = zoh(TF_BASE, N_BASE)
    print(f"base: r0 = {BASE.r0}, v0 = {BASE.v0}, rho = [{BASE.rho1}, {BASE.rho2}]"
          f" m/s^2, g = {BASE.gmag}, gs = 55 deg, tf = {TF_BASE}, N = {N_BASE}")
    print("=" * 77)
    check1_analytic_1d()
    check2_lcvx_census()
    tfmin3d = check3_boundary()
    check4_limits(tfmin3d)
    check5_glideslope()
    check6_costate()
    check7_cross_solver()
    check8_refinement()
    print("=" * 77)
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
