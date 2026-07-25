# =============================================================================
# tracking_oracle.py -- self-contained validation for E14b (fuel + tracking)
#
# E14's chain with a quadratic state-tracking term (dense-Q habitat):
#
#   min  Sum_t sigma_t dt  +  (gamma/2) Sum_{t} (x_t - xbar_t)' W (x_t - xbar_t)
#   s.t. E14's constraint set unchanged (dynamics, pins, thrust SOC + bounds,
#        glideslope).
#
# W = inv(P_disp), P_disp the Van Loan dispersion of white accelerometer
# noise through the double integrator over a characteristic time tau_c --
# dense 6x6 with physically real r-v cross terms (velocity error integrates
# into position error). Normalized to unit spectral norm so gamma is
# interpretable. Reference xbar: straight-line descent corridor (constant
# velocity -r0/tf), deliberately endpoint-inconsistent with the pins.
#
# Purpose (besides the missing dense-Q habitat cell): the STALL-CURE
# experiment. E14's affine-IPM stall is a degenerate optimal face from
# bang-off-bang bounds (e14/stall.md). Tracking makes the optimum unique in
# (x, u) (the u -> state map is injective), so strict complementarity should
# hold generically for gamma above a threshold. The oracle measures the
# degeneracy directly (check [4]) as solver-independent evidence.
#
# Checks:
#   [1] INVARIANCE: with xbar := the E14 optimal trajectory, x(gamma) == xbar
#       and fuel(gamma) == fuel_E14 for all gamma (simultaneous minimizer).
#   [2] LIMITS: gamma -> 0 recovers the E14 solution at measured rate.
#   [3] FRONTIER: fuel(gamma) nondecreasing, tracking(gamma) nonincreasing.
#   [4] DEGENERACY vs gamma: strict-complementarity margin min(a + lam_a)
#       over rho1-active arcs, and the count of arcs with BOTH a ~ 0 and
#       lam_a ~ 0, at gamma in {0, base, large}. Prediction: margin grows,
#       doubly-degenerate count -> 0.
#   [5] LCvx census vs gamma (measured claim; tracking may pull ||u|| below
#       rho1, growing the census -- the relaxation-tightness caveat).
#   [6] CROSS-SOLVER: Clarabel vs SCS.
#   [7] BOUNDARY INVARIANCE: feasibility is objective-independent; statuses
#       at 0.98 / 1.02 tf_min match E14.
#   [8] COSTATE: adjoint recursion gains the tracking forcing,
#       eta_{t-1} = A' eta_t + gamma W (x_t - xbar_t) on cone-inactive nodes.
#
# Usage: python3 tracking_oracle.py [--quick]
# =============================================================================

import sys
import numpy as np
import cvxpy as cp
import warnings
warnings.filterwarnings("ignore")

from soft_landing_oracle import (BASE, TF_BASE, zoh, build, solve, solve_fuel,
                                 FEASIBLE, CLARABEL_TIGHT, bisect_tfmin_status)

QUICK = "--quick" in sys.argv
TFMIN40 = 35.5376

# -----------------------------------------------------------------------------
# W (Van Loan dispersion metric) and the reference corridor
# -----------------------------------------------------------------------------

def dispersion_W(tau_c=8.0, q=1.0, eps=1e-3):
    """W = inv(Van Loan P + eps*I), unit spectral norm. Dense by nature."""
    Ac = np.block([[np.zeros((3, 3)), np.eye(3)], [np.zeros((3, 6))]])
    L = np.vstack([np.zeros((3, 3)), np.eye(3)])
    # Van Loan: F = expm([[-Ac, L q L'], [0, Ac']] tau) ; P = F22' F12
    import scipy.linalg as sla
    M = np.block([[-Ac, q * (L @ L.T)], [np.zeros((6, 6)), Ac.T]])
    F = sla.expm(M * tau_c)
    P = F[6:, 6:].T @ F[:6, 6:]
    P = 0.5 * (P + P.T) + eps * np.eye(6)
    W = np.linalg.inv(P)
    return 0.5 * (W + W.T) / np.linalg.norm(W, 2)

W_BASE = dispersion_W()

def corridor(inst, tf, N):
    """Straight-line descent: r linear r0 -> 0, v = -r0/tf constant."""
    xb = np.zeros((N + 1, 6))
    vbar = -inst.r0 / tf
    for t in range(N + 1):
        xb[t, :3] = inst.r0 * (1 - t / N)
        xb[t, 3:] = vbar
    return xb

GAMMA_BASE = 2e-3   # tuned below: tracking cost ~ fuel/3 at base

# -----------------------------------------------------------------------------
# Solves
# -----------------------------------------------------------------------------

def build_track(inst, tf, N, gamma, xbar, W=W_BASE, glideslope=True):
    dt, A, Bd = zoh(tf, N)
    X = cp.Variable((N + 1, 6))
    U = cp.Variable((N, 3))
    s = cp.Variable(N)
    cons = [X[0] == np.hstack([inst.r0, inst.v0])]
    for t in range(N):
        cons.append(X[t + 1] == A @ X[t] + Bd @ (U[t] + inst.gvec))
    cons += [cp.norm(U, 2, axis=1) <= s, s <= inst.rho2]
    lower = (s >= inst.rho1) if inst.rho1 > 0 else (s >= 0)
    cons.append(lower)
    if glideslope:
        cons.append(cp.norm(X[:, :2], 2, axis=1) * inst.tan_gs <= X[:, 2])
    cons.append(X[N] == np.zeros(6))
    track = cp.sum([cp.quad_form(X[t] - xbar[t], cp.psd_wrap(W))
                    for t in range(N + 1)])
    prob = cp.Problem(cp.Minimize(cp.sum(s) * dt + 0.5 * gamma * track), cons)
    return prob, X, U, s, lower

def solve_track(inst, tf, N, gamma, xbar, **kw):
    prob, X, U, s, lower = build_track(inst, tf, N, gamma, xbar, **kw)
    try:
        prob.solve(**CLARABEL_TIGHT)
    except Exception:
        prob.solve(solver=cp.CLARABEL)
    assert prob.status in FEASIBLE, (gamma, prob.status)
    dt = tf / N
    fuel = dt * float(np.sum(s.value))
    tr = float(sum((X.value[t] - xbar[t]) @ W_BASE @ (X.value[t] - xbar[t])
                   for t in range(N + 1)))
    return prob, X, U, s, lower, fuel, tr

# -----------------------------------------------------------------------------
# Checks
# -----------------------------------------------------------------------------

def check1_invariance(X14, fuel14):
    print("[1] INVARIANCE (xbar := E14 optimum -> simultaneous minimizer)")
    for gamma in (0.1, 10.0):
        _, X, U, s, _, fuel, tr = solve_track(BASE, TF_BASE, 60, gamma, X14)
        dx = np.max(np.abs(X.value - X14))
        df = abs(fuel - fuel14) / fuel14
        print(f"    gamma = {gamma:5.1f}: ||x - xbar||_inf = {dx:.2e} m, "
              f"rel dfuel = {df:.2e}, tracking = {tr:.2e}")
        assert dx < 1e-3 and df < 1e-6
    return True

def check2_gamma0_limit(X14, xbar):
    print("[2] LIMIT gamma -> 0 (recovers E14)")
    errs, gs = [], (1e-4, 1e-5, 1e-6)
    for g in gs:
        _, X, *_ = solve_track(BASE, TF_BASE, 60, g, xbar)
        errs.append(np.max(np.abs(X.value - X14)))
        print(f"    gamma = {g:.0e}: ||x - x_E14||_inf = {errs[-1]:.2e} m")
    rate = np.log(errs[0] / errs[-1]) / np.log(gs[0] / gs[-1])
    print(f"    measured rate: gamma^{rate:.2f}")
    assert errs[0] > errs[1] > errs[2], "not monotone in gamma"
    assert 0.8 < rate < 1.2, f"expected ~linear rate, got {rate}" 

def check3_frontier(xbar):
    print("[3] FRONTIER (fuel up, tracking down in gamma)")
    gs = np.array([2e-4, 6e-4, 2e-3, 6e-3, 2e-2])
    fuels, trs = [], []
    for g in gs:
        _, X, U, s, _, fuel, tr = solve_track(BASE, TF_BASE, 60, g, xbar)
        fuels.append(fuel); trs.append(tr)
        print(f"    gamma = {g:.0e}: fuel = {fuel:8.3f}  tracking = {tr:10.2f}"
              f"  gamma*tr/2 = {0.5 * g * tr:7.3f}")
    assert all(fuels[i] <= fuels[i + 1] + 1e-6 for i in range(len(gs) - 1))
    assert all(trs[i] >= trs[i + 1] - 1e-6 for i in range(len(gs) - 1))

def check4_degeneracy(xbar):
    print("[4] DEGENERACY vs gamma (the stall-cure evidence)")
    out = {}
    for g in (0.0, GAMMA_BASE, 10 * GAMMA_BASE):
        prob, X, U, s, lower, fuel, tr = solve_track(BASE, TF_BASE, 60, g, xbar)
        lam = np.asarray(lower.dual_value).ravel()
        a = s.value - BASE.rho1
        act = a < 1e-6 * (1 + BASE.rho1)
        both = int(np.sum(act & (lam < 1e-6)))
        marg = np.min(a + lam) if act.any() else np.inf
        n_act = int(act.sum())
        out[g] = (n_act, both, marg)
        print(f"    gamma = {g:.0e}: rho1-active arcs = {n_act:3d};"
              f" doubly-degenerate (a ~ 0 AND lam ~ 0) = {both:3d};"
              f" strict-compl margin min(a + lam) = {marg:.2e}")
    m0, mb = out[0.0][2], out[GAMMA_BASE][2]
    assert out[GAMMA_BASE][1] <= out[0.0][1]
    assert (mb > 10 * m0) or out[GAMMA_BASE][1] == 0
    # [4b] the ONSET SCALING LAW: at gamma = 0 the min strict-compl margin
    # collapses as ~1/N^2 (measured 1.41e-2 -> 2.20e-4 over N = 60 -> 480,
    # exactly x4 per doubling); at N = 480 it numerically matches the
    # affine-IPM limit-cycle amplitude (mu ~ 2.1-2.5e-4, e14/stall.md).
    # gamma_base lifts it ~30x. Exact degeneracy never occurs (doubly-deg = 0
    # at every N): the face is NEAR-degenerate, margin -> 0 as 1/N^2.
    print("    [4b] margin scaling vs N (onset law; stall at N >= 240):")
    Ns = (60, 120, 240) if QUICK else (60, 120, 240, 480)
    margins0, marginsb = [], []
    for N in Ns:
        xbN = corridor(BASE, TF_BASE, N)
        for g, acc in ((0.0, margins0), (GAMMA_BASE, marginsb)):
            _, X, U, s, lower, fuel, tr = solve_track(BASE, TF_BASE, N, g, xbN)
            lam = np.asarray(lower.dual_value).ravel()
            a = s.value - BASE.rho1
            act = a < 1e-6 * (1 + BASE.rho1)
            acc.append(float(np.min((a + lam)[act])))
        print(f"      N = {N:4d}: min(a + lam) gamma=0: {margins0[-1]:.2e}"
              f"   gamma_base: {marginsb[-1]:.2e}"
              f"   lift x{marginsb[-1] / margins0[-1]:.0f}")
    slope = (np.log(margins0[-1] / margins0[0])
             / np.log(Ns[-1] / Ns[0]))
    print(f"      gamma = 0 scaling: margin ~ N^{slope:.2f} (law: -2)")
    assert -2.4 < slope < -1.6, f"onset law broken: N^{slope}"
    assert all(mb > 5 * m0 for m0, mb in zip(margins0, marginsb))
    return out

def check5_lcvx(xbar):
    print("[5] LCvx CENSUS vs gamma (measured; the tightness caveat)")
    for g in (0.0, GAMMA_BASE, 10 * GAMMA_BASE):
        _, X, U, s, *_ = solve_track(BASE, TF_BASE, 60, g, xbar)
        unorm = np.linalg.norm(U.value, axis=1)
        bad = int(np.sum(unorm < BASE.rho1 * (1 - 1e-6)))
        print(f"    gamma = {g:.0e}: nodes with ||u|| < rho1: {bad}")

def check6_cross(xbar):
    print("[6] CROSS-SOLVER (Clarabel vs SCS at base gamma)")
    p1, X1, *_ = solve_track(BASE, TF_BASE, 60, GAMMA_BASE, xbar)
    p2, X2, U2, s2, lower2 = build_track(BASE, TF_BASE, 60, GAMMA_BASE, xbar)
    p2.solve(solver=cp.SCS, eps=1e-9, max_iters=200_000)
    assert p2.status in FEASIBLE
    dx = np.max(np.abs(X1.value - X2.value))
    print(f"    ||dX||_inf = {dx:.2e} m")
    assert dx < 1e-3

def check7_boundary(xbar):
    print("[7] BOUNDARY INVARIANCE (objective-independent feasibility)")
    for f, want in ((0.98, False), (1.02, True)):
        prob, *_ = build_track(BASE, f * TFMIN40, 40, GAMMA_BASE,
                               corridor(BASE, f * TFMIN40, 40))
        try:
            prob.solve(solver=cp.CLARABEL)
        except Exception:
            pass
        ok = prob.status in FEASIBLE
        print(f"    tf = {f:.2f} tf_min: {prob.status}")
        assert ok == want

def check8_costate(xbar):
    print("[8] COSTATE (adjoint recursion with tracking forcing)")
    prob, X, U, s, lower, fuel, tr = solve_track(BASE, TF_BASE, 60, GAMMA_BASE, xbar)
    dt, A, Bd = zoh(TF_BASE, 60)
    margin = X.value[:, 2] - BASE.tan_gs * np.linalg.norm(X.value[:, :2], axis=1)
    inactive = margin > 1e-3
    etas = [prob.constraints[1 + t].dual_value for t in range(60)]
    res, den = [], []
    for t in range(1, 60):
        if inactive[t]:
            forcing = GAMMA_BASE * W_BASE @ (X.value[t] - xbar[t])
            # cvxpy dual sign: eta_cvxpy = -eta_pmp; recursion picks up -forcing
            res.append(np.linalg.norm(etas[t - 1] - A.T @ etas[t] + forcing))
            den.append(np.linalg.norm(etas[t - 1]))
    rel = max(res) / max(den)
    print(f"    rel residual on {len(res)} cone-inactive nodes: {rel:.2e}")
    assert rel < 1e-5

def main():
    print("=" * 77)
    print("E14b tracking oracle" + ("  [--quick]" if QUICK else ""))
    print(f"W: Van Loan dispersion metric, tau_c = 8 s (dense 6x6, r-v cross "
          f"terms, ||W|| = 1); gamma_base = {GAMMA_BASE}")
    print("=" * 77)
    st, p14, X14v, U14, s14 = solve_fuel(BASE, TF_BASE, 60)
    assert st in FEASIBLE
    X14, fuel14 = X14v.value, p14.value
    xbar = corridor(BASE, TF_BASE, 60)
    check1_invariance(X14, fuel14)
    check2_gamma0_limit(X14, xbar)
    check3_frontier(xbar)
    check4_degeneracy(xbar)
    check5_lcvx(xbar)
    check6_cross(xbar)
    check7_boundary(xbar)
    check8_costate(xbar)
    print("=" * 77)
    print("ALL CHECKS PASSED")

if __name__ == "__main__":
    main()
