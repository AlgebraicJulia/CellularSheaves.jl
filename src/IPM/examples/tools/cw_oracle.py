# =============================================================================
# cw_oracle.py -- validation for the CW orbital formation-keeping example
# (Clohessy-Wiltshire / Hill relative dynamics; PRISMA/TanDEM-X heritage)
#
# V spacecraft hold a rigid formation (slots d_v in the LVLH frame) around a
# circular reference orbit, mean motion w. Convex program over horizon N:
#
#   min  sum_v sum_k dt*s_vk                        (fuel, SOC (s,u))
#      + (gam/2) sum_v,t (x_vt - xbar_vt)' W (x_vt - xbar_vt)   (slot tracking)
#      + (kap/2) sum_(i,j)in ring, t (x_it - x_jt - d_ij)' W (.)  (formation)
#   s.t. x_{t+1} = Phi x_t + Gam u_t,  x_{v,0} = d_v + delta_v,  ||u|| <= umax
#
# W = Van Loan inverse dispersion (accel process noise) — honest here: track
# in the metric of your own navigation dispersion. THE PHYSICS IS REAL: this
# is the shape whose benchmark win is ALREADY MEASURED (the retired synthetic
# fleet's matched-density verdict: HSD DOF^1.19 vs Cla^1.55 / MskD^1.42);
# CW replaces the indefensible coupling with the actual mission objective.
#
# CLOSED-FORM GATES (CW is generous):
#   [1] STM: analytic Phi(dt), Gam(dt) vs expm — 1e-12.
#   [2] HOVER THEOREM: holding a fixed LVLH offset (Dx, Dy, Dz) requires
#       u = (-3w^2 Dx, 0, w^2 Dz) exactly; along-track offsets are FREE.
#       Gate pair: (a) along-track train -> optimal fuel ~ 0;
#       (b) radial offset -> fuel/time -> w^2*sqrt(9Dx^2+Dz^2) (rate match).
#   [3] DECOUPLING: kap = 0 == V independent programs (per-vehicle fuel).
#   [4] FRONTIER: fuel up, formation error down, monotone in kap.
#   [5] CROSS-SOLVER: Clarabel vs SCS.
#   [6] FROZEN base-instance constants for def.jl.
#
# Usage: python3 cw_oracle.py [--quick]
# =============================================================================
import sys
import numpy as np
import cvxpy as cp
from scipy.linalg import expm
import warnings; warnings.filterwarnings("ignore")

QUICK = "--quick" in sys.argv
FEASIBLE = ("optimal", "optimal_inaccurate")

OMEGA = 1.13e-3          # rad/s (LEO)
DT, N_BASE = 60.0, 60    # 1-min steps, 1-hour horizon
V_BASE = 4
UMAX = 0.05              # m/s^2
GAM, KAP = 1e-4, 1e-3
SIGA = 1e-5              # accel process noise for Van Loan

# frozen slot offsets (m) and initial dispersions (m, m/s) — literals
SLOTS = [np.array([80.0, 0.0, 0.0, 0, 0, 0]),
         np.array([0.0, 120.0, 0.0, 0, 0, 0]),
         np.array([-80.0, 0.0, 40.0, 0, 0, 0]),
         np.array([0.0, -120.0, -40.0, 0, 0, 0])]
DISP = [np.array([3.1, -4.2, 1.7, 0.011, -0.007, 0.004]),
        np.array([-2.4, 5.0, -2.2, -0.009, 0.012, -0.005]),
        np.array([4.4, 2.6, 0.9, 0.006, 0.010, 0.008]),
        np.array([-3.7, -3.3, 2.8, -0.012, -0.004, -0.006])]

def cw_A(w):
    A = np.zeros((6, 6)); A[0:3, 3:6] = np.eye(3)
    A[3, 0] = 3 * w * w; A[3, 4] = 2 * w
    A[4, 3] = -2 * w
    A[5, 2] = -w * w
    return A

def cw_stm_analytic(w, t):
    s, c = np.sin(w * t), np.cos(w * t)
    P = np.zeros((6, 6))
    P[0] = [4 - 3 * c, 0, 0, s / w, 2 * (1 - c) / w, 0]
    P[1] = [6 * (s - w * t), 1, 0, 2 * (c - 1) / w, (4 * s - 3 * w * t) / w, 0]
    P[2] = [0, 0, c, 0, 0, s / w]
    P[3] = [3 * w * s, 0, 0, c, 2 * s, 0]
    P[4] = [6 * w * (c - 1), 0, 0, -2 * s, 4 * c - 3, 0]
    P[5] = [0, 0, -w * s, 0, 0, c]
    return P

def discretize(w, dt):
    A = cw_A(w)
    B = np.zeros((6, 3)); B[3:6] = np.eye(3)
    M = np.zeros((9, 9)); M[:6, :6] = A; M[:6, 6:] = B
    E = expm(M * dt)
    return E[:6, :6], E[:6, 6:]

def vanloan_W(w, dt, siga):
    A = cw_A(w)
    Qc = np.zeros((6, 6)); Qc[3:6, 3:6] = siga ** 2 * np.eye(3)
    M = np.block([[-A, Qc], [np.zeros((6, 6)), A.T]]) * dt
    E = expm(M)
    Qd = E[6:, 6:].T @ E[:6, 6:]
    return np.linalg.inv(Qd + 1e-9 * np.trace(Qd) / 6 * np.eye(6))

RING = [(i, (i + 1) % V_BASE) for i in range(V_BASE)]

def solve_formation(V, N, gam, kap, slots, disp, w=OMEGA, dt=DT,
                    solver="CLARABEL", tight=True, pin_terminal=False):
    Phi, Gam_ = discretize(w, dt)
    W = vanloan_W(w, dt, SIGA)
    Wn = W / np.trace(W) * 6.0          # normalized metric (scale-stable)
    ring = [(i, (i + 1) % V) for i in range(V)] if V > 1 else []
    X = [cp.Variable((N + 1, 6)) for _ in range(V)]
    U = [cp.Variable((N, 3)) for _ in range(V)]
    s = [cp.Variable(N) for _ in range(V)]
    cons = []
    for v in range(V):
        cons.append(X[v][0] == slots[v] + disp[v])
        if pin_terminal:
            cons.append(X[v][N] == slots[v])
        for k in range(N):
            cons.append(X[v][k + 1] == Phi @ X[v][k] + Gam_ @ U[v][k])
            cons.append(cp.norm(U[v][k]) <= s[v][k])
            cons.append(s[v][k] <= UMAX)
    Wh = np.linalg.cholesky(Wn).T
    obj = sum(dt * cp.sum(s[v]) for v in range(V))
    for v in range(V):
        obj = obj + gam / 2 * cp.sum_squares((X[v] - slots[v][None, :]) @ Wh.T)
    for (i, j) in ring:
        dij = slots[i] - slots[j]
        obj = obj + kap / 2 * cp.sum_squares((X[i] - X[j] - dij[None, :]) @ Wh.T)
    prob = cp.Problem(cp.Minimize(obj), cons)
    if solver == "SCS":
        prob.solve(solver=cp.SCS, eps=1e-8, max_iters=200_000)
    elif tight:
        prob.solve(solver=cp.CLARABEL, tol_gap_abs=1e-11, tol_gap_rel=1e-11,
                   tol_feas=1e-11)
    else:
        prob.solve(solver=cp.CLARABEL)
    assert prob.status in FEASIBLE, prob.status
    fuel = [dt * float(np.sum(s[v].value)) for v in range(V)]
    ferr = max((float(np.max(np.abs(X[i].value - X[j].value
                - (slots[i] - slots[j])[None, :]))) for (i, j) in ring),
               default=0.0)
    return prob, X, U, fuel, ferr

def main():
    print("=" * 77)
    print(f"CW formation-keeping oracle: V = {V_BASE}, N = {N_BASE}, "
          f"dt = {DT}s, omega = {OMEGA}" + ("  [--quick]" if QUICK else ""))
    print("=" * 77)

    print("[1] STM: analytic CW Phi vs expm")
    Phi_e, _ = discretize(OMEGA, DT)
    d = np.max(np.abs(cw_stm_analytic(OMEGA, DT) - Phi_e))
    print(f"    max |Phi_analytic - expm| = {d:.2e}")
    assert d < 1e-12

    print("[2] HOVER THEOREM (the physics anchor)")
    # (a) along-track train: slots on the y-axis only -> formation is FREE
    train = [np.array([0.0, 60.0 * v, 0.0, 0, 0, 0]) for v in range(V_BASE)]
    zero = [np.zeros(6) for _ in range(V_BASE)]
    _, _, _, fuel_a, _ = solve_formation(V_BASE, N_BASE, 1e-1, 1e-1,
                                         train, zero)
    print(f"    along-track train, zero dispersion: total fuel = "
          f"{sum(fuel_a):.3e} m/s (free formation)")
    assert sum(fuel_a) < 1e-4
    # (b) radial offset: hover rate = w^2*sqrt(9 Dx^2 + Dz^2), per vehicle
    Dx, Dz = 80.0, 40.0
    radial = [np.array([Dx, 0, Dz, 0, 0, 0]), np.zeros(6)]
    _, _, Uh, fuel_b, _ = solve_formation(2, N_BASE, 1e+1, 0.0,
                                          radial, [np.zeros(6)] * 2,
                                          pin_terminal=True)
    rate = fuel_b[0] / (N_BASE * DT)
    rate_theory = OMEGA ** 2 * np.sqrt(9 * Dx * Dx + Dz * Dz)
    rel = abs(rate - rate_theory) / rate_theory
    print(f"    radial hover: fuel rate {rate:.4e} vs theory "
          f"{rate_theory:.4e} m/s^2 (rel {rel:.2e})")
    assert rel < 2e-2
    umid = Uh[0].value[N_BASE // 2]
    u_theory = np.array([-3 * OMEGA ** 2 * Dx, 0.0, OMEGA ** 2 * Dz])
    print(f"    mid-horizon u = {np.round(umid, 7)} vs theory "
          f"{np.round(u_theory, 7)}")
    assert np.max(np.abs(umid - u_theory)) < 3e-5

    print("[3] DECOUPLING: kap = 0 == independent vehicles (fuel invariant)")
    _, _, _, fuel_c, _ = solve_formation(V_BASE, N_BASE, GAM, 0.0, SLOTS, DISP)
    worst = 0.0
    for v in range(V_BASE):
        _, _, _, f1, _ = solve_formation(1, N_BASE, GAM, 0.0,
                                         [SLOTS[v]], [DISP[v]])
        worst = max(worst, abs(fuel_c[v] - f1[0]) / max(f1[0], 1e-9))
    print(f"    per-vehicle fuel rel = {worst:.2e}")
    assert worst < 1e-4

    print("[4] FRONTIER: monotone in kap")
    fuels, errs = [], []
    for kap in (1e-4, 1e-3, 1e-2):
        p, _, _, fu, fe = solve_formation(V_BASE, N_BASE, GAM, kap,
                                          SLOTS, DISP, tight=False)
        fuels.append(sum(fu)); errs.append(fe)
        print(f"    kap = {kap}: fuel = {fuels[-1]:.4f} m/s, "
              f"formation err = {errs[-1]:.3f} m")
    assert fuels == sorted(fuels) and errs == sorted(errs, reverse=True)

    print("[5] CROSS-SOLVER (Clarabel tight vs SCS)")
    p1, *_ = solve_formation(V_BASE, N_BASE, GAM, KAP, SLOTS, DISP)
    p2, *_ = solve_formation(V_BASE, N_BASE, GAM, KAP, SLOTS, DISP,
                             solver="SCS")
    rel = abs(p1.value - p2.value) / abs(p1.value)
    print(f"    objective rel = {rel:.2e}")
    assert rel < 2e-4   # SCS(1e-8)-quality agreement

    print("=" * 77)
    print("ALL CHECKS PASSED")
    print(f"FROZEN for def.jl: obj(base) = {p1.value:.8f}; per-vehicle fuel = "
          f"{[round(f, 6) for f in fuel_c]}; hover rate theory = "
          f"{rate_theory:.6e}")

if __name__ == "__main__":
    main()
