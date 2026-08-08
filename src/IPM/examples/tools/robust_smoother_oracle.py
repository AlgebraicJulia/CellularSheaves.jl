# =============================================================================
# robust_smoother_oracle.py — self-contained validation for E11 (robust
# fixed-interval smoothing under burst-corrupted measurements)
#
# Linear-Gaussian state-space model over a horizon T:
#
#     x_{t+1} = A x_t + w_t,  w_t ~ N(0, W),   x_0 ~ N(mu0, P0)
#     y_t     = C x_t + v_t,  v_t ~ N(0, V),   t = 0..T,
#
# with A = expm(Ac dt) for a damped spring-mass chain (the matrix exponential
# is dense BY NATURE, even though Ac is tridiagonal), W from the Van Loan
# (1978) block-exponential (dense by nature for the same reason), C dense
# sensing rows, and V = sig^2 (F F' + I) a dense factor model (common-mode
# sensor drift, cf. E9's factor Sigma). Whitened MAP smoothing with the
# SQUARE-ROOT measurement loss (E7's loss, dynamic edition):
#
#     min_x  (1/2)||L0^{-1}(x_0 - mu0)||^2
#          + (1/2) sum_t ||Lw^{-1}(x_{t+1} - A x_t)||^2
#          + lam   sum_t ||Lv^{-1}(y_t - C x_t)||        (norms, NOT squares)
#
# The sqrt loss bounds each timestep's pull on the trajectory at lam (the
# gradient of ||r|| is a unit vector), so measurement BURSTS — contiguous
# windows where every sensor is corrupted (telemetry dropouts, EMI) — cannot
# drag the estimate: the chain coasts across the outage on dynamics alone.
# With the quadratic loss (lam -> RTS limit) the same program IS the Kalman
# smoother, which is what makes this family an analytic oracle. Model and
# robust-loss lineage: Rauch-Tung-Striebel (1965); Aravkin, Burke, Ljung,
# Lozano & Pillonetto, "Generalized Kalman smoothing" (Automatica 2017).
#
# Conic form (the e11.jl layout, dense-Q + SOC — a habitat cell no prior
# example occupies: E9 is dense-Q + pow, E7 is SOC with Q = 0). Vertices
# t = 1..T are OVERLAPPING PAIR STALKS z_t = (x_{t-1}, x_t) in R^{2n}
# (CofreeCone) carrying the dense process Gram
#     G = [A'Om A, -A'Om; -Om A, Om],  Om = W^{-1}
# (+ prior P0^{-1} on the first block of z_1) — B1's overlapping-patch
# discipline for coupled quadratics, with n-row agreement edges between
# consecutive pairs. Each timestep hangs one SOC stalk (u_t, w_t, s_t) of
# dim m+2 off the chain: m dense tie rows w_t = M x_t - b_t with
# M = Lv^{-1} C (dense by nature, nonzero rhs — E7's pattern) and a one-row
# pin s_t = delta; the objective adds lam * u_t. delta = 0 is the flagship
# sqrt loss; delta > 0 with lam = delta is the pseudo-Huber deformation
#     psi_delta(r) = delta (||(r, delta)|| - delta) -> ||r||^2/2  (delta -> inf)
# whose limit is EXACTLY the RTS smoother.
#
# Checks:
#   [1] DENSE BY NATURE: density of A, W, V, M at 1e-8 of the max entry.
#   [2] ANALYTIC (RTS): three independent algorithms agree — RTS recursion,
#       banded information-form solve, conic quadratic MAP. Machine to
#       solver precision. (The suite's third analytic oracle.)
#   [3] FORMULATION CERTIFICATE: the pair-stalk + agreement + SOC-tie
#       assembly (exactly what e11.jl builds) equals the native norm-sum
#       formulation; plus a random-point algebra identity for the pair-Gram
#       split (machine precision).
#   [4] DEFORMATION: pseudo-Huber x(delta) -> RTS at measured rate ~delta^-2.
#   [5] ESCALATION: burst corruption x1 -> x128; robust clean-region state
#       RMSE saturates while RTS blows up; clean-data price is a few percent.
#   [6] SUPPORT RECOVERY + INFLUENCE CAP: the largest whitened residuals
#       identify the corrupted steps exactly, with a wide margin; the
#       quadratic loss lets corrupted steps pull with ||r|| >> lam.
#   [7] KKT: subgradient-selection stationarity residual at the optimum
#       (interpolation-active steps ||w_t|| = 0 get a least-squares
#       subgradient projected to the unit ball, E6's discipline).
#   [8] CLASSICAL (DARE): the influence of a boundary data perturbation on
#       interior states decays geometrically at the closed-loop rate
#       rho((I - K_inf C) A) predicted by the discrete algebraic Riccati
#       equation — measured vs predicted slope.
#   [9] CROSS-SOLVER: Clarabel vs SCS on the flagship.
#
# Usage: python3 robust_smoother_oracle.py [--quick]
# =============================================================================

import sys
import numpy as np
import cvxpy as cp
from scipy.linalg import expm, cholesky, solve_discrete_are

QUICK = "--quick" in sys.argv
SEED_SYS, SEED_SIM = 7, 3
LAM = 1.5
TIGHT = dict(solver=cp.CLARABEL, tol_gap_abs=1e-10, tol_gap_rel=1e-10,
             tol_feas=1e-10, max_iter=2000)


def conic_solve(pr, **kw):
    try:
        pr.solve(**(TIGHT | kw))
        assert pr.status in ("optimal", "optimal_inaccurate")
    except Exception:
        try:
            pr.solve(solver=cp.CLARABEL, **kw)
            assert pr.status in ("optimal", "optimal_inaccurate")
        except Exception:
            pr.solve(solver=cp.SCS, eps=1e-9, max_iters=200_000)
            assert pr.status in ("optimal", "optimal_inaccurate"), pr.status
    return pr.value


# -----------------------------------------------------------------------------
# Instance: damped spring-mass chain, Van Loan discretization, factor-noise V
# -----------------------------------------------------------------------------

def make_instance(K=6, m=4, dt=0.5, k=4.0, kg=1.0, damp=0.4, q=0.5,
                  sig=0.12, sigf=0.7, seed=SEED_SYS):
    n = 2 * K
    rng = np.random.default_rng(seed)
    Kmat = (np.diag((2 * k + kg) * np.ones(K))
            - np.diag(k * np.ones(K - 1), 1) - np.diag(k * np.ones(K - 1), -1))
    Ac = np.block([[np.zeros((K, K)), np.eye(K)], [-Kmat, -damp * np.eye(K)]])
    Wc = np.zeros((n, n)); Wc[K:, K:] = q * np.eye(K)
    M = np.block([[-Ac, Wc], [np.zeros((n, n)), Ac.T]]) * dt  # Van Loan (1978)
    F = expm(M)
    A = F[n:, n:].T
    W = A @ F[:n, n:]; W = 0.5 * (W + W.T)
    C = rng.normal(size=(m, n)) / np.sqrt(n)                  # dense sensing
    Fn = sigf * rng.normal(size=(m, 2))
    V = sig ** 2 * (Fn @ Fn.T + np.eye(m))                    # dense factor V
    return dict(n=n, m=m, A=A, W=W, C=C, V=V,
                P0=0.5 * np.eye(n), mu0=np.zeros(n))


def whiten(inst):
    Lw = cholesky(inst["W"], lower=True)
    Lv = cholesky(inst["V"], lower=True)
    L0 = cholesky(inst["P0"], lower=True)
    return Lw, Lv, L0, np.linalg.solve(Lv, inst["C"])


def simulate(inst, T, nburst=3, burst_len=8, corrupt_scale=32.0, seed=SEED_SIM):
    n, m = inst["n"], inst["m"]
    A, C, mu0 = inst["A"], inst["C"], inst["mu0"]
    Lw, Lv, L0, _ = whiten(inst)
    rng = np.random.default_rng(seed)
    x = np.zeros((T + 1, n)); y = np.zeros((T + 1, m))
    x[0] = mu0 + L0 @ rng.normal(size=n)
    for t in range(T):
        x[t + 1] = A @ x[t] + Lw @ rng.normal(size=n)
    corrupted = np.zeros(T + 1, dtype=bool)
    starts = np.sort(rng.choice(np.arange(5, T - burst_len - 5),
                                size=nburst, replace=False))
    for i in range(1, nburst):
        starts[i] = max(starts[i], starts[i - 1] + burst_len + 4)
    for s in starts:
        corrupted[s:s + burst_len] = True
    for t in range(T + 1):
        s = corrupt_scale if corrupted[t] else 1.0
        y[t] = C @ x[t] + s * (Lv @ rng.normal(size=m))
    return x, y, corrupted


# -----------------------------------------------------------------------------
# Reference algorithms: RTS recursion and banded information form
# -----------------------------------------------------------------------------

def rts(inst, y):
    n = inst["n"]; T = y.shape[0] - 1
    A, W, C, V, mu0, P0 = (inst[k] for k in ("A", "W", "C", "V", "mu0", "P0"))
    xf = np.zeros((T + 1, n)); Pf = np.zeros((T + 1, n, n))
    xp = np.zeros((T + 1, n)); Pp = np.zeros((T + 1, n, n))
    xpred, Ppred = mu0.copy(), P0.copy()
    for t in range(T + 1):
        xp[t], Pp[t] = xpred, Ppred
        S = C @ Ppred @ C.T + V
        K = np.linalg.solve(S, C @ Ppred).T
        xf[t] = xpred + K @ (y[t] - C @ xpred)
        Pf[t] = (np.eye(n) - K @ C) @ Ppred; Pf[t] = 0.5 * (Pf[t] + Pf[t].T)
        xpred = A @ xf[t]; Ppred = A @ Pf[t] @ A.T + W
    xs = xf.copy()
    for t in range(T - 1, -1, -1):
        G = np.linalg.solve(Pp[t + 1], A @ Pf[t]).T
        xs[t] = xf[t] + G @ (xs[t + 1] - xp[t + 1])
    return xs


def info_form(inst, y):
    n = inst["n"]; T = y.shape[0] - 1; N = (T + 1) * n
    A, W, C, V, mu0, P0 = (inst[k] for k in ("A", "W", "C", "V", "mu0", "P0"))
    Om, Ov, O0 = np.linalg.inv(W), np.linalg.inv(V), np.linalg.inv(P0)
    J = np.zeros((N, N)); h = np.zeros(N)
    sl = lambda t: slice(t * n, (t + 1) * n)
    J[sl(0), sl(0)] += O0; h[sl(0)] += O0 @ mu0
    for t in range(T):
        J[sl(t), sl(t)] += A.T @ Om @ A
        J[sl(t + 1), sl(t + 1)] += Om
        J[sl(t), sl(t + 1)] += -A.T @ Om
        J[sl(t + 1), sl(t)] += -Om @ A
    for t in range(T + 1):
        J[sl(t), sl(t)] += C.T @ Ov @ C
        h[sl(t)] += C.T @ Ov @ y[t]
    return np.linalg.solve(J, h).reshape(T + 1, n), J, h


# -----------------------------------------------------------------------------
# Conic solvers
# -----------------------------------------------------------------------------

def smooth_grad(inst, X):
    """Gradient of prior + process quadratics at X (for the KKT check)."""
    T = X.shape[0] - 1
    A, mu0 = inst["A"], inst["mu0"]
    Om, O0 = np.linalg.inv(inst["W"]), np.linalg.inv(inst["P0"])
    g = np.zeros_like(X)
    g[0] += O0 @ (X[0] - mu0)
    for t in range(T):
        r = Om @ (X[t + 1] - A @ X[t])
        g[t] += -A.T @ r; g[t + 1] += r
    return g


def solve_native(inst, y, lam=None, delta=None, solver_kw=None):
    """Direct formulation: quadratic (both None), sqrt loss (lam), or
    pseudo-Huber (delta)."""
    T = y.shape[0] - 1
    Lw, Lv, L0, M = whiten(inst)
    A, mu0 = inst["A"], inst["mu0"]
    X = cp.Variable((T + 1, inst["n"]))
    Li0, Liw = np.linalg.inv(L0), np.linalg.inv(Lw)
    obj = 0.5 * cp.sum_squares(Li0 @ (X[0] - mu0))
    for t in range(T):
        obj = obj + 0.5 * cp.sum_squares(Liw @ (X[t + 1] - A @ X[t]))
    for t in range(T + 1):
        r = M @ X[t] - np.linalg.solve(Lv, y[t])
        if delta is not None:
            obj = obj + delta * (cp.norm(cp.hstack([r, delta])) - delta)
        elif lam is not None:
            obj = obj + lam * cp.norm(r, 2)
        else:
            obj = obj + 0.5 * cp.sum_squares(r)
    pr = cp.Problem(cp.Minimize(obj))
    conic_solve(pr, **(solver_kw or {}))
    return X.value, pr.value


def pair_gram(inst):
    A = inst["A"]; Om = np.linalg.inv(inst["W"])
    return np.block([[A.T @ Om @ A, -A.T @ Om], [-Om @ A, Om]])


def solve_sheaf_assembly(inst, y, lam, delta=0.0, solver_kw=None):
    """EXACTLY the e11.jl layout: T pair stalks (CofreeCone, 2n) with dense
    Gram Q blocks, n-row agreement edges, T+1 SOC stalks (u, w, s) with m
    dense tie rows and a one-row pin s = delta."""
    n, m, T = inst["n"], inst["m"], y.shape[0] - 1
    Lw, Lv, L0, M = whiten(inst)
    mu0 = inst["mu0"]
    O0 = np.linalg.inv(inst["P0"])
    G = pair_gram(inst)
    Z = [cp.Variable(2 * n) for _ in range(T)]
    U = [cp.Variable(1 + m + 1) for _ in range(T + 1)]
    cons, obj = [], 0
    for t in range(T):
        Qt = G + (np.pad(O0, ((0, n), (0, n))) if t == 0 else 0)
        obj = obj + 0.5 * cp.quad_form(Z[t], cp.psd_wrap(Qt))
    obj = obj - (O0 @ mu0) @ Z[0][:n]
    for t in range(T + 1):
        obj = obj + lam * U[t][0]
        cons.append(cp.SOC(U[t][0], U[t][1:]))
    for t in range(T - 1):
        cons.append(Z[t][n:] == Z[t + 1][:n])                 # agreement
    for t in range(T + 1):
        xt = Z[t][:n] if t < T else Z[T - 1][n:]
        cons.append(U[t][1:1 + m] == M @ xt - np.linalg.solve(Lv, y[t]))
        cons.append(U[t][1 + m] == delta)                     # pin
    pr = cp.Problem(cp.Minimize(obj), cons)
    conic_solve(pr, **(solver_kw or {}))
    X = np.zeros((T + 1, n))
    for t in range(T):
        X[t] = Z[t].value[:n]
    X[T] = Z[T - 1].value[n:]
    return X, pr.value


# -----------------------------------------------------------------------------
# Checks
# -----------------------------------------------------------------------------

def density(Mx, tol=1e-8):
    return float(np.mean(np.abs(Mx) > tol * np.abs(Mx).max()))


def check_dense_by_nature(inst):
    _, _, _, M = whiten(inst)
    d = [density(inst["A"]), density(inst["W"]), density(inst["V"]), density(M)]
    print(f"[1] dense by nature: density@1e-8  A {d[0]:.3f}  W {d[1]:.3f}  "
          f"V {d[2]:.3f}  M=Lv^-1 C {d[3]:.3f}")
    assert min(d) > 0.95, d


def check_analytic(inst, T):
    x, y, _ = simulate(inst, T, corrupt_scale=1.0)
    xs = rts(inst, y)
    xi, _, _ = info_form(inst, y)
    e1 = np.abs(xs - xi).max()
    Xq, _ = solve_native(inst, y)
    e2 = np.abs(Xq - xs).max()
    print(f"[2] ANALYTIC (RTS): recursion vs information form {e1:.2e}; "
          f"conic quadratic MAP vs RTS {e2:.2e}")
    assert e1 < 1e-10 and e2 < 1e-8, (e1, e2)
    return y, xs


def check_formulation(inst, y):
    # random-point algebra identity for the pair-Gram split
    T = y.shape[0] - 1
    rng = np.random.default_rng(0)
    Lw, Lv, L0, M = whiten(inst)
    Li0, Liw = np.linalg.inv(L0), np.linalg.inv(Lw)
    G = pair_gram(inst); O0 = np.linalg.inv(inst["P0"])
    worst = 0.0
    for _ in range(10):
        X = rng.normal(size=(T + 1, inst["n"]))
        direct = 0.5 * np.sum((Li0 @ (X[0] - inst["mu0"])) ** 2) + 0.5 * sum(
            np.sum((Liw @ (X[t + 1] - inst["A"] @ X[t])) ** 2) for t in range(T))
        z = np.hstack([X[:-1], X[1:]])
        split = 0.5 * sum(z[t] @ G @ z[t] for t in range(T)) \
            + 0.5 * X[0] @ O0 @ X[0] - (O0 @ inst["mu0"]) @ X[0] \
            + 0.5 * inst["mu0"] @ O0 @ inst["mu0"]
        worst = max(worst, abs(direct - split) / max(abs(direct), 1.0))
    Xr, vr = solve_native(inst, y, lam=LAM)
    Xs, vs = solve_sheaf_assembly(inst, y, lam=LAM)
    dx, dv = np.abs(Xs - Xr).max(), abs(vs - vr) / abs(vr)
    print(f"[3] formulation certificate: pair-Gram identity {worst:.2e}; "
          f"sheaf assembly vs native ||dx||_inf {dx:.2e} (rel dv {dv:.2e})")
    assert worst < 1e-12 and dx < 1e-4, (worst, dx)
    return Xr


def check_deformation(inst, T):
    x, y, _ = simulate(inst, T, corrupt_scale=1.0)     # clean: residuals O(1)
    xs = rts(inst, y)
    deltas = [2.0, 4.0, 8.0, 16.0]
    errs = []
    for d in deltas:
        Xd, _ = solve_native(inst, y, delta=d)
        errs.append(np.abs(Xd - xs).max())
    rate = np.polyfit(np.log(deltas), np.log(errs), 1)[0]
    print(f"[4] deformation: ||x(delta) - RTS||_inf = "
          f"{errs[0]:.2e} -> {errs[-1]:.2e} over delta {deltas[0]:.0f} -> "
          f"{deltas[-1]:.0f}; measured rate delta^{rate:.2f} (predicted -2)")
    assert rate < -1.6 and errs[-1] < 1e-2, (rate, errs)


def check_escalation(inst, T):
    scales = [1.0, 32.0] if QUICK else [1.0, 8.0, 32.0, 128.0]
    rows = []
    for s in scales:
        x, y, corr = simulate(inst, T, corrupt_scale=s)
        xs = rts(inst, y)
        Xr, _ = solve_native(inst, y, lam=LAM)
        cl = ~corr
        rmse = lambda Z, msk: np.sqrt(np.mean((Z - x)[msk] ** 2))
        rows.append((s, rmse(xs, cl), rmse(Xr, cl),
                     np.sqrt(np.mean((xs - x) ** 2)), np.sqrt(np.mean((Xr - x) ** 2))))
        print(f"[5]   x{s:<4.0f} clean-region RMSE: RTS {rows[-1][1]:.4f} "
              f"robust {rows[-1][2]:.4f}   overall: RTS {rows[-1][3]:.4f} "
              f"robust {rows[-1][4]:.4f}")
    s1, sN = rows[0], rows[-1]
    price = s1[2] / s1[1] - 1
    print(f"[5] escalation: robust clean-region saturates "
          f"{s1[2]:.3f} -> {sN[2]:.3f}; RTS blows up {s1[1]:.3f} -> {sN[1]:.3f}; "
          f"clean price {price:+.1%}")
    assert sN[2] < 0.5 * sN[1], (sN[2], sN[1])       # robust beats RTS off-burst
    assert sN[2] < 1.6 * s1[2], (sN[2], s1[2])       # saturation
    assert price < 0.10, price                       # modest clean price
    assert sN[4] < 0.6 * sN[3], (sN[4], sN[3])       # headline: overall RMSE


def check_support_influence(inst, y, corr, Xr):
    _, Lv, _, M = whiten(inst)
    T = y.shape[0] - 1
    rn = np.array([np.linalg.norm(M @ Xr[t] - np.linalg.solve(Lv, y[t]))
                   for t in range(T + 1)])
    k = int(corr.sum())
    top = set(np.argsort(rn)[-k:])
    tp = len(top & set(np.where(corr)[0]))
    margin = rn[corr].min() / rn[~corr].max()
    print(f"[6] support recovery: {tp}/{k} of the top-{k} residuals are "
          f"corrupted steps (margin x{margin:.1f}); influence cap: quadratic "
          f"pull up to ||r|| = {rn.max():.1f} vs lam = {LAM} "
          f"({rn.max() / LAM:.0f}x)")
    assert tp == k and margin > 1.5, (tp, k, margin)


def check_kkt(inst, y, Xr):
    _, Lv, _, M = whiten(inst)
    T = y.shape[0] - 1
    g = smooth_grad(inst, Xr)
    nact = 0
    for t in range(T + 1):
        w = M @ Xr[t] - np.linalg.solve(Lv, y[t])
        nw = np.linalg.norm(w)
        if nw > 1e-7:
            g[t] += LAM * (M.T @ (w / nw))
        else:                                # subgradient selection (E6)
            nact += 1
            s, *_ = np.linalg.lstsq(LAM * M.T, -g[t], rcond=None)
            s = s / max(1.0, np.linalg.norm(s))
            g[t] += LAM * (M.T @ s)
    res = np.abs(g).max()
    print(f"[7] KKT: stationarity residual {res:.2e} "
          f"({nact} interpolation-active steps, subgradient-selected)")
    assert res < 1e-3, res


def check_dare_decay(inst, y):
    n, m = inst["n"], inst["m"]
    T = y.shape[0] - 1
    _, J, _ = info_form(inst, y)
    Ov = np.linalg.inv(inst["V"])
    dh = np.zeros((T + 1) * n); dh[:n] = inst["C"].T @ Ov @ np.ones(m)
    dx = np.linalg.solve(J, dh).reshape(T + 1, n)
    mag = np.linalg.norm(dx, axis=1)
    t1, t2 = 5, min(60, T - 5)
    slope = np.polyfit(np.arange(t1, t2), np.log(mag[t1:t2]), 1)[0]
    P = solve_discrete_are(inst["A"].T, inst["C"].T, inst["W"], inst["V"])
    K = P @ inst["C"].T @ np.linalg.inv(inst["C"] @ P @ inst["C"].T + inst["V"])
    rho = max(abs(np.linalg.eigvals((np.eye(n) - K @ inst["C"]) @ inst["A"])))
    ratio = slope / np.log(rho)
    print(f"[8] CLASSICAL (DARE): boundary-perturbation decay e^{slope:.4f}"
          f"/step vs closed-loop rho((I-KC)A) = {rho:.4f} "
          f"(log {np.log(rho):.4f}); ratio {ratio:.3f}")
    assert 0.9 < ratio < 1.1, ratio


def check_cross_solver(inst, y):
    Xc, vc = solve_native(inst, y, lam=LAM)
    T = y.shape[0] - 1
    Xs2, vs2 = solve_native(inst, y, lam=LAM,
                            solver_kw=dict(solver=cp.SCS, eps=1e-9,
                                           max_iters=200_000))
    dx = np.abs(Xc - Xs2).max()
    print(f"[9] cross-solver: Clarabel vs SCS ||dx||_inf {dx:.2e} "
          f"(rel dv {abs(vc - vs2) / abs(vc):.2e})")
    assert dx < 1e-4, dx


# -----------------------------------------------------------------------------

if __name__ == "__main__":
    T = 120 if QUICK else 200
    inst = make_instance()
    print(f"instance: n = {inst['n']} (spring-mass chain), m = {inst['m']} "
          f"dense sensors, T = {T}, lam = {LAM}")
    check_dense_by_nature(inst)
    check_analytic(inst, T)
    x, y, corr = simulate(inst, T, corrupt_scale=32.0)
    Xr = check_formulation(inst, y)
    check_deformation(inst, 120)
    check_escalation(inst, T)
    check_support_influence(inst, y, corr, Xr)
    check_kkt(inst, y, Xr)
    check_dare_decay(inst, y)
    check_cross_solver(inst, y)
    print("ALL CHECKS PASSED")
