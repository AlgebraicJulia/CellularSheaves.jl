# =============================================================================
# almgren_pow_oracle.py — self-contained validation for E9 (optimal execution)
#
# Multi-period, multi-asset optimal execution with 3/2-power market impact
# (power-cone habitat with genuinely dense per-period Q):
#
#     min_{x_1..x_{T-1}}  γ Σ_{t=1}^{T-1} x_t' Σ x_t
#                         + Σ_{t=1}^{T} Σ_i η_i |v_{t,i}|^{3/2}
#     with  v_t = x_{t-1} - x_t,   x_0 = X given,   x_T = 0,
#
# Σ = B F B' + D a dense factor-model covariance. The 3/2 exponent is the
# empirically established temporary-impact power law (Almgren, Thum,
# Hauptmann & Li 2005); the multi-period convex-optimization treatment with
# this exact term is Boyd et al. (2017), and the continuous-time power-law
# analysis is Almgren (2003). Conic form: one PowerCone(2/3) triple
# (s, 1, v) per asset per trade, s >= |v|^{3/2}; risk enters through dense
# per-period Q = 2γΣ — the dense-Q habitat (cf. E3/E4), which E6-E8 lack.
#
# Checks:
#   [1] ANALYTIC: continuous-time single-asset reference. The Euler-Lagrange
#       problem min ∫ η|ẋ|^{3/2} + κ x² dt, x(0)=X, x(T)=0, conserves
#       E = (η/2)|ẋ|^{3/2} - κ x², reducing the BVP to a scalar bisection;
#       the discrete conic solution (η_dt = η/√Δt, γ_dt = κΔt) must converge
#       to it under Δt-refinement at a measured rate.
#   [2] LIMITS: η/γ -> ∞ recovers TWAP (equal trades, by convexity of
#       |·|^{3/2}); η -> 0 recovers immediate liquidation. Both analytic.
#   [3] DECOUPLING: with diagonal Σ the multi-asset problem separates into
#       n single-asset problems exactly; with factor coupling it must not.
#   [4] CROSS-FAMILY: conic solve (Clarabel power cones) vs smooth
#       first-order solve (L-BFGS on the C^1 objective) vs SCS.
#   [5] KKT: interior stationarity residual
#       2γΣx_t + (3/2)η∘(|v_{t+1}|^{1/2}sgn v_{t+1} - |v_t|^{1/2}sgn v_t) = 0.
#   [6] FRONTIER: γ-sweep gives monotone risk/impact trade-off.
#
# Usage: python3 almgren_pow_oracle.py [--quick]
# =============================================================================

import sys
import numpy as np
import scipy.integrate as si
import scipy.optimize as so
import cvxpy as cp

QUICK = "--quick" in sys.argv
SEED = 7
TIGHT = dict(solver=cp.CLARABEL, tol_gap_abs=1e-10, tol_gap_rel=1e-10,
             tol_feas=1e-10, max_iter=2000)


def conic_solve(prob):
    try:
        prob.solve(**TIGHT)
        assert prob.status in ("optimal", "optimal_inaccurate")
    except Exception:
        prob.solve(solver=cp.SCS, eps=1e-8, max_iters=100_000)
        assert prob.status in ("optimal", "optimal_inaccurate"), prob.status
    return prob.value


# -----------------------------------------------------------------------------
# Instance
# -----------------------------------------------------------------------------

def make_instance(n=8, kf=3, seed=SEED):
    rng = np.random.default_rng(seed)
    B = rng.normal(size=(n, kf)) / np.sqrt(kf)
    F = np.diag(rng.uniform(0.5, 1.5, kf))
    D = np.diag(rng.uniform(0.1, 0.3, n))
    Sigma = B @ F @ B.T + D
    Sigma = 0.5 * (Sigma + Sigma.T)
    X = rng.uniform(0.5, 2.0, n)
    eta = 60.0 * rng.uniform(0.5, 2.0, n)   # mid-regime: dev-from-TWAP ~ 0.3
    return Sigma, X, eta


# -----------------------------------------------------------------------------
# Objective, gradient, solvers
# -----------------------------------------------------------------------------

def trades(Xm, X0):
    """v_t = x_{t-1} - x_t for t=1..T, with x_0 = X0 and x_T = 0."""
    path = np.vstack([X0, Xm, np.zeros_like(X0)])
    return -np.diff(path, axis=0)


def objective(Xm, X0, Sigma, eta, gamma):
    v = trades(Xm, X0)
    risk = gamma * np.einsum("ti,ij,tj->", Xm, Sigma, Xm)
    impact = (np.abs(v) ** 1.5 @ eta).sum()
    return risk + impact


def grad(Xm, X0, Sigma, eta, gamma):
    v = trades(Xm, X0)
    dv = 1.5 * np.sqrt(np.abs(v)) * np.sign(v) * eta        # d|v|^{3/2}/dv
    return 2.0 * gamma * (Xm @ Sigma) + dv[1:] - dv[:-1]    # x_t hits v_t, v_{t+1}


def solve_smooth(X0, Sigma, eta, gamma, T, x0=None):
    n = len(X0)
    if x0 is None:                                          # TWAP warm start
        x0 = np.outer(1.0 - np.arange(1, T) / T, X0)
    r = so.minimize(lambda z: objective(z.reshape(T - 1, n), X0, Sigma, eta, gamma),
                    x0.ravel(), method="L-BFGS-B",
                    jac=lambda z: grad(z.reshape(T - 1, n), X0, Sigma, eta, gamma).ravel(),
                    options=dict(maxiter=50_000, ftol=1e-15, gtol=1e-12))
    return r.x.reshape(T - 1, n), r.fun


def solve_conic(X0, Sigma, eta, gamma, T, solver=None):
    n = len(X0)
    Xm = cp.Variable((T - 1, n))
    path = cp.vstack([X0.reshape(1, n), Xm, np.zeros((1, n))])
    v = path[:-1] - path[1:]
    s = cp.Variable((T, n))
    cons = [cp.constraints.PowCone3D(cp.vec(s, order="C"),
                                     np.ones(T * n),
                                     cp.vec(v, order="C"),
                                     2.0 / 3.0)]
    L = np.linalg.cholesky(Sigma)
    obj = gamma * cp.sum_squares(Xm @ L) + cp.sum(s @ eta)
    prob = cp.Problem(cp.Minimize(obj), cons)
    if solver == "SCS":
        prob.solve(solver=cp.SCS, eps=1e-8, max_iters=100_000)
        assert prob.status in ("optimal", "optimal_inaccurate")
        f = prob.value
    else:
        f = conic_solve(prob)
    return Xm.value, f


# -----------------------------------------------------------------------------
# Continuous-time analytic reference (single asset)
# -----------------------------------------------------------------------------
# min ∫_0^T η|ẋ|^{3/2} + κ x² dt, x(0)=X, x(T)=0. Autonomous Lagrangian =>
# conserved E = ẋ L_ẋ − L = (η/2)|ẋ|^{3/2} − κx², so along the (monotone)
# liquidation path |ẋ| = ((2/η)(κx² + E))^{2/3}, and the horizon condition
#     T(E) = ∫_0^X ((2/η)(κ x² + E))^{-2/3} dx
# is strictly decreasing in E with T(0+) = ∞, T(∞) = 0: bisect for E.

def horizon_of_E(E, X, eta, kappa):
    # x = sqrt(E/k) sinh(u) removes the boundary layer at x=0 analytically:
    # T(E) = (eta/2E)^{2/3} sqrt(E/k) * int_0^U cosh(u)^{-1/3} du,
    # U = asinh(X sqrt(k/E)); the integrand is smooth, <= 1, ~ e^{-u/3}.
    U = np.arcsinh(X * np.sqrt(kappa / E))
    val, _ = si.quad(lambda u: np.cosh(u) ** (-1.0 / 3.0), 0.0, U, limit=200)
    return (eta / (2.0 * E)) ** (2.0 / 3.0) * np.sqrt(E / kappa) * val


def continuous_reference(X, eta, kappa, T):
    lo, hi = 1e-14, 1.0
    while horizon_of_E(hi, X, eta, kappa) > T:
        hi *= 4.0
    while horizon_of_E(lo, X, eta, kappa) < T:
        lo /= 4.0
    for _ in range(90):
        mid = np.sqrt(lo * hi)
        (lo, hi) = (mid, hi) if horizon_of_E(mid, X, eta, kappa) > T else (lo, mid)
    E = np.sqrt(lo * hi)

    def rhs(t, x):
        return -(((2.0 / eta) * (kappa * x[0] ** 2 + E)) ** (2.0 / 3.0))

    sol = si.solve_ivp(rhs, (0.0, T), [X], dense_output=True,
                       rtol=1e-12, atol=1e-14, max_step=T / 200)
    xfun = lambda t: float(np.clip(sol.sol(t)[0], 0.0, None))
    speed = lambda t: ((2.0 / eta) * (kappa * xfun(t) ** 2 + E)) ** (2.0 / 3.0)
    cost, _ = si.quad(lambda t: eta * speed(t) ** 1.5 + kappa * xfun(t) ** 2,
                      0.0, T, limit=400)
    return E, xfun, cost


# -----------------------------------------------------------------------------
# Checks
# -----------------------------------------------------------------------------

def main():
    n, T = 8, 64
    gamma = 1.0
    Sigma, X0, eta = make_instance(n)
    print(f"n={n} T={T} gamma={gamma}: tr(Sigma)={np.trace(Sigma):.2f} "
          f"(dense factor model), X in [{X0.min():.2f},{X0.max():.2f}]")

    # [1] analytic continuous-time reference, Δt-refinement -------------------
    Xs, etas, kap, Th = 1.5, 1.0, 1.0, 4.0        # single-asset scenario
    E, xfun, c_cont = continuous_reference(Xs, etas, kap, Th)
    print(f"[1] continuous reference: E={E:.6e}, cost={c_cont:.8f}")
    errs, costs, Ns = [], [], ([8, 16, 32] if QUICK else [8, 16, 32, 64, 128])
    for Np in Ns:
        dt = Th / Np
        Xm, f = solve_conic(np.array([Xs]), np.array([[kap / gamma]]) * dt / gamma,
                            np.array([etas / np.sqrt(dt)]), gamma, Np)
        tgrid = dt * np.arange(1, Np)
        errs.append(float(np.max(np.abs(Xm.ravel() - [xfun(t) for t in tgrid]))))
        costs.append(f)
    slope = np.polyfit(np.log(Ns), np.log(errs), 1)[0]
    print("    sup|x_dt - x(t)|:",
          "  ".join(f"N={N}: {e:.2e}" for N, e in zip(Ns, errs)),
          f" (rate Δt^{-slope:.2f})")
    print("    cost:", "  ".join(f"{c:.6f}" for c in costs),
          f"-> continuous {c_cont:.6f}")
    assert errs[-1] < errs[0] / 4 and -slope > 0.5
    assert abs(costs[-1] - c_cont) < 3 * abs(costs[0] - c_cont) / (Ns[-1] / Ns[0])

    # [2] analytic limits -------------------------------------------------------
    twap = np.outer(1.0 - np.arange(1, T) / T, X0)
    devs = []
    for scale in (1e0, 1e2, 1e4):
        Xm, _ = solve_conic(X0, Sigma, scale * eta, gamma, T)
        devs.append(float(np.max(np.abs(Xm - twap)) / X0.max()))
    Xm0, _ = solve_conic(X0, Sigma, 1e-8 * eta, gamma, T)
    dump = float(np.max(np.abs(Xm0)) / X0.max())
    print(f"[2] limits: |x - TWAP|/|X| at eta x1/x100/x1e4: "
          f"{devs[0]:.3f} / {devs[1]:.4f} / {devs[2]:.5f}; "
          f"immediate-dump residual at eta->0: {dump:.2e}")
    assert devs[2] < devs[1] < devs[0] and devs[2] < 5e-3 and dump < 1e-2
    assert 0.1 < devs[0] < 0.8, "reference regime not mid: retune eta/gamma"

    # [3] decoupling -------------------------------------------------------------
    Ddiag = np.diag(np.diag(Sigma))
    Xm_d, _ = solve_conic(X0, Ddiag, eta, gamma, T)
    Xm_stack = np.column_stack([
        solve_conic(X0[i:i + 1], Ddiag[i:i + 1, i:i + 1], eta[i:i + 1],
                    gamma, T)[0].ravel() for i in range(n)])
    dd = float(np.max(np.abs(Xm_d - Xm_stack)))
    Xm_full, f_full = solve_conic(X0, Sigma, eta, gamma, T)
    coupling = float(np.max(np.abs(Xm_full - Xm_d)))
    print(f"[3] decoupling: diagonal-Sigma multi-asset == stacked 1-asset "
          f"to {dd:.2e}; factor coupling moves the solution by {coupling:.3f}")
    assert dd < 1e-6 and coupling > 1e-2

    # [4] cross-family agreement --------------------------------------------------
    Xm_s, f_s = solve_smooth(X0, Sigma, eta, gamma, T)
    Xm_scs, f_scs = solve_conic(X0, Sigma, eta, gamma, T, solver="SCS")
    print(f"[4] cross-family: conic vs L-BFGS ||dx||inf="
          f"{np.max(np.abs(Xm_full - Xm_s)):.2e} (dF={abs(f_full - f_s):.2e}); "
          f"Clarabel vs SCS ||dx||inf={np.max(np.abs(Xm_full - Xm_scs)):.2e}")
    assert np.max(np.abs(Xm_full - Xm_s)) < 1e-4
    assert np.max(np.abs(Xm_full - Xm_scs)) < 1e-5

    # [5] KKT stationarity ---------------------------------------------------------
    g = grad(Xm_full, X0, Sigma, eta, gamma)
    # |v|^{3/2} is C^1 (derivative 0 at v = 0), so the gradient formula is
    # valid everywhere -- no interiority requirement.
    res = float(np.max(np.abs(g)) / max(1.0, np.max(np.abs(2 * gamma * Xm_full @ Sigma))))
    print(f"[5] KKT: stationarity residual {res:.2e}")
    assert res < 1e-6

    # [6] frontier monotonicity ------------------------------------------------------
    risks, impacts = [], []
    for gm in (0.25, 0.5, 1.0, 2.0, 4.0):
        Xm, _ = solve_conic(X0, Sigma, eta, gm, T)
        v = trades(Xm, X0)
        risks.append(float(np.einsum("ti,ij,tj->", Xm, Sigma, Xm)))
        impacts.append(float((np.abs(v) ** 1.5 @ eta).sum()))
    print(f"[6] frontier (gamma 0.25->4): risk {risks[0]:.3f} -> {risks[-1]:.3f} "
          f"(monotone down), impact {impacts[0]:.3f} -> {impacts[-1]:.3f} "
          f"(monotone up)")
    assert all(a > b for a, b in zip(risks, risks[1:]))
    assert all(a < b for a, b in zip(impacts, impacts[1:]))

    print("\nAll checks passed.")


if __name__ == "__main__":
    main()
