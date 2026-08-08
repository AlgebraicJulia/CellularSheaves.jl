# =============================================================================
# poisson_tv_oracle.py — self-contained validation for E8 (Poisson-TV)
#
# Photon-limited deconvolution with nonlocal TV, exponential-cone habitat:
#
#     min_{u >= 0}  sum_g [ mu_g - y_g log mu_g ]  +  lam * sum_i ||z_i||_2
#     with          mu = tau * (G u) + b,   z_i = sqrt(w_i) (u_i - u_{nbrs(i)})
#
# y_g ~ Poisson(mu_g(u_true)); G a reflect-padded Gaussian PSF; b > 0 the
# dark/background rate (keeps mu interior even where u = 0). Scene: dark
# plateaus + a compact source + a smooth textured emission band, so zero
# counts occur and nonnegativity is live. Conic form: one exponential-cone
# triple per bin with y_g > 0 (t_g <= log mu_g; bins with y_g = 0 contribute
# only the linear term sum mu and need no cone), one SOC per node (TV),
# nonnegative orthant on u.
#
# Checks:
#   [1] MLEM (EM with additive background; Shepp-Vardi 1982) vs conic MLE at
#       lam = 0: EM's NLL converges to the conic optimum from above. (The
#       unregularized deconvolution MLE is ill-posed, so u itself is not
#       compared — only the optimal value, which is what EM certifies.)
#   [2] nodal Poisson-TV reference: cross-solver agreement (Clarabel vs SCS)
#       and a subgradient-selection stationarity certificate
#   [3] model value at low counts: exact Poisson likelihood clearly beats
#       both Gaussianized workflows -- the weighted-Gaussian approximation
#       (same TV, same blur) and no-deblur Anscombe-TV -- each at its own
#       best lambda
#   [4] Gaussian limit: Poisson and weighted-Gaussian solutions coalesce as
#       exposure tau grows (likelihoods agree to O(1/sqrt(counts)))
#   [5] positivity is load-bearing: sub-graph-resolution spike dipoles are
#       TV-invisible and numerically blur-null, so dropping u >= 0 does not
#       merely dip negative -- it destroys the reconstruction (the "nearly
#       black object" role of positivity in deconvolution)
#   [6] tile/coefficient (Chebyshev + C^k jets) formulation:
#       (a) completeness: full basis, no jets == nodal optimum
#       (b) compression: m = 32 of Tsz = 128 preserves reconstruction quality
#
# Usage: python3 poisson_tv_oracle.py [--quick]
# =============================================================================

import sys
import numpy as np
import scipy.sparse as sp
import cvxpy as cp

QUICK = "--quick" in sys.argv
SEED = 1
TIGHT = dict(solver=cp.CLARABEL, tol_gap_abs=1e-10, tol_gap_rel=1e-10,
             tol_feas=1e-10, max_iter=2000)


def conic_solve(prob):
    try:
        prob.solve(**TIGHT)
        assert prob.status in ("optimal", "optimal_inaccurate")
    except Exception:
        prob.solve(solver=cp.SCS, eps=1e-9, max_iters=200_000)
        assert prob.status in ("optimal", "optimal_inaccurate"), prob.status
    return prob.value


# -----------------------------------------------------------------------------
# Scene, PSF, counts
# -----------------------------------------------------------------------------

def make_scene(N):
    """Dark plateaus + compact source + smooth textured band (all C^inf)."""
    x = np.arange(N, dtype=float)
    ramp = lambda s: 0.5 - 0.5 * np.cos(np.pi * np.clip(s, 0.0, 1.0))
    src = 2.2 * np.exp(-0.5 * ((x - N / 4) / (N / 48)) ** 2)
    win = ramp((x - 0.52 * N) / (0.06 * N)) * ramp((0.92 * N - x) / (0.06 * N))
    band = win * (0.9 + 0.5 * np.sin(2 * np.pi * x / 16))
    return 0.05 + src + band


def psf_matrix(N, sigma=1.5, radius=4):
    taps = np.exp(-0.5 * (np.arange(-radius, radius + 1) / sigma) ** 2)
    taps /= taps.sum()
    rows, cols, vals = [], [], []
    for i in range(N):
        for k, t in zip(range(-radius, radius + 1), taps):
            j = i + k
            j = -j if j < 0 else (2 * (N - 1) - j if j >= N else j)  # reflect
            rows.append(i); cols.append(j); vals.append(t)
    return sp.csr_matrix((vals, (rows, cols)), shape=(N, N))


def make_counts(N, tau, b, seed=SEED):
    u_true = make_scene(N)
    G = psf_matrix(N)
    rng = np.random.default_rng(seed)
    return u_true, G, rng.poisson(tau * (G @ u_true) + b).astype(float)


# -----------------------------------------------------------------------------
# Nonlocal graph on Anscombe-stabilized counts (E6's construction)
# -----------------------------------------------------------------------------

def patch_matrix(v, q):
    N = len(v)
    vp = np.concatenate([v[q:0:-1], v, v[N - 2:N - q - 2:-1]])
    return np.stack([vp[i:i + 2 * q + 1] for i in range(N)])


def nonlocal_graph(y, K=12, q=4, R=64):
    ya = 2.0 * np.sqrt(y + 0.375)                     # Anscombe (1948)
    N = len(ya)
    P = patch_matrix(ya, q)
    h = np.median(np.abs(np.diff(ya)))
    nbrs = np.zeros((N, K), dtype=int)
    wts = np.zeros((N, K))
    for i in range(N):
        js = np.array([j for j in range(max(0, i - R), min(N, i + R + 1))
                       if abs(j - i) > 2])
        d2 = ((P[js] - P[i]) ** 2).sum(axis=1) / (2 * q + 1)
        order = np.argsort(d2)[:K]
        nbrs[i] = js[order]
        w = np.exp(-d2[order] / (2 * h ** 2))
        wts[i] = w * (K / w.sum())
    return nbrs, wts


def tv_ops(nbrs, wts):
    N, K = nbrs.shape
    rows, cols, vals = [], [], []
    for i in range(N):
        for k in range(K):
            sw = np.sqrt(wts[i, k])
            r = i * K + k
            rows += [r, r]; cols += [i, nbrs[i, k]]; vals += [sw, -sw]
    return sp.csr_matrix((vals, (rows, cols)), shape=(N * K, N))


def nll(y, mu):
    pos = y > 0
    return mu.sum() - y[pos] @ np.log(mu[pos])


# -----------------------------------------------------------------------------
# MLEM (EM with additive background) and conic solves
# -----------------------------------------------------------------------------

def mlem(y, G, tau, b, iters):
    N = G.shape[1]
    u = np.full(N, max(y.mean(), 1.0) / tau)
    sens = tau * np.asarray(G.sum(axis=0)).ravel()
    A_T = (tau * G).T.tocsr()
    for _ in range(iters):
        mu = tau * (G @ u) + b
        u = u * (A_T @ (y / mu)) / sens
    return u


def poisson_objective(u, y, G, tau, b, D, N, K, lam):
    z = (D @ u).reshape(N, K)
    return nll(y, tau * (G @ u) + b) + lam * np.linalg.norm(z, axis=1).sum()


def solve_poisson_tv(y, G, tau, b, D, N, K, lam, nonneg=True, solver=None):
    u = cp.Variable(N)
    mu = tau * (G @ u) + b
    pos = y > 0
    obj = cp.sum(mu) - y[pos] @ cp.log(mu[pos])
    if lam > 0:
        obj = obj + lam * cp.sum(
            cp.norm(cp.reshape(D @ u, (N, K), order="C"), 2, axis=1))
    prob = cp.Problem(cp.Minimize(obj), [u >= 0] if nonneg else [])
    if solver == "SCS":
        prob.solve(solver=cp.SCS, eps=1e-9, max_iters=200_000)
        assert prob.status in ("optimal", "optimal_inaccurate")
        f = prob.value
    else:
        f = conic_solve(prob)
    return u.value, f


def solve_gauss_tv(yhat, wobs, D, N, K, lam, nonneg=True):
    """Weighted-LS data fit + the same nonlocal TV (Gaussian surrogate)."""
    u = cp.Variable(N)
    obj = 0.5 * cp.sum(cp.multiply(wobs, cp.square(u - yhat))) + \
        lam * cp.sum(cp.norm(cp.reshape(D @ u, (N, K), order="C"), 2, axis=1))
    prob = cp.Problem(cp.Minimize(obj), [u >= 0] if nonneg else [])
    conic_solve(prob)
    return u.value


# -----------------------------------------------------------------------------
# Chebyshev tile machinery (numpy port of E6's)
# -----------------------------------------------------------------------------

def chebvander(us, m):
    V = np.zeros((len(us), m))
    V[:, 0] = 1.0
    if m > 1:
        V[:, 1] = us
    for j in range(2, m):
        V[:, j] = 2 * us * V[:, j - 1] - V[:, j - 2]
    return V


def chebder_mat(m):
    Dm = np.zeros((m, m))
    for j in range(1, m):
        c = np.zeros(m); c[j] = 1.0
        d = np.zeros(m + 2)
        for k in range(j, 0, -1):
            d[k - 1] = d[k + 1] + 2 * k * c[k]
        d[0] /= 2
        Dm[:, j] = d[:m]
    return Dm


def tile_eval(N, Tsz, m):
    P = N // Tsz
    us = 2.0 * np.arange(Tsz) / (Tsz - 1) - 1.0
    return sp.block_diag([chebvander(us, m)] * P).tocsr(), P


def jet_constraints(P, Tsz, m, k):
    Dm = chebder_mat(m)
    pairs = []
    for r in range(k + 1):
        Dr = np.linalg.matrix_power(Dm, r) * (2.0 / Tsz) ** r
        ep = (chebvander(np.array([1.0]), m) @ Dr).ravel()
        em = (chebvander(np.array([-1.0]), m) @ Dr).ravel()
        s = np.linalg.norm(np.concatenate([ep, em]))
        pairs.append((ep / s, em / s))
    A = []
    for t in range(P - 1):
        for ep, em in pairs:
            row = np.zeros(P * m)
            row[t * m:(t + 1) * m] = ep
            row[(t + 1) * m:(t + 2) * m] = -em
            A.append(row)
    return np.array(A) if A else np.zeros((0, P * m))


def solve_coeff_poisson_tv(y, G, tau, b, D, N, K, lam, Tsz, m, k):
    E, P = tile_eval(N, Tsz, m)
    J = jet_constraints(P, Tsz, m, k)
    th = cp.Variable(P * m)
    u = E @ th
    mu = tau * (G @ u) + b
    pos = y > 0
    obj = cp.sum(mu) - y[pos] @ cp.log(mu[pos]) + \
        lam * cp.sum(cp.norm(cp.reshape(D @ u, (N, K), order="C"), 2, axis=1))
    cons = [u >= 0]
    if J.shape[0] > 0:
        cons.append(J @ th == 0)
    prob = cp.Problem(cp.Minimize(obj), cons)
    f = conic_solve(prob)
    return np.asarray(E @ th.value).ravel(), f


# -----------------------------------------------------------------------------
# Stationarity certificate (subgradient selection, following the E6 oracle)
# -----------------------------------------------------------------------------

def kkt_residual(u, y, G, tau, b, D, N, K, lam):
    mu = tau * (G @ u) + b
    g = tau * np.asarray(G.T @ (1.0 - y / mu)).ravel()
    z = (D @ u).reshape(N, K)
    nz = np.linalg.norm(z, axis=1)
    act = nz <= 1e-3 * max(nz.max(), 1e-12)           # groups at the kink
    shat = np.zeros_like(z)
    shat[~act] = z[~act] / nz[~act, None]
    r0 = g + lam * (D.T @ shat.ravel())
    if act.any():
        idx = np.where(act)[0]
        Dact = sp.vstack([D[i * K:(i + 1) * K] for i in idx]).T.tocsc()
        s = sp.linalg.lsqr(Dact * lam, -r0, atol=1e-12, btol=1e-12)[0]
        sn = np.linalg.norm(s.reshape(-1, K), axis=1)
        s = (s.reshape(-1, K) / np.maximum(sn, 1.0)[:, None]).ravel()
        r = r0 + lam * (Dact @ s)
    else:
        r = r0
    scale = max(1.0, np.abs(g).max())
    free = u > 1e-6
    viol_free = np.abs(r[free]).max() / scale
    viol_bnd = max(0.0, -(r[~free].min() if (~free).any() else 0.0)) / scale
    return max(viol_free, viol_bnd)


# -----------------------------------------------------------------------------
# Checks
# -----------------------------------------------------------------------------

def main():
    N = 256 if QUICK else 512
    tau, b, lam = 15.0, 0.2, 0.8
    u_true, G, y = make_counts(N, tau, b)
    nbrs, wts = nonlocal_graph(y)
    K = nbrs.shape[1]
    D = tv_ops(nbrs, wts)
    mse = lambda u: float(((u - u_true) ** 2).mean())

    print(f"N={N} tau={tau} b={b}: counts min={y.min():.0f} "
          f"med={np.median(y):.0f} max={y.max():.0f}, zeros={int((y == 0).sum())}")

    # [1] MLEM certifies the conic MLE optimum --------------------------------
    u_ml, f_ml = solve_poisson_tv(y, G, tau, b, D, N, K, lam=0.0)
    f_em = nll(y, tau * (G @ mlem(y, G, tau, b, 1000 if QUICK else 4000)) + b)
    gap = (f_em - f_ml) / abs(f_ml)
    print(f"[1] MLEM vs conic MLE: NLL {f_em:.4f} -> {f_ml:.4f} "
          f"(rel gap {gap:.1e}, one-sided)")
    assert -1e-8 < gap < 5e-4, gap

    # [2] Poisson-TV reference: cross-solver + stationarity -------------------
    u_tv, f_tv = solve_poisson_tv(y, G, tau, b, D, N, K, lam=lam)
    u_scs, f_scs = solve_poisson_tv(y, G, tau, b, D, N, K, lam=lam, solver="SCS")
    res = kkt_residual(u_tv, y, G, tau, b, D, N, K, lam)
    print(f"[2] Poisson-TV: F={f_tv:.4f}; Clarabel vs SCS "
          f"||du||inf={np.abs(u_tv - u_scs).max():.1e}, "
          f"dF={abs(f_tv - f_scs):.1e}; KKT residual={res:.1e}")
    assert np.abs(u_tv - u_scs).max() < 1e-4 and res < 2e-3

    # [3] model value at low counts --------------------------------------------
    grid = (0.2, 0.4, 0.8, 1.6, 3.2)
    m_pois = min(mse(solve_poisson_tv(y, G, tau, b, D, N, K, lg)[0])
                 for lg in grid)
    wobs = tau ** 2 / np.maximum(y, 1.0)
    m_gauss = min(mse(solve_gauss_tv((y - b) / tau, wobs, D, N, K, lg))
                  for lg in grid)
    ya = 2.0 * np.sqrt(y + 0.375)
    m_ans = np.inf
    for lg in grid:
        ua = solve_gauss_tv(ya, np.ones(N), D, N, K, 2 * lg, nonneg=False)
        m_ans = min(m_ans, mse(np.maximum((ua / 2) ** 2 - 0.375 - b, 0) / tau))
    print(f"[3] model value (MSE, best lam each): Poisson-TV {m_pois:.4f} vs "
          f"weighted-Gauss-TV {m_gauss:.4f}, Anscombe-TV(no deblur) {m_ans:.4f}, "
          f"MLE {mse(np.maximum(u_ml, 0)):.4f}")
    assert m_pois < 0.75 * min(m_gauss, m_ans)   # exact likelihood wins clearly

    # [4] Gaussian limit: median |u_pois - u_wgauss| ~ 1/sqrt(tau) ------------
    taus, gaps = (15.0, 60.0, 240.0), []
    print("[4] Gaussian limit, median |u_pois - u_wgauss|:", end=" ")
    for tt in taus:
        ut2, G2, y2 = make_counts(N, tt, b, seed=2)
        nb2, wt2 = nonlocal_graph(y2)
        D2 = tv_ops(nb2, wt2)
        lt = lam * tt / tau
        up, _ = solve_poisson_tv(y2, G2, tt, b, D2, N, K, lam=lt)
        uw = solve_gauss_tv((y2 - b) / tt, tt ** 2 / np.maximum(y2, 1.0),
                            D2, N, K, lt)
        gaps.append(float(np.median(np.abs(up - uw))))
        print(f"tau={tt:.0f}: {gaps[-1]:.4f}", end="  ")
    slope = np.polyfit(np.log(taus), np.log(gaps), 1)[0]
    print(f"(slope tau^{slope:.2f})")
    assert gaps[2] < gaps[1] < gaps[0] and slope < -0.25

    # [5] positivity is load-bearing (pinned to the N = 512 reference
    # instance; the blur-null dipole blow-up is a property of that geometry).
    # The nonlocal graph excludes offsets <= 2, so spike dipoles at
    # sub-graph resolution are TV-invisible AND numerically blur-null
    # (Gaussian transfer ~1.5e-5 at Nyquist): without the cone they blow up.
    # A tiny ridge keeps the diagnostic solve well-posed. Cf. Donoho,
    # Johnstone, Hoch & Stern (1992), "nearly black" objects.
    N5 = 512
    if N == N5:
        ut5, G5, y5, D5, utv5 = u_true, G, y, D, u_tv
    else:
        ut5, G5, y5 = make_counts(N5, tau, b)
        nb5, wt5 = nonlocal_graph(y5)
        D5 = tv_ops(nb5, wt5)
        utv5, _ = solve_poisson_tv(y5, G5, tau, b, D5, N5, K, lam=lam)
    u2 = cp.Variable(N5)
    mu2 = tau * (G5 @ u2) + b
    pos5 = y5 > 0
    obj2 = cp.sum(mu2) - y5[pos5] @ cp.log(mu2[pos5]) + \
        lam * cp.sum(cp.norm(cp.reshape(D5 @ u2, (N5, K), order="C"), 2, axis=1)) + \
        0.5e-3 * tau * cp.sum_squares(u2)
    conic_solve(cp.Problem(cp.Minimize(obj2)))
    u_free = u2.value
    mse5 = lambda u: float(((u - ut5) ** 2).mean())
    print(f"[5] positivity load-bearing (N=512): min u unconstrained "
          f"{u_free.min():+.1f} (MSE {mse5(u_free):.1f}) vs constrained "
          f"{utv5.min():.1e} (MSE {mse5(utv5):.4f})")
    assert u_free.min() < -1.0 and mse5(u_free) > 100 * mse5(utv5)

    # [6] tile/coefficient formulation -------------------------------------------
    Nc, Tc = 128, 16
    utc, Gc, yc = make_counts(Nc, tau, b, seed=3)
    nbc, wtc = nonlocal_graph(yc, K=6, q=3, R=12)
    Dc = tv_ops(nbc, wtc)
    _, f_nod = solve_poisson_tv(yc, Gc, tau, b, Dc, Nc, 6, lam=lam)
    _, f_cof = solve_coeff_poisson_tv(yc, Gc, tau, b, Dc, Nc, 6, lam,
                                      Tsz=Tc, m=Tc, k=-1)
    print(f"[6a] completeness (Tsz=m={Tc}, no jets): F {f_cof:.6f} vs nodal "
          f"{f_nod:.6f} (excess {f_cof - f_nod:+.1e})")
    assert abs(f_cof - f_nod) < 1e-4 * abs(f_nod)

    u32, _ = solve_coeff_poisson_tv(y, G, tau, b, D, N, K, lam,
                                    Tsz=128, m=32, k=2)
    print(f"[6b] compression (Tsz=128, m=32, C^2, u-DOF/4): "
          f"MSE {mse(np.maximum(u32, 0)):.4f} vs nodal {mse(u_tv):.4f}")

    print("\nAll checks passed.")


if __name__ == "__main__":
    main()
