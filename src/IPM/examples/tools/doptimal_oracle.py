# =============================================================================
# doptimal_oracle.py — self-contained validation for E10 (D-optimal design)
#
# Distributed D-optimal sensor design on a corridor of R regions
# (PSD + exp + NOC habitat: symmetric and nonsymmetric cones coupled
# within every region):
#
#     max_w  Σ_r log det M_r(w)
#     s.t.   M_r(w) = Σ_{i: footprint ∋ r} w_i a_i^{(r)} a_i^{(r)T}
#            Σ_{i ∈ region r} w_i = b_r,   w ≥ 0
#
# Each region estimates an m-dim local parameter; interior sensors inform
# one region, boundary sensors inform BOTH adjacent regions (the chain
# coupling), and every sensor's weight sits in its home region's budget.
# Classical continuous experiment design: Kiefer–Wolfowitz (1960),
# Fedorov (1972), Titterington (1976), Pukelsheim (1993), Boyd &
# Vandenberghe §7.5, Joshi & Boyd (2009).
#
# Conic form (what the Julia driver builds): per region, the log-det
# triangular trick — a 2m×2m PSD block [[M, Z],[Zᵀ, diag(Z)]] ⪰ 0 with Z
# lower triangular gives det M ≥ Π_j Z_jj, plus m exponential cones
# t_j ≤ log Z_jj; maximize Σ t. PSD and exp cones tied within one region.
#
# Checks:
#   [1] CLASSICAL SINGLE REGION: the multiplicative algorithm
#       w ← w ∘ leverage/m (Titterington 1976) ascends to the conic
#       optimum (one-sided); the Kiefer–Wolfowitz certificate holds at the
#       solution: max_i a_iᵀM⁻¹a_i = m, with equality on the support; the
#       support obeys Carathéodory: |supp| ≤ m(m+1)/2 (Fedorov).
#   [2] FORMULATION CERTIFICATE: the explicit triangular-PSD + exp-cone
#       program (exactly what e10.jl assembles) equals cvxpy's native
#       log_det formulation in value and weights.
#   [3] MULTI-REGION KW: generalized equivalence certificate — leverage
#       lev_i = Σ_{r ∋ i} a_i^{(r)T} M_r⁻¹ a_i^{(r)} is flat (= λ_r) on
#       each region's support and dominated off it; the exact identity
#       Σ_r b_r λ_r = Σ_r m_r holds; Clarabel vs SCS cross-check.
#   [4] DECOUPLING: deleting boundary sensors separates the problem into
#       R independent classical designs (multi-region == stacked
#       single-region exactly, λ_r = m/b_r); restoring them moves the
#       solution by a measurable amount.
#   [5] SCALING IDENTITY (exact): M(κw) = κ M(w), so
#       v(κ·b) = v(b) + (Σ_r m) log κ to machine precision.
#
# Usage: python3 doptimal_oracle.py [--quick]
# =============================================================================

import sys
import numpy as np
import cvxpy as cp

QUICK = "--quick" in sys.argv
SEED = 11
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
# Instance: R regions, m parameters each, nI interior sensors per region,
# nB boundary sensors between each adjacent pair (footprint in both).
# -----------------------------------------------------------------------------

def make_corridor(R=4, m=6, nI=40, nB=10, seed=SEED):
    rng = np.random.default_rng(seed)
    vec = lambda: rng.normal(size=m) * rng.uniform(0.5, 2.0)
    sensors = []          # (home region, {region: vector})
    for r in range(R):
        for _ in range(nI):
            sensors.append((r, {r: vec()}))
    for r in range(R - 1):
        for _ in range(nB):
            sensors.append((r, {r: vec(), r + 1: vec()}))
    b = np.ones(R)
    return sensors, b, R, m


def info_matrices(w, sensors, R, m):
    Ms = [np.zeros((m, m)) for _ in range(R)]
    for wi, (_, foot) in zip(w, sensors):
        for r, a in foot.items():
            Ms[r] += wi * np.outer(a, a)
    return Ms


def leverages(w, sensors, R, m):
    Ms = info_matrices(w, sensors, R, m)
    Minv = [np.linalg.inv(M) for M in Ms]
    return np.array([sum(a @ Minv[r] @ a for r, a in foot.items())
                     for _, foot in sensors])


def value_of(w, sensors, R, m):
    return sum(np.linalg.slogdet(M)[1] for M in info_matrices(w, sensors, R, m))


# -----------------------------------------------------------------------------
# Solvers
# -----------------------------------------------------------------------------

def solve_logdet(sensors, b, R, m):
    """Native cvxpy log_det formulation (reference)."""
    n = len(sensors)
    w = cp.Variable(n, nonneg=True)
    Ms = [0 for _ in range(R)]
    for i, (_, foot) in enumerate(sensors):
        for r, a in foot.items():
            Ms[r] = Ms[r] + w[i] * np.outer(a, a)
    cons = [cp.sum(w[[i for i, (h, _) in enumerate(sensors) if h == r]]) == b[r]
            for r in range(R)]
    prob = cp.Problem(cp.Maximize(cp.sum([cp.log_det(M) for M in Ms])), cons)
    conic_solve(prob)
    return w.value, prob.value


def solve_triangular(sensors, b, R, m):
    """Explicit triangular-PSD + exp-cone form — what e10.jl builds."""
    n = len(sensors)
    w = cp.Variable(n, nonneg=True)
    obj = 0
    cons = [cp.sum(w[[i for i, (h, _) in enumerate(sensors) if h == r]]) == b[r]
            for r in range(R)]
    for r in range(R):
        M = 0
        for i, (_, foot) in enumerate(sensors):
            if r in foot:
                a = foot[r]
                M = M + w[i] * np.outer(a, a)
        Z = cp.Variable((m, m))                    # lower triangular
        t = cp.Variable(m)
        cons += [Z[j, k] == 0 for j in range(m) for k in range(j + 1, m)]
        U = cp.bmat([[M, Z], [Z.T, cp.diag(cp.diag(Z))]])
        cons.append(U >> 0)
        cons += [cp.constraints.ExpCone(t[j], np.ones(1)[0] * 1.0, Z[j, j])
                 for j in range(m)]                # t_j <= log Z_jj
        obj = obj + cp.sum(t)
    prob = cp.Problem(cp.Maximize(obj), cons)
    conic_solve(prob)
    return w.value, prob.value


def solve_logdet_scs(sensors, b, R, m):
    n = len(sensors)
    w = cp.Variable(n, nonneg=True)
    Ms = [0 for _ in range(R)]
    for i, (_, foot) in enumerate(sensors):
        for r, a in foot.items():
            Ms[r] = Ms[r] + w[i] * np.outer(a, a)
    cons = [cp.sum(w[[i for i, (h, _) in enumerate(sensors) if h == r]]) == b[r]
            for r in range(R)]
    prob = cp.Problem(cp.Maximize(cp.sum([cp.log_det(M) for M in Ms])), cons)
    prob.solve(solver=cp.SCS, eps=1e-9, max_iters=200_000)
    assert prob.status in ("optimal", "optimal_inaccurate")
    return w.value, prob.value


def multiplicative(sensors, b, R, m, iters=5000):
    """Titterington's multiplicative algorithm, per-region normalized.

    w_i <- w_i * lev_i / lambda_r(w), lambda_r = (sum of on-region
    w*lev)/b_r, which preserves the budgets and is the classical
    w <- w * lev/m update when R = 1, b = 1.
    """
    n = len(sensors)
    home = np.array([h for h, _ in sensors])
    w = np.zeros(n)
    for r in range(R):
        idx = home == r
        w[idx] = b[r] / idx.sum()
    vals = []
    for _ in range(iters):
        lev = leverages(w, sensors, R, m)
        for r in range(R):
            idx = home == r
            lam = (w[idx] * lev[idx]).sum() / b[r]
            w[idx] *= lev[idx] / lam
        vals.append(value_of(w, sensors, R, m))
    return w, np.array(vals)


# -----------------------------------------------------------------------------
# Certificates
# -----------------------------------------------------------------------------

def kw_certificate(w, sensors, b, R, m, wtol=1e-5):
    """Generalized Kiefer-Wolfowitz residuals.

    Returns (flatness, dominance, euler):
      flatness  = max_r max_{i in supp_r} |lev_i - lambda_r| / lambda_r
      dominance = max_r max_{i off supp} (lev_i - lambda_r) / lambda_r
      euler     = |sum_r b_r lambda_r - R*m| / (R*m)
    with lambda_r = max on-support leverage in region r.
    """
    home = np.array([h for h, _ in sensors])
    lev = leverages(w, sensors, R, m)
    supp = w > wtol * w.max()
    flat, dom, blam = 0.0, 0.0, 0.0
    for r in range(R):
        on = supp & (home == r)
        off = (~supp) & (home == r)
        lam = lev[on].max()
        flat = max(flat, float(np.abs(lev[on] - lam).max() / lam))
        if off.any():
            dom = max(dom, float((lev[off] - lam).max() / lam))
        blam += b[r] * lam
    euler = abs(blam - R * m) / (R * m)
    return flat, dom, euler


# -----------------------------------------------------------------------------
# Checks
# -----------------------------------------------------------------------------

def main():
    R, m, nI, nB = (3, 5, 30, 8) if QUICK else (4, 6, 40, 10)
    sensors, b, R, m = make_corridor(R=R, m=m, nI=nI, nB=nB)
    n = len(sensors)
    print(f"R={R} regions, m={m}, {n} sensors "
          f"({nI}/region interior, {nB}/interface boundary)")

    # [1] classical single region ---------------------------------------------
    s1, b1, _, _ = make_corridor(R=1, m=m, nI=nI, nB=0, seed=SEED + 1)
    w_c, f_c = solve_logdet(s1, b1, 1, m)
    w_mult, vals = multiplicative(s1, b1, 1, m,
                                  iters=2000 if QUICK else 8000)
    mono = np.all(np.diff(vals) > -1e-11)
    gap = (f_c - vals[-1]) / abs(f_c)
    flat, dom, euler = kw_certificate(w_c, s1, b1, 1, m)
    supp = int((w_c > 1e-5 * w_c.max()).sum())
    cara = m * (m + 1) // 2
    print(f"[1] single region: multiplicative ascends monotonically ({mono}), "
          f"one-sided gap {gap:.1e};")
    print(f"    KW: flatness {flat:.1e}, dominance {dom:.1e} "
          f"(max leverage = m = {m}); support {supp} <= Caratheodory "
          f"{cara}")
    assert mono and -1e-9 < gap < 5e-4
    assert flat < 1e-3 and dom < 1e-6 and supp <= cara

    # [2] formulation certificate ------------------------------------------------
    w_ld, f_ld = solve_logdet(sensors, b, R, m)
    w_tr, f_tr = solve_triangular(sensors, b, R, m)
    dv = abs(f_ld - f_tr) / abs(f_ld)
    dw = float(np.max(np.abs(w_ld - w_tr)))
    print(f"[2] formulation: triangular-PSD+exp == log_det to rel {dv:.1e} "
          f"in value, {dw:.1e} in weights")
    assert dv < 1e-6 and dw < 1e-4

    # [3] multi-region KW + cross-solver ------------------------------------------
    flat, dom, euler = kw_certificate(w_ld, sensors, b, R, m)
    w_scs, f_scs = solve_logdet_scs(sensors, b, R, m)
    print(f"[3] multi-region KW: flatness {flat:.1e}, dominance {dom:.1e}, "
          f"Euler |Σ b_r λ_r - Rm|/Rm = {euler:.1e}; Clarabel vs SCS "
          f"dF {abs(f_ld - f_scs):.1e}, ||dw||inf "
          f"{np.max(np.abs(w_ld - w_scs)):.1e}")
    assert flat < 1e-3 and dom < 1e-6 and euler < 1e-4
    assert abs(f_ld - f_scs) < 1e-5 * abs(f_ld)

    # [4] decoupling ---------------------------------------------------------------
    s0 = [(h, f) for (h, f) in sensors if len(f) == 1]
    w_0, f_0 = solve_logdet(s0, b, R, m)
    f_stack, off = 0.0, 0
    w_stack = np.zeros(len(s0))
    for r in range(R):
        sr = [(0, {0: f[r]}) for (h, f) in s0 if h == r]
        wr, fr = solve_logdet(sr, b[r:r + 1], 1, m)
        w_stack[off:off + len(sr)] = wr
        f_stack += fr
        off += len(sr)
    dd = max(abs(f_0 - f_stack) / abs(f_0), float(np.max(np.abs(w_0 - w_stack))))
    # coupling moves the interior weights: compare on shared (interior) sensors
    int_idx = [i for i, (_, f) in enumerate(sensors) if len(f) == 1]
    coupling = float(np.max(np.abs(w_ld[int_idx] - w_0)))
    wB = w_ld[[i for i, (_, f) in enumerate(sensors) if len(f) == 2]]
    print(f"[4] decoupling: no-boundary multi-region == stacked to {dd:.1e}; "
          f"boundary sensors take weight {wB.sum():.3f} and move interior "
          f"weights by {coupling:.3f}")
    assert dd < 1e-5 and coupling > 1e-3 and wB.sum() > 0.05

    # [5] exact scaling identity ------------------------------------------------------
    kappa = 3.0
    _, f_k = solve_logdet(sensors, kappa * b, R, m)
    pred = f_ld + R * m * np.log(kappa)
    print(f"[5] scaling identity: v(κb) = {f_k:.8f} vs v(b) + Rm·logκ = "
          f"{pred:.8f} (diff {abs(f_k - pred):.1e})")
    assert abs(f_k - pred) < 1e-5 * abs(pred)

    print("\nAll checks passed.")


if __name__ == "__main__":
    main()
