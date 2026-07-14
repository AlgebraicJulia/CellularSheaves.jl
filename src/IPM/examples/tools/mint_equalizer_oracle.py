# =============================================================================
# mint_equalizer_oracle.py — self-contained validation for E14 (short
# multichannel FIR equalizer-bank design under colored noise: regularized
# MINT, i.e. the D = 1 synthesis filter bank)
#
# M channels with impulse responses h_m (length Lh, exponentially decaying
# reverberant FIR), one equalizer f_m (length Lf) per channel, composite
# response c = sum_m h_m * f_m in R^P, P = Lh + Lf - 1, target g = e_d
# (pure delay d). Stacked convolution map T = [Toep(h_1) ... Toep(h_M)],
# c = T f. Design problem:
#
#     min_f  (1/2) sum_m f_m' (Theta + lam I) f_m    (amplified noise power)
#     s.t.   || T f - g || <= eps                    (near-perfect equalization)
#
# Theta is the autocorrelation Gram of the AR(2) ambient noise field (dense
# BY NATURE — the noise is colored), lam = sensor self-noise power (white
# floor; makes the objective strictly PD, E13's uniqueness discipline —
# physics, not decoration). Lineage: MINT, Miyoshi & Kaneda (IEEE TASSP
# 1988); regularized/short multichannel equalization (Kodrasi & Doclo lineage);
# channel shortening (Melsa-Younce-Rohrs 1996); Bezout/coprimality (Kailath).
#
# SHORT-EQUALIZER (TALL) REGIME. The suite formulation lives BELOW the MINT
# length threshold Lf* = ceil((Lh-1)/(M-1)): T is tall with full column rank,
# exact inversion is impossible (eps_LS > 0), and the design trades residual
# ISI against noise amplification — the practically relevant regime (exact
# MINT is notoriously noise-fragile). eps_LS = ||(I - T G^{-1} T') g|| is the
# short-equalizer floor, G = T'T.
#
# WHITENING (E7's trick, relocated from objective to CONSTRAINT). With
# G = L L' (Cholesky) and b = L^{-1} T' g:
#
#     || T f - g ||^2 = || L' f - b ||^2 + eps_LS^2          (joint identity)
#
# so eps = hypot(eps_LS, eps_x) turns the ball into || L' f - b || <= eps_x:
# the excess budget eps_x IS the whitened radius. Two closed-form endpoints:
# eps_x -> 0 collapses the feasible set to the pinv bank f_hat = G^{-1} T' g;
# eps >= ||g|| = 1 admits f = 0.
#
# SPLIT (E13's Agler move, relocated from the constraint slack to the
# reconstruction budget). Center the target split at the pinv bank,
# g_m = T_m f_hat_m + r_hat / M with r_hat = g - T f_hat (orthogonal to
# range(T), hence to every range(T_m)); then EXACTLY
#
#     || T_m f_m - g_m ||^2 = || L_m'(f_m - f_hat_m) ||^2 + (eps_LS/M)^2
#
# with G_m = T_m' T_m = L_m L_m' (per-channel autocorrelation Gram, dense by
# nature for reverberant channels with Lf <= Lh). One SOC stalk
# (s_m, w_m, sigma_m) of dim Lf + 2 per channel: dense ties
# w_m = L_m' f_m - b_m (b_m = L_m' f_hat_m, nonzero rhs — E7's pattern), a
# one-row pin sigma_m = eps_LS / M, and ONE budget row sum_m s_m = eps. By
# the triangle inequality sum_m ||T_m f_m - g_m|| >= ||T f - g||, feasibility
# of the split program CERTIFIES ||T f - g|| <= eps — a compositional
# certificate whose conservatism is measured, not hidden (checks [5]), and
# whose eps_x -> 0 endpoint is EXACTLY the joint one (both pin f = f_hat).
#
# EXACT SLICE SPLIT (zero conservatism, O(M^2) coupling mass). The whitened
# residual partitions EXACTLY by Cholesky row blocks,
#
#     || L' f - b ||^2 = sum_i || (L' f - b)_i ||^2,
#
# so M slice stalks (s_i, w_i) of dim Lf + 1 with dense triangular ties
# w_i = sum_{j >= i} (L')_{ij} f_j - b_i, plus one (M+1)-dim hub SOC
# enforcing ||s|| <= eps_x, reproduce the joint ball with NO relaxation —
# the gluing maps ARE the Cholesky blocks of the joint convolution Gram.
# The price is coupling mass ~ M^2 Lf^2 / 2 (slice i touches filters
# i..M) vs the budget split's M Lf^2: a pre-registered tradeoff pair —
# conservatism vs coupling mass — that the benchmark measures on both
# sides. The joint ball as a single SOC stalk of dim M*Lf + 2 is the
# monolithic comparator (cf. E13's monolithic LMI).
#
# WHY D = 1 (the polyphase exclusion). For a decimated bank (decimation
# D > 1) the per-filter constraint Gram is EXACTLY polyphase-block-diagonal:
# (T_m' T_m)[j, j'] = 0 whenever j != j' (mod D) — structural zeros no
# whitening can remove, i.e. structural 1/D block density: sheaf shape
# without sheaf density, in B rather than Q. The decimated bank is therefore
# out of habitat BY STRUCTURE (a boundary-suite control, pre-registered:
# bloat ~ D), and E14 is the D = 1 member of the filter-bank family, where
# every Gram is honestly dense. Check [1] measures both facts.
#
# Checks:
#   [1] DENSE BY NATURE + EXCLUSION: densities of Theta, G_m, L_m', joint G
#       at 1e-8 of the max entry; the D = 4 decimated control has off-class
#       Gram entries EXACTLY zero.
#   [2] MINT / BEZOUT (classical certificate): eps_LS(Lf) drops to machine
#       zero exactly at Lf* = ceil((Lh-1)/(M-1)) (Miyoshi-Kaneda length
#       condition); channels sharing a common zero keep eps_LS bounded away
#       from zero ABOVE the threshold (coprimality is necessary — Bezout).
#   [3] FORMULATION CERTIFICATE: whitening identities (joint and split) at
#       random points to machine precision; the stalk assembly (Q, B, g and
#       cone layout exactly as e14.jl builds) equals the native cvxpy
#       formulation.
#   [4] TRS (classical certificate): the joint problem is a trust-region
#       subproblem; the secular-equation solve (Gay 1981; More-Sorensen 1983)
#       agrees with the conic solver, and the More-Sorensen optimality
#       conditions hold at the solution.
#   [5] EXACTNESS + CONSERVATISM: v_exact == v_joint (zero conservatism);
#       the budget split satisfies its certificate ||T f - g|| <= eps but
#       uses only ~1/3 of the whitened ball, with v_budget/v_joint -> 1 as
#       eps_x -> 0 — the measured price of O(M) coupling mass.
#   [6] ALLOCATION LOAD-BEARING (budget form): the optimized budget beats
#       the uniform allocation s_m = eps/M; allocation spread reported.
#   [7] DEFORMATION: ||f(eps_x) - f_hat||_inf -> 0 at measured rate ~eps_x^1.
#   [8] NOISE-SHAPING CONTROL (diagonal-Q ablation, exact form): Theta ->
#       white moves the solution by >10% and pays a large excess — the
#       D6/smooth-bary control, pre-registering the e14.jl --white row.
#   [9] EQUALIZATION (exact form): simulated source through channels + AR
#       noise; the designed bank cuts output noise power ~2x vs the pinv
#       bank at a certified ISI budget; measured ratio matches the
#       quadratic-form prediction.
#  [10] CROSS-SOLVER: Clarabel vs SCS on the exact program.
#
# FLAGSHIP = EXACT SLICE SPLIT. Measured at the default operating point:
# the exact form uses 100% of the whitened ball and cuts amplified noise
# power to ~0.45x the pinv bank; the budget split's triangle-inequality
# conservatism throttles ball use to ~0.33 and the noise gain to ~0.93x.
# e14.jl therefore builds the exact form by default, the budget form under
# --budget, and the oracle's TRS solve is the joint comparator for both.
#
# Figures quoted in e14.jl are from this file executed at M = 8, Lh = 1024,
# Lf = 64, trev = 256, AR(2) r = 0.97 phi = 0.55 pi, lam = 1e-3,
# eps_x = 0.05, d = (Lh + Lf) // 2.
# =============================================================================

import numpy as np
from scipy.signal import lfilter

np.random.seed(0)
RNG = np.random.default_rng(7)

CLARABEL = dict(solver="CLARABEL")
SCS = dict(solver="SCS", eps=1e-9, max_iters=200_000)


# -----------------------------------------------------------------------------
# Instance
# -----------------------------------------------------------------------------

def make_channels(M=8, Lh=1024, trev=256, seed=7):
    rng = np.random.default_rng(seed)
    env = np.exp(-np.arange(Lh) / trev)
    chans = []
    for _ in range(M):
        h = rng.standard_normal(Lh) * env
        chans.append(h / np.linalg.norm(h))
    return chans


def ar2_gram(Lf, r=0.97, phi=0.55 * np.pi, nimp=8192):
    """Autocorrelation Toeplitz of an AR(2) noise field, normalized diag = 1."""
    a = np.array([1.0, -2 * r * np.cos(phi), r * r])
    psi = lfilter([1.0], a, np.eye(1, nimp).ravel())
    rho = np.array([psi[:nimp - k] @ psi[k:] for k in range(Lf)])
    rho = rho / rho[0]
    idx = np.abs(np.subtract.outer(np.arange(Lf), np.arange(Lf)))
    return rho[idx], a


def toep(h, Lf):
    Lh = len(h)
    T = np.zeros((Lh + Lf - 1, Lf))
    for j in range(Lf):
        T[j:j + Lh, j] = h
    return T


class Instance:
    def __init__(self, M=8, Lh=1024, Lf=64, trev=256, r=0.97,
                 phi=0.55 * np.pi,
                 lam=1e-3, epsx=0.05, seed=7):
        self.M, self.Lh, self.Lf, self.lam, self.epsx = M, Lh, Lf, lam, epsx
        self.d = (Lh + Lf) // 2
        self.chans = make_channels(M, Lh, trev, seed)
        self.Theta, self.ar_a = ar2_gram(Lf, r, phi)
        # heterogeneous sensor-noise scales (mics of different quality /
        # distance) — this is what makes the budget allocation load-bearing
        self.sigmas = np.geomspace(0.5, 2.0, M)
        self.Qms = [sg ** 2 * self.Theta + lam * np.eye(Lf)
                    for sg in self.sigmas]
        self.Ts = [toep(h, Lf) for h in self.chans]
        self.T = np.hstack(self.Ts)
        self.P = Lh + Lf - 1
        self.g = np.zeros(self.P)
        self.g[self.d] = 1.0
        # joint whitening
        self.G = self.T.T @ self.T
        self.L = np.linalg.cholesky(self.G)              # G = L L'
        self.fhat = np.linalg.solve(self.G, self.T.T @ self.g)
        self.rhat = self.g - self.T @ self.fhat
        self.eLS = np.linalg.norm(self.rhat)
        self.b = np.linalg.solve(self.L, self.T.T @ self.g)   # = L' fhat
        self.eps = np.hypot(self.eLS, epsx)
        # per-channel whitening
        self.Gms = [Tm.T @ Tm for Tm in self.Ts]
        self.Lms = [np.linalg.cholesky(Gm) for Gm in self.Gms]
        self.fhm = self.fhat.reshape(M, Lf)
        self.bms = [self.Lms[m].T @ self.fhm[m] for m in range(M)]

    def Qm(self, m, white=False):
        """White control replaces Theta by an equal-power diagonal, KEEPING
        the per-channel scales — it isolates coloredness, not heterogeneity."""
        if white:
            return (self.sigmas[m] ** 2 * np.trace(self.Theta) / self.Lf
                    + self.lam) * np.eye(self.Lf)
        return self.Qms[m]

    def Qbar(self, white=False):
        from scipy.linalg import block_diag
        return block_diag(*[self.Qm(m, white) for m in range(self.M)])


# -----------------------------------------------------------------------------
# Solvers
# -----------------------------------------------------------------------------

def solve_split(inst, epsx=None, white=False, uniform=False, solver_kw=None):
    """Native cvxpy form of the split program e14.jl assembles."""
    import cvxpy as cp
    M, Lf = inst.M, inst.Lf
    epsx = inst.epsx if epsx is None else epsx
    eps = np.hypot(inst.eLS, epsx)
    F = cp.Variable((M, Lf))
    s = cp.Variable(M)
    cons = []
    for m in range(M):
        resid = cp.hstack([inst.Lms[m].T @ F[m] - inst.bms[m],
                           np.array([inst.eLS / M])])
        cons.append(cp.SOC(s[m], resid))
    if uniform:
        cons.append(s == eps / M)
    else:
        cons.append(cp.sum(s) == eps)
    obj = 0.5 * cp.sum([cp.quad_form(F[m], cp.psd_wrap(inst.Qm(m, white)))
                        for m in range(M)])
    pr = cp.Problem(cp.Minimize(obj), cons)
    pr.solve(**(solver_kw or CLARABEL))
    assert pr.status in ("optimal", "optimal_inaccurate"), pr.status
    return F.value.ravel(), pr.value


def solve_joint(inst, epsx=None, solver_kw=None):
    import cvxpy as cp
    epsx = inst.epsx if epsx is None else epsx
    f = cp.Variable(inst.M * inst.Lf)
    resid = cp.hstack([inst.L.T @ f - inst.b, np.array([inst.eLS])])
    cons = [cp.SOC(np.hypot(inst.eLS, epsx), resid)]
    Qb = inst.Qbar()
    pr = cp.Problem(cp.Minimize(0.5 * cp.quad_form(f, cp.psd_wrap(Qb))), cons)
    pr.solve(**(solver_kw or CLARABEL))
    assert pr.status in ("optimal", "optimal_inaccurate"), pr.status
    return f.value, pr.value


def solve_joint_trs(inst, epsx=None):
    """Secular-equation TRS solve of the joint problem (Gay 1981; MS 1983).

    In whitened coordinates y = L' f - b the joint program is
        min (1/2)(y + b)' H (y + b),  ||y|| <= eps_x,   H = L^{-1} Qbar L^{-T}.
    H is PD, the unconstrained minimizer is y = -b (i.e. f = 0), and for
    eps_x < ||b|| the constraint is active: y(nu) = -(H + nu I)^{-1} H b with
    the multiplier nu > 0 solving ||y(nu)|| = eps_x.
    """
    epsx = inst.epsx if epsx is None else epsx
    n = inst.M * inst.Lf
    Qb = inst.Qbar()
    Linv = np.linalg.inv(inst.L)
    H = Linv @ Qb @ Linv.T
    H = 0.5 * (H + H.T)
    lamH, U = np.linalg.eigh(H)
    beta = U.T @ (H @ inst.b)
    assert np.linalg.norm(inst.b) > epsx, "constraint inactive: f = 0 optimal"

    def ynorm(nu):
        return np.sqrt(np.sum((beta / (lamH + nu)) ** 2))

    lo, hi = 0.0, 1.0
    while ynorm(hi) > epsx:
        hi *= 4.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if ynorm(mid) > epsx:
            lo = mid
        else:
            hi = mid
    nu = 0.5 * (lo + hi)
    y = U @ (-beta / (lamH + nu))
    f = np.linalg.solve(inst.L.T, y + inst.b)
    val = 0.5 * f @ (Qb @ f)
    ms = dict(nu=nu,
              mineig=lamH.min() + nu,                      # H + nu I > 0
              resid=abs(np.linalg.norm(y) - epsx),         # boundary
              stat=np.linalg.norm(H @ (y + inst.b) + nu * y))
    return f, val, ms


def solve_exact(inst, epsx=None, solver_kw=None):
    """Native cvxpy form of the exact slice split (zero conservatism)."""
    import cvxpy as cp
    M, Lf = inst.M, inst.Lf
    epsx = inst.epsx if epsx is None else epsx
    f = cp.Variable(M * Lf)
    s = cp.Variable(M)
    cons = []
    for i in range(M):
        sl = slice(i * Lf, (i + 1) * Lf)
        cons.append(cp.SOC(s[i], inst.L.T[sl, :] @ f - inst.b[sl]))
    cons.append(cp.SOC(epsx, s))
    Qb = inst.Qbar()
    pr = cp.Problem(cp.Minimize(0.5 * cp.quad_form(f, cp.psd_wrap(Qb))), cons)
    pr.solve(**(solver_kw or CLARABEL))
    assert pr.status in ("optimal", "optimal_inaccurate"), pr.status
    return f.value, pr.value


def solve_exact_white(inst, epsx=None, solver_kw=None):
    import cvxpy as cp
    M, Lf = inst.M, inst.Lf
    epsx = inst.epsx if epsx is None else epsx
    f = cp.Variable(M * Lf)
    s = cp.Variable(M)
    cons = []
    for i in range(M):
        sl = slice(i * Lf, (i + 1) * Lf)
        cons.append(cp.SOC(s[i], inst.L.T[sl, :] @ f - inst.b[sl]))
    cons.append(cp.SOC(epsx, s))
    Qw = inst.Qbar(white=True)
    pr = cp.Problem(cp.Minimize(0.5 * cp.quad_form(f, cp.psd_wrap(Qw))), cons)
    pr.solve(**(solver_kw or CLARABEL))
    assert pr.status in ("optimal", "optimal_inaccurate"), pr.status
    return f.value, pr.value


def assemble_exact(inst, epsx=None):
    """The exact-slice (Q, B, g, cones) layout e14.jl builds by default. Stalk
    order: filter stalks 1..M (dim Lf, free), slice SOC stalks (s_i, w_i)
    (dim Lf + 1), hub SOC stalk (h0, shat) (dim M + 1). Edges: slice ties
    (Lf rows, w_i - sum_{j>=i} (L')_{ij} f_j = -b_i), s-links
    (shat_i - s_i = 0), hub pin (h0 = eps_x)."""
    M, Lf = inst.M, inst.Lf
    epsx = inst.epsx if epsx is None else epsx
    dsoc = Lf + 1
    n0 = M * Lf + M * dsoc + (M + 1)
    fcol = lambda m: slice(m * Lf, (m + 1) * Lf)
    scol = lambda i: slice(M * Lf + i * dsoc, M * Lf + (i + 1) * dsoc)
    hcol = slice(M * Lf + M * dsoc, n0)
    n1 = M * Lf + M + 1
    B = np.zeros((n1, n0))
    gv = np.zeros(n1)
    row = 0
    for i in range(M):                                    # slice ties
        B[row:row + Lf, scol(i)][:, 1:] = np.eye(Lf)
        for j in range(i, M):
            B[row:row + Lf, fcol(j)] = -inst.L.T[i * Lf:(i + 1) * Lf,
                                                 j * Lf:(j + 1) * Lf]
        gv[row:row + Lf] = -inst.b[i * Lf:(i + 1) * Lf]
        row += Lf
    for i in range(M):                                    # s-links
        B[row, hcol.start + 1 + i] = 1.0
        B[row, scol(i).start] = -1.0
        row += 1
    B[row, hcol.start] = 1.0                              # hub pin
    gv[row] = epsx
    Q = np.zeros((n0, n0))
    for m in range(M):
        Q[fcol(m), fcol(m)] = inst.Qm(m)
    return Q, B, gv, n0, [scol(i) for i in range(M)] + [hcol]


def solve_assembly_exact(inst, epsx=None, solver_kw=None):
    import cvxpy as cp
    Q, B, gv, n0, socs = assemble_exact(inst, epsx)
    x = cp.Variable(n0)
    cons = [B @ x == gv]
    for sl in socs:
        cons.append(cp.SOC(x[sl.start], x[sl.start + 1:sl.stop]))
    pr = cp.Problem(cp.Minimize(0.5 * cp.quad_form(x, cp.psd_wrap(Q))), cons)
    pr.solve(**(solver_kw or CLARABEL))
    assert pr.status in ("optimal", "optimal_inaccurate"), pr.status
    f = np.concatenate([x.value[m * inst.Lf:(m + 1) * inst.Lf]
                        for m in range(inst.M)])
    return f, pr.value


def assemble_split(inst, epsx=None, white=False):
    """The budget-split (Q, B, g, cones) layout of e14.jl --budget. Stalk order:
    filter stalks 1..M (dim Lf, free) then SOC stalks (s_m, w_m, sigma_m)
    (dim Lf + 2). Edges: per-channel ties (Lf rows, w_m - L_m' f_m = -b_m),
    per-channel pins (sigma_m = eps_LS/M), one budget row (sum s_m = eps)."""
    M, Lf = inst.M, inst.Lf
    epsx = inst.epsx if epsx is None else epsx
    eps = np.hypot(inst.eLS, epsx)
    dsoc = Lf + 2
    n0 = M * Lf + M * dsoc
    fcol = lambda m: slice(m * Lf, (m + 1) * Lf)
    scol = lambda m: slice(M * Lf + m * dsoc, M * Lf + (m + 1) * dsoc)
    n1 = M * (Lf + 1) + 1
    B = np.zeros((n1, n0))
    gv = np.zeros(n1)
    row = 0
    for m in range(M):                                    # ties
        B[row:row + Lf, scol(m)][:, 1:Lf + 1] = np.eye(Lf)
        B[row:row + Lf, fcol(m)] = -inst.Lms[m].T
        gv[row:row + Lf] = -inst.bms[m]
        row += Lf
    for m in range(M):                                    # pins
        B[row, scol(m).start + Lf + 1] = 1.0
        gv[row] = inst.eLS / M
        row += 1
    for m in range(M):                                    # budget
        B[row, scol(m).start] = 1.0
    gv[row] = eps
    Q = np.zeros((n0, n0))
    for m in range(M):
        Q[fcol(m), fcol(m)] = inst.Qm(m, white)
    return Q, B, gv, n0, [scol(m) for m in range(M)]


def solve_assembly(inst, epsx=None, solver_kw=None):
    import cvxpy as cp
    Q, B, gv, n0, socs = assemble_split(inst, epsx)
    x = cp.Variable(n0)
    cons = [B @ x == gv]
    for sl in socs:
        cons.append(cp.SOC(x[sl.start], x[sl.start + 1:sl.stop]))
    pr = cp.Problem(cp.Minimize(0.5 * cp.quad_form(x, cp.psd_wrap(Q))), cons)
    pr.solve(**(solver_kw or CLARABEL))
    assert pr.status in ("optimal", "optimal_inaccurate"), pr.status
    f = np.concatenate([x.value[m * inst.Lf:(m + 1) * inst.Lf]
                        for m in range(inst.M)])
    return f, pr.value


# -----------------------------------------------------------------------------
# Checks
# -----------------------------------------------------------------------------

def density(A, tol=1e-8):
    return (np.abs(A) > tol * np.abs(A).max()).mean()


def tri_density(L, tol=1e-8):
    n = L.shape[0]
    m = np.abs(L[np.tril_indices(n)])
    return (m > tol * m.max()).mean()


def check_dense_and_exclusion(inst):
    dT = density(inst.Theta)
    dG = min(density(Gm) for Gm in inst.Gms)
    dL = min(tri_density(Lm) for Lm in inst.Lms)
    dJ = density(inst.G)
    assert min(dT, dG, dL, dJ) > 0.95, (dT, dG, dL, dJ)
    print(f"  [PASS] dense by nature: density@1e-8 Theta {dT:.3f} "
          f"G_m >= {dG:.3f} L_m' >= {dL:.3f} joint G {dJ:.3f}")
    # decimated control: D = 4, small
    D, Lh, Lf = 4, 64, 32
    h = RNG.standard_normal(Lh)
    P = Lh + Lf - 1
    Td = np.zeros((D * P, Lf))
    for r in range(D):
        for tau in range(P):
            for j in range(Lf):
                u = tau - j
                if 0 <= u < Lh and (u + r) % D == 0:
                    Td[r * P + tau, j] = h[u]
    Gd = Td.T @ Td
    off = np.abs(Gd[(np.subtract.outer(np.arange(Lf), np.arange(Lf)) % D) != 0])
    inc = Gd[(np.subtract.outer(np.arange(Lf), np.arange(Lf)) % D) == 0]
    assert off.max() == 0.0
    print(f"  [PASS] polyphase exclusion (D = 4 control): off-class Gram "
          f"entries EXACTLY 0 (max {off.max():.1f}); in-class density "
          f"{(np.abs(inc) > 1e-12).mean():.2f} — structural 1/D, whitening "
          f"cannot remove it")


def check_mint(seed=3):
    M, Lh, trev = 4, 121, 40
    rng = np.random.default_rng(seed)
    ch = [rng.standard_normal(Lh) * np.exp(-np.arange(Lh) / trev)
          for _ in range(M)]
    Lstar = int(np.ceil((Lh - 1) / (M - 1)))

    def els(chans, Lf):
        T = np.hstack([toep(h, Lf) for h in chans])
        g = np.zeros(T.shape[0])
        g[(len(chans[0]) + Lf) // 2] = 1.0
        return np.linalg.norm(g - T @ np.linalg.lstsq(T, g, rcond=None)[0])

    e_below, e_at = els(ch, Lstar - 1), els(ch, Lstar)
    assert e_below > 1e-3 and e_at < 1e-9, (e_below, e_at)
    print(f"  [PASS] MINT length threshold: eps_LS = {e_below:.2e} at "
          f"Lf* - 1 = {Lstar - 1} -> {e_at:.2e} at Lf* = {Lstar} "
          f"(Miyoshi-Kaneda)")
    # Common zeros must sit ON the unit circle to obstruct a delayed target:
    # for |z0| > 1 the obstruction functional decays like z0^{-d} and eps_LS
    # collapses despite non-coprimality. A shared spectral null at angle th
    # gives the ANALYTIC lower bound eps_LS >= ||proj_{span(cos th t, sin
    # th t)} e_d|| ~ sqrt(2/P), since range(T) lies in the codim-2 subspace
    # of polynomials vanishing at e^{+-i th}.
    th = 0.7 * np.pi
    chz = [np.convolve(h, [1.0, -2 * np.cos(th), 1.0]) for h in ch]
    Lfz = Lstar + 8
    Pz = len(chz[0]) + Lfz - 1
    tt = np.arange(Pz)
    V = np.stack([np.cos(th * tt), np.sin(th * tt)], axis=1)
    ed = np.zeros(Pz)
    ed[(len(chz[0]) + Lfz) // 2] = 1.0
    bound = np.linalg.norm(V @ np.linalg.lstsq(V, ed, rcond=None)[0])
    e_z = els(chz, Lfz)
    assert e_z >= bound * (1 - 1e-10) and e_z > 1e-2, (e_z, bound)
    print(f"  [PASS] Bezout falsification: shared spectral null at theta = "
          f"0.7 pi keeps eps_LS = {e_z:.4f} >= analytic bound "
          f"sqrt-projection {bound:.4f} ~ sqrt(2/P) = {np.sqrt(2/Pz):.4f} "
          f"above the length threshold (coprimality is necessary)")


def check_formulation(inst):
    # random-point identities
    f = RNG.standard_normal(inst.M * inst.Lf)
    lhs = np.linalg.norm(inst.T @ f - inst.g) ** 2
    rhs = np.linalg.norm(inst.L.T @ f - inst.b) ** 2 + inst.eLS ** 2
    ej = abs(lhs - rhs) / lhs
    m = inst.M // 2
    fm = RNG.standard_normal(inst.Lf)
    gm = inst.Ts[m] @ inst.fhm[m] + inst.rhat / inst.M
    lhs2 = np.linalg.norm(inst.Ts[m] @ fm - gm) ** 2
    rhs2 = (np.linalg.norm(inst.Lms[m].T @ fm - inst.bms[m]) ** 2
            + (inst.eLS / inst.M) ** 2)
    es = abs(lhs2 - rhs2) / lhs2
    assert max(ej, es) < 1e-12, (ej, es)
    print(f"  [PASS] whitening identities at random points: joint "
          f"{ej:.1e}, split {es:.1e} (T' rhat max "
          f"{np.abs(inst.T.T @ inst.rhat).max():.1e})")
    f_nat, v_nat = solve_split(inst)
    f_asm, v_asm = solve_assembly(inst)
    df = np.abs(f_nat - f_asm).max()
    dv = abs(v_nat - v_asm) / abs(v_nat)
    assert df < 1e-5 and dv < 1e-7, (df, dv)
    print(f"  [PASS] stalk assembly (the e14.jl --budget layout) == native "
          f"formulation: ||Df||_inf = {df:.1e}, rel Dv = {dv:.1e}")
    f_en, v_en = solve_exact(inst)
    f_ea, v_ea = solve_assembly_exact(inst)
    dfe = np.abs(f_en - f_ea).max()
    dve = abs(v_en - v_ea) / abs(v_en)
    assert dfe < 1e-4 and dve < 1e-6, (dfe, dve)
    print(f"  [PASS] stalk assembly (the e14.jl default (exact) layout) == native "
          f"formulation: ||Df||_inf = {dfe:.1e}, rel Dv = {dve:.1e}")
    return f_nat, v_nat, f_en, v_en


def check_trs(inst):
    f_c, v_c = solve_joint(inst)
    f_t, v_t, ms = solve_joint_trs(inst)
    df = np.abs(f_c - f_t).max()
    dv = abs(v_c - v_t) / abs(v_t)
    assert dv < 1e-6 and df < 1e-4, (dv, df)
    assert ms["mineig"] > 0 and ms["resid"] < 1e-10 and ms["stat"] < 1e-8
    print(f"  [PASS] TRS classical certificate: secular solve == conic "
          f"(||Df||_inf {df:.1e}, rel Dv {dv:.1e}); More-Sorensen: "
          f"nu = {ms['nu']:.4f}, boundary resid {ms['resid']:.1e}, "
          f"stationarity {ms['stat']:.1e}")
    return v_t


def check_split(inst, f_split, v_split, f_exact, v_exact, v_joint):
    dve = abs(v_exact - v_joint) / abs(v_joint)
    assert dve < 1e-6, dve
    use_e = np.linalg.norm(inst.L.T @ f_exact - inst.b) / inst.epsx
    use_b = np.linalg.norm(inst.L.T @ f_split - inst.b) / inst.epsx
    assert use_e > 0.99
    print(f"  [PASS] exact slice split == joint ball: rel Dv = {dve:.1e}, "
          f"ball use {use_e:.3f} (zero conservatism, coupling mass "
          f"~ M^2 Lf^2 / 2)")
    isi = np.linalg.norm(inst.T @ f_split - inst.g)
    assert isi <= inst.eps * (1 + 1e-6), (isi, inst.eps)
    assert v_split >= v_joint * (1 - 1e-7)
    print(f"  [PASS] budget-split certificate: ||T f - g|| = {isi:.4f} <= "
          f"eps = {inst.eps:.4f}, but ball use only {use_b:.3f} — the "
          f"triangle inequality throttles the design")
    ratios = []
    exs = (4 * inst.epsx, inst.epsx, inst.epsx / 4, inst.epsx / 16)
    for ex in exs:
        _, vs = solve_split(inst, epsx=ex)
        _, vj, _ = solve_joint_trs(inst, epsx=ex)
        ratios.append(vs / vj)
    assert all(a > b for a, b in zip(ratios, ratios[1:])), ratios
    assert ratios[-1] < 1.10, ratios
    print(f"  [PASS] budget-split conservatism (v_budget/v_joint): "
          + " -> ".join(f"{r:.3f}" for r in ratios)
          + f" over eps_x = {exs[0]} -> {exs[-1]} (-> 1 at the pinv "
            f"endpoint; the certified budget buys O(M) coupling mass)")


def check_allocation(inst, v_split):
    _, v_u = solve_split(inst, uniform=True)
    gain = (v_u - v_split) / v_u
    f_o, _ = solve_split(inst)
    s = [np.hypot(np.linalg.norm(inst.Lms[m].T @ f_o.reshape(inst.M, -1)[m]
                                 - inst.bms[m]), inst.eLS / inst.M)
         for m in range(inst.M)]
    spread = max(s) / min(s)
    assert gain > 0.005, gain
    print(f"  [PASS] allocation load-bearing: optimized budget beats uniform "
          f"by {100*gain:.1f}% (allocation spread x{spread:.2f})")


def check_deformation(inst):
    errs = []
    exs = (2e-2, 5e-3)
    for ex in exs:
        f, _ = solve_split(inst, epsx=ex)
        errs.append(np.abs(f - inst.fhat).max())
    rate = np.log(errs[1] / errs[0]) / np.log(exs[1] / exs[0])
    assert errs[1] < errs[0] and 0.7 < rate < 1.3, (errs, rate)
    print(f"  [PASS] deformation to pinv bank: ||f(eps_x) - fhat||_inf "
          f"{errs[0]:.2e} -> {errs[1]:.2e} over eps_x = {exs[0]} -> {exs[1]} "
          f"(rate eps_x^{rate:.2f})")


def check_noise_shaping(inst, f_col):
    f_w, _ = solve_exact_white(inst)
    move = np.linalg.norm(f_w - f_col) / np.linalg.norm(f_col)
    Qb = inst.Qbar()
    excess = (f_w @ Qb @ f_w) / (f_col @ Qb @ f_col)
    assert move > 0.05 and excess > 1.2, (move, excess)
    print(f"  [PASS] noise shaping load-bearing (exact form): white-noise "
          f"(diagonal-Q) control moves the solution by {100*move:.1f}% and "
          f"pays x{excess:.2f} amplified colored-noise power — the e14.jl "
          f"--white ablation row, pre-registered")


def check_equalization(inst, f_exact, sigman=0.05, N=200_000, seed=11):
    rng = np.random.default_rng(seed)
    x = rng.standard_normal(N)
    F = f_exact.reshape(inst.M, inst.Lf)
    Fh = inst.fhm

    def run(Fb):
        noise = 0.0
        for m in range(inst.M):
            n_m = sigman * inst.sigmas[m] * lfilter(
                [1.0], inst.ar_a, rng.standard_normal(N))
            noise = noise + np.convolve(Fb[m], n_m)[:N]
        c = sum(np.convolve(inst.chans[m], Fb[m]) for m in range(inst.M))
        isi = np.linalg.norm(c - inst.g)
        return isi, noise.var()

    rng = np.random.default_rng(seed)
    isi_d, npow_d = run(F)
    rng = np.random.default_rng(seed)
    isi_p, npow_p = run(Fh)
    Th = inst.Theta
    pred = (sum(inst.sigmas[m] ** 2 * F[m] @ Th @ F[m] for m in range(inst.M))
            / sum(inst.sigmas[m] ** 2 * Fh[m] @ Th @ Fh[m]
                  for m in range(inst.M)))
    meas = npow_d / npow_p
    assert meas < 0.7 and abs(meas - pred) / pred < 0.15
    print(f"  [PASS] equalization: designed bank cuts output noise power to "
          f"{meas:.3f}x pinv (predicted {pred:.3f}x) at ISI {isi_d:.4f} vs "
          f"{isi_p:.4f} (budget {inst.eps:.4f})")


def check_cross_solver(inst):
    f_c, _ = solve_exact(inst, solver_kw=CLARABEL)
    f_s, _ = solve_exact(inst, solver_kw=SCS)
    df = np.abs(f_c - f_s).max()
    assert df < 1e-4, df
    print(f"  [PASS] cross-solver (exact form): Clarabel vs SCS "
          f"||Df||_inf = {df:.1e}")


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

if __name__ == "__main__":
    inst = Instance()
    print(f"E14 oracle: M = {inst.M}, Lh = {inst.Lh}, Lf = {inst.Lf}, "
          f"d = {inst.d}, lam = {inst.lam}, eps_x = {inst.epsx}")
    print(f"  short-equalizer floor eps_LS = {inst.eLS:.4f} "
          f"(threshold Lf* = {int(np.ceil((inst.Lh - 1) / (inst.M - 1)))}, "
          f"||g|| = 1), eps = {inst.eps:.4f}")
    print("\n[1] dense by nature + polyphase exclusion")
    check_dense_and_exclusion(inst)
    print("\n[2] MINT / Bezout classical certificates")
    check_mint()
    print("\n[3] formulation certificate")
    f_split, v_split, f_exact, v_exact = check_formulation(inst)
    print("\n[4] TRS classical certificate (joint comparator)")
    v_joint = check_trs(inst)
    print("\n[5] exactness + conservatism")
    check_split(inst, f_split, v_split, f_exact, v_exact, v_joint)
    print("\n[6] allocation load-bearing (budget form)")
    check_allocation(inst, v_split)
    print("\n[7] deformation to the pinv endpoint")
    check_deformation(inst)
    print("\n[8] noise-shaping (diagonal-Q) control")
    check_noise_shaping(inst, f_exact)
    print("\n[9] equalization on simulated data")
    check_equalization(inst, f_exact)
    print("\n[10] cross-solver")
    check_cross_solver(inst)
    print("\nAll checks passed.")
