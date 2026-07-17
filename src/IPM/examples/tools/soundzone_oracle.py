# =============================================================================
# soundzone_oracle.py — self-contained validation for E13 (distributed
# sound-zone control: cooperative multizone rendering with local coupling —
# E11 tiled over a zone graph)
#
# Z zones on a path, each with an array of S loudspeakers and C control
# points. Acoustic paths are reverberant FIRs (length Lh, decay tau) whose
# inter-zone amplitude falls off as gamma^{|Dz|} with propagation onset
# d0*|Dz|. Every array within radius rho of a program COOPERATES in
# rendering it (pressure matching drives all nearby speakers per program —
# Betlehem-Zhang-Poulsen-Abhayapala, IEEE SPM 2015; Coleman et al.; car-
# cabin personal audio: Cheer & Elliott), so the stalks are f_{z',p} in
# R^{S*Lf}, one per (array z', program p) with |z' - p| <= rho. Per (zone z,
# program p), |p - z| <= rho, the composite response SUMS across arrays,
#
#     r_{z,p} = sum_{z' in W(p)} A_{z,z'} f_{z',p} - g_{z,p},
#
# with REFERENCE-FIELD bright targets g_{z,z} = A_{z,z} f_ref (the field of
# the own array under nominal driving — virtual-source pressure matching;
# delta targets would demand full dereverberation and push the floors to
# ~0.9 of the signal, admitting f = 0) and dark targets 0, constrained per
# ball ||r_{z,p}|| <= eps_{z,p} with SIGNAL-SCALE radii
#
#     eps_{z,p} = hypot(floor_{z,p}, eta * s_z),   s_z = ||g_{z,z}||,
#
# floor_{z,p} = ||r_{z,p}(f_hat)|| at the global stacked LS point f_hat.
# Objective: (1/2) sum f' Q_p f, Q_p the RADIATED-POWER Gram (mutual
# radiation resistance sinc(w*delta*|ds|), Morse-Ingard, integrated against
# program p's AR(2) PSD; heterogeneous peak frequencies across programs)
# plus the driver floor lam*I. PER-ARRAY ACOUSTIC-OUTPUT BUDGETS (SPL /
# noise-ordinance caps) couple the programs at each array:
#
#     sum_{p in P_a} f_{a,p}' Q_p f_{a,p} <= beta_a^2,
#     beta_a = theta * sqrt(power_a(f_hat)),  theta = 0.85,
#
# with the SAME radiated-power Grams as the objective — dense across
# speakers by the sinc radiation coupling. (A signal-power budget,
# I_S (x) Toeplitz, would have EXACTLY block-diagonal Cholesky ties,
# tri-density 1/S — a third structural-density observation, avoided here
# by the physics of acoustic rather than electrical caps.)
#
# feasible because the effort-minimal solution runs at 0.70-0.79 of f_hat's
# per-array powers (measured), and structurally load-bearing: WITHOUT them
# the problem decomposes EXACTLY by program (balls never mix programs, the
# objective is stalk-separable) — Z independent subproblems, fatal for a
# benchmark; WITH them the stalk-constraint graph is one component
# (check [6], a machine-checkable connectivity certificate).
#
# TWO RECORDED EXCLUSIONS (design-process negatives, in the polyphase
# tradition):
#   (i) ONE ARRAY PER PROGRAM makes the per-zone gluing EXACTLY block-
#       diagonal across neighborhood stalks — Cholesky tri-density
#       1/(2rho+1), structural, whitening-proof. Cooperative rendering
#       restores density 1.000 (check [1] measures both).
#  (ii) LONG-RANGE INFLUENCE IS NOT REQUIRED and should not be engineered:
#       the physics of local acoustic coupling is SCREENING — influence
#       decays geometrically at a gamma-tracked rate (check [7], the E09-
#       DARE species: decay as certificate). What a benchmark requires is
#       NON-DECOMPOSABILITY, a different property, certified by [6].
#
# LOCAL WHITENING + EXACT SLICE SPLIT (E11's machinery as the LOCAL
# structure at every ball). Per reproduction ball, G = M'M = LL' (tall, PD),
# ||r||^2 = ||L'phi - b||^2 + c^2, rows partition by window-stalk slices:
# k slice SOCs (s_j, w_j) with dense triangular ties, a (k+2)-dim hub
# (h, s_hat, sigma) with pins h = eps, sigma = c. Budget balls are already
# program-partitioned: |P_a| slices w_{a,p} = Rhat_p' f_{a,p} (dense
# Toeplitz-Cholesky ties) + a (|P_a|+1)-dim hub pinned at beta_a. Zero
# conservatism (check [5]); bounded arity — coupling mass O(Z), the
# pre-registered answer to E11's M-dial reading (slope 2.29 from O(M^2)
# global coupling).
#
# ENDPOINT. Radii at the global floors: any point in the intersection is a
# global LS minimizer, the stacked map has full column rank, so eta -> 0
# pins f_hat EXACTLY (budgets excluded: f_hat sits above the theta < 1
# caps) — measured rate eta^~1 (check [2]).
#
# Checks:
#   [1] DENSE BY NATURE + EXCLUSION (i).
#   [2] ENDPOINT: ||f(eta) - f_hat||_inf -> 0 at rate ~eta^1 (no-budget).
#   [3] FORMULATION CERTIFICATE: whitening identities (reproduction and
#       budget partitions) at random points; the stalk assembly (exactly
#       what e13/def.jl builds) == native cvxpy formulation.
#   [4] TRS (classical certificate): Z = 1, rho = 0, no-budget reduces to a
#       single-ball trust-region subproblem — secular solve (Gay 1981;
#       More-Sorensen 1983) == conic, MS conditions at the solution.
#   [5] EXACTNESS: slice-split assembly == joint per-ball native (zero
#       conservatism, values and states, budgets included on both sides).
#   [6] STRUCTURE: connectivity certificate (Z components without budgets —
#       the recorded decomposability exclusion — vs 1 with); gamma = 0
#       separates per ARRAY (stacked == joint); the non-cooperative
#       baseline violates the dark balls measurably.
#   [7] SCREENING (classical certificate, E09-DARE species): a bright-
#       target perturbation in zone 0 decays geometrically across arrays,
#       hop ratio tracking gamma.
#   [8] TRUNCATION REACH: uncontrolled leakage beyond the constrained
#       radius decays one power of gamma per zone of separation (measured
#       slope increment ~1.0), from a base STEEPER than the naive first-
#       neighbor gamma^1 — the dark balls at distance rho shadow the
#       emissions that would leak farther.
#   [9] DELIVERABLE: simulated AR programs through the true physics — per-
#       speaker signal power matches the spectral-Gram prediction; the
#       cooperative design improves RESPONSE-DOMAIN acoustic contrast (the
#       certified quantity) by +2.3 to +4.4 dB per zone over the non-
#       cooperative baseline; budget utilization reported.
#  [10] CROSS-SOLVER: Clarabel vs SCS on the assembled program.
#
# Figures quoted in e13/def.jl are from this file executed at Z = 6, S = 4,
# C = 3, Lf = 16, Lh = 96, rho = 1, gamma = 0.25, tau = 32, d0 = 4,
# lam = 1e-3, eta = 0.1, theta = 0.85.
# =============================================================================

import numpy as np
from scipy.signal import lfilter
from scipy.linalg import block_diag

RNG = np.random.default_rng(7)

CLARABEL = dict(solver="CLARABEL", tol_feas=1e-10, tol_gap_abs=1e-10,
                tol_gap_rel=1e-10)
SCS = dict(solver="SCS", eps=1e-9, max_iters=200_000)


# -----------------------------------------------------------------------------
# Instance
# -----------------------------------------------------------------------------

class Instance:
    def __init__(self, Z=6, S=4, C=3, Lf=16, Lh=96, rho=1, gam=0.25,
                 tau=32.0, d0=4, delta=2.0, lam=1e-3, eta=0.1, theta=0.85,
                 dr=4, r_ar=0.9, seed=7, whiten=True, budget=True):
        self.__dict__.update(Z=Z, S=S, C=C, Lf=Lf, Lh=Lh, rho=rho, gam=gam,
                             tau=tau, d0=d0, delta=delta, lam=lam, eta=eta,
                             theta=theta, dr=dr, r_ar=r_ar, seed=seed,
                             budget=budget)
        self.nf = S * Lf
        self.P = Lh + Lf - 1
        rng = np.random.default_rng(seed)
        self.raw = {(z, zp): [[self._rir0(z - zp, rng) for s in range(S)]
                              for c in range(C)]
                    for z in range(Z) for zp in range(Z)}
        self.phis = np.linspace(0.3 * np.pi, 0.7 * np.pi, max(Z, 2))[:Z]
        self.Qs = [self._effort_gram(p) for p in range(Z)]
        self.Rhat = [np.linalg.cholesky(Q) for Q in self.Qs]
        self._build(gam, whiten)

    # ---- physics ----
    def _rir0(self, dz, rng):
        h = np.zeros(self.Lh)
        on = min(self.d0 * abs(dz), self.Lh - 8)
        n = self.Lh - on
        h[on:] = rng.standard_normal(n) * np.exp(-np.arange(n) / self.tau)
        return h

    def rir(self, z, zp, c, s, gam=None):
        gam = self.gam if gam is None else gam
        return gam ** abs(z - zp) * self.raw[z, zp][c][s]

    def _toep(self, h):
        T = np.zeros((self.P, self.Lf))
        for j in range(self.Lf):
            T[j:j + self.Lh, j] = h
        return T

    def Ablock(self, z, zp, gam=None):
        return np.block([[self._toep(self.rir(z, zp, c, s, gam))
                          for s in range(self.S)] for c in range(self.C)])

    def ar_coeffs(self, p):
        return np.array([1.0, -2 * self.r_ar * np.cos(self.phis[p]),
                         self.r_ar ** 2])

    def spec_rho(self, p, K):
        psi = lfilter([1.0], self.ar_coeffs(p), np.eye(1, 8192).ravel())
        rho = np.array([psi[:8192 - k] @ psi[k:] for k in range(K)])
        return rho / rho[0]

    def _effort_gram(self, p):
        S, Lf, nf = self.S, self.Lf, self.nf
        w = np.linspace(1e-4, np.pi, 2048)
        Sw = 1.0 / np.abs(np.polyval(self.ar_coeffs(p)[::-1],
                                     np.exp(-1j * w))) ** 2
        Sw /= Sw.max()
        Q = np.zeros((nf, nf))
        for s in range(S):
            for sp in range(s, S):
                kern = np.sinc(w * self.delta * (sp - s) / np.pi)
                for lag in range(Lf):
                    v = np.trapezoid(Sw * kern * np.cos(w * lag), w) / np.pi
                    for i in range(Lf - lag):
                        Q[s * Lf + i, sp * Lf + i + lag] = v
                        Q[sp * Lf + i + lag, s * Lf + i] = v
                        Q[s * Lf + i + lag, sp * Lf + i] = v
                        Q[sp * Lf + i, s * Lf + i + lag] = v
        return Q + self.lam * np.eye(nf)

    # ---- structure, targets, whitening, budgets ----
    def window(self, p):
        return [z for z in range(self.Z) if abs(z - p) <= self.rho]

    def programs_at(self, a):
        return [p for p in range(self.Z) if a in self.window(p)]

    def _build(self, gam, whiten):
        Z, C, P, nf = self.Z, self.C, self.P, self.nf
        self.stalks = [(zp, p) for p in range(Z) for zp in self.window(p)]
        self.sidx = {sk: i for i, sk in enumerate(self.stalks)}
        self.balls = [(z, p) for p in range(Z) for z in range(Z)
                      if abs(p - z) <= self.rho]
        f_ref = np.zeros(nf)
        for s in range(self.S):
            f_ref[s * self.Lf + self.dr] = 1.0 / self.S
        self.f_ref = f_ref
        self.Mb, self.gb, self.colsb = {}, {}, {}
        for (z, p) in self.balls:
            cols = [(zp, p) for zp in self.window(p)]
            self.Mb[z, p] = np.hstack([self.Ablock(z, zp, gam)
                                       for zp, _ in cols])
            self.gb[z, p] = (self.Ablock(z, z, gam) @ f_ref if p == z
                             else np.zeros(C * P))
            self.colsb[z, p] = cols
        self.scale = {z: np.linalg.norm(self.gb[z, z]) for z in range(Z)}
        # global stacked LS point (normal equations accumulated blockwise)
        nstk = len(self.stalks)
        Gg = np.zeros((nstk * nf, nstk * nf))
        hg = np.zeros(nstk * nf)
        for ball in self.balls:
            M, g, cols = self.Mb[ball], self.gb[ball], self.colsb[ball]
            for i, ci in enumerate(cols):
                ii = self.sidx[ci]
                hg[ii * nf:(ii + 1) * nf] += \
                    M[:, i * nf:(i + 1) * nf].T @ g
                for j, cj in enumerate(cols):
                    jj = self.sidx[cj]
                    Gg[ii * nf:(ii + 1) * nf, jj * nf:(jj + 1) * nf] += \
                        M[:, i * nf:(i + 1) * nf].T @ M[:, j * nf:(j + 1) * nf]
        self.mineig_global = np.linalg.eigvalsh(Gg).min()
        self.fhat = np.linalg.solve(Gg, hg)
        self.floors = {}
        for ball in self.balls:
            phi = self.gather(self.fhat, ball)
            self.floors[ball] = np.linalg.norm(
                self.Mb[ball] @ phi - self.gb[ball])
        # amplifier budgets from f_hat's per-array signal powers
        pw = self.powers(self.fhat)
        self.betas = self.theta * np.sqrt(pw)
        # per-ball whitening (gamma = 0 makes G structurally singular)
        self.Lb, self.bb, self.cb = {}, {}, {}
        if not whiten:
            return
        for ball in self.balls:
            M, g = self.Mb[ball], self.gb[ball]
            L = np.linalg.cholesky(M.T @ M)
            b = np.linalg.solve(L, M.T @ g)
            self.Lb[ball], self.bb[ball] = L, b
            self.cb[ball] = np.sqrt(max(g @ g - b @ b, 0.0))

    def gather(self, f, ball):
        return np.concatenate([f[self.sidx[c] * self.nf:
                                 (self.sidx[c] + 1) * self.nf]
                               for c in self.colsb[ball]])

    def eps(self, ball, eta=None):
        eta = self.eta if eta is None else eta
        return np.hypot(self.floors[ball], eta * self.scale[ball[0]])

    def powers(self, f):
        pw = np.zeros(self.Z)
        for i, (zp, p) in enumerate(self.stalks):
            x = f[i * self.nf:(i + 1) * self.nf]
            pw[zp] += x @ self.Qs[p] @ x
        return pw

    def Qbar(self):
        return block_diag(*[self.Qs[p] for (_, p) in self.stalks])


# -----------------------------------------------------------------------------
# Solvers
# -----------------------------------------------------------------------------

def solve_native(inst, eta=None, budget=None, only_own=False, gpert=None,
                 gam_maps=None, solver_kw=None):
    """Raw per-ball SOC balls + per-array budget SOCs (the native form).
    only_own drops darks AND budgets (the non-cooperative baseline).
    gpert = (ball, g_new, eps_add) perturbs one ball's target."""
    import cvxpy as cp
    budget = inst.budget if budget is None else budget
    nf = inst.nf
    f = cp.Variable(len(inst.stalks) * nf)
    cons = []
    for ball in inst.balls:
        z, p = ball
        if only_own and p != z:
            continue
        g, eps = inst.gb[ball], inst.eps(ball, eta)
        if gpert is not None and ball == gpert[0]:
            g, eps = gpert[1], eps + gpert[2]
        M = inst.Mb[ball]
        expr = sum(M[:, j * nf:(j + 1) * nf]
                   @ f[inst.sidx[c] * nf:(inst.sidx[c] + 1) * nf]
                   for j, c in enumerate(inst.colsb[ball])) - g
        cons.append(cp.SOC(eps, expr))
    if budget and not only_own:
        for a in range(inst.Z):
            terms = [inst.Rhat[p].T @ f[inst.sidx[a, p] * nf:
                                        (inst.sidx[a, p] + 1) * nf]
                     for p in inst.programs_at(a)]
            cons.append(cp.SOC(inst.betas[a], cp.hstack(terms)))
    obj = 0.5 * cp.sum([cp.quad_form(f[i * nf:(i + 1) * nf],
                                     cp.psd_wrap(inst.Qs[p]))
                        for i, (_, p) in enumerate(inst.stalks)])
    pr = cp.Problem(cp.Minimize(obj), cons)
    pr.solve(**(solver_kw or CLARABEL))
    assert pr.status in ("optimal", "optimal_inaccurate"), pr.status
    return f.value, pr.value


def assemble(inst, eta=None):
    """The exact (Q, B, g, cones) layout e13/def.jl builds. Stalk order:
    filter stalks (array, program) | per reproduction ball: k slice SOCs
    (s_j, w_j), dim nf+1, then hub (h, s_hat_1..k, sigma), dim k+2 | per
    array: |P_a| budget slices (u_p, v_p), dim nf+1, then hub (h, u_hat),
    dim |P_a|+1. Edges: slice ties (dense triangular L' blocks, nonzero
    rhs), budget ties (dense Rhat' blocks, zero rhs), s-links, pins
    (h = eps, sigma = c; h = beta)."""
    nf = inst.nf
    nstk = len(inst.stalks)
    off = nstk * nf
    rsl, rhub, bsl, bhub = {}, {}, {}, {}
    for ball in inst.balls:
        k = len(inst.colsb[ball])
        for j in range(k):
            rsl[ball, j] = slice(off, off + nf + 1)
            off += nf + 1
        rhub[ball] = slice(off, off + k + 2)
        off += k + 2
    for a in range(inst.Z):
        Pa = inst.programs_at(a)
        for p in Pa:
            bsl[a, p] = slice(off, off + nf + 1)
            off += nf + 1
        bhub[a] = slice(off, off + len(Pa) + 1)
        off += len(Pa) + 1
    n0 = off
    nrows = (sum(len(inst.colsb[b]) * (nf + 2) + 2 for b in inst.balls)
             + sum(len(inst.programs_at(a)) * (nf + 1) + 1
                   for a in range(inst.Z)))
    B = np.zeros((nrows, n0))
    gv = np.zeros(nrows)
    r = 0
    for ball in inst.balls:
        k = len(inst.colsb[ball])
        LT = inst.Lb[ball].T
        for j in range(k):                        # reproduction slice ties
            sl = rsl[ball, j]
            B[r:r + nf, sl.start + 1:sl.stop] = np.eye(nf)
            for m in range(j, k):
                cm = inst.sidx[inst.colsb[ball][m]]
                B[r:r + nf, cm * nf:(cm + 1) * nf] = \
                    -LT[j * nf:(j + 1) * nf, m * nf:(m + 1) * nf]
            gv[r:r + nf] = -inst.bb[ball][j * nf:(j + 1) * nf]
            r += nf
        for j in range(k):                        # s-links
            B[r, rhub[ball].start + 1 + j] = 1.0
            B[r, rsl[ball, j].start] = -1.0
            r += 1
        B[r, rhub[ball].start] = 1.0              # hub pin h = eps
        gv[r] = inst.eps(ball, eta)
        r += 1
        B[r, rhub[ball].start + k + 1] = 1.0      # sigma pin = c
        gv[r] = inst.cb[ball]
        r += 1
    for a in range(inst.Z):
        Pa = inst.programs_at(a)
        for jj, p in enumerate(Pa):               # budget slice ties
            sl = bsl[a, p]
            B[r:r + nf, sl.start + 1:sl.stop] = np.eye(nf)
            ia = inst.sidx[a, p]
            B[r:r + nf, ia * nf:(ia + 1) * nf] = -inst.Rhat[p].T
            r += nf
        for jj, p in enumerate(Pa):               # u-links
            B[r, bhub[a].start + 1 + jj] = 1.0
            B[r, bsl[a, p].start] = -1.0
            r += 1
        B[r, bhub[a].start] = 1.0                 # hub pin h = beta
        gv[r] = inst.betas[a]
        r += 1
    Q = np.zeros((n0, n0))
    for i, (_, p) in enumerate(inst.stalks):
        Q[i * nf:(i + 1) * nf, i * nf:(i + 1) * nf] = inst.Qs[p]
    socs = ([rsl[b, j] for b in inst.balls
             for j in range(len(inst.colsb[b]))]
            + [rhub[b] for b in inst.balls]
            + [bsl[a, p] for a in range(inst.Z)
               for p in inst.programs_at(a)]
            + [bhub[a] for a in range(inst.Z)])
    return Q, B, gv, n0, socs


def solve_assembly(inst, eta=None, solver_kw=None):
    import cvxpy as cp
    Q, B, gv, n0, socs = assemble(inst, eta)
    x = cp.Variable(n0)
    cons = [B @ x == gv]
    for sl in socs:
        cons.append(cp.SOC(x[sl.start], x[sl.start + 1:sl.stop]))
    pr = cp.Problem(cp.Minimize(0.5 * cp.quad_form(x, cp.psd_wrap(Q))), cons)
    pr.solve(**(solver_kw or CLARABEL))
    assert pr.status in ("optimal", "optimal_inaccurate"), pr.status
    return x.value[:len(inst.stalks) * inst.nf], pr.value


def solve_trs_single(inst1):
    """Z = 1, rho = 0, no budget: single-ball trust-region subproblem by the
    secular equation (Gay 1981; More-Sorensen 1983)."""
    ball = inst1.balls[0]
    L, b = inst1.Lb[ball], inst1.bb[ball]
    Q = inst1.Qs[0]
    ex = np.sqrt(max(inst1.eps(ball) ** 2 - inst1.cb[ball] ** 2, 0.0))
    Linv = np.linalg.inv(L)
    H = Linv @ Q @ Linv.T
    H = 0.5 * (H + H.T)
    lamH, U = np.linalg.eigh(H)
    beta = U.T @ (H @ b)
    assert np.linalg.norm(b) > ex

    def ynorm(nu):
        return np.sqrt(np.sum((beta / (lamH + nu)) ** 2))

    lo, hi = 0.0, 1.0
    while ynorm(hi) > ex:
        hi *= 4
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if ynorm(mid) > ex else (lo, mid)
    nu = 0.5 * (lo + hi)
    y = U @ (-beta / (lamH + nu))
    f = np.linalg.solve(L.T, y + b)
    ms = dict(nu=nu, resid=abs(np.linalg.norm(y) - ex),
              stat=np.linalg.norm(H @ (y + b) + nu * y))
    return f, 0.5 * f @ Q @ f, ms


# -----------------------------------------------------------------------------
# Checks
# -----------------------------------------------------------------------------

def density(A, tol=1e-8):
    return (np.abs(A) > tol * np.abs(A).max()).mean()


def tridensity(L, tol=1e-8):
    v = np.abs(L[np.tril_indices(L.shape[0])])
    return (v > tol * v.max()).mean()


def check_dense_and_exclusion(inst):
    dQ = min(density(Q) for Q in inst.Qs)
    dG = min(density(inst.Mb[b].T @ inst.Mb[b]) for b in inst.balls)
    dL = min(tridensity(inst.Lb[b]) for b in inst.balls)
    dR = min(tridensity(R) for R in inst.Rhat)
    ev = min(np.linalg.eigvalsh(inst.Mb[b].T @ inst.Mb[b]).min()
             for b in inst.balls)
    assert min(dQ, dG, dL, dR) > 0.95 and ev > 1e-8, (dQ, dG, dL, dR, ev)
    print(f"  [PASS] dense by nature: density@1e-8 Q_p >= {dQ:.3f}, "
          f"per-ball G >= {dG:.3f}, reproduction ties >= {dL:.3f}, "
          f"budget ties >= {dR:.3f} (min eig {ev:.1e}; tall "
          f"{inst.C * inst.P} rows vs <= {(2 * inst.rho + 1) * inst.nf} "
          f"cols)")
    z = inst.Z // 2
    k = 2 * inst.rho + 1
    Ms = []
    for i, zp in enumerate(inst.window(z)):
        Ms.append(np.hstack([inst.Ablock(z, zq) if q == i else
                             np.zeros((inst.C * inst.P, inst.nf))
                             for q, zq in enumerate(inst.window(z))]))
    Gx = (V := np.vstack(Ms)).T @ V
    td = tridensity(np.linalg.cholesky(Gx + 1e-12 * np.eye(Gx.shape[0])))
    off = max(np.abs(Gx[i * inst.nf:(i + 1) * inst.nf,
                        j * inst.nf:(j + 1) * inst.nf]).max()
              for i in range(k) for j in range(k) if i != j)
    assert td < 1.2 / k and off == 0.0, (td, off)
    print(f"  [PASS] exclusion (i), one array per program: cross-stalk Gram "
          f"blocks EXACTLY 0, Cholesky tri-density {td:.3f} ~ 1/(2rho+1) = "
          f"{1 / k:.3f} — structural, whitening-proof; cooperative "
          f"rendering restores {dL:.3f}")


def check_endpoint(inst):
    errs, etas = [], (2e-2, 5e-3)
    for eta in etas:
        f, _ = solve_native(inst, eta=eta, budget=False)
        errs.append(np.abs(f - inst.fhat).max())
    rate = np.log(errs[1] / errs[0]) / np.log(etas[1] / etas[0])
    assert errs[1] < errs[0] and 0.7 < rate < 1.3, (errs, rate)
    print(f"  [PASS] endpoint (no-budget): radii at the GLOBAL floors pin "
          f"f_hat (stacked min eig {inst.mineig_global:.2f}); "
          f"||f - f_hat||_inf {errs[0]:.2e} -> {errs[1]:.2e} "
          f"(rate eta^{rate:.2f})")


def check_formulation(inst):
    rng = np.random.default_rng(3)
    worst, cfl = 0.0, True
    for ball in inst.balls:
        phi = rng.standard_normal(inst.Mb[ball].shape[1])
        lhs = np.linalg.norm(inst.Mb[ball] @ phi - inst.gb[ball]) ** 2
        rhs = (np.linalg.norm(inst.Lb[ball].T @ phi - inst.bb[ball]) ** 2
               + inst.cb[ball] ** 2)
        worst = max(worst, abs(lhs - rhs) / max(lhs, 1e-30))
        cfl &= inst.cb[ball] <= inst.floors[ball] + 1e-10
    a = inst.Z // 2
    x = {p: rng.standard_normal(inst.nf) for p in inst.programs_at(a)}
    lhs = sum(x[p] @ inst.Qs[p] @ x[p] for p in x)
    rhs = sum(np.linalg.norm(inst.Rhat[p].T @ x[p]) ** 2 for p in x)
    wb = abs(lhs - rhs) / lhs
    assert worst < 1e-11 and wb < 1e-12 and cfl, (worst, wb, cfl)
    print(f"  [PASS] whitening identities at random points: reproduction "
          f"{worst:.1e}, budget partition {wb:.1e}; local floors "
          f"c <= global floors on every ball")
    f_nat, v_nat = solve_native(inst)
    f_asm, v_asm = solve_assembly(inst)
    df = np.abs(f_nat - f_asm).max()
    dv = abs(v_nat - v_asm) / abs(v_nat)
    assert df < 1e-4 and dv < 1e-5, (df, dv)
    print(f"  [PASS] stalk assembly (the e13/def.jl layout) == native "
          f"formulation: ||Df||_inf = {df:.1e}, rel Dv = {dv:.1e}")
    return f_nat, v_nat, (df, dv)


def check_trs():
    inst1 = Instance(Z=1, rho=0, budget=False)
    f_c, v_c = solve_native(inst1, budget=False)
    f_t, v_t, ms = solve_trs_single(inst1)
    dv = abs(v_c - v_t) / abs(v_t)
    assert dv < 1e-6 and ms["resid"] < 1e-10 and ms["stat"] < 1e-8, (dv, ms)
    print(f"  [PASS] TRS classical certificate (Z = 1, rho = 0, no-budget "
          f"reduction): secular == conic, rel Dv {dv:.1e}; More-Sorensen: "
          f"nu = {ms['nu']:.3f}, boundary {ms['resid']:.1e}, stationarity "
          f"{ms['stat']:.1e}")


def check_exactness(dfdv):
    print(f"  [PASS] exactness: slice-split assembly == joint per-ball "
          f"native, budgets included (zero conservatism): ||Df||_inf = "
          f"{dfdv[0]:.1e}, rel Dv = {dfdv[1]:.1e}")


def components(inst, budget):
    parent = list(range(len(inst.stalks)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i, j):
        parent[find(i)] = find(j)

    for ball in inst.balls:
        ids = [inst.sidx[c] for c in inst.colsb[ball]]
        for i in ids[1:]:
            union(ids[0], i)
    if budget:
        for a in range(inst.Z):
            ids = [inst.sidx[a, p] for p in inst.programs_at(a)]
            for i in ids[1:]:
                union(ids[0], i)
    return len({find(i) for i in range(len(inst.stalks))})


def check_structure(inst):
    c0, c1 = components(inst, False), components(inst, True)
    assert c0 == inst.Z and c1 == 1, (c0, c1)
    print(f"  [PASS] connectivity certificate (exclusion (ii) recorded): "
          f"{c0} components without budgets — the problem decomposes "
          f"EXACTLY by program — vs {c1} with the amplifier budgets")
    inst0 = Instance(Z=inst.Z, gam=0.0, seed=inst.seed, whiten=False)
    f_j, v_j = solve_native(inst0)
    import cvxpy as cp
    f_s = np.zeros_like(f_j)
    v_s = 0.0
    for a in range(inst0.Z):                      # per-array subproblems
        Pa = inst0.programs_at(a)
        fs = {p: cp.Variable(inst0.nf) for p in Pa}
        cons = []
        for p in Pa:                              # single-stalk balls at a
            for z in inst0.window(p):
                if z != a:                        # A[z,zp]=0 unless z=zp
                    continue
                j = inst0.colsb[z, p].index((a, p))
                M = inst0.Mb[z, p][:, j * inst0.nf:(j + 1) * inst0.nf]
                cons.append(cp.SOC(inst0.eps((z, p)),
                                   M @ fs[p] - inst0.gb[z, p]))
        cons.append(cp.SOC(inst0.betas[a],
                           cp.hstack([inst0.Rhat[p].T @ fs[p] for p in Pa])))
        pr = cp.Problem(cp.Minimize(0.5 * cp.sum(
            [cp.quad_form(fs[p], cp.psd_wrap(inst0.Qs[p])) for p in Pa])),
            cons)
        pr.solve(**CLARABEL)
        v_s += pr.value
        for p in Pa:
            i = inst0.sidx[a, p]
            f_s[i * inst0.nf:(i + 1) * inst0.nf] = fs[p].value
    dv = abs(v_j - v_s) / abs(v_j)
    df = np.abs(f_j - f_s).max()
    assert dv < 1e-6 and df < 1e-4, (dv, df)
    print(f"  [PASS] gamma = 0 separation: stacked PER-ARRAY solves == "
          f"joint (rel Dv {dv:.1e}, ||Df||_inf {df:.1e}) — the budgets set "
          f"the gamma = 0 granularity")
    f_nc, _ = solve_native(inst, only_own=True)
    viol = max(np.linalg.norm(inst.Mb[b] @ inst.gather(f_nc, b)
                              - inst.gb[b]) / inst.eps(b)
               for b in inst.balls if b[0] != b[1])
    assert viol > 1.5, viol
    print(f"  [PASS] coupling load-bearing: the non-cooperative baseline "
          f"violates the dark balls by up to x{viol:.1f}")
    return f_nc


def check_screening(inst):
    ratios = {}
    for gam in (0.25, 0.4):
        k = Instance(Z=inst.Z, gam=gam, seed=inst.seed)
        f0, _ = solve_native(k)
        gpert = k.gb[0, 0].copy()
        gpert[k.C * k.P // (2 * k.C)] += 0.5 * k.scale[0]
        f1, _ = solve_native(k, gpert=((0, 0), gpert, 0.5 * k.scale[0]))
        mv = [max(np.abs(f1[i * k.nf:(i + 1) * k.nf]
                         - f0[i * k.nf:(i + 1) * k.nf]).max()
                  for i, (zp, _) in enumerate(k.stalks) if zp == zt)
              for zt in range(k.Z)]
        ratios[gam] = (mv[1] / mv[0], mv)
    r25, r40 = ratios[0.25][0], ratios[0.4][0]
    assert r25 < r40 < 0.6, (r25, r40)
    print(f"  [PASS] screening (E09-DARE species): hop-1 influence ratio "
          f"{r25:.2f} at gamma 0.25 vs {r40:.2f} at gamma 0.4 — geometric "
          f"decay at a gamma-tracked rate; long-range influence is "
          f"screening-suppressed BY PHYSICS, non-decomposability is "
          f"certified by connectivity, not by reach")


def check_truncation(inst):
    gams = (0.15, 0.25, 0.35)
    leak = {1: [], 2: []}
    for gam in gams:
        k = Instance(Z=inst.Z, gam=gam, seed=inst.seed)
        f, _ = solve_native(k)
        for dist in (1, 2):
            worst = 0.0
            for p in range(k.Z):
                for z in range(k.Z):
                    if abs(z - p) != k.rho + dist:
                        continue
                    resid = np.zeros(k.C * k.P)
                    for zp in k.window(p):
                        i = k.sidx[zp, p]
                        resid += k.Ablock(z, zp, gam) @ \
                            f[i * k.nf:(i + 1) * k.nf]
                    worst = max(worst, np.linalg.norm(resid) / k.scale[z])
            leak[dist].append(worst)
    s1 = np.polyfit(np.log(gams), np.log(leak[1]), 1)[0]
    s2 = np.polyfit(np.log(gams), np.log(leak[2]), 1)[0]
    assert s1 > 0.8 and 0.7 < s2 - s1 < 1.3, (s1, s2)
    print(f"  [PASS] truncation reach: uncontrolled leakage scales as "
          f"gamma^{s1:.2f} at distance rho+1 and gamma^{s2:.2f} at rho+2 — "
          f"increment {s2 - s1:.2f} (each zone of separation costs one "
          f"power of gamma); the base exceeds the naive first-neighbor "
          f"gamma^1 because the dark balls at distance rho already SHADOW "
          f"the emissions that would leak farther")


def check_deliverable(inst, f, f_nc, N=100_000, seed=11):
    rng = np.random.default_rng(seed)
    progs = [lfilter([1.0], inst.ar_coeffs(p), rng.standard_normal(N))
             for p in range(inst.Z)]
    p = inst.Z // 2
    i0 = inst.sidx[p, p]
    fblk = f[i0 * inst.nf:(i0 + 1) * inst.nf]
    rho_sp = inst.spec_rho(p, inst.Lf)
    Th = rho_sp[np.abs(np.subtract.outer(np.arange(inst.Lf),
                                         np.arange(inst.Lf)))]
    x = progs[p] / progs[p].std()
    meas = sum(np.convolve(fblk[s * inst.Lf:(s + 1) * inst.Lf], x)[:N].var()
               for s in range(inst.S))
    pred = sum(fblk[s * inst.Lf:(s + 1) * inst.Lf] @ Th
               @ fblk[s * inst.Lf:(s + 1) * inst.Lf]
               for s in range(inst.S))
    rel = abs(meas - pred) / pred
    assert rel < 0.1, (meas, pred)
    print(f"  [PASS] deliverable (i): per-speaker signal power matches the "
          f"spectral-Gram prediction (rel {rel:.3f})")

    def resp_contrast(fv):
        vals = []
        for z in range(inst.Z):
            bright = np.linalg.norm(inst.Mb[z, z] @ inst.gather(fv, (z, z)))
            dark2 = sum(np.linalg.norm(inst.Mb[z, p]
                                       @ inst.gather(fv, (z, p))) ** 2
                        for p in range(inst.Z)
                        if p != z and abs(p - z) <= inst.rho)
            vals.append(10 * np.log10(bright ** 2 / dark2))
        return np.array(vals)

    def sig_contrast(fv):
        z, c = inst.Z // 2, 0
        bright = dark = 0.0
        for p2 in range(inst.Z):
            if abs(p2 - z) > inst.rho and p2 != z:
                continue
            y = np.zeros(N)
            for zp in inst.window(p2):
                i = inst.sidx[zp, p2]
                for s in range(inst.S):
                    fs = fv[i * inst.nf + s * inst.Lf:
                            i * inst.nf + (s + 1) * inst.Lf]
                    y += np.convolve(np.convolve(inst.rir(z, zp, c, s), fs),
                                     progs[p2])[:N]
            if p2 == z:
                bright = y.var()
            else:
                dark += y.var()
        return 10 * np.log10(bright / dark)
    gain = resp_contrast(f) - resp_contrast(f_nc)
    util = np.sqrt(inst.powers(f)) / inst.betas
    assert gain.min() > 1.5 and gain.mean() > 2.5, gain
    assert util.max() < 1.0 + 1e-6, util
    print(f"  [PASS] deliverable (ii): response-domain acoustic contrast "
          f"(the certified quantity) improves by "
          f"[+{gain.min():.1f}, +{gain.max():.1f}] dB per zone (mean "
          f"+{gain.mean():.1f}) over the non-cooperative baseline; budget "
          f"utilization in [{util.min():.2f}, {util.max():.2f}]. "
          f"Signal-domain contrast at the center point: "
          f"{sig_contrast(f):.1f} vs {sig_contrast(f_nc):.1f} dB — the "
          f"program spectra concentrate power away from the suppressed "
          f"band, diluting the single-point figure (reported, not gated)")


def check_cross_solver(inst):
    f_c, _ = solve_assembly(inst, solver_kw=CLARABEL)
    f_s, _ = solve_assembly(inst, solver_kw=SCS)
    df = np.abs(f_c - f_s).max()
    assert df < 1e-4, df
    print(f"  [PASS] cross-solver: Clarabel vs SCS ||Df||_inf = {df:.1e}")


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

if __name__ == "__main__":
    inst = Instance()
    print(f"E13 oracle: Z = {inst.Z}, S = {inst.S}, C = {inst.C}, "
          f"Lf = {inst.Lf}, Lh = {inst.Lh}, rho = {inst.rho}, gamma = "
          f"{inst.gam}, eta = {inst.eta}, theta = {inst.theta}")
    fb = [inst.floors[z, z] / inst.scale[z] for z in range(inst.Z)]
    fd = [inst.floors[b] / inst.scale[b[0]] for b in inst.balls
          if b[0] != b[1]]
    print(f"  {len(inst.stalks)} stalks of dim {inst.nf}, "
          f"{len(inst.balls)} reproduction balls + {inst.Z} budgets; "
          f"bright floors/scale [{min(fb):.3f}, {max(fb):.3f}], "
          f"dark [{min(fd):.3f}, {max(fd):.3f}]")
    print("\n[1] dense by nature + exclusion (i)")
    check_dense_and_exclusion(inst)
    print("\n[2] endpoint at the global floors")
    check_endpoint(inst)
    print("\n[3] formulation certificate")
    f_nat, v_nat, dfdv = check_formulation(inst)
    print("\n[4] TRS classical certificate (Z = 1 reduction)")
    check_trs()
    print("\n[5] exactness of the slice split")
    check_exactness(dfdv)
    print("\n[6] structure: connectivity, gamma = 0 separation, "
          "load-bearing")
    f_nc = check_structure(inst)
    print("\n[7] screening (exclusion (ii) resolved)")
    check_screening(inst)
    print("\n[8] truncation reach")
    check_truncation(inst)
    print("\n[9] deliverable: simulated programs through the true physics")
    check_deliverable(inst, f_nat, f_nc)
    print("\n[10] cross-solver")
    check_cross_solver(inst)
    print("\nAll checks passed.")
