# =============================================================================
# torus_lyapunov_oracle.py — self-contained validation for E13 (compositional
# Lyapunov certification on a K x K torus of subsystems; Agler-split
# dense SDP with Gramian-weighted dense Q)
#
# The suite's control-theory pillar. A K x K torus of coupled subsystems,
# each with DENSE internal dynamics A_loc in R^{d x d} (skew + dense SPD
# damping) and dense nearest-neighbor couplings eps*{Cp, Cm, Dp, Dm}, is
# certified stable compositionally:
#
#   min  sum_v [ tr(P_v) + (gamma/2) || W^{1/2} P_v W^{1/2} ||_F^2 ]
#   s.t. A'P + PA + I + S = 0,   P = diag(P_v),  P_v >= 0,
#        S = sum_edges embed(S_e),  S_e >= 0  (2d x 2d, per torus edge),
#
# where W solves A_loc W + W A_loc' + I = 0 (the local controllability
# Gramian — a canonical, dynamics-derived, fully dense weight). The
# quadratic term is the flagship (gamma = 0.1): it makes the program
# strictly convex in P, so the optimal certificate is UNIQUE (Tikhonov
# selection on the linear SDP's optimal face) and cross-solver gates can
# compare STATES, not just values — a first for the suite's SDP examples.
# gamma -> 0 recovers the pure linear certificate at measured rate
# gamma^1.00 (the textbook first-order perturbation rate).
#
# WHY THIS EXAMPLE EXISTS — the conjunction of the suite's two measured
# win-conditions, on the topology that E13's predecessor proved is
# survivable only when blocks are small. (1) DENSE Q: the Gramian-weighted
# quadratic supplies dense 10-dim svec Q blocks on every P_v stalk —
# dense Q is free for the sheaf solver (it factors b^3 per block
# regardless) and a pure tax on baselines (dense cliques in their KKT
# fill, or dense factor rows after conic reformulation). (2) SMALL BLOCKS,
# BALANCED SIDES: vertex stalks svec 10 (P_v) and 36 (S_e) against row
# blocks of 10-16 — b_v ~ b_e, so no side is thin and nobody inherits a
# cubed-constant advantage (the failure mode of the quantum-marginal
# torus, whose 136-dim stalks against 10-dim rows handed Mosek a
# ~2500x smaller cubic constant). And coker B = 0 BY THE PRIVATE-SLOT
# THEOREM: every row group pins a fresh S_e block, so the sum-split
# never manufactures redundancy — stated in E12's docs, verified here
# by numerical rank. Two new species: sum-decomposition edges (vs
# agreement / partial-trace / jets) and a spectral (2D-DFT) oracle.
#
# The Agler split is CONSERVATIVE here by design: the cyclic band pattern
# of the torus is non-chordal, so requiring S = sum of clique-supported
# PSD pieces is a strict restriction of S >= 0 — that is what a
# compositional certificate is, and the conservatism is measured, not
# hidden: the monolithic comparator (P block-diag, one big LMI) has an
# exact spectral oracle at ANY K (group averaging -> per-wavenumber
# d-dim Hermitian LMIs), and the gap curve is reported (8.1% at K = 3
# shrinking to 3.4% at K = 4 at gamma = 0). E1's non-chordal argument,
# relocated to the constraint side.
#
# Checks:
#   [1]  setup identities: Gramian solves the local Lyapunov equation to
#        machine precision, is dense and SPD; global A is Hurwitz.
#   [2]  COKER = 0 (the private-slot theorem, verified): the full svec
#        coboundary (exactly the e13.jl layout) has FULL ROW RANK at
#        K = 3, 4 — the design converse of E12/quantum-torus's
#        rank-deficiency gates. Aspect reported (wide, 42/82 ~ 0.51).
#   [3]  FEASIBILITY of the Agler split at the flagship coupling
#        (conservatism can mean "no certificate exists"; here one does).
#   [4]  SYMMETRY ORACLE: v(K) = K^2 * v_cell to ~1e-10 for K = 3, 4
#        (K = 2 is degenerate — wraparound doubles couplings — and is
#        documented, not gated), and STATE UNIQUENESS: every P_v from the
#        full solve equals the cell solution to ~1e-8.
#   [5]  MONOLITHIC SPECTRAL ORACLE + CONSERVATISM CURVE: per-wavenumber
#        DFT LMIs equal the direct monolithic solve at K = 3 (~3e-10);
#        the per-cell gap (split vs monolithic) is reported for
#        K = 3..8 and shrinks toward its continuum limit.
#   [6]  DYNAMICAL CERTIFICATE: simulate xdot = Ax; V = x'Px decays
#        monotonically and V(t) <= V(0) exp(-t / lambda_max(P)), the rate
#        the certificate promises, over random initial conditions.
#   [7]  INSTABILITY LOCALIZATION: destabilize ONE subsystem; the global
#        system goes unstable, the certificate program goes infeasible,
#        and the minimal per-vertex relaxation (slack t_v on the -I
#        margin) CONCENTRATES on the broken vertex. Control: the healthy
#        system needs zero relaxation.
#   [8]  DEFORMATION: |v_gamma - v_0| -> 0 at measured rate gamma^1.00.
#   [9]  FORMULATION CERTIFICATE: the svec sheaf assembly (B, Q, c —
#        exactly what e13.jl builds) equals the matrix-native program,
#        in VALUE and STATE.
#   [10] CROSS-SOLVER: Clarabel vs SCS, on STATES (uniqueness makes this
#        meaningful).
#
# Usage: python3 torus_lyapunov_oracle.py [--quick]
# =============================================================================

import sys
import numpy as np
import cvxpy as cp
from itertools import product
from scipy.linalg import solve_lyapunov, sqrtm, expm

QUICK = "--quick" in sys.argv
D = 4                       # subsystem dimension
EPS = 0.15                  # coupling strength
GAMMA = 0.1                 # flagship Gramian-weighted regularization
SEED = 11
SP = D*(D+1)//2             # svec dim of P_v (10)
SS = 2*D*(2*D+1)//2         # svec dim of S_e (36)

TIGHT = dict(solver=cp.CLARABEL, tol_gap_abs=1e-11, tol_gap_rel=1e-11,
             tol_feas=1e-11, max_iter=2000)


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
# Instance
# -----------------------------------------------------------------------------

def make_instance(seed=SEED, d=D, eps=EPS):
    rng = np.random.default_rng(seed)
    G = rng.normal(size=(d, d))
    S0 = rng.normal(size=(d, d)); S0 = S0 - S0.T
    Aloc = S0 - (np.eye(d) + 0.3*(G @ G.T))            # Hurwitz, dense
    Cp = rng.normal(size=(d, d))/np.sqrt(d)
    Cm = rng.normal(size=(d, d))/np.sqrt(d)
    Dp = rng.normal(size=(d, d))/np.sqrt(d)
    Dm = rng.normal(size=(d, d))/np.sqrt(d)
    W = solve_lyapunov(Aloc, -np.eye(d))               # controllability Gramian
    Wh = np.real(sqrtm(W))
    return dict(d=d, eps=eps, Aloc=Aloc, Cp=Cp, Cm=Cm, Dp=Dp, Dm=Dm,
                W=W, Wh=Wh)


def global_A(inst, K, destabilize=None, shift=0.0):
    d, eps = inst["d"], inst["eps"]
    idx = {(x, y): (x*K + y) for x, y in product(range(K), range(K))}
    A = np.zeros((K*K*d, K*K*d))
    for x, y in product(range(K), range(K)):
        v = idx[(x, y)]
        Av = inst["Aloc"].copy()
        if destabilize == (x, y):
            Av += shift*np.eye(d)
        A[v*d:(v+1)*d, v*d:(v+1)*d] = Av
        for (C, w) in ((inst["Cp"], idx[((x+1) % K, y)]),
                       (inst["Cm"], idx[((x-1) % K, y)]),
                       (inst["Dp"], idx[(x, (y+1) % K)]),
                       (inst["Dm"], idx[(x, (y-1) % K)])):
            A[v*d:(v+1)*d, w*d:(w+1)*d] += eps*C
    return A, idx


# -----------------------------------------------------------------------------
# svec machinery and the linear maps of the assembly
# -----------------------------------------------------------------------------

def svec(M):
    Dm_ = M.shape[0]
    v = np.zeros(Dm_*(Dm_+1)//2); k = 0
    for c in range(Dm_):
        for r in range(c, Dm_):
            v[k] = M[c, c] if r == c else np.sqrt(2.0)*M[r, c]
            k += 1
    return v


def smat(v, Dm_):
    M = np.zeros((Dm_, Dm_)); k = 0
    for c in range(Dm_):
        for r in range(c, Dm_):
            if r == c:
                M[c, c] = v[k]
            else:
                M[r, c] = M[c, r] = v[k]/np.sqrt(2.0)
            k += 1
    return M


def lin_map(f, din, dout_shape):
    """Matrix of a linear map given by f: svec-coords -> array."""
    cols = []
    for k in range(din):
        e = np.zeros(din); e[k] = 1.0
        cols.append(f(e).ravel())
    return np.array(cols).T


def build_maps(inst):
    d = inst["d"]
    A, Wh = inst["Aloc"], inst["Wh"]
    LY = lin_map(lambda e: svec(A.T @ smat(e, d) + smat(e, d) @ A), SP, (SP,))
    GG = lin_map(lambda e: svec(Wh @ smat(e, d) @ Wh), SP, (SP,))
    E11 = lin_map(lambda e: svec(smat(e, 2*d)[:d, :d]), SS, (SP,))
    E22 = lin_map(lambda e: svec(smat(e, 2*d)[d:, d:]), SS, (SP,))
    EOD = lin_map(lambda e: smat(e, 2*d)[:d, d:], SS, (d, d))     # d^2 x SS
    # off-diagonal block maps: vec(eps*(C'_in P)) and vec(eps*(P C_out))
    odL = lambda C: lin_map(lambda e: inst["eps"]*(C.T @ smat(e, d)), SP, (d, d))
    odR = lambda C: lin_map(lambda e: inst["eps"]*(smat(e, d) @ C), SP, (d, d))
    return dict(LY=LY, GG=GG, E11=E11, E22=E22, EOD=EOD,
                H_v=odR(inst["Cp"]), H_w=odL(inst["Cm"]),
                V_v=odR(inst["Dp"]), V_w=odL(inst["Dm"]))


def assemble_B(inst, K):
    """Rows of the coboundary in svec coordinates — the exact e13.jl layout.
    Stalk order: P_1..P_{K^2} | Sh_1..Sh_{K^2} | Sv_1..Sv_{K^2}."""
    W = K*K
    mp = build_maps(inst)
    ncol = W*SP + 2*W*SS
    pc = lambda v: v*SP
    hc = lambda v: W*SP + v*SS
    vc = lambda v: W*SP + W*SS + v*SS
    idx = {(x, y): (x*K + y) for x, y in product(range(K), range(K))}
    rows, rhs = [], []
    for x, y in product(range(K), range(K)):
        v = idx[(x, y)]
        vr, vu = idx[((x+1) % K, y)], idx[(x, (y+1) % K)]
        vl, vd = idx[((x-1) % K, y)], idx[(x, (y-1) % K)]
        R = np.zeros((SP, ncol))                     # vertex equation
        R[:, pc(v):pc(v)+SP] += mp["LY"]
        R[:, hc(v):hc(v)+SS] += mp["E11"]
        R[:, vc(v):vc(v)+SS] += mp["E11"]
        R[:, hc(vl):hc(vl)+SS] += mp["E22"]
        R[:, vc(vd):vc(vd)+SS] += mp["E22"]
        rows.append(R); rhs.append(svec(-np.eye(inst["d"])))
        Rh = np.zeros((inst["d"]**2, ncol))          # horizontal edge equation
        Rh[:, pc(v):pc(v)+SP] += mp["H_v"]
        Rh[:, pc(vr):pc(vr)+SP] += mp["H_w"]
        Rh[:, hc(v):hc(v)+SS] += mp["EOD"]
        rows.append(Rh); rhs.append(np.zeros(inst["d"]**2))
        Rv = np.zeros((inst["d"]**2, ncol))          # vertical edge equation
        Rv[:, pc(v):pc(v)+SP] += mp["V_v"]
        Rv[:, pc(vu):pc(vu)+SP] += mp["V_w"]
        Rv[:, vc(v):vc(v)+SS] += mp["EOD"]
        rows.append(Rv); rhs.append(np.zeros(inst["d"]**2))
    return np.vstack(rows), np.concatenate(rhs)


# -----------------------------------------------------------------------------
# Solvers
# -----------------------------------------------------------------------------

def solve_full(inst, K, gamma=GAMMA, A=None, slack=False, solver_kw=None):
    d = inst["d"]
    if A is None:
        A, _ = global_A(inst, K)
    idx = {(x, y): (x*K + y) for x, y in product(range(K), range(K))}
    P = {v: cp.Variable((d, d), symmetric=True) for v in range(K*K)}
    Sh = {v: cp.Variable((2*d, 2*d), symmetric=True) for v in range(K*K)}
    Sv = {v: cp.Variable((2*d, 2*d), symmetric=True) for v in range(K*K)}
    t = cp.Variable(K*K, nonneg=True) if slack else None
    cons = ([p >> 0 for p in P.values()] + [s >> 0 for s in Sh.values()]
            + [s >> 0 for s in Sv.values()])

    def blk(v, w):
        return (A[w*d:(w+1)*d, v*d:(v+1)*d].T @ P[w]
                + P[v] @ A[v*d:(v+1)*d, w*d:(w+1)*d])

    for x, y in product(range(K), range(K)):
        v = idx[(x, y)]
        vr, vu = idx[((x+1) % K, y)], idx[(x, (y+1) % K)]
        vl, vd = idx[((x-1) % K, y)], idx[(x, (y-1) % K)]
        margin = (1 - t[v])*np.eye(d) if slack else np.eye(d)
        cons.append(blk(v, v) + margin + Sh[v][:d, :d] + Sv[v][:d, :d]
                    + Sh[vl][d:, d:] + Sv[vd][d:, d:] == 0)
        cons.append(blk(v, vr) + Sh[v][:d, d:] == 0)
        cons.append(blk(v, vu) + Sv[v][:d, d:] == 0)
    if slack:
        obj = cp.sum(t) + 1e-4*cp.sum([cp.trace(p) for p in P.values()])
    else:
        obj = cp.sum([cp.trace(p)
                      + 0.5*gamma*cp.sum_squares(inst["Wh"] @ p @ inst["Wh"])
                      for p in P.values()])
    pr = cp.Problem(cp.Minimize(obj), cons)
    conic_solve(pr, **(solver_kw or {}))
    out = [P[v].value for v in range(K*K)]
    return (pr.value, out, t.value) if slack else (pr.value, out)


def solve_cell(inst, gamma=GAMMA, solver_kw=None):
    d, eps = inst["d"], inst["eps"]
    Aloc = inst["Aloc"]
    P = cp.Variable((d, d), symmetric=True)
    Sh = cp.Variable((2*d, 2*d), symmetric=True)
    Sv = cp.Variable((2*d, 2*d), symmetric=True)
    cons = [P >> 0, Sh >> 0, Sv >> 0,
            Aloc.T @ P + P @ Aloc + np.eye(d)
            + Sh[:d, :d] + Sv[:d, :d] + Sh[d:, d:] + Sv[d:, d:] == 0,
            eps*(inst["Cm"].T @ P + P @ inst["Cp"]) + Sh[:d, d:] == 0,
            eps*(inst["Dm"].T @ P + P @ inst["Dp"]) + Sv[:d, d:] == 0]
    obj = (cp.trace(P)
           + 0.5*gamma*cp.sum_squares(inst["Wh"] @ P @ inst["Wh"]))
    pr = cp.Problem(cp.Minimize(obj), cons)
    conic_solve(pr, **(solver_kw or {}))
    return pr.value, P.value


def solve_monolithic_dft(inst, K, gamma=GAMMA):
    """Monolithic comparator (P block-diag, ONE big LMI) via group averaging:
    per-wavenumber d-dim Hermitian LMIs, exact at any K."""
    d, eps = inst["d"], inst["eps"]
    P = cp.Variable((d, d), symmetric=True)
    cons = [P >> 0]
    for kx, ky in product(range(K), range(K)):
        tx, ty = 2*np.pi*kx/K, 2*np.pi*ky/K
        Ath = (inst["Aloc"]
               + eps*(inst["Cp"]*np.exp(1j*tx) + inst["Cm"]*np.exp(-1j*tx)
                      + inst["Dp"]*np.exp(1j*ty) + inst["Dm"]*np.exp(-1j*ty)))
        M = Ath.conj().T @ P + P @ Ath + np.eye(d)
        cons.append(cp.bmat([[cp.real(M), -cp.imag(M)],
                             [cp.imag(M), cp.real(M)]]) << 0)
    obj = K*K*(cp.trace(P)
               + 0.5*gamma*cp.sum_squares(inst["Wh"] @ P @ inst["Wh"]))
    pr = cp.Problem(cp.Minimize(obj), cons)
    conic_solve(pr)
    return pr.value, P.value


def solve_monolithic_direct(inst, K, gamma=GAMMA):
    d = inst["d"]
    A, _ = global_A(inst, K)
    n = K*K*d
    P = {v: cp.Variable((d, d), symmetric=True) for v in range(K*K)}
    Pbig = cp.bmat([[P[v] if v == w else np.zeros((d, d))
                     for w in range(K*K)] for v in range(K*K)])
    cons = [p >> 0 for p in P.values()]
    cons.append(A.T @ Pbig + Pbig @ A + np.eye(n) << 0)
    obj = cp.sum([cp.trace(p)
                  + 0.5*gamma*cp.sum_squares(inst["Wh"] @ p @ inst["Wh"])
                  for p in P.values()])
    pr = cp.Problem(cp.Minimize(obj), cons)
    conic_solve(pr)
    return pr.value


def solve_sheaf_assembly(inst, K, gamma=GAMMA, solver_kw=None):
    """EXACTLY the e13.jl program: z over svec stalks, B z = g, per-stalk
    PSD, Q = Gramian Gram on P stalks, c = svec(I) on P stalks."""
    W = K*K
    B, g = assemble_B(inst, K)
    mp = build_maps(inst)
    z = cp.Variable(B.shape[1])
    cons = [B @ z == g]
    for v in range(W):
        Zp = cp.Variable((inst["d"], inst["d"]), symmetric=True)
        cons.append(Zp >> 0)
        k = 0
        for c in range(inst["d"]):
            cons.append(Zp[c, c] == z[v*SP + k]); k += 1
            for r in range(c+1, inst["d"]):
                cons.append(Zp[r, c] == z[v*SP + k]/np.sqrt(2.0)); k += 1
    for j in range(2*W):
        base = W*SP + j*SS
        Zs = cp.Variable((2*inst["d"], 2*inst["d"]), symmetric=True)
        cons.append(Zs >> 0)
        k = 0
        for c in range(2*inst["d"]):
            cons.append(Zs[c, c] == z[base + k]); k += 1
            for r in range(c+1, 2*inst["d"]):
                cons.append(Zs[r, c] == z[base + k]/np.sqrt(2.0)); k += 1
    ci = svec(np.eye(inst["d"]))
    Qg = gamma*(mp["GG"].T @ mp["GG"])
    obj = sum(ci @ z[v*SP:(v+1)*SP]
              + 0.5*cp.quad_form(z[v*SP:(v+1)*SP], cp.psd_wrap(Qg))
              for v in range(W))
    pr = cp.Problem(cp.Minimize(obj), cons)
    conic_solve(pr, **(solver_kw or {}))
    Ps = [smat(z.value[v*SP:(v+1)*SP], inst["d"]) for v in range(W)]
    return pr.value, Ps


# -----------------------------------------------------------------------------
# Checks
# -----------------------------------------------------------------------------

def check_setup(inst):
    d = inst["d"]
    res = np.abs(inst["Aloc"] @ inst["W"] + inst["W"] @ inst["Aloc"].T
                 + np.eye(d)).max()
    ev = np.linalg.eigvalsh(inst["W"]).min()
    A3, _ = global_A(inst, 3)
    stab = np.linalg.eigvals(A3).real.max()
    print(f"[1] setup: Gramian residual {res:.2e}; min eig(W) = {ev:.3f} "
          f"(dense, min |W_ij| = {np.abs(inst['W']).min():.3f}); "
          f"max Re eig(A, K=3) = {stab:.3f} (Hurwitz)")
    assert res < 1e-12 and ev > 0 and stab < 0


def check_coker(inst):
    print("[2] COKER = 0 (private-slot theorem, verified by rank):")
    for K in (3, 4):
        B, _ = assemble_B(inst, K)
        n1, n0 = B.shape
        r = np.linalg.matrix_rank(B, tol=1e-9)
        print(f"      K={K}: B is {n1} x {n0} (aspect {n1/n0:.3f}, wide); "
              f"rank = {r} = N1 -> coker = {n1 - r}")
        assert r == n1


def check_feasibility_symmetry(inst):
    v_cell, P_cell = solve_cell(inst)
    print(f"[3] Agler-split feasibility: cell value v_bar = {v_cell:.8f} "
          f"(optimal; conservatism did not preclude a certificate)")
    Ks = (3,) if QUICK else (3, 4)
    vals = {}
    for K in Ks:
        v, Ps = solve_full(inst, K)
        vals[K] = (v, Ps)
        dev = max(np.abs(p - P_cell).max() for p in Ps)
        print(f"[4]   K={K}: v = {v:.8f} vs K^2 v_bar = {K*K*v_cell:.8f} "
              f"(|diff| {abs(v - K*K*v_cell):.2e}); state uniqueness "
              f"max|P_v - P_bar| = {dev:.2e}")
        assert abs(v - K*K*v_cell) < 1e-5 and dev < 1e-5
    print("[4] symmetry oracle + state uniqueness PASS "
          "(K = 2 degenerate by doubled wraparound couplings; documented)")
    return vals, v_cell, P_cell


def check_monolithic_gap(inst):
    v_dft, _ = solve_monolithic_dft(inst, 3)
    v_dir = solve_monolithic_direct(inst, 3)
    print(f"[5] spectral oracle: DFT {v_dft:.8f} vs direct monolithic "
          f"{v_dir:.8f} (|diff| {abs(v_dft - v_dir):.2e})")
    assert abs(v_dft - v_dir) < 1e-5
    v_cell0, _ = solve_cell(inst, gamma=0.0)
    Ks = (3, 4, 6) if QUICK else (3, 4, 6, 8)
    print("[5] conservatism curve (per-cell, gamma = 0): "
          "split - monolithic =")
    for K in Ks:
        vm, _ = solve_monolithic_dft(inst, K, gamma=0.0)
        gap = v_cell0 - vm/(K*K)
        print(f"      K={K}: {gap:.5f} ({100*gap*(K*K)/vm:.2f}% of monolithic)")
        assert gap > -1e-7      # split is a restriction: never below


def check_dynamical(inst, vals):
    K = 3
    _, Ps = vals[K]
    d = inst["d"]
    A, _ = global_A(inst, K)
    Pg = np.zeros((K*K*d, K*K*d))
    for v in range(K*K):
        Pg[v*d:(v+1)*d, v*d:(v+1)*d] = Ps[v]
    lam = np.linalg.eigvalsh(Pg).max()
    rng = np.random.default_rng(0)
    worst = 0.0
    ts = np.linspace(0.0, 6.0, 25)
    for _ in range(5):
        x0 = rng.normal(size=K*K*d)
        V = [ (expm(A*t) @ x0) @ Pg @ (expm(A*t) @ x0) for t in ts]
        V = np.array(V)
        assert np.all(np.diff(V) < 1e-10), "V not monotone decreasing"
        worst = max(worst, (V/(V[0]*np.exp(-ts/lam))).max())
    print(f"[6] dynamical certificate: V = x'Px monotone decreasing; "
          f"V(t) <= V(0) exp(-t/lambda_max(P)) with max ratio "
          f"{worst:.3f} <= 1")
    assert worst <= 1.0 + 1e-9


def check_localization(inst):
    K = 3
    vstar = (1, 1)
    shift = 3.0
    Abad, idx = global_A(inst, K, destabilize=vstar, shift=shift)
    stab = np.linalg.eigvals(Abad).real.max()
    assert stab > 0, "destabilization failed"
    _, _, t = solve_full(inst, K, A=Abad, slack=True)
    _, _, t0 = solve_full(inst, K, slack=True)
    vs = idx[vstar]
    others = np.delete(t, vs)
    print(f"[7] instability localization: broken subsystem at {vstar} "
          f"(max Re eig = {stab:.2f}); minimal relaxation t_broken = "
          f"{t[vs]:.4f}, max elsewhere = {others.max():.2e}; healthy "
          f"control: max t = {t0.max():.2e}")
    assert np.argmax(t) == vs and t[vs] > 0.5
    assert others.max() < 0.05*t[vs] and t0.max() < 1e-5


def check_deformation(inst, v_cell):
    v0, _ = solve_cell(inst, gamma=0.0)
    gs = [0.4, 0.2, 0.1, 0.05]
    errs = []
    for g in gs:
        vg, _ = solve_cell(inst, gamma=g)
        errs.append(abs(vg - v0))
    rate = np.polyfit(np.log(gs), np.log(errs), 1)[0]
    print(f"[8] deformation gamma -> 0: |v_gamma - v_0| {errs[0]:.2e} -> "
          f"{errs[-1]:.2e}; rate gamma^{rate:.2f} (predicted 1)")
    assert 0.8 < rate < 1.2


def check_formulation(inst, vals):
    K = 3
    v_mat, Ps_mat = vals[K]
    v_sh, Ps_sh = solve_sheaf_assembly(inst, K)
    dv = abs(v_sh - v_mat)
    dx = max(np.abs(a - b).max() for a, b in zip(Ps_sh, Ps_mat))
    print(f"[9] formulation certificate: svec assembly vs matrix-native "
          f"|dv| = {dv:.2e}, state |dx| = {dx:.2e}")
    assert dv < 1e-6 and dx < 1e-5


def check_cross_solver(inst):
    v_c, P_c = solve_cell(inst)
    P = cp.Variable((inst["d"], inst["d"]), symmetric=True)
    Sh = cp.Variable((2*inst["d"], 2*inst["d"]), symmetric=True)
    Sv = cp.Variable((2*inst["d"], 2*inst["d"]), symmetric=True)
    d, eps = inst["d"], inst["eps"]
    cons = [P >> 0, Sh >> 0, Sv >> 0,
            inst["Aloc"].T @ P + P @ inst["Aloc"] + np.eye(d)
            + Sh[:d, :d] + Sv[:d, :d] + Sh[d:, d:] + Sv[d:, d:] == 0,
            eps*(inst["Cm"].T @ P + P @ inst["Cp"]) + Sh[:d, d:] == 0,
            eps*(inst["Dm"].T @ P + P @ inst["Dp"]) + Sv[:d, d:] == 0]
    obj = (cp.trace(P)
           + 0.5*GAMMA*cp.sum_squares(inst["Wh"] @ P @ inst["Wh"]))
    pr = cp.Problem(cp.Minimize(obj), cons)
    pr.solve(solver=cp.SCS, eps=1e-9, max_iters=200_000)
    dx = np.abs(P.value - P_c).max()
    print(f"[10] cross-solver ON STATES (unique optimum): Clarabel vs SCS "
          f"|dv| = {abs(pr.value - v_c):.2e}, |dP| = {dx:.2e}")
    assert abs(pr.value - v_c) < 1e-6 and dx < 1e-5


# -----------------------------------------------------------------------------

if __name__ == "__main__":
    inst = make_instance()
    print(f"instance: K x K torus of d = {D} subsystems, eps = {EPS}, "
          f"gamma = {GAMMA} (Gramian-weighted dense Q); stalks svec "
          f"{SP} (P_v) + 2 x {SS} (S_e)")
    check_setup(inst)
    check_coker(inst)
    vals, v_cell, P_cell = check_feasibility_symmetry(inst)
    check_monolithic_gap(inst)
    check_dynamical(inst, vals)
    check_localization(inst)
    check_deformation(inst, v_cell)
    check_formulation(inst, vals)
    check_cross_solver(inst)
    print("ALL CHECKS PASSED")
