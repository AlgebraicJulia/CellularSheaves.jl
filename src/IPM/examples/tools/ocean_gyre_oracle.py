"""
ocean_gyre_oracle.py -- self-contained numpy/scipy validation for E14:
gyre-aware robust flow recovery on a punctured K x K grid.

Complex: planar grid, open outer boundary; ISLANDS = removed faces (holes).
  d0: m x n     vertex potentials -> edge flows        (coboundary)
  d1: f x m     edge flows -> face circulations        (islands' rows absent)
  eta: m x nI   orthonormal harmonic basis, ker d1 cap ker d0'; nI = #islands.

Model. Each edge carries a T-sample flow time series g_e. Components live in a
smooth temporal basis Phi (T x p, orthonormal; constant + first Fourier pair):

    g_e = (d0 x)_e Phi' + (d1' z)_e Phi' + (eta h)_e Phi' + u_e + v_e

  x (n x p)  potential coefficients        -- smooth by construction
  z (f x p)  vorticity coefficients        -- smooth by construction
  h (nI x p) GYRE coefficients (estimand)  -- circulation around each island
  u (m x T)  inlier residual, cost 1/2 u_e' W_e u_e with W_e DENSE SPD
             (inverse squared-exponential GP kernel, per-edge length scale)
  v (m x T)  outlier bursts, group-lasso  lam * sum_e ||v_e||_2  (SOC epigraph)

The smooth-basis restriction is the identifiability mechanism: white bursts
are OUTSIDE the model span at any coefficient magnitude, so burst/clean
separation is structural, not tuned. (Penalty-based discrimination provably
fails here: shrinkage misfit on true signal grows faster than burst residual
-- measured during design, see conversation log.)

Gates (this file, fully executed; numbers quoted in e14/def.jl):
  G1 structure     d1 d0 = 0; dim H1 = nI; Hodge orthogonality;
                   SPAN: rank[d0 | d1' | eta] = m -- the combined design
                   covers C^1 although coker d0 alone has dim m - rank d0.
  G2 LS recovery   gyre coefficients at the noise floor; component orthogonality.
  G3 island-blind  drop eta -> residual grows with gyre strength and
                   CONCENTRATES on high-harmonic-amplitude edges (dose-response).
  G4 localization  screening margin ~11x; group-lasso support recovery EXACT;
                   debiased refit == mask-oracle to machine precision.
  G5 deformation   lam sweep: active set 6 -> 1, distance to V=0 fit shrinks.
"""
import numpy as np
from scipy.linalg import cho_factor, cho_solve

np.random.seed(7)

# ---------------------------------------------------------------- complex
def punctured_grid(K, islands):
    nv = (K + 1) ** 2
    vid = lambda i, j: i * (K + 1) + j
    E, eid = [], {}
    for i in range(K):
        for j in range(K + 1):
            eid[("h", i, j)] = len(E); E.append((vid(i, j), vid(i + 1, j)))
    for i in range(K + 1):
        for j in range(K):
            eid[("v", i, j)] = len(E); E.append((vid(i, j), vid(i, j + 1)))
    m = len(E)
    d0 = np.zeros((m, nv))
    for e, (a, b) in enumerate(E):
        d0[e, b] += 1.0; d0[e, a] -= 1.0
    faces = [(i, j) for i in range(K) for j in range(K) if (i, j) not in islands]
    d1 = np.zeros((len(faces), m))
    for pfc, (i, j) in enumerate(faces):
        d1[pfc, eid[("h", i, j)]] += 1.0
        d1[pfc, eid[("v", i + 1, j)]] += 1.0
        d1[pfc, eid[("h", i, j + 1)]] -= 1.0
        d1[pfc, eid[("v", i, j)]] -= 1.0
    return d0, d1

K = 8
islands = [(K // 4, K // 4), (3 * K // 4, 3 * K // 4)]
d0, d1 = punctured_grid(K, islands)
m, n = d0.shape
f = d1.shape[0]
_, sv, Vt = np.linalg.svd(np.vstack([d0.T, d1]))
eta = Vt[np.sum(sv > 1e-10):].T
nI = eta.shape[1]
D = np.hstack([d0, d1.T, eta])
q = D.shape[1]

T, p = 6, 3
tt = np.arange(T)
Phi, _ = np.linalg.qr(np.vstack([np.ones(T), np.sin(2 * np.pi * tt / T),
                                 np.cos(2 * np.pi * tt / T)]).T)

def dense_W(scale):
    Kmat = np.exp(-0.5 * ((tt[:, None] - tt[None, :]) / scale) ** 2)
    return np.linalg.inv(Kmat + 0.05 * np.eye(T))

Ws = [dense_W(1.0 + 1.5 * np.random.rand()) for _ in range(m)]
RIDGE = 1e-6
nz = [np.nonzero(D[e])[0] for e in range(m)]

def make_chol(qq, mask=None):
    H = np.zeros((qq * p, qq * p))
    for e in range(m):
        if mask and e in mask:
            continue
        idx = [i for i in nz[e] if i < qq]
        de = D[e, idx]
        ii = np.concatenate([np.arange(i * p, (i + 1) * p) for i in idx])
        H[np.ix_(ii, ii)] += np.kron(np.outer(de, de), Phi.T @ Ws[e] @ Phi)
    H += RIDGE * np.eye(qq * p)
    return cho_factor(H)

def fit(g, chol, qq, mask=None):
    GW = np.stack([Phi.T @ (Ws[e] @ g[e]) if not (mask and e in mask)
                   else np.zeros(p) for e in range(m)])
    return cho_solve(chol, (D[:, :qq].T @ GW).reshape(qq * p)).reshape(qq, p)

pred = lambda C, qq: (D[:, :qq] @ C) @ Phi.T
chol_full = make_chol(q)
chol_blind = make_chol(n + f)

# ---------------------------------------------------------------- synthesis
xc = 1.0 * np.random.randn(n, p)
zc = np.zeros((f, p))
vort = np.random.choice(f, 5, replace=False)
zc[vort] = 2.0 * np.random.randn(5, p)          # smooth-in-time vortices
hc = 3.0 * np.random.randn(nI, p)               # gyre circulations (estimand)
sig = 0.05
g_noisy = (D @ np.vstack([xc, zc, hc])) @ Phi.T + sig * np.random.randn(m, T)
burst = set(int(i) for i in np.random.choice(m, 6, replace=False))
g_data = g_noisy.copy()
for e in burst:
    g_data[e] += 4.0 * np.random.randn(T)       # white bursts: outside span(Phi)

# ---------------------------------------------------------------- G1
print(f"complex: K={K}, n={n} vertices, m={m} edges, f={f} faces, "
      f"T={T}, p={p}, islands={nI}")
print(f"[G1] ||d1 d0||_max = {np.abs(d1 @ d0).max():.1e};  dim H1 = {nI}")
print(f"[G1] Hodge orthogonality: |d0'd1'| = {np.abs(d0.T @ d1.T).max():.1e}, "
      f"|d0'eta| = {np.abs(d0.T @ eta).max():.1e}, |d1 eta| = {np.abs(d1 @ eta).max():.1e}")
print(f"[G1] SPAN: rank[d0 | d1' | eta] = {np.linalg.matrix_rank(D)} of m = {m} "
      f"(coker d0 alone: {m - np.linalg.matrix_rank(d0)})")

# ---------------------------------------------------------------- G2
C = fit(g_noisy, chol_full, q)
r = g_noisy - pred(C, q)
print(f"[G2] gyre coeff err = {np.abs(C[n + f:] - hc).max():.3e} "
      f"(noise floor ~{2.5 * sig:.3f}); resid rms = {np.sqrt(np.mean(r ** 2)):.3e} "
      f"(sig = {sig}); <d0 x, eta h> = "
      f"{abs(np.sum((d0 @ C[:n] @ Phi.T) * (eta @ C[n + f:] @ Phi.T))):.1e}")

# ---------------------------------------------------------------- G3
amp = np.linalg.norm(eta, axis=1)
top = np.argsort(-amp)[:m // 10]
for gy in (0.0, 1.0, 3.0):
    gg = (np.hstack([d0, d1.T]) @ np.vstack([xc, zc])) @ Phi.T \
         + gy * (eta @ hc @ Phi.T) / 3.0 + sig * np.random.randn(m, T)
    Cb = fit(gg, chol_blind, n + f)
    rb = gg - pred(Cb, n + f)
    conc = np.sum(rb[top] ** 2) / np.sum(rb ** 2)
    print(f"[G3] gyre {gy:.0f}: blind resid rms = {np.sqrt(np.mean(rb ** 2)):.3e}; "
          f"top-10% harmonic-amplitude edges carry {100 * conc:5.1f}% of residual "
          f"(uniform = 10%)")

# ---------------------------------------------------------------- G4
C0 = fit(g_data, chol_full, q)
R = g_data - pred(C0, q)
s = np.array([np.linalg.norm(Ws[e] @ R[e]) for e in range(m)])
clean = np.setdiff1d(np.arange(m), sorted(burst))
bmin, cmax = s[sorted(burst)].min(), s[clean].max()
print(f"[G4a] screening: bursts min = {bmin:.2f}, clean max = {cmax:.2f}, "
      f"margin = {bmin / cmax:.1f}x")
lam = float(np.sqrt(bmin * cmax))

def robust(g, lam, iters=400, tol=1e-9):
    V = np.zeros((m, T))
    L = max(np.linalg.eigvalsh(W).max() for W in Ws)
    st = 1.0 / L
    for it in range(iters):
        Cc = fit(g - V, chol_full, q)
        Rr = g - pred(Cc, q)
        Vp = V.copy()
        for e in range(m):
            w = V[e] - st * (Ws[e] @ (V[e] - Rr[e]))
            nw = np.linalg.norm(w)
            V[e] = np.zeros(T) if nw <= st * lam else (1 - st * lam / nw) * w
        if it > 3 and np.abs(V - Vp).max() < tol:
            break
    return Cc, V, it + 1

Cr, V, it = robust(g_data, lam)
found = set(int(i) for i in np.where(np.linalg.norm(V, axis=1) > 1e-2)[0])
print(f"[G4b] group-lasso (lam = {lam:.2f}, {it} its): found = {sorted(found)}, "
      f"true = {sorted(burst)}, exact = {found == burst}")
chol_db = make_chol(q, mask=found)
Cdb = fit(g_data, chol_db, q, mask=found)     # debiased robust estimator
Cor = fit(g_noisy, chol_db, q, mask=found)    # mask-oracle: clean data, same mask
print(f"[G4c] gyre err: naive LS = {np.abs(C0[n + f:] - hc).max():.3e}, "
      f"lasso = {np.abs(Cr[n + f:] - hc).max():.3e}, "
      f"debiased = {np.abs(Cdb[n + f:] - hc).max():.3e}, "
      f"clean-data LS = {np.abs(C[n + f:] - hc).max():.3e}")
print(f"[G4d] debiased vs mask-oracle: max|dh| = "
      f"{np.abs(Cdb[n + f:] - Cor[n + f:]).max():.2e} "
      f"(gap to clean-data LS is the measured identifiability price of the mask)")

# ---------------------------------------------------------------- G5
for mult in (1.0, 2.0, 4.0, 8.0, 16.0):
    Cd, Vd, _ = robust(g_data, mult * lam, iters=200)
    nact = int(np.sum(np.linalg.norm(Vd, axis=1) > 1e-2))
    print(f"[G5] lam = {mult * lam:7.2f}: active = {nact:3d}, "
          f"max|C - C_V0| = {np.abs(Cd - C0).max():.3e}")
print("done.")
