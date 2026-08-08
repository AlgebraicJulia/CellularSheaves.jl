# Cross-validate e14b's build_tracking convention: solve the RAW solver-form
# program min 0.5 p'Qp - c'p (Q_state = gamma*W, c_state = +gamma*W*xbar_t,
# c_sigma = -dt) over the E14 standard-form variables and compare against
# tracking_oracle's natural-form build.
import numpy as np, cvxpy as cp, warnings
warnings.filterwarnings("ignore")
from tracking_oracle import W_BASE, GAMMA_BASE, corridor, solve_track
from soft_landing_oracle import BASE, TF_BASE, zoh
from xcheck import defjl_assembly  # reuse the validated E14 B, g layout pieces

# rebuild the standard-form matrices exactly as xcheck.defjl_assembly does,
# but return them instead of solving (duplicate the layout inline):
def standard_form(inst, tf, N):
    dt, A, Bd = zoh(tf, N)
    ns, nt, nb = 6*(N+1), 4*N, 2*N
    ngs = 3*(N-1)
    ncol = ns + nt + nb + ngs
    def cs(t): return slice(6*(t-1), 6*t)
    def ct(k): return slice(ns+4*(k-1), ns+4*k)
    def cb(k): return slice(ns+nt+2*(k-1), ns+nt+2*k)
    def cg(t): return slice(ns+nt+nb+3*(t-2), ns+nt+nb+3*(t-1))
    rows, g = [], []
    for t in range(1, N+1):
        R = np.zeros((6, ncol)); R[:, cs(t)] = -A; R[:, cs(t+1)] = np.eye(6)
        R[:, ct(t)][:, 1:4] = -Bd
        rows.append(R); g.append(Bd @ inst.gvec)
    R = np.zeros((6, ncol)); R[:, cs(1)] = np.eye(6)
    rows.append(R); g.append(np.hstack([inst.r0, inst.v0]))
    R = np.zeros((6, ncol)); R[:, cs(N+1)] = np.eye(6)
    rows.append(R); g.append(np.zeros(6))
    for k in range(1, N+1):
        R = np.zeros((2, ncol)); R[0, ct(k).start] = 1; R[1, ct(k).start] = 1
        R[0, cb(k).start] = -1; R[1, cb(k).start+1] = 1
        rows.append(R); g.append(np.array([inst.rho1, inst.rho2]))
    for t in range(2, N+1):
        R = np.zeros((3, ncol)); R[:, cg(t)] = np.eye(3)
        R[0, cs(t).start+2] = -1/inst.tan_gs
        R[1, cs(t).start] = -1; R[2, cs(t).start+1] = -1
        rows.append(R); g.append(np.zeros(3))
    return np.vstack(rows), np.hstack(g), ncol, cs, ct, cb, cg, dt

N = 60
xbar = corridor(BASE, TF_BASE, N)
Bmat, gv, ncol, cs, ct, cb, cg, dt = standard_form(BASE, TF_BASE, N)

# solver-form Q, c per e14b/def.jl build_tracking
Q = np.zeros((ncol, ncol)); c = np.zeros(ncol)
for t in range(1, N+2):
    Q[cs(t), cs(t)] = GAMMA_BASE * W_BASE
    c[cs(t)] = GAMMA_BASE * (W_BASE @ xbar[t-1])
for k in range(1, N+1):
    c[ct(k).start] = -dt

x = cp.Variable(ncol)
cons = [Bmat @ x == gv]
for k in range(1, N+1):
    cons += [cp.norm(x[ct(k)][1:4]) <= x[ct(k)][0], x[cb(k)] >= 0]
for t in range(2, N+1):
    cons.append(cp.norm(x[cg(t)][1:3]) <= x[cg(t)][0])
prob = cp.Problem(cp.Minimize(0.5*cp.quad_form(x, cp.psd_wrap(Q)) - c @ x), cons)
prob.solve(solver=cp.CLARABEL)
assert prob.status in ("optimal", "optimal_inaccurate"), prob.status
Xjl = np.vstack([x.value[cs(t)] for t in range(1, N+2)])

# oracle natural form
po, Xo, Uo, so, lo, fuel_o, tr_o = solve_track(BASE, TF_BASE, N, GAMMA_BASE, xbar)

const = 0.5 * GAMMA_BASE * sum(xbar[t] @ W_BASE @ xbar[t] for t in range(N+1))
print("solver-form value + const = %.6f vs oracle objective %.6f  |D| = %.2e"
      % (prob.value + const, po.value, abs(prob.value + const - po.value)))
print("||DX||_inf = %.2e m" % np.max(np.abs(Xjl - Xo.value)))
fuel_jl = dt * sum(x.value[ct(k).start] for k in range(1, N+1))
print("fuel: %.4f vs %.4f" % (fuel_jl, fuel_o))
