# # Spring Oscillator: TrajectorySheaf vs Euler Colocation
#
# This example demonstrates the [`TrajectorySheaf`](@ref) API on a classic
# two-point boundary-value problem (BVP): the simple harmonic oscillator.
#
# We compare two approaches to solving the BVP:
# 1. **Sheaf colocation** — build a `TrajectorySheaf` with the *exact*
#    discrete-time dynamics (matrix exponential) and call
#    [`colocation_trajectory`](@ref), which uses harmonic extension.
# 2. **Euler colocation** — discretize the same ODE with a forward-Euler
#    step, set up the analogous block-bidiagonal linear system by hand, and
#    solve with Julia's built-in `\` operator.
#
# Comparing both solutions highlights two things:
# * the sheaf formulation provides an elegant, modular interface to what is
#   otherwise a bespoke linear-algebra construction;
# * the exact matrix-exponential dynamics incur no discretization error,
#   whereas the Euler approach accumulates O(h) error per step.

using CellularSheaves
using LinearAlgebra
using BlockArrays

# ## Spring oscillator in phase space
#
# The equation of motion is
# ```math
#   x''(t) = -\omega^2 x(t),
# ```
# with angular frequency ``\omega``.  Introducing the phase-space state
# ``q(t) = [x(t),\; \dot{x}(t)]^\top`` gives the first-order system
# ```math
#   \dot{q} = J q, \quad J = \begin{pmatrix} 0 & 1 \\ -\omega^2 & 0 \end{pmatrix}.
# ```

ω = 2.0                          # angular frequency
J = [0.0  1.0; -ω^2  0.0]       # continuous-time dynamics matrix

# We discretize time into `N` equal steps per period.  The period is
# ``T = 2\pi/\omega``.

N        = 20
T_period = 2π / ω
h        = T_period / N          # step size

# ### Exact discrete-time transition matrix
#
# The exact solution of ``\dot{q} = J q`` over one step of length ``h`` is
# ``q(t+h) = e^{hJ} q(t)``.  For the harmonic oscillator this has the
# closed form
# ```math
#   A = \begin{pmatrix}
#         \cos(\omega h) & \sin(\omega h)/\omega \\
#        -\omega\sin(\omega h) & \cos(\omega h)
#       \end{pmatrix}.
# ```

A_exact = [cos(ω*h)       sin(ω*h)/ω;
           -ω*sin(ω*h)    cos(ω*h)  ]

# Sanity-check against Julia's matrix exponential:
@assert norm(A_exact - exp(h * J)) < 1e-12

# ### Forward-Euler approximation
#
# The simplest first-order discretization replaces ``e^{hJ}`` by the
# first-order Taylor approximation ``I + h J``:

A_euler = I + h * J

# ## Boundary conditions
#
# We solve the BVP over half a period (``k = N/2`` steps), starting at
# maximum displacement and ending at minimum displacement.

k  = N ÷ 2                      # 10 steps = half period
x0 = [1.0, 0.0]                 # initial: max displacement, zero velocity
xk = A_exact^k * x0             # exact final state after half period ≈ [-1, 0]

println("x0 = ", round.(x0; digits=6))
println("xk = ", round.(xk; digits=6))

# ## Approach 1: Sheaf-based colocation (exact dynamics)
#
# A [`TrajectorySheaf`](@ref) encodes the dynamics ``q_{t+1} = A\,q_t`` as a
# cellular sheaf on a path graph with `k+1` vertices.  Its global sections are
# exactly the valid `k`-step trajectories.
#
# 1. We start from a base sheaf `F` with a single ``d``-dimensional vertex
#    stalk (``d = 2``).  The total 0-cochain dimension of `F` is `d`, which
#    becomes the state dimension of the trajectory sheaf.
# 2. We construct the `TrajectorySheaf` from `F`, `A_exact`, and `k`.
# 3. [`colocation_trajectory`](@ref) calls `harmonic_extension` to solve
#    ``L_{II}\,x_I = -L_{IB}\,x_B`` where ``B = \{1,\, k+1\}`` are the
#    boundary vertices carrying `x0` and `xk`.

d = 2                                        # state dimension
F = EuclideanSheaf{Float64}(fill(d, 1))      # base sheaf: one d-dim stalk

tsheaf = TrajectorySheaf(F, A_exact, k)

println("Inner sheaf: ", k+1, " vertices, ", k, " edges")

traj_sheaf = colocation_trajectory(tsheaf, x0, xk)

# Extract positions and velocities at each time step:
xs_sheaf = [Array(traj_sheaf[Block(t)])[1] for t in 1:k+1]
vs_sheaf = [Array(traj_sheaf[Block(t)])[2] for t in 1:k+1]

# Dynamics residual: ``\max_t \|q_{t+1} - A_{\text{exact}}\,q_t\|``
dyn_err_sheaf = maximum(
    norm(Array(traj_sheaf[Block(t+1)]) - A_exact * Array(traj_sheaf[Block(t)]))
    for t in 1:k)

println("Sheaf (exact A) — max dynamics residual: ", dyn_err_sheaf)

# ## Approach 2: Euler colocation (direct linear system)
#
# We now set up the analogous BVP directly using the forward-Euler matrix
# ``A_{\text{euler}} = I + h J``.
#
# Unknowns: interior states ``q_1, \ldots, q_{k-1}`` (a ``d(k-1)``-vector).
# The dynamics constraints ``q_{t+1} = A_{\text{euler}}\,q_t`` for
# ``t = 0, \ldots, k-1`` yield ``k`` block equations; substituting the
# boundary values ``q_0 = x_0`` and ``q_k = x_k`` produces a
# ``k \times (k-1)`` block over-determined system that we solve in the
# least-squares sense (equivalently, the harmonic extension of the same
# path-graph sheaf built with ``A_{\text{euler}}``):
#
# ```math
#   \underbrace{%
#     \begin{pmatrix}
#       I     &       &        \\
#       -A_e  & I     &        \\
#             & \ddots & \ddots \\
#             &       & -A_e   & I     \\
#             &       &        & -A_e
#     \end{pmatrix}}_{M \in \mathbb{R}^{kd \times (k-1)d}}
#   \begin{pmatrix} q_1 \\ \vdots \\ q_{k-1} \end{pmatrix}
#   =
#   \underbrace{%
#     \begin{pmatrix}
#       A_e x_0 \\ 0 \\ \vdots \\ 0 \\ x_k
#     \end{pmatrix}}_{b \in \mathbb{R}^{kd}}
# ```

Ae = Matrix{Float64}(A_euler)   # ensure concrete Matrix type
n_int = k - 1                   # number of interior time steps
M_sys = zeros(k * d, n_int * d)
b_sys = zeros(k * d)

for t in 0:k-1
    rows = t*d+1 : (t+1)*d

    # Contribution from q_{t+1}
    if t + 1 <= k - 1
        cols_next = t*d+1 : (t+1)*d
        M_sys[rows, cols_next] += I(d)
    else                          # q_k is a boundary value
        b_sys[rows] += xk
    end

    # Contribution from -A_euler * q_t
    if t >= 1
        cols_curr = (t-1)*d+1 : t*d
        M_sys[rows, cols_curr] -= Ae
    else                          # q_0 is a boundary value
        b_sys[rows] += Ae * x0
    end
end

# Solve the over-determined system in the least-squares sense:
q_int_euler = M_sys \ b_sys

# Reconstruct the full trajectory (including boundaries):
xs_euler = zeros(k + 1)
vs_euler = zeros(k + 1)
xs_euler[1], vs_euler[1] = x0
xs_euler[k+1], vs_euler[k+1] = xk
for t in 1:n_int
    idx = (t-1)*d+1 : t*d
    xs_euler[t+1] = q_int_euler[idx[1]]
    vs_euler[t+1] = q_int_euler[idx[2]]
end

# Dynamics residual with respect to the *exact* matrix:
dyn_err_euler = maximum(
    norm([xs_euler[t+1], vs_euler[t+1]] - A_exact * [xs_euler[t], vs_euler[t]])
    for t in 1:k)

println("Euler direct  — max dynamics residual: ", dyn_err_euler)

# ## Comparison
#
# The sheaf approach with the exact matrix exponential achieves near-machine-
# epsilon dynamics residual, while the Euler approach accumulates ``O(h)``
# discretization error per step.

println()
println("Step size h = ", round(h; digits=6), " (N=$N steps/period)")
println("Sheaf (exact A)  max residual: ", round(dyn_err_sheaf; sigdigits=3))
println("Euler direct     max residual: ", round(dyn_err_euler; sigdigits=3))
println("Ratio (Euler/Sheaf): ", round(dyn_err_euler / max(dyn_err_sheaf, eps()); sigdigits=3))

# As expected, the Euler residual is ``O(h \omega) \approx O(0.63)`` per step,
# consistent with the first-order truncation error of forward Euler.
#
# Both boundary conditions are satisfied exactly by construction, because
# harmonic extension pins the boundary vertices regardless of the dynamics
# matrix used.

println("Sheaf x(0) = ", round.(Array(traj_sheaf[Block(1)]);   digits=10))
println("Sheaf x(T/2) = ", round.(Array(traj_sheaf[Block(k+1)]); digits=10))
println("Euler x(0) = ", round.([xs_euler[1], vs_euler[1]];   digits=10))
println("Euler x(T/2) = ", round.([xs_euler[k+1], vs_euler[k+1]]; digits=10))
