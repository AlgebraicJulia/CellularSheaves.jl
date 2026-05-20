# # Controlled Trajectory Examples: 1 — Double Integrator
#
# This is the **first example** in a four-part progression of controlled-trajectory
# benchmarks built with `CellularSheaves.jl`:
#
# 1. **Double integrator** (this page) — the canonical first LQR system.
# 2. **Vehicle platoon** — same control story in a multi-agent / path-graph setting.
# 3. **Planar quadrotor** — a more realistic multi-input robotic system.
# 4. **Mass-spring-damper chain** — graph-coupled mechanical system bridging
#    toward later consensus and network examples.
#
# Each example follows the same conceptual pipeline:
#
# ```math
# (A_c, B_c)
# \;\longrightarrow\;
# (A_d, B_d)
# \;\longrightarrow\;
# \mathcal{T}(x_1, x_{k+1})
# \;\longrightarrow\;
# \min_{z \in \mathcal{T}(x_1, x_{k+1})} J(z).
# ```
#
# ## Physical system
#
# A unit point mass moves along a line.  The state is position and velocity;
# the control input is acceleration:
#
# ```math
# x = \begin{bmatrix} p \\ \dot{p} \end{bmatrix}, \qquad u = a,
# \qquad
# A_c = \begin{bmatrix} 0 & 1 \\ 0 & 0 \end{bmatrix},
# \qquad
# B_c = \begin{bmatrix} 0 \\ 1 \end{bmatrix}.
# ```
#
# Despite its simplicity this system exhibits every feature of the
# controlled-trajectory workflow and serves as the natural reference point for
# the rest of the progression.
#
# **References:**
# - Anderson, B. D. O. & Moore, J. B. (1990). *Optimal Control: Linear Quadratic
#   Methods*. Prentice-Hall.
# - Lewis, F. L., Vrabie, D. L., & Syrmos, V. L. (2012). *Optimal Control*, 3rd ed.
#   John Wiley & Sons. ISBN 978-0-470-63349-6.

using CellularSheaves
using LinearAlgebra
using BlockArrays
using Plots
using SparseArrays
using CellularSheaves.TrajectorySheaves
# Import the API extension that provides `nullspace_trajectory_family`.
using CellularSheaves.TrajectorySheaves: nullspace_trajectory_family

# ## Step 1: Define the continuous-time model
#
# The double integrator in standard state-space form.

Ac = [0.0  1.0;
      0.0  0.0]
Bc = reshape([0.0, 1.0], 2, 1)

# ## Step 2: Discretize with zero-order hold
#
# A [`ControlledTrajectorySheaf`](@ref) wraps the ZOH discretization
# ``x_{t+1} = A_d x_t + B_d u_t`` as a path-graph sheaf.  Its global
# sections are exactly the feasible sampled state-control trajectories.

h = 0.25   # sample period (seconds)
k = 8      # number of time steps

F  = EuclideanSheaf{Float64}(fill(2, 1))   # base sheaf: one vertex, 2D stalk
ts = ControlledTrajectorySheaf(F, Ac, Bc, h, k)

println("State dimension  n = ", ts.state_dim)
println("Control dimension m = ", ts.control_dim)
println("Steps             k = ", ts.k)

# ## Step 3: Fix endpoint states
#
# We start at rest and aim to reach a target position with zero velocity.

x1  = [0.0, 0.0]   # initial state:  position 0, velocity 0
xk1 = [1.0, 0.0]   # terminal state: position 1, velocity 0

# ## Step 4: Build the feasible trajectory affine space
#
# [`feasible_control_trajectory_basis`](@ref) returns:
# - `z_p` — a particular feasible trajectory (the harmonic extension of the
#   boundary data);
# - `N`   — columns spanning all endpoint-preserving perturbations.
#
# Any feasible trajectory satisfies ``z = z_p + N \alpha`` for some
# ``\alpha \in \mathbb{R}^r``.

z_p_raw, N = feasible_control_trajectory_basis(ts, x1, xk1)

n = ts.state_dim
m = ts.control_dim
r = size(N, 2)
println("Trajectory dimension p = ", (k+1)*n + k*m)
println("Free parameters      r = ", r)

# ## Step 5: Assemble the LQR objective
#
# The conditioning diagnostic above helps distinguish two different effects:
# the full sheaf Laplacian is singular for structural reasons, whereas the
# restricted interior block `L_II` is the matrix actually inverted in the
# harmonic extension step.
#
# As in the scalar integrator example, a nonzero running state penalty `Q`
# pulls the optimizer toward smaller intermediate states and can visibly squash
# the trajectory. To emphasize the minimum-control-effort solution between fixed
# endpoints, we set `Q = 0` here.
#
# We therefore minimize the quadratic running cost
#
# ```math
# J(z) = \tfrac{1}{2} \sum_{t=1}^{k}
#   \bigl( x_t^\top Q\, x_t + u_t^\top R_u\, u_t \bigr)
#   + \tfrac{1}{2}\, x_{k+1}^\top Q_f\, x_{k+1}.
# ```
#
# The assembled cost Hessian ``H`` is block-diagonal: ``k`` copies of ``Q``
# on the running state blocks, ``Q_f`` on the terminal block, and ``k`` copies
# of ``R_u`` on the control blocks. Because `x_{k+1}` is fixed as a hard
# boundary condition, the `Q_f` term is constant over the feasible set and does
# not affect the optimizer.

Q  = zeros(n, n)
Ru = Matrix{Float64}(I, m, m)
Qf = 0.0 * Matrix{Float64}(I, n, n)   # constant over the feasible set here

H, f, _ = lqr_objective(ts, Q, Ru; Qf=Qf)
# Rows/columns are ordered: [x₁, x₂, …, x_{k+1}, u₁, …, u_k].
# The state blocks (left upper) carry Q and Qf; the control blocks (lower right) carry Ru.

p_H = heatmap(Matrix(H);
    color=:viridis, colorbar=true,
    title="Cost Hessian H ($(size(H,1))×$(size(H,2)))",
    xlabel="trajectory index", ylabel="trajectory index",
    yflip=true, aspect_ratio=:equal, size=(500, 450))
p_H

# ## Step 6: Solve for the optimal feasible trajectory
#
# [`optimal_control_trajectory`](@ref) restricts the quadratic objective to
# the feasible affine space and solves the reduced quadratic programme
#
# ```math
# \min_{\alpha} \tfrac{1}{2}\, \alpha^\top (N^\top H N)\, \alpha
#   + \alpha^\top N^\top (H z_p + f).
# ```

z_opt, α_opt, z_p_block, null_basis =
    optimal_control_trajectory(ts, x1, xk1, H, f)

# Verify first-order optimality: the reduced gradient should vanish at α*.
Rred = null_basis' * H * null_basis
rred = null_basis' * (H * Array(z_p_block) + f)
println("Optimality residual ‖Rred α* + rred‖ = ", norm(Rred * α_opt + rred))

# Verify that z_opt lies in the feasible affine space.
println("Affine-space error ‖z_opt − (z_p + N α*)‖ = ",
        norm(Array(z_opt) - (Array(z_p_block) + null_basis * α_opt)))


# ## Step 7: Plot only the optimal trajectory (reference view)
#
# State blocks: `z_opt[Block(t)]` for `t = 1, …, k+1`.
# Control blocks: `z_opt[Block(k+1+t)]` for `t = 1, …, k`.

times_state   = h .* (0:k)
times_control = h .* (0:k-1)

positions  = [Array(z_opt[Block(t)])[1] for t in 1:k+1]
velocities = [Array(z_opt[Block(t)])[2] for t in 1:k+1]
controls   = [Array(z_opt[Block(k+1+t)])[1] for t in 1:k]

p_pos = plot(times_state, positions;
    lw=2, marker=:circle, label="position p(t)",
    xlabel="time (s)", ylabel="p, ṗ",
    title="Double integrator: optimal state trajectory")
plot!(p_pos, times_state, velocities;
    lw=2, marker=:square, linestyle=:dash, label="velocity ṗ(t)")

p_ctrl = plot(times_control, controls;
    lw=2, marker=:diamond, color=:red, label="control u(t)",
    xlabel="time (s)", ylabel="u",
    title="Double integrator: optimal control input")

di_plot = plot(p_pos, p_ctrl; layout=(2, 1), size=(700, 500))
di_plot

# ## Step 8: Spy plot of the nullspace basis matrix
#
# The nullspace basis `null_basis` (size p × r) spans the feasible
# perturbations around the optimal trajectory. A spy plot gives a quick visual
# cue of its sparsity pattern, which can be useful for debugging or
# understanding how many degrees of freedom are active.
#
# The `Plots.spy` recipe works directly with this returned matrix.
spy(null_basis, markersize=2, title="Sparsity pattern of null_basis")

# ## Verification (converted to warnings)
#
# Endpoint checks – issue warnings instead of failing the build.
if !(Array(z_opt[Block(1)]) ≈ x1)
    @warn "Initial state not satisfied"
end
if !(Array(z_opt[Block(k + 1)]) ≈ xk1)
    @warn "Terminal state not satisfied"
end

# Dynamics consistency – report any violations as warnings.
for t in 1:k
    xt  = Array(z_opt[Block(t)])
    xt1 = Array(z_opt[Block(t + 1)])
    ut  = Array(z_opt[Block(k + 1 + t)])
    resid = norm(ts.Ad * xt + ts.Bd * ut - xt1)
    if resid ≥ 1e-10
        @warn "Dynamics violated at step $t (residual = $resid)"
    end
end
println("All endpoint and dynamics checks completed (warnings, if any, were shown above).")

# ## Relationship to the sheaf boundary-value problem
#
# The particular solution `z_p` is the **harmonic extension** of the boundary
# data ``(x_1, x_{k+1})`` to the full trajectory space.  It minimises the
# sheaf Laplacian energy (i.e., the squared coboundary ``\|dz\|^2``).
#
# The null-basis columns ``N_j`` are the *global sections of the controlled
# sheaf* after zeroing the two boundary vertices — exactly the
# endpoint-preserving degrees of freedom.
#
# The optimal trajectory lives in the affine space ``z_p + N\alpha^*``, where
# ``\alpha^*`` minimises the LQR cost restricted to that affine space.  The two
# viewpoints — sheaf / harmonic-extension and quadratic / LQR — combine cleanly
# because the feasible set is a finite-dimensional affine subspace.
#
# ## Step 9: Visualize the nullspace as a trajectory family
#
# The spy plot in Step 8 revealed which rows of `null_basis` are non-zero,
# showing the sparsity structure of the feasible perturbations.  Step 9 brings
# those abstract columns to life: each one becomes a concrete state-control
# trajectory that starts and ends at the prescribed boundary conditions.
#
# Each nullspace basis vector gives one endpoint-preserving control mode.
# We visualize the family members
#
# ```math
# z^{(j)} = z_p + n_j
# ```
#
# where `n_j` is the `j`-th basis column scaled by `amplitude`.  These
# trajectories are feasible but not generally optimal; the optimizer chooses a
# linear combination of these directions that minimises the LQR cost.

family = nullspace_trajectory_family(ts, Array(z_p_block), null_basis; amplitude=3.0)

# Plot the trajectories as perturbations of the optimal trajectory.

max_family = 8
if length(family) > max_family
    println("Plotting first $max_family of $(length(family)) basis trajectories.")
    family = family[1:max_family]
end


p_family_state = plot(
    xlabel="time (s)", ylabel="state",
    title="Nullspace family: state trajectories")
for (j, zf) in enumerate(family)
    pj = [Array(zf[Block(t)])[1] for t in 1:k+1]
    vj = [Array(zf[Block(t)])[2] for t in 1:k+1]
    plot!(p_family_state, times_state, pj; lw=1.5, alpha=0.6, label="n$(j): p")
    plot!(p_family_state, times_state, vj; lw=1.5, alpha=0.6, ls=:dash, label="n$(j): ṗ")
end
plot!(p_family_state, times_state, positions; lw=3, color=:black, label="optimal p")
plot!(p_family_state, times_state, velocities; lw=3, color=:black, ls=:dash, label="optimal ṗ")

p_family_ctrl = plot(
    xlabel="time (s)", ylabel="u",
    title="Nullspace family: control trajectories")
for (j, zf) in enumerate(family)
    uj = [Array(zf[Block(k+1+t)])[1] for t in 1:k]
    plot!(p_family_ctrl, times_control, uj; lw=1.5, alpha=0.6, label="n$(j)")
end
plot!(p_family_ctrl, times_control, controls; lw=3, color=:black, label="optimal")

family_plot = plot(p_family_state, p_family_ctrl; layout=(2, 1), size=(800, 560))
family_plot

# ## What this example establishes
#
# The double integrator illustrates the complete pipeline in its simplest form.
# Steps 8 and 9 together reveal two complementary views of the feasible
# subspace: the spy plot (Step 8) shows the algebraic sparsity structure of the
# null-basis matrix, while the trajectory family (Step 9) shows what those
# basis directions look like as physical state-control trajectories.
#
# The next example, **Vehicle Platoon**, keeps the same control story but
# scales to a multi-agent setting by stacking two double integrators on a
# two-vertex base sheaf, making the connection to network-structured problems
# explicit for the first time.