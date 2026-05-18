# # Multi-Quadrotor Target Tracking via LQR on a Trajectory Sheaf
#
# This example extends the single-agent planar-quadrotor example to a
# **two-agent target-tracking** scenario.  Each quadrotor must fly from a
# hover equilibrium to a distinct 2-D target position — a prototypical
# multi-target assignment problem.
#
# The controlled-trajectory-sheaf framework encodes both agents and both
# targets in a single, jointly-feasible LQR problem.  The feasible set is
# the affine subspace
#
# ```math
# \mathcal{T}(x_1, x_{k+1}) = \{ z_p + N\alpha : \alpha \in \mathbb{R}^r \},
# ```
#
# and the LQR cost is minimised by solving a **reduced quadratic programme**
# whose Hessian is factored with the sparse LDL⊤ solver
# `ChordalLDLt` from `CliqueTrees.Multifrontal`.
#
# ## Physical model
#
# Each quadrotor is linearised around hover following the same model used in
# [Controlled Trajectory Examples: 3 — Planar Quadrotor](@ref):
#
# ```math
# x_i = \begin{bmatrix} y_i \\ z_i \\ \phi_i \\ \dot{y}_i \\ \dot{z}_i \\ \dot{\phi}_i
#         \end{bmatrix},
# \qquad u_i = \begin{bmatrix} u_{i,1} \\ u_{i,2} \end{bmatrix},
# \qquad
# \dot{x}_i = A_c x_i + B_c u_i.
# ```
#
# The two agents are **decoupled** (no aerodynamic interaction at this level
# of modelling), so the composite state-space is block-diagonal:
#
# ```math
# A_c^{\text{full}} = \begin{bmatrix} A_c & 0 \\ 0 & A_c \end{bmatrix},
# \qquad
# B_c^{\text{full}} = \begin{bmatrix} B_c & 0 \\ 0 & B_c \end{bmatrix}.
# ```
#
# Encoding this as an `EuclideanSheaf` with **two** vertices, each carrying a
# 6-dimensional stalk (one per quadrotor), makes the multi-agent structure
# explicit.
#
# **References:**
# - Hansen, J. & Ghrist, R. (2011). "Toward a Spectral Theory of Cellular Sheaves."
#   *Journal of Applied and Computational Topology*.
# - Mueller & D'Andrea (2013). "Stability and control of a quadrocopter despite
#   the complete loss of one, two, or three propellers." *ICRA 2014*.
# - Fairbanks, J. (2024). *CellularSheaves.jl* — trajectory sheaf and LQR docs.

using CellularSheaves
using CellularSheaves.TrajectorySheaves: nullspace_trajectory_family
using LinearAlgebra
using BlockArrays
using Plots
using SparseArrays

# ## Step 1: Physical parameters

g      = 9.81
m_veh  = 0.5
I_quad = 0.01
ℓ      = 0.25

# Single-quadrotor hover linearization (same as in controlled_planar_quadrotor.jl).

Ac_single = [0.0  0.0   0.0   1.0  0.0  0.0;
             0.0  0.0   0.0   0.0  1.0  0.0;
             0.0  0.0   0.0   0.0  0.0  1.0;
             0.0  0.0  -g     0.0  0.0  0.0;
             0.0  0.0   0.0   0.0  0.0  0.0;
             0.0  0.0   0.0   0.0  0.0  0.0]

Bc_single = [0.0              0.0;
             0.0              0.0;
             0.0              0.0;
             0.0              0.0;
             1.0/m_veh        1.0/m_veh;
             ℓ/(2.0*I_quad)  -ℓ/(2.0*I_quad)]

n_single = size(Ac_single, 1)   # 6
m_single = size(Bc_single, 2)   # 2

# ## Step 2: Stack into a two-quadrotor composite system

n = 2 * n_single   # composite state dim = 12
m = 2 * m_single   # composite control dim = 4

Ac = zeros(n, n)
Bc = zeros(n, m)
Ac[1:n_single, 1:n_single]               = Ac_single
Ac[n_single+1:end, n_single+1:end]       = Ac_single
Bc[1:n_single, 1:m_single]               = Bc_single
Bc[n_single+1:end, m_single+1:end]       = Bc_single

println("Composite state dim  n = ", n)
println("Composite control dim m = ", m)

# ## Step 3: Build the ControlledTrajectorySheaf
#
# The base sheaf has **two vertices**, each with a 6-dimensional stalk
# (one per quadrotor).  The path-graph structure of the trajectory sheaf
# uses this as the per-step fiber.

h = 0.05    # sample period (s)
k = 16      # time steps (0.8 s horizon)

F  = EuclideanSheaf{Float64}([n_single, n_single])
ts = ControlledTrajectorySheaf(F, Ac, Bc, h, k)

println("Steps k = ", ts.k, "  (horizon = $(k*h) s)")

# ## Step 4: Endpoint conditions (two targets)
#
# Quadrotor 1 begins at hover and must reach a point 1 m to the right
# and 0.5 m up, returning to zero roll and velocity.
# Quadrotor 2 begins at the same height but 2 m to the right and must
# reach a point 3 m to the right at the same altitude.
#
# State layout per agent: [y, z, φ, ẏ, ż, φ̇].

x1  = zeros(n)
xk1 = zeros(n)

# Quadrotor 1: y = 0 → 1.0 m,  z = 0 → 0.5 m
xk1[1] = 1.0    # y₁ target
xk1[2] = 0.5    # z₁ target

# Quadrotor 2: y = 2.0 → 3.0 m,  z = 0 → 0.0 m
x1[7]  = 2.0    # y₂ initial
xk1[7] = 3.0    # y₂ target

println("Initial state:  ", x1)
println("Terminal state: ", xk1)

# ## Step 5: LQR objective
#
# We use a diagonal state weight that penalises position and roll more
# heavily than velocity, and a stronger terminal weight to drive the agents
# to rest at the targets.

Q_single  = Diagonal([1.0, 1.0, 2.0, 0.5, 0.5, 0.5])
Qf_single = Diagonal([10.0, 10.0, 20.0, 1.0, 1.0, 1.0])
Ru_single = Matrix{Float64}(I, m_single, m_single)

Q  = Matrix(cat(Q_single,  Q_single;  dims=(1,2)))
Qf = Matrix(cat(Qf_single, Qf_single; dims=(1,2)))
Ru = Matrix(cat(Ru_single, Ru_single; dims=(1,2)))

H, f, _ = lqr_objective(ts, Q, Ru; Qf=Qf)

# The cost Hessian H is block-diagonal: state cost blocks (Q or Qf) alternate
# with control cost blocks (Ru).  Visualising H as a heatmap makes this
# band structure visible.

p_H = heatmap(Matrix(H);
    color=:viridis, colorbar=true,
    title="Cost Hessian H ($(size(H,1))×$(size(H,2)))",
    xlabel="trajectory index", ylabel="trajectory index",
    yflip=true, aspect_ratio=:equal, size=(560, 500))
p_H

# ## Step 6: Solve the LQR problem via the trajectory affine space
#
# `optimal_control_trajectory` uses the following pipeline internally:
#
# 1. Compute the particular feasible solution `z_p` (harmonic extension).
# 2. Extract the nullspace basis `N` of the coboundary map.
# 3. Form the **reduced Hessian** `R_red = Nᵀ H N` and **reduced gradient**
#    `r_red = Nᵀ(H z_p + f)`.
# 4. Factor `R_red` using `ChordalLDLt` (LDL⊤ factorisation with
#    `RowMaximum` fill-reducing pivot) and solve `R_red α* = -r_red`.
#
# The sparse LDL⊤ path avoids densifying H and exploits the banded
# structure of the reduced Hessian.

z_opt, α_opt, z_p_block, null_basis =
    optimal_control_trajectory(ts, x1, xk1, H, f)

r = size(null_basis, 2)
println("Nullspace dimension r = ", r, "  (free parameters)")
println("All entries finite: ", all(isfinite, Array(z_opt)))

# First-order optimality check.
q_p  = Array(z_p_block)
Rred = Symmetric(null_basis' * H * null_basis)
rred = null_basis' * (H * q_p + f)
println("Optimality residual ‖R_red α* + r_red‖ = ",
        norm(Rred * α_opt + rred))

# ## Step 7: Extract per-quadrotor trajectories
#
# State ordering in each block: [y₁, z₁, φ₁, ẏ₁, ż₁, φ̇₁, y₂, z₂, φ₂, ẏ₂, ż₂, φ̇₂].

times_state   = h .* (0:k)
times_control = h .* (0:k-1)

y1   = [Array(z_opt[Block(t)])[1] for t in 1:k+1]
z1   = [Array(z_opt[Block(t)])[2] for t in 1:k+1]
phi1 = [Array(z_opt[Block(t)])[3] for t in 1:k+1]
y2   = [Array(z_opt[Block(t)])[7] for t in 1:k+1]
z2   = [Array(z_opt[Block(t)])[8] for t in 1:k+1]
phi2 = [Array(z_opt[Block(t)])[9] for t in 1:k+1]

u11  = [Array(z_opt[Block(k+1+t)])[1] for t in 1:k]
u12  = [Array(z_opt[Block(k+1+t)])[2] for t in 1:k]
u21  = [Array(z_opt[Block(k+1+t)])[3] for t in 1:k]
u22  = [Array(z_opt[Block(k+1+t)])[4] for t in 1:k]

# ## Step 8: Time-series plots
#
# Position and roll for both quadrotors, plus the four control inputs.

p_pos = plot(times_state, y1;
    lw=2, marker=:circle, label="Q1: y",
    xlabel="time (s)", ylabel="position (m) / angle (rad)",
    title="Multi-quadrotor: position and roll")
plot!(p_pos, times_state, z1; lw=2, marker=:square,  linestyle=:dash, label="Q1: z")
plot!(p_pos, times_state, phi1; lw=2, marker=:diamond, linestyle=:dot, label="Q1: φ")
plot!(p_pos, times_state, y2;  lw=2, marker=:circle,  color=:darkorange, label="Q2: y")
plot!(p_pos, times_state, z2;  lw=2, marker=:square,  color=:darkorange,
      linestyle=:dash, label="Q2: z")
plot!(p_pos, times_state, phi2; lw=2, marker=:diamond, color=:darkorange,
      linestyle=:dot, label="Q2: φ")

p_ctrl = plot(times_control, u11;
    lw=2, marker=:circle, label="Q1: u₁",
    xlabel="time (s)", ylabel="thrust deviation (N)",
    title="Multi-quadrotor: control inputs")
plot!(p_ctrl, times_control, u12; lw=2, marker=:square, linestyle=:dash, label="Q1: u₂")
plot!(p_ctrl, times_control, u21; lw=2, marker=:circle, color=:darkorange, label="Q2: u₁")
plot!(p_ctrl, times_control, u22; lw=2, marker=:square, color=:darkorange,
      linestyle=:dash, label="Q2: u₂")

multi_ts_plot = plot(p_pos, p_ctrl; layout=(2, 1), size=(800, 560))
multi_ts_plot

# ## Step 9: Trajectories in the physical (y, z) plane
#
# Both quadrotors trace arcs in their respective vertical planes.  The roll
# coupling drives a simultaneous pitch-up / pitch-down manoeuvre as each
# vehicle accelerates horizontally.

t_norm = range(0.0, 1.0; length=k+1)

p_plane1 = scatter(y1, z1;
    marker_z=t_norm, color=:plasma, colorbar=false,
    label="", xlabel="y (m)", ylabel="z (m)",
    title="Q1 trajectory  (y,z)-plane",
    markerstrokewidth=0, markersize=6)
plot!(p_plane1, y1, z1; lw=1.5, color=:gray, label="")
scatter!(p_plane1, [y1[1]],   [z1[1]];   color=:green, markersize=9,
    markershape=:star5, label="start")
scatter!(p_plane1, [y1[end]], [z1[end]]; color=:red,   markersize=9,
    markershape=:star5, label="end")

p_plane2 = scatter(y2, z2;
    marker_z=t_norm, color=:viridis, colorbar=false,
    label="", xlabel="y (m)", ylabel="z (m)",
    title="Q2 trajectory  (y,z)-plane",
    markerstrokewidth=0, markersize=6)
plot!(p_plane2, y2, z2; lw=1.5, color=:gray, label="")
scatter!(p_plane2, [y2[1]],   [z2[1]];   color=:green, markersize=9,
    markershape=:star5, label="start")
scatter!(p_plane2, [y2[end]], [z2[end]]; color=:red,   markersize=9,
    markershape=:star5, label="end")

plane_plot = plot(p_plane1, p_plane2; layout=(1, 2), size=(900, 380))
plane_plot

# ## Step 10: Nullspace basis as a trajectory family
#
# The columns of `null_basis` are the free degrees of freedom around the
# optimal trajectory: endpoint-preserving perturbations that the LQR optimiser
# weighted against each other to minimise cost.
#
# We visualise the first few basis directions as concrete state trajectories.
# Each member of the family shares the same boundary conditions but differs
# in the interior.  This reveals the "shape" of the feasible subspace
# — columns that perturb mostly in position versus mostly in control.
#
# The spy plot shows the sparsity structure: non-zero rows correspond to
# time steps where the basis direction is active.

p_spy = spy(null_basis;
    markersize=2,
    title="Null-basis sparsity  ($(size(null_basis,1)) × $r)")
p_spy

# Now construct the perturbed trajectory family.

family = nullspace_trajectory_family(ts, q_p, null_basis; amplitude=0.5)

max_show = min(length(family), 6)

p_fam_y1 = plot(xlabel="time (s)", ylabel="y₁ (m)",
    title="Nullspace family: Q1 lateral position")
p_fam_y2 = plot(xlabel="time (s)", ylabel="y₂ (m)",
    title="Nullspace family: Q2 lateral position")

for j in 1:max_show
    zf  = family[j]
    y1j = [Array(zf[Block(t)])[1] for t in 1:k+1]
    y2j = [Array(zf[Block(t)])[7] for t in 1:k+1]
    plot!(p_fam_y1, times_state, y1j; lw=1.5, alpha=0.65, label="n$j")
    plot!(p_fam_y2, times_state, y2j; lw=1.5, alpha=0.65, label="n$j")
end

plot!(p_fam_y1, times_state, y1; lw=3, color=:black, label="optimal")
plot!(p_fam_y2, times_state, y2; lw=3, color=:black, label="optimal")

family_plot = plot(p_fam_y1, p_fam_y2; layout=(1, 2), size=(900, 360))
family_plot

# ## Verification

if !(Array(z_opt[Block(1)]) ≈ x1)
    @warn "Initial state not satisfied"
end
if !(isapprox(Array(z_opt[Block(k+1)]), xk1; atol=1e-8))
    @warn "Terminal state not satisfied"
end
if !all(isfinite, Array(z_opt))
    @warn "Trajectory contains non-finite values"
end

for t in 1:k
    xt  = Array(z_opt[Block(t)])
    xt1 = Array(z_opt[Block(t+1)])
    ut  = Array(z_opt[Block(k+1+t)])
    resid = norm(ts.Ad * xt + ts.Bd * ut - xt1)
    if resid >= 1e-9
        @warn "Dynamics violated at step $t (residual = $resid)"
    end
end
println("Endpoint and dynamics checks complete (any violations shown as warnings above).")

# ## Key takeaways
#
# - **Multi-agent LQR via a trajectory sheaf.** Stacking two planar-quadrotor
#   models on a two-vertex base sheaf yields a composite system whose feasible
#   trajectories are jointly optimised.  The boundary conditions impose
#   *independent* targets for each agent; the shared affine structure and cost
#   matrix couple them through the optimiser.
#
# - **Sparse LDL⊤ solver.** The reduced Hessian `N⊤ H N` inherits the banded
#   structure of the path-graph coboundary.  `ChordalLDLt` factorises it
#   without densification, keeping memory proportional to the trajectory length
#   rather than its square.
#
# - **Nullspace trajectory family.** The `r` columns of `null_basis` span the
#   interior degrees of freedom.  Visualising them as perturbed trajectories
#   gives geometric intuition for the LQR optimiser's trade-off: directions
#   that produce large state excursions incur high cost and are suppressed;
#   directions that spread the effort smoothly across the horizon are preferred.
