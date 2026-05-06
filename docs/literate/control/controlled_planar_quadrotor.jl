# # Controlled Trajectory Examples: 3 — Planar Quadrotor
#
# This is the **third example** in the four-part controlled-trajectory
# progression.  The [second example](@ref "Controlled Trajectory Examples: 2 — Vehicle Platoon")
# introduced a multi-agent stacked system.  Here we stay with a single agent
# but add the physical coupling that arises in a real robotic platform: the
# **planar quadrotor**.
#
# ## Physical system
#
# A planar (2-D) quadrotor has position ``(y, z)`` in the vertical plane,
# roll angle ``\phi``, and the corresponding velocities.  The two motor
# thrusts ``u_1`` and ``u_2`` provide total thrust and differential torque.
#
# ```math
# x = \begin{bmatrix} y \\ z \\ \phi \\ \dot{y} \\ \dot{z} \\ \dot{\phi}
#     \end{bmatrix}, \qquad
# u = \begin{bmatrix} u_1 \\ u_2 \end{bmatrix}.
# ```
#
# Linearizing around the hover equilibrium
# ``(y, z, \phi, \dot{y}, \dot{z}, \dot{\phi}) = 0``,
# ``u_1 = u_2 = mg/2`` and using small-angle approximations, gives the
# linear time-invariant model
#
# ```math
# A_c = \begin{bmatrix}
#   0 & 0 & 0 & 1 & 0 & 0 \\
#   0 & 0 & 0 & 0 & 1 & 0 \\
#   0 & 0 & 0 & 0 & 0 & 1 \\
#   0 & 0 & -g & 0 & 0 & 0 \\
#   0 & 0 & 0 & 0 & 0 & 0 \\
#   0 & 0 & 0 & 0 & 0 & 0
# \end{bmatrix},
# \qquad
# B_c = \begin{bmatrix}
#   0 & 0 \\ 0 & 0 \\ 0 & 0 \\ 0 & 0 \\
#   1/m & 1/m \\
#   \ell/(2I) & -\ell/(2I)
# \end{bmatrix},
# ```
#
# where ``g = 9.81\,\text{m/s}^2``, ``m`` is the total mass, ``I`` the moment
# of inertia about the roll axis, and ``\ell`` the half-arm length.  The
# coupling term ``-g\,\phi`` in the ``\ddot{y}`` row captures the lateral
# acceleration due to roll.
#
# **References:**
# - Mahony, Müller & Ganguli (2012), "A non-linear observer for attitude estimation
#   of a fixed-wing unmanned aerial vehicle without GPS measurements",
#   *Control Engineering Practice*.
# - Mueller & D'Andrea (2013), "Stability and control of a quadrocopter despite
#   the complete loss of one, two, or three propellers",
#   *ICRA 2014*.

using CellularSheaves
using LinearAlgebra
using BlockArrays
using Plots

# ## Step 1: Physical parameters and hover linearization

g  = 9.81    # gravitational acceleration (m/s²)
m  = 0.5     # vehicle mass (kg)
I_quad = 0.01    # moment of inertia about roll axis (kg·m²)
ℓ  = 0.25    # half arm length (m)

# The coupling ``-g`` in the (4,3) entry of ``A_c`` is the hallmark of
# the hover linearization: lateral acceleration is proportional to roll angle.

Ac = [0.0  0.0   0.0   1.0  0.0  0.0;
      0.0  0.0   0.0   0.0  1.0  0.0;
      0.0  0.0   0.0   0.0  0.0  1.0;
      0.0  0.0  -g     0.0  0.0  0.0;
      0.0  0.0   0.0   0.0  0.0  0.0;
      0.0  0.0   0.0   0.0  0.0  0.0]

Bc = [0.0           0.0;
      0.0           0.0;
      0.0           0.0;
      0.0           0.0;
      1.0/m         1.0/m;
      ℓ/(2.0*I_quad)  -ℓ/(2.0*I_quad)]

println("State dim n = ", size(Ac, 1))
println("Control dim m = ", size(Bc, 2))

# ## Step 2: Build the ControlledTrajectorySheaf

h = 0.05    # sample period (seconds) — short step for fast quadrotor dynamics
k = 12      # number of time steps (0.6 s total horizon)

F  = EuclideanSheaf{Float64}(fill(6, 1))   # one vertex, 6D stalk
ts = ControlledTrajectorySheaf(F, Ac, Bc, h, k)

println("Ad =\n", round.(ts.Ad; digits=4))
println("Bd =\n", round.(ts.Bd; digits=4))

# ## Step 3: Fix endpoint states
#
# We execute a short lateral manoeuvre: move 0.5 m in the ``y``-direction while
# returning to hover at the end.  The ``z``-position and roll angle start and
# end at zero.

x1  = zeros(6)                    # start at hover
xk1 = [0.5, 0.0, 0.0, 0.0, 0.0, 0.0]   # 0.5 m lateral displacement, back to hover

# ## Step 4: Assemble the LQR objective
#
# Penalize lateral position and roll angle more heavily in the terminal cost
# to discourage residual roll at landing.

n = ts.state_dim
m_ctrl = ts.control_dim

Q  = Diagonal([1.0, 1.0, 2.0, 0.5, 0.5, 0.5])   # state running cost
Ru = Matrix{Float64}(I, m_ctrl, m_ctrl)            # control cost
Qf = Diagonal([10.0, 10.0, 20.0, 1.0, 1.0, 1.0]) # heavier terminal penalty

H, f, _ = lqr_objective(ts, Matrix(Q), Ru; Qf=Matrix(Qf))

# ## Step 5: Solve for the optimal trajectory

z_opt, α_opt, z_p, null_basis = optimal_control_trajectory(ts, x1, xk1, H, f)

println("Free parameters r = ", size(null_basis, 2))
println("All entries finite: ", all(isfinite, Array(z_opt)))

# ## Step 6: Plot the state and control trajectories
#
# The six state components and two control inputs are plotted separately.

times_state   = h .* (0:k)
times_control = h .* (0:k-1)

y_traj   = [Array(z_opt[Block(t)])[1] for t in 1:k+1]
z_traj   = [Array(z_opt[Block(t)])[2] for t in 1:k+1]
phi_traj = [Array(z_opt[Block(t)])[3] for t in 1:k+1]
u1_traj  = [Array(z_opt[Block(k+1+t)])[1] for t in 1:k]
u2_traj  = [Array(z_opt[Block(k+1+t)])[2] for t in 1:k]

p_state = plot(times_state, y_traj;
    lw=2, marker=:circle, label="y (lateral)",
    xlabel="time (s)", ylabel="position (m) / angle (rad)",
    title="Planar quadrotor: state trajectory")
plot!(p_state, times_state, z_traj;
    lw=2, marker=:square, linestyle=:dash, label="z (vertical)")
plot!(p_state, times_state, phi_traj;
    lw=2, marker=:diamond, linestyle=:dot, label="φ (roll)")

p_ctrl = plot(times_control, u1_traj;
    lw=2, marker=:circle, label="u₁(t)",
    xlabel="time (s)", ylabel="thrust deviation (N)",
    title="Planar quadrotor: control inputs")
plot!(p_ctrl, times_control, u2_traj;
    lw=2, marker=:square, linestyle=:dash, label="u₂(t)")

quadrotor_plot = plot(p_state, p_ctrl; layout=(2, 1), size=(700, 500))
quadrotor_plot

# ## Verification

@assert Array(z_opt[Block(1)])     ≈ x1  "Initial state not satisfied"
@assert Array(z_opt[Block(k + 1)]) ≈ xk1 atol=1e-8 "Terminal state not satisfied"
@assert all(isfinite, Array(z_opt))     "Trajectory contains non-finite values"

for t in 1:k
    xt  = Array(z_opt[Block(t)])
    xt1 = Array(z_opt[Block(t + 1)])
    ut  = Array(z_opt[Block(k + 1 + t)])
    @assert norm(ts.Ad * xt + ts.Bd * ut - xt1) < 1e-9 "Dynamics violated at step $t"
end
println("All endpoint and dynamics constraints satisfied.")

# ## What this example adds
#
# Compared with the vehicle platoon, the planar quadrotor introduces:
#
# - A **physically coupled** state space: the lateral acceleration is driven by
#   the roll angle through the ``-g\,\phi`` term in ``A_c``.
# - **Two control inputs** (``u_1``, ``u_2``) whose sum and difference
#   independently govern vertical thrust and roll torque.
# - A **6-dimensional state**, requiring a more structured choice of ``Q`` and
#   ``Q_f`` to balance position and angular degrees of freedom.
#
# The next and final example, **Mass-Spring-Damper Chain**, returns to a
# graph-coupled mechanical system that bridges toward sheaf-theoretic
# network analysis.
