# # Feedforward Layered Control for 6-Agent Escort Formation ($SE(3)$ Homogeneous Sheaf)
# 
# This example demonstrates eliminating tracking lag in the layered control architecture by 
# combining **joint Tikhonov reference-and-velocity filtering** with a **feedforward LQR control law**.
# 
# ## Theoretical & Mathematical Background
# 
# ### 1. Harmonic Reference & Velocity Solves
# High-level spatial coordination computes the harmonic extension $q^*(t)$ by solving the 
# restricted Laplacian linear system:
# 
# $$H q^*(t) = -L_{AB} p(t)$$
# 
# where $p(t)$ is the vector of target positions. For time-invariant sheaf restrictions $H = L_{AA}$,
# differentiating this identity with respect to time yields the exact harmonic reference velocity $\dot{q}^*(t)$:
# 
# $$H \dot{q}^*(t) = -L_{AB} \dot{p}(t)$$
# 
# A second distributed tree solve on worker processes computes $\dot{q}^*(t)$ using the 
# pre-factorized restricted Laplacian $H$.
# 
# ### 2. Joint Tikhonov Filtering
# On each agent's flight computer, a `JointTikhonovFilter` smooths both $q^*(t)$ and $\dot{q}^*(t)$:
# 
# $$\epsilon \dot{x}_{\text{ref}} = -x_{\text{ref}} + q^*(t), \quad \epsilon \dot{v}_{\text{ref}} = -v_{\text{ref}} + \dot{q}^*(t)$$
# 
# ### 3. Feedforward + Feedback Control Law
# Standard LQR uses feedback $u_{\text{fb}} = -K(x - x_{\text{ref}})$, which incurs tracking lag when following moving references.
# The feedforward controller adds $u_{\text{ff}} = B^\dagger (v_{\text{ref}} - A_c x)$ where $B^\dagger = (B_c^\top B_c)^{-1} B_c^\top$ is the pseudoinverse of the control input matrix, yielding total control input:
# 
# $$u = -K (x - x_{\text{ref}}) + B^\dagger (v_{\text{ref}} - A_c x)$$

using CellularSheaves
using CellularSheaves.Formations
using CellularSheaves.AgentControllers
using CellularSheaves.DistributedLayeredControl
using LinearAlgebra
using Statistics
using Distributed
using Plots
using Printf

# A single house style for every figure below
default(framestyle = :box, grid = true, gridalpha = 0.18, gridstyle = :dot,
    titlefontsize = 10, guidefontsize = 9, legendfontsize = 8, tickfontsize = 8,
    markerstrokewidth = 0, size = (720, 380))

const RING_COLOR = :steelblue
const TARGET_COLOR = :black

# ## Setup 10D Quadrotor Dynamics & Controllers

dyn = QuadrotorDynamics()
DT = 0.05
nx = 10

# Compute Optimal LQR Cost Matrices
Q_diag = [500.0, 500.0, 500.0, 150.0, 150.0, 100.0, 100.0, 100.0, 5.0, 5.0]
Q_lqr = Matrix(Diagonal(Q_diag))
R_lqr = Matrix(Diagonal([0.005, 0.005, 0.005]))

# Construct Standard Feedback LQR Controller
lqr_ctrl = LQRController(dyn, DT, Q_lqr, R_lqr)

# ## SE(3) Homogeneous Coordination Sheaf Construction (D = 4)

const NA = 6
const NT = 1
const TV1 = NA + 1
r_ring = 0.3

# Build escort ring with Agent 1 as observer
sheaf = build_escort_ring(NA, TV1, r_ring; observers=[1])

# Target position p(t) and target velocity p_dot(t)
target1_pos(node, t) = [0.5cos(0.5*t), 0.5sin(0.5*t), 1.5 + 0.1sin(1.0*t), 1.0]
target1_vel(node, t) = [-0.25sin(0.5*t), 0.25cos(0.5*t), 0.1cos(1.0*t), 0.0]

# ## Provision Worker Processes

workers_pids = addprocs(NA; exeflags = "--project=$(Base.active_project())")

@everywhere workers_pids begin
    using CellularSheaves
end

# ## Simulation Framework

STEPS = 200
epsilon = 0.02
init_states = [[r_ring*i/NA, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] for i in 1:NA]
dyns = [QuadrotorDynamics(m=0.5 + 0.05*i, Ixx=0.01 + 0.002*i, Iyy=0.01 + 0.002*i) for i in 1:NA]

# Standard LQR agent configurations used for both position-only and velocity-enhanced tracking
lqr_configs = [(init_states[i], dyns[i], LQRController(dyns[i], DT, Q_lqr, R_lqr).K) for i in 1:NA]

prob_lqr = LayeredControlProblem(;
    sheaf=sheaf,
    target_nodes=[TV1],
    target_trajectory_func=target1_pos,
    agent_configs=lqr_configs,
    dt=DT,
    steps=STEPS,
    r_ring=r_ring
)

prob_vel = LayeredControlProblem(;
    sheaf=sheaf,
    target_nodes=[TV1],
    target_trajectory_func=target1_pos,
    target_velocity_func=target1_vel,
    agent_configs=lqr_configs,
    dt=DT,
    steps=STEPS,
    r_ring=r_ring
)

# Run Position-Only Reference Tracking Simulation
init_distributed_agents!(workers_pids, lqr_configs, DT, epsilon; use_velocity=false)
sim_lqr_res = run_layered_simulation(prob_lqr, workers_pids; mode=:distributed)

# Run Joint Position & Velocity Reference Tracking Simulation
init_distributed_agents!(workers_pids, lqr_configs, DT, epsilon; use_velocity=true)
sim_vel_res = run_layered_simulation(prob_vel, workers_pids; mode=:distributed)

# Clean up worker processes
rmprocs(workers_pids)

# ## Visualizations

# ### 1. Standard LQR Escort Animation
animate_layered_escort(sim_lqr_res; frame_step = 2, filename = "layered_escort_lqr.gif", fps = 15)
nothing # hide

# ![Standard LQR Escort Tracking](layered_escort_lqr.gif)

# ### 2. Joint Velocity-Enhanced LQR Escort Animation
animate_layered_escort(sim_vel_res; frame_step = 2, filename = "layered_escort_feedforward.gif", fps = 15)
nothing # hide

# ![Velocity-Enhanced LQR Escort Tracking](layered_escort_feedforward.gif)

# ### 3. Tracking Lag Elimination Comparison
#
# We compare the per-agent tracking lag $||x_i[1:3] - q^*_i[1:3]||$ over time between 
# position-only reference tracking and joint position + velocity reference tracking.

ts = (1:STEPS) .* DT

err_lqr = zeros(STEPS, NA)
err_vel = zeros(STEPS, NA)

for step in 1:STEPS
    for i in 1:NA
        qstar_3d = sim_lqr_res.qstar_history[step, i, 1:3]
        err_lqr[step, i] = norm(sim_lqr_res.sim_data[step, i, 1:3] - qstar_3d)
        err_vel[step, i] = norm(sim_vel_res.sim_data[step, i, 1:3] - qstar_3d)
    end
end

mean_err_lqr = mean(err_lqr, dims=2)[:]
mean_err_vel = mean(err_vel, dims=2)[:]

p_comp = plot(;
    title="Tracking Lag Elimination (Mean Agent Distance to Harmonic Extension)",
    xlabel="time [s]",
    ylabel="tracking lag [m]",
    yscale=:log10,
    xlims=(0, ts[end]),
    ylims=(1e-4, 1.0),
    legend=:topright,
    bottom_margin=5Plots.mm,
    left_margin=5Plots.mm
)

plot!(p_comp, ts, mean_err_lqr, lw=2.0, color=:crimson, label="Position-Only Reference Tracking")
plot!(p_comp, ts, mean_err_vel, lw=2.0, color=:forestgreen, label="Joint Position & Velocity Tracking")

savefig(p_comp, "tracking_lag_comparison.png")

# ![Tracking Lag Comparison](tracking_lag_comparison.png)
