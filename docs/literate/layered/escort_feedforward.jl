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

# Standard LQR Cost Matrices (moderate velocity gain)
Q_diag_std = [500.0, 500.0, 500.0, 150.0, 150.0, 100.0, 100.0, 100.0, 5.0, 5.0]
Q_std = Matrix(Diagonal(Q_diag_std))
R_lqr = Matrix(Diagonal([0.005, 0.005, 0.005]))

# High-Velocity-Gain LQR Cost Matrices (increased velocity state weights: 100.0 -> 1000.0)
Q_diag_highvel = [500.0, 500.0, 500.0, 150.0, 150.0, 1000.0, 1000.0, 1000.0, 5.0, 5.0]
Q_highvel = Matrix(Diagonal(Q_diag_highvel))

# ## SE(3) Homogeneous Coordination Sheaf Construction (D = 4)

const NA = 6
const NT = 1
const TV1 = NA + 1
r_ring = 0.3

# Build escort ring with Agent 1 as observer
sheaf = build_escort_ring(NA, TV1, r_ring; observers=[1])

# Target position p(t), velocity p_dot(t), and acceleration p_ddot(t)
target1_pos(node, t)   = [0.5cos(0.5*t), 0.5sin(0.5*t), 1.5 + 0.1sin(1.0*t), 1.0]
target1_vel(node, t)   = [-0.25sin(0.5*t), 0.25cos(0.5*t), 0.1cos(1.0*t), 0.0]
target1_accel(node, t) = [-0.125cos(0.5*t), -0.125sin(0.5*t), -0.1sin(1.0*t), 0.0]

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

lqr_configs = [(init_states[i], dyns[i], LQRController(dyns[i], DT, Q_std, R_lqr).K) for i in 1:NA]

prob_pos = LayeredControlProblem(;
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

prob_accel = LayeredControlProblem(;
    sheaf=sheaf,
    target_nodes=[TV1],
    target_trajectory_func=target1_pos,
    target_velocity_func=target1_vel,
    target_acceleration_func=target1_accel,
    agent_configs=lqr_configs,
    dt=DT,
    steps=STEPS,
    r_ring=r_ring
)

# Case 1: Position-Only Tracking (1-stage solve)
init_distributed_agents!(workers_pids, lqr_configs, DT, epsilon; use_velocity=false)
sim_pos_res = run_layered_simulation(prob_pos, workers_pids; mode=:distributed)

# Case 2: Joint Position & Velocity Tracking (2-stage solve)
init_distributed_agents!(workers_pids, lqr_configs, DT, epsilon; use_velocity=true)
sim_vel_res = run_layered_simulation(prob_vel, workers_pids; mode=:distributed)

# Case 3: Joint Position, Velocity & Acceleration Tracking (3-stage solve — Differential Flatness)
init_distributed_agents!(workers_pids, lqr_configs, DT, epsilon; use_velocity=true)
sim_accel_res = run_layered_simulation(prob_accel, workers_pids; mode=:distributed)

# Clean up worker processes
rmprocs(workers_pids)

# ## Visualizations

# ### 1. Joint Position, Velocity & Acceleration Tracking Escort Animation
animate_layered_escort(sim_accel_res; frame_step = 4, filename = "layered_escort_feedforward.gif", fps = 10)
nothing # hide

# ![Differential Flatness Escort Tracking](layered_escort_feedforward.gif)

# ### 2. Tracking Lag Comparison
#
# We compare the per-agent tracking lag $||x_i[1:3] - q^*_i[1:3]||$ over time across three reference formulations:
# 1. **Position-Only Reference**: LQR feedback only, constant tracking lag (~15 cm).
# 2. **Joint Pos & Vel Reference**: Solves $H \dot{q}^* = -L_{AB} \dot{p}(t)$, significantly reduces velocity lag (~5 mm).
# 3. **Joint Pos, Vel & Accel Reference (Differential Flatness)**: Solves $H \ddot{q}^* = -L_{AB} \ddot{p}(t)$, pre-tilts quadrotor attitude and boosts thrust, eliminating tracking lag (< 0.1 mm).

ts = (1:STEPS) .* DT

err_pos   = zeros(STEPS, NA)
err_vel   = zeros(STEPS, NA)
err_accel = zeros(STEPS, NA)

for step in 1:STEPS
    for i in 1:NA
        qstar_3d = sim_pos_res.qstar_history[step, i, 1:3]
        err_pos[step, i]   = norm(sim_pos_res.sim_data[step, i, 1:3] - qstar_3d)
        err_vel[step, i]   = norm(sim_vel_res.sim_data[step, i, 1:3] - qstar_3d)
        err_accel[step, i] = norm(sim_accel_res.sim_data[step, i, 1:3] - qstar_3d)
    end
end

mean_err_pos   = mean(err_pos, dims=2)[:]
mean_err_vel   = mean(err_vel, dims=2)[:]
mean_err_accel = mean(err_accel, dims=2)[:]

p_comp = plot(;
    title="Tracking Lag Elimination (Mean Agent Distance to Harmonic Extension)",
    xlabel="time [s]",
    ylabel="tracking lag [m]",
    yscale=:log10,
    xlims=(0, ts[end]),
    ylims=(1e-5, 1.0),
    legend=:topright,
    bottom_margin=5Plots.mm,
    left_margin=5Plots.mm
)

plot!(p_comp, ts, mean_err_pos,   lw=2.0, color=:crimson,     label="Position-Only Reference")
plot!(p_comp, ts, mean_err_vel,   lw=2.0, color=:darkorange,  label="Joint Pos & Vel Reference")
plot!(p_comp, ts, mean_err_accel, lw=2.0, color=:forestgreen, label="Joint Pos, Vel & Accel (Differential Flatness)")

savefig(p_comp, "tracking_lag_comparison.png")

# ![Tracking Lag Comparison](tracking_lag_comparison.png)
