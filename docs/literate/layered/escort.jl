# # Layered Control Architecture for 6-Agent Escort Formation ($SE(3)$ Homogeneous Sheaf)
# 
# This example demonstrates scaling the layered control architecture to a 6-agent formation 
# escorting a slow-moving target in 3D space using an **$SE(3)$ Homogeneous Affine Cellular Sheaf**.
# 
# ## High-Level API
# 
# We use the new `Formations` and `DistributedLayeredControl` modules to significantly 
# simplify the setup of the escort ring and the execution of the distributed local controllers.

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

# ## Setup 10D Quadrotor Dynamics & DARE Solver

dyn = QuadrotorDynamics()
DT = 0.05
nx = 10

# Compute Optimal LQR Gain via Discrete Algebraic Riccati Equation (DARE)
Q_diag = [500.0, 500.0, 500.0, 150.0, 150.0, 100.0, 100.0, 100.0, 5.0, 5.0]
Q_lqr = Matrix(Diagonal(Q_diag))
R_lqr = Matrix(Diagonal([0.005, 0.005, 0.005]))

lqr_controller = LQRController(dyn, DT, Q_lqr, R_lqr)
K_lqr = lqr_controller.K

# ## SE(3) Homogeneous Coordination Sheaf Construction (D = 4)

# 6 agents in a single escort ring.
const NA = 6
const NT = 1
const TV1 = NA + 1
r_ring = 0.3

# Build the escort ring, with Agent 1 pinned as the observer
sheaf = build_escort_ring(NA, TV1, r_ring; observers=[1])

# Fast-moving target trajectory
target1_pos(node, t) = [0.5cos(0.5*t), 0.5sin(0.5*t), 1.5 + 0.1sin(1.0*t), 1.0]

# ## Provision Worker Processes

# Provision exactly NA worker processes (one for each agent's flight computer)
workers_pids = addprocs(NA; exeflags = "--project=$(Base.active_project())")

# Provide the cellular sheaves environment to all workers
@everywhere workers_pids begin
    using CellularSheaves
end

# ## Simulation Framework

STEPS = 200
epsilon = 0.02
# Start agents in a line along the x-axis
init_states = [[r_ring*i/NA, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] for i in 1:NA]

# Create heterogeneous agent properties (varying mass and inertia)
dyns = [QuadrotorDynamics(m=0.5 + 0.05*i, Ixx=0.01 + 0.002*i, Iyy=0.01 + 0.002*i) for i in 1:NA]
K_lqrs = [LQRController(d, DT, Q_lqr, R_lqr).K for d in dyns]
agent_configs = [(init_states[i], dyns[i], K_lqrs[i]) for i in 1:NA]

prob = LayeredControlProblem(sheaf, [TV1], target1_pos, agent_configs, DT, STEPS, r_ring)

# Run distributed simulation
init_distributed_agents!(workers_pids, agent_configs, DT, epsilon)
sim_d_res = run_layered_simulation(prob, workers_pids; mode=:distributed)

# Run centralized simulation
init_distributed_agents!(workers_pids, agent_configs, DT, epsilon)
sim_c_res = run_layered_simulation(prob, workers_pids; mode=:centralised)

divergence = maximum(abs.(sim_d_res.sim_data .- sim_c_res.sim_data))
@printf("Max divergence between centralized and distributed 6-agent simulation: %.3e\n", divergence)

# Clean up worker processes
rmprocs(workers_pids)

# ## Multi-Projection Trajectory & Attitude Dynamics Visualization

# We can now easily animate the simulation using the Plots.jl recipe framework.
# The `animate_layered_escort` function abstracts away the complex multi-panel visualization.
# (Note: this relies on the `CellularSheavesPlots` package extension)

animate_layered_escort(sim_d_res; frame_step = 2, filename = "layered_escort_tracking.gif", fps = 15)
nothing # hide

# ![6-Agent Escort Tracking Projections](layered_escort_tracking.gif)
