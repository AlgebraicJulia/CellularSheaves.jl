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
Q_diag = [150.0, 150.0, 150.0, 50.0, 50.0, 30.0, 30.0, 30.0, 1.0, 1.0]
Q_lqr = Matrix(Diagonal(Q_diag))
R_lqr = Matrix(Diagonal([0.01, 0.01, 0.01]))

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

# Run distributed simulation
init_distributed_agents!(workers_pids, agent_configs, DT, epsilon)
sim_d, q_d = run_layered_simulation(sheaf, workers_pids, [TV1], target1_pos, DT, STEPS; mode=:distributed)

# Run centralized simulation
init_distributed_agents!(workers_pids, agent_configs, DT, epsilon)
sim_c, q_c = run_layered_simulation(sheaf, workers_pids, [TV1], target1_pos, DT, STEPS; mode=:centralised)

divergence = maximum(abs.(sim_d .- sim_c))
@printf("Max divergence between centralized and distributed 6-agent simulation: %.3e\n", divergence)

# Clean up worker processes
rmprocs(workers_pids)

# ## Multi-Projection Trajectory & Attitude Dynamics Visualization

lims_xy = (-1.8, 1.8)
lims_z = (0.0, 1.8)
lims_rel = (-1.5, 1.5)
ts = (1:STEPS) .* DT

anim = @animate for k in 1:2:STEPS
    t_curr = k * DT
    
    p1 = plot(; aspect_ratio = 1, xlims = lims_xy, ylims = lims_xy,
              xlabel = "x position (m)", ylabel = "y position (m)",
              title = "World Top-Down View (x-y Plane)", legend = false)
    target_orbit_t = range(0, STEPS*DT; length=100) # slow moving target trajectory path
    plot!(p1, [target1_pos(TV1, t)[1] for t in target_orbit_t], [target1_pos(TV1, t)[2] for t in target_orbit_t]; color = :gray80, linestyle = :dot, linewidth = 1)
    scatter!(p1, [target1_pos(TV1, t_curr)[1]], [target1_pos(TV1, t_curr)[2]]; marker = :star5, markersize = 10, color = TARGET_COLOR)
    for i in 1:NA
        plot!(p1, sim_d[1:k, i, 1], sim_d[1:k, i, 2];
              seriestype = :path, marker = :circle, markersize = 3, alpha = 0.6,
              linewidth = 1.4, color = RING_COLOR)
    end
    
    ## Draw communication topology (ring)
    ring_x = [sim_d[k, i, 1] for i in 1:NA]
    push!(ring_x, sim_d[k, 1, 1])
    ring_y = [sim_d[k, i, 2] for i in 1:NA]
    push!(ring_y, sim_d[k, 1, 2])
    plot!(p1, ring_x, ring_y; color = :gray80, linestyle = :dot, linewidth = 1)

    p2 = plot(; aspect_ratio = 1, xlims = lims_rel, ylims = lims_rel,
              xlabel = "rel x to target (m)", ylabel = "rel y to target (m)",
              title = "Target-Centered Escort Ring (1.2m)", legend = false)
    circ_ang = range(0, 2π; length = 100) # reference 1.2m circle
    plot!(p2, r_ring .* cos.(circ_ang), r_ring .* sin.(circ_ang); color = :gray80, linestyle = :dash, linewidth = 1)
    scatter!(p2, [0.0], [0.0]; marker = :star5, markersize = 10, color = TARGET_COLOR)
    ## Historical trajectory in target-centered frame
    for i in 1:NA
        rel_x_hist = [sim_d[step, i, 1] - target1_pos(TV1, ts[step])[1] for step in 1:k]
        rel_y_hist = [sim_d[step, i, 2] - target1_pos(TV1, ts[step])[2] for step in 1:k]
        plot!(p2, rel_x_hist, rel_y_hist; 
              seriestype = :path, marker = :circle, markersize = 3, alpha = 0.6,
              linewidth = 1.4, color = RING_COLOR)
    end
    
    ## Draw communication topology (ring) in target-centered frame
    ring_rel_x = [sim_d[k, i, 1] - target1_pos(TV1, t_curr)[1] for i in 1:NA]
    push!(ring_rel_x, sim_d[k, 1, 1] - target1_pos(TV1, t_curr)[1])
    ring_rel_y = [sim_d[k, i, 2] - target1_pos(TV1, t_curr)[2] for i in 1:NA]
    push!(ring_rel_y, sim_d[k, 1, 2] - target1_pos(TV1, t_curr)[2])
    plot!(p2, ring_rel_x, ring_rel_y; color = :gray80, linestyle = :dot, linewidth = 1)

    p3 = plot(; xlabel = "time (s)", ylabel = "roll angle ϕ (deg)",
              title = "Roll Tilt Dynamics ϕ(t)", legend = false, xlims = (0, 10.0), ylims = (-10.0, 10.0))
    for i in 1:NA
        plot!(p3, ts[1:k], rad2deg.(sim_d[1:k, i, 4]); linewidth = 1.2, color = RING_COLOR)
    end

    p4 = plot(; xlabel = "time (s)", ylabel = "pitch angle θ (deg)",
              title = "Pitch Tilt Dynamics θ(t)", legend = false, xlims = (0, 10.0), ylims = (-10.0, 10.0))
    for i in 1:NA
        plot!(p4, ts[1:k], rad2deg.(sim_d[1:k, i, 5]); linewidth = 1.2, color = RING_COLOR)
    end

    plot(p1, p2, p3, p4; layout = (2, 2), size = (900, 700),
         plot_title = @sprintf("6-Agent SE(3) Moving Escort Ring (t = %.2f s)", t_curr))
end
gif(anim, "layered_escort_tracking.gif"; fps = 15)
nothing # hide

# ![6-Agent Escort Tracking Projections](layered_escort_tracking.gif)
