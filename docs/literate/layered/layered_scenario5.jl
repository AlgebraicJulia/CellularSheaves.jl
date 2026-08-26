# # Layered Control Architecture for Multi-Quadrotor Tracking (Scenario 5)
# 
# This example demonstrates a layered control architecture designed for multi-agent target tracking 
# where high-level spatial coordination is handled by a static cellular sheaf and low-level tracking is 
# managed by independent LQR controllers running on distributed Julia worker processes.
# 
# ## Architecture Overview
# 
# The control pipeline is structured as follows:
# 
# 1. **2D Coordination Sheaf**: Operates on spatial positions ($y, z$), resolving the conflict 
#    between target tracking and agent consensus by computing a harmonic extension $\mathbf{q}^*$.
# 2. **Distributed Tree Solve**: Solves for $\mathbf{q}^*$ across parallel worker processes 
#    using tree-parallel message passing on the clique tree of the restricted Laplacian.
# 3. **Local Stalk Embedding & Tikhonov Filter**: Each worker embeds its 2D spatial reference into its 
#    full 6D state space ($[y^*, z^*, 0, 0, 0, 0]^\top$) and smooths it using a Tikhonov reference filter.
# 4. **LQR Control**: Drives each agent to the filtered reference on its dedicated worker process using 
#    an optimal state-feedback gain computed via the Discrete Algebraic Riccati Equation (DARE).
# 
# This separation of concerns ensures that the coordination Laplacian $H$ remains strictly positive-definite 
# (full rank) and allows high-level spatial coordination to run efficiently across worker processes.

using CellularSheaves
using CellularSheaves.ControlSheaves.AgentControllers
using CellularSheaves.ControlSheaves.DistributedLayeredControl
using LinearAlgebra
using Plots
using Printf
using Distributed

# A single house style for every figure below — matches the escort.jl style
default(framestyle = :box, grid = true, gridalpha = 0.18, gridstyle = :dot,
    titlefontsize = 10, guidefontsize = 9, legendfontsize = 8, tickfontsize = 8,
    markerstrokewidth = 0, size = (1200, 380))

# ## Setup Dynamics & Constants

h = 0.05
T_end = 2.0
steps = Int(T_end / h) + 1

# We use the PlanarQuadrotorDynamics built into AgentControllers.jl
dyn = PlanarQuadrotorDynamics()
Ac, Bc = AgentControllers.continuous_matrices(dyn)

nx = size(Ac, 1)
nu = size(Bc, 2)

# Compute Optimal LQR Gain (Discrete-time)
Q_diag = zeros(nx)
Q_diag[1:2] .= 10000.0
Q_diag[3] = 50.0
Q_diag[4] = 500.0
Q_diag[5:6] .= 10.0
Q_lqr = Matrix{Float64}(Diagonal(Q_diag))
R_lqr = Matrix{Float64}(I, nu, nu) * 0.0001

ctrl = LQRController(dyn, h, Q_lqr, R_lqr)

# ## 2D Coordination Sheaf Construction

# Scenario 5 Configuration:
# - Agent 1 tracks Target 1 in z
# - Agent 2 tracks Target 2 in yz
# - Agents agree in y (Consensus)
# We set D = 2 for spatial positions [y, z].
D = 2
NA = 2
NT = 2
TotalV = NA + NT

sheaf = EuclideanSheaf{Float64}(fill(D, TotalV))

# Coordinate projection matrices
R_y  = [1.0 0.0]
R_z  = [0.0 1.0]
R_yz = Matrix{Float64}(I, 2, 2)

# Setup coordination edges
add_sheaf_edge!(sheaf, 1, 2, R_y, R_y)
add_sheaf_edge!(sheaf, 1, 3, R_z, R_z)
add_sheaf_edge!(sheaf, 2, 4, R_yz, R_yz)

# ## Simulation Framework

# Ensure we have enough worker processes
workers_pids = workers()
if length(workers_pids) < NA
    addprocs(NA - length(workers_pids);
        exeflags = ["--project=$(Base.active_project())",
            "--sysimage=$(unsafe_string(Base.JLOptions().image_file))"])
    workers_pids = workers()
    @eval @everywhere using CellularSheaves
    @eval @everywhere using CellularSheaves.ControlSheaves.AgentControllers
    @eval @everywhere using CellularSheaves.ControlSheaves.DistributedLayeredControl
end

function run_scenario(target_type=:bobbing, mode=:distributed)
    local t1_pos, t2_pos
    if target_type == :bobbing
        omega = 2π * 2 / (40 * h)
        t1_pos = t -> [0.0, 1.0 + 0.3sin(omega*t)]
        t2_pos = t -> [1.5, 2.0 + 0.3sin(omega*t)]
    else
        t1_pos = t -> [0.0, 1.0]
        t2_pos = t -> [1.5, 2.0]
    end

    target_trajectory_func = (v, t) -> v == 3 ? t1_pos(t) : t2_pos(t)
    
    init_states = [[-2.0, 0.5, 0.0, 0.0, 0.0, 0.0], 
                   [2.0, 1.0, 0.0, 0.0, 0.0, 0.0]]

    agent_configs = [ (init_states[i], dyn, ctrl.K) for i in 1:NA ] # Initialize agent configs for the workers

    init_distributed_agents!(workers_pids, agent_configs, h, 0.02) # Initialize workers

    prob = LayeredControlProblem(
        sheaf=sheaf,
        target_nodes=collect((NA+1):(NA+NT)),
        target_trajectory_func=target_trajectory_func,
        agent_configs=agent_configs,
        dt=h,
        steps=steps,
        pos_dim=D  # pass both y and z to the physical controller (no homogeneous coordinate)
    )

    result = run_layered_simulation(prob, workers_pids; mode=mode, nx=nx)
    return result
end


# ## Execution and Visualization

# Run distributed simulation scenarios
res_fixed = run_scenario(:fixed, :distributed)
res_bob   = run_scenario(:bobbing, :distributed)

# Compare distributed against centralized simulation to verify precision
res_fixed_c = run_scenario(:fixed, :centralised)
divergence = maximum(abs.(res_fixed.sim_data .- res_fixed_c.sim_data))
@printf("Max divergence between centralized and distributed simulation: %.3e\n", divergence)

# ## Animated Multi-panel Visualization
#
# We animate each scenario so the bobbing motion is clearly visible.
# For each frame k the three panels show:
#   1. y-z position plane with agent trails, target positions, and harmonic extensions
#   2. Tilt angle θ for both agents vs time (up to current frame)
#   3. Per-agent tracking error ||x_i[1:2] - q*_i|| vs time
#
# The visual style matches the escort example: steelblue agents, purple harmonic
# extension squares (alpha 0.3), grey star targets.

animate_scenario5(res_fixed;
    frame_step = 2,
    filename   = "scenario5_fixed.gif",
    fps        = 15,
    label_suffix = "fixed targets")
nothing # hide

animate_scenario5(res_bob;
    frame_step = 2,
    filename   = "scenario5_bobbing.gif",
    fps        = 15,
    label_suffix = "bobbing targets")
nothing # hide

# ![Scenario 5 [fixed targets]](scenario5_fixed.gif)

# ![Scenario 5 [bobbing targets]](scenario5_bobbing.gif)

# Clean up worker processes
rmprocs(workers_pids)
