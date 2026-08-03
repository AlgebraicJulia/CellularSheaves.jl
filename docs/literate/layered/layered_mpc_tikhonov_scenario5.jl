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
    addprocs(NA - length(workers_pids))
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
        steps=steps
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

# ## Multi-panel Visualization
#
# For each scenario (fixed / bobbing) we build three panels:
#   1. y-z position plane with agent trajectories, targets, and harmonic extensions
#   2. Roll tilt angle θ (state index 3) for both agents vs time
#   3. Per-agent tracking error:  ||x_i[1:2] - q*_i|| over time
#
# The six panels are arranged in a 2×3 grid: rows = scenario, columns = panel type.

function build_scenario_panels(res, scenario_label)
    prob = res.problem
    sim  = res.sim_data        # [steps, NA, nx]
    qh   = res.qstar_history   # [steps, NA, D]
    ts   = (1:prob.steps) .* prob.dt
    NA   = 2                   # two agents
    NT   = 2                   # two targets
    T1, T2 = prob.target_nodes

    ## ---------- panel 1: y-z position plane ----------

    all_y = vcat(sim[:, :, 1]...)
    all_z = vcat(sim[:, :, 2]...)
    t1_y = [prob.target_trajectory_func(T1, t)[1] for t in ts]
    t1_z = [prob.target_trajectory_func(T1, t)[2] for t in ts]
    t2_y = [prob.target_trajectory_func(T2, t)[1] for t in ts]
    t2_z = [prob.target_trajectory_func(T2, t)[2] for t in ts]

    min_y = minimum(vcat(all_y, t1_y, t2_y))
    max_y = maximum(vcat(all_y, t1_y, t2_y))
    min_z = minimum(vcat(all_z, t1_z, t2_z))
    max_z = maximum(vcat(all_z, t1_z, t2_z))
    cy = (min_y + max_y) / 2;  cz = (min_z + max_z) / 2
    span = max(max_y - min_y, max_z - min_z) / 2 + 0.4

    pyz = plot(;
        xlabel = "y position [m]",
        ylabel = "z position [m]",
        title  = "y-z Plane [$scenario_label]",
        aspect_ratio = 1,
        xlims = (cy - span, cy + span),
        ylims = (cz - span, cz + span),
        legend = :topright,
    )
    plot!(pyz, t1_y, t1_z; ls=:dot, lw=1, color=:gray60, label="T1 orbit")
    plot!(pyz, t2_y, t2_z; ls=:dot, lw=1, color=:gray30, label="T2 orbit")
    scatter!(pyz, [t1_y[end]], [t1_z[end]]; marker=:star5, ms=10, color=:gray60, label="Target 1")
    scatter!(pyz, [t2_y[end]], [t2_z[end]]; marker=:star5, ms=10, color=:gray30, label="Target 2")
    for i in 1:NA
        plot!(pyz, sim[:, i, 1], sim[:, i, 2];
              lw=1.2, color=:steelblue, alpha=0.7, label=(i==1 ? "Agents" : false))
        scatter!(pyz, [sim[end, i, 1]], [sim[end, i, 2]];
                 marker=:circle, ms=7, color=:steelblue, label=false)
        scatter!(pyz, qh[:, i, 1], qh[:, i, 2];
                 marker=:square, ms=3, color=:purple, alpha=0.3, label=(i==1 ? "Harmonic ext." : false))
    end

    ## ---------- panel 2: tilt angle θ vs time ----------

    all_theta = sim[:, :, 3]
    max_th = maximum(abs.(all_theta))
    pad_th = max(max_th * 0.2, 0.05)

    pth = plot(;
        xlabel = "time [s]",
        ylabel = "tilt angle θ [rad]",
        title  = "Tilt Angle [$scenario_label]",
        ylims  = (-max_th - pad_th, max_th + pad_th),
        xlims  = (0, ts[end]),
        legend = :topright,
    )
    for i in 1:NA
        plot!(pth, ts, sim[:, i, 3]; lw=1.2, label="Agent $i")
    end

    ## ---------- panel 3: per-agent tracking error ||pos_i - q*_i|| ----------

    terr = [norm.(eachrow(sim[:, i, 1:2] .- qh[:, i, :])) for i in 1:NA]
    max_terr = maximum(vcat(terr...))

    pte = plot(;
        xlabel = "time [s]",
        ylabel = "position error [m]",
        title  = "Tracking Error [$scenario_label]",
        ylims  = (0, max_terr * 1.2 + 0.01),
        xlims  = (0, ts[end]),
        legend = :topright,
    )
    for i in 1:NA
        plot!(pte, ts, terr[i]; lw=1.2, label="Agent $i")
    end

    return pyz, pth, pte
end

p1, p2, p3 = build_scenario_panels(res_fixed, "fixed targets")
p4, p5, p6 = build_scenario_panels(res_bob,   "bobbing targets")

fig = plot(p1, p2, p3, p4, p5, p6;
    layout = (2, 3),
    size   = (1400, 900),
    plot_title = "Layered Control Scenario 5"
)
savefig(fig, "layered_control_combined.png")

# ![Layered Control Scenario 5](layered_control_combined.png)
