# # Multilayer Hierarchical Control for Team-of-Teams Escort Formation
#
# This example demonstrates multi-agent escort formation control comparing a **3-Layer Pushforward Control Architecture** ($f_\star F$)
# against a **Direct 2-Layer Solve** on the full sheaf $F$ using the generic `Layered` architecture.
#
# ## Mathematical & Categorical Overview
# - **Generic Topology Specification (`LayeredEscortSpec`):**
#   - **Ring 1 (6 agents):** Escorts Target 1 (node 16) at radius $r_1 = 0.3\,\text{m}$.
#   - **Ring 2 (6 agents):** Escorts Target 2 (node 17) at radius $r_2 = 0.36\,\text{m}$.
#   - **Support Pool (3 agents):** Bridges Ring 1 and Ring 2, observing Target 1 and Target 2 directly.
# - **Pushforward Sheaf $f_\star F$ (3 coarse nodes):**
#   - Coarsens Fiber 1 $\to$ Node 1, Fiber 2 $\to$ Node 2, Support Pool $\to$ Node 3.
#   - Fiber bases are structured and managed via `LayeredFiberBases`.
#   - Pseudoinverse lifting $T^+$ projects high-level midpoint references back onto the **subspace of exact subteam consensus**.

using CellularSheaves
using CellularSheaves.Formations
using CellularSheaves.AgentControllers
using CellularSheaves.ControlSheaves.Layered
using CellularSheaves.Pushforwards
using CellularSheaves.GraphHomomorphisms
using LinearAlgebra
using Statistics
using Plots
using Printf

# --- Topology Specification ---
const r_ring = 0.3

rings = [
    RingSpec(1, 6, r_ring),
    RingSpec(2, 6, r_ring * 1.2)
]
supports = [
    SupportSpec(1, 2, 3)
]

spec = LayeredEscortSpec(rings, supports)

# --- Sheaf Topology & Homomorphism ---
F = build_layered_escort_sheaf(spec)
f = build_layered_homomorphism(spec)

PfF = pushforward_sheaf(f, F)
bases = build_layered_fiber_bases(f, F, spec)

# --- Simulation Parameters ---
DT = 0.05
STEPS = 200

# Setup LQR gain for tests
dyn_test = QuadrotorDynamics()
Ad, Bd = CellularSheaves.AgentControllers.discrete_matrices(dyn_test, DT)
# Standard LQR Cost Matrices (moderate velocity gain)
Q_diag = [500.0, 500.0, 500.0, 150.0, 150.0, 100.0, 100.0, 100.0, 5.0, 5.0]
Q_lqr = Matrix(Diagonal(Q_diag))
R_lqr = Matrix(Diagonal([0.005, 0.005, 0.005]))

K_test = CellularSheaves.AgentControllers.solve_dare(Ad, Bd, 10*Q_lqr, R_lqr)

# Target position, velocity, and acceleration trajectories
target1_pos(t) = [1.0 + 0.5 * cos(0.5 * t), 0.5 * sin(0.5 * t), 1.5 + 0.01 * sin(1.0 * t), 1.0]
target2_pos(t) = [-1.0 - 0.5 * cos(0.5 * t), -0.5 * sin(0.5 * t), 1.5 + 0.01 * cos(1.0 * t), 1.0]

target1_vel(t) = [-0.25 * sin(0.5 * t), 0.25 * cos(0.5 * t), 0.01 * cos(1.0 * t), 0.0]
target2_vel(t) = [0.25 * sin(0.5 * t), -0.25 * cos(0.5 * t), -0.01 * sin(1.0 * t), 0.0]

target1_acc(t) = [-0.125 * cos(0.5 * t), -0.125 * sin(0.5 * t), -0.01 * sin(1.0 * t), 0.0]
target2_acc(t) = [0.125 * cos(0.5 * t), 0.125 * sin(0.5 * t), -0.01 * cos(1.0 * t), 0.0]

target_trajs = [target1_pos, target2_pos]
target_vels  = [target1_vel, target2_vel]
target_accels = [target1_acc, target2_acc]

prob_pf = LayeredEscortProblem(
    spec, F, f, PfF, bases,
    HomogeneousDynamics(dyn_test, K_test),
    target_trajs, target_vels, target_accels,
    DT, STEPS
)

#  ## Simulation 1: 3-Layer Pushforward Architecture (f_★F) ---
@time res_pf = run_layered_escort_simulation(prob_pf; use_feedforward=true)
@time res_pf = run_layered_escort_simulation(prob_pf; use_feedforward=true)

# ## Simulation 2: Direct 2-Layer Solve on Sheaf F ---

function direct_simulation()
    sim_data_dir = Vector{Vector{Vector{Float64}}}()
    qstar_history_dir = Vector{Matrix{Float64}}()
    target_history_dir = Vector{Vector{Vector{Float64}}}()

    time_grid = 0:DT:(STEPS*DT)
    ## Initialize agent states for direct simulation using AgentState and JointTikhonovFilter
    agent_states_dir = [AgentState(zeros(10), dyn_test, DT, K_test, 0.02; use_velocity=true) for _ in 1:spec.n_agents]
    for state in agent_states_dir
        state.x[3] = 1.5
    end

    for step in 1:STEPS
        t = time_grid[step]
        p1 = target1_pos(t)
        p2 = target2_pos(t)
        p3_center = (p1 + p2) ./ 2.0
        push!(target_history_dir, [p1, p2, p3_center])

        ## Direct Harmonic Solves on F for position, velocity, and acceleration
        qstar_agents = solve_direct_harmonic(F, spec.target_nodes, [p1, p2])
        qstar_dot_agents = solve_direct_harmonic(F, spec.target_nodes, [target1_vel(t), target2_vel(t)])
        qstar_ddot_agents = solve_direct_harmonic(F, spec.target_nodes, [target1_acc(t), target2_acc(t)])
        push!(qstar_history_dir, qstar_agents)

        ## Agent Dynamics Integration with Feedforward Velocity & Acceleration Support
        step_states = Vector{Vector{Float64}}()
        for i in 1:spec.n_agents
            q_ref_i = qstar_agents[i, 1:3]
            q_dot_i = qstar_dot_agents[i, 1:3]
            q_ddot_i = qstar_ddot_agents[i, 1:3]
            
            step_agent!(agent_states_dir[i], q_ref_i, q_dot_i, q_ddot_i, DT)
            push!(step_states, copy(agent_states_dir[i].x))
        end
        push!(sim_data_dir, step_states)
    end
    return (sim_data_dir, qstar_history_dir, target_history_dir, time_grid)
end

@time direct_simulation()
@time sim_data_dir, qstar_history_dir, target_history_dir, time_grid = direct_simulation()

# --- 3. Render Animations & Visualizations ---
r1, r2, r3 = r_ring, r_ring * 1.2, r_ring * 0.8

# Render 3-Layer Pushforward Animation
animate_comprehensive_escort(res_pf.sim_data, res_pf.qstar_history, res_pf.target_history, time_grid, r1, r2, r3, "multilayer_comprehensive_animation.gif"; frame_step=4, fps=10, title_prefix="3-Layer Pushforward Architecture (f_★F)")

# Render Direct 2-Layer Solve Animation
animate_comprehensive_escort(sim_data_dir, qstar_history_dir, target_history_dir, time_grid, r1, r2, r3, "direct_sheaf_comprehensive_animation.gif"; frame_step=4, fps=10, title_prefix="Direct 2-Layer Solve (F)")

println("Multilayer Escort simulation with generic architecture complete.")

# ## Visual Comparison Animations
#
# ![3-Layer Pushforward Architecture Animation](multilayer_comprehensive_animation.gif)
#
# ![Direct 2-Layer Solve Animation](direct_sheaf_comprehensive_animation.gif)
