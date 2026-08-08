# # Multilayer Hierarchical Control for Team-of-Teams Escort Formation
#
# This example demonstrates multi-agent escort formation control comparing a **3-Layer Pushforward Control Architecture** ($f_\star F$)
# against a **Direct 2-Layer Solve** on the full sheaf $F$.
#
# ## Mathematical & Categorical Overview
# - **Full Sheaf $F$ (17 nodes, 15 agents + 2 targets):**
#   - **Ring 1 (6 agents):** Escorts Target 1 (node 16) at radius $r_1 = 0.3\,\text{m}$.
#   - **Ring 2 (6 agents):** Escorts Target 2 (node 17) at radius $r_2 = 0.36\,\text{m}$.
#   - **Ring 3 (3 support agents):** Tracks the midpoint between Target 1 and Target 2 with identity consensus restriction maps ($I_{4\times 4}, I_{4\times 4}$).
# - **Pushforward Sheaf $f_\star F$ (3 coarse nodes):**
#   - Coarsens Fiber 1 $\to$ Node 1, Fiber 2 $\to$ Node 2, Fiber 3 $\to$ Node 3.
#   - Under $f_\star F$, Fiber 3 coarsens into a 4D consensus stalk.
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

# --- Simulation Parameters ---
const NA1 = 6
const NA2 = 6
const NA3 = 3
const NA = NA1 + NA2 + NA3
const r_ring = 0.3

DT = 0.05
STEPS = 200
T_END = STEPS * DT
time_grid = 0:DT:T_END

# Target trajectories
target1_pos(t) = [1.0 + 0.05 * cos(0.5 * t), 0.05 * sin(0.5 * t), 1.5 + 0.01 * sin(1.0 * t), 1.0]
target2_pos(t) = [-1.0 - 0.05 * cos(0.5 * t), -0.05 * sin(0.5 * t), 1.5 + 0.01 * cos(1.0 * t), 1.0]

# --- 1. Build Sheaf Topology & Homomorphism ---
F = build_multilayer_sheaf(NA1, NA2, NA3, r_ring, r_ring * 1.2, r_ring * 0.8)
f = build_homomorphism(NA1, NA2, NA3)

PfF = pushforward_sheaf(f, F)
fiber_bases = all_fiber_bases(f, F)

B1_T1 = fiber_bases[1][25:28, :]
B2_T2 = fiber_bases[2][25:28, :]

# --- 2. Simulation 1: 3-Layer Pushforward Architecture (f_★F) ---
dyns_pf = [QuadrotorDynamics(m=0.5 + 0.01 * i, Ixx=0.01, Iyy=0.01) for i in 1:NA]
init_states_pf = [[0.0, 0.0, 1.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] for _ in 1:NA]

sim_data_pf = []
qstar_history_pf = []
target_history_pf = []

current_states_pf = copy(init_states_pf)
for step in 1:STEPS
    t = time_grid[step]
    p1 = target1_pos(t)
    p2 = target2_pos(t)
    
    # 1. High-Level Planner: Solve harmonic extension on coarse pushforward sheaf PfF
    q_H = solve_high_level_harmonic(PfF, p1, p2, B1_T1, B2_T2)
    p3_center = fiber_bases[3][1:4, :] * q_H[3]
    push!(target_history_pf, [p1, p2, p3_center])

    # 2. Mid-Level Planner: Lift high-level section to agent space G using fiber bases
    qstar_agents = solve_mid_level_harmonic(q_H, fiber_bases)
    push!(qstar_history_pf, qstar_agents)

    # 3. Agent Dynamics Integration
    step_states = []
    for i in 1:NA
        x_actual = current_states_pf[i][1:3]
        x_ref = qstar_agents[i, 1:3]
        error = x_actual - x_ref
        current_states_pf[i][1:3] .-= 0.3 * error * DT
        push!(step_states, copy(current_states_pf[i]))
    end
    push!(sim_data_pf, step_states)
end

# --- 3. Simulation 2: Direct 2-Layer Solve on Sheaf F ---
dyns_dir = [QuadrotorDynamics(m=0.5 + 0.01 * i, Ixx=0.01, Iyy=0.01) for i in 1:NA]
init_states_dir = [[0.0, 0.0, 1.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] for _ in 1:NA]

sim_data_dir = []
qstar_history_dir = []
target_history_dir = []

current_states_dir = copy(init_states_dir)
for step in 1:STEPS
    t = time_grid[step]
    p1 = target1_pos(t)
    p2 = target2_pos(t)
    p3_center = (p1 + p2) ./ 2.0
    push!(target_history_dir, [p1, p2, p3_center])

    # Direct Harmonic Solve on F
    qstar_agents = solve_direct_harmonic(F, p1, p2)
    push!(qstar_history_dir, qstar_agents)

    # Agent Dynamics Integration
    step_states = []
    for i in 1:NA
        x_actual = current_states_dir[i][1:3]
        x_ref = qstar_agents[i, 1:3]
        error = x_actual - x_ref
        current_states_dir[i][1:3] .-= 0.3 * error * DT
        push!(step_states, copy(current_states_dir[i]))
    end
    push!(sim_data_dir, step_states)
end

# --- 4. Render Animations & Visualizations ---
r1, r2, r3 = r_ring, r_ring * 1.2, r_ring * 0.8

# Render 3-Layer Pushforward Animation
animate_comprehensive_escort(sim_data_pf, qstar_history_pf, target_history_pf, time_grid, r1, r2, r3, "multilayer_comprehensive_animation.gif"; frame_step=4, fps=10, title_prefix="3-Layer Pushforward Architecture (f_★F)")

# Render Direct 2-Layer Solve Animation
animate_comprehensive_escort(sim_data_dir, qstar_history_dir, target_history_dir, time_grid, r1, r2, r3, "direct_sheaf_comprehensive_animation.gif"; frame_step=4, fps=10, title_prefix="Direct 2-Layer Solve (F)")

println("Multilayer & Direct Escort comparison simulation complete.")
