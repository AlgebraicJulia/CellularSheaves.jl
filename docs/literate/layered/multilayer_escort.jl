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

# Target trajectories
target1_pos(t) = [1.0 + 0.05 * cos(0.5 * t), 0.05 * sin(0.5 * t), 1.5 + 0.01 * sin(1.0 * t), 1.0]
target2_pos(t) = [-1.0 - 0.05 * cos(0.5 * t), -0.05 * sin(0.5 * t), 1.5 + 0.01 * cos(1.0 * t), 1.0]

target_trajs = [target1_pos, target2_pos]

prob_pf = LayeredEscortProblem(
    spec, F, f, PfF, bases,
    HomogeneousDynamics(QuadrotorDynamics()),
    target_trajs, target_trajs, target_trajs,
    DT, STEPS
)

# --- 1. Simulation 1: 3-Layer Pushforward Architecture (f_★F) ---
res_pf = run_layered_escort_simulation(prob_pf)

# --- 2. Simulation 2: Direct 2-Layer Solve on Sheaf F ---
sim_data_dir = []
qstar_history_dir = []
target_history_dir = []

time_grid = 0:DT:(STEPS*DT)
init_states_dir = [[0.0, 0.0, 1.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] for _ in 1:spec.n_agents]
current_states_dir = copy(init_states_dir)

for step in 1:STEPS
    t = time_grid[step]
    p1 = target1_pos(t)
    p2 = target2_pos(t)
    p3_center = (p1 + p2) ./ 2.0
    push!(target_history_dir, [p1, p2, p3_center])

    ## Direct Harmonic Solve on F
    qstar_agents = solve_direct_harmonic(F, spec.target_nodes, [p1, p2])
    push!(qstar_history_dir, qstar_agents)

    ## Agent Dynamics Integration
    step_states = []
    for i in 1:spec.n_agents
        x_actual = current_states_dir[i][1:3]
        x_ref = qstar_agents[i, 1:3]
        error = x_actual - x_ref
        current_states_dir[i][1:3] .-= 0.3 * error * DT
        push!(step_states, copy(current_states_dir[i]))
    end
    push!(sim_data_dir, step_states)
end

# --- 3. Render Animations & Visualizations ---
r1, r2, r3 = r_ring, r_ring * 1.2, r_ring * 0.8

# Render 3-Layer Pushforward Animation
animate_comprehensive_escort(res_pf.sim_data, res_pf.qstar_history, res_pf.target_history, time_grid, r1, r2, r3, "multilayer_comprehensive_animation.gif"; frame_step=4, fps=10, title_prefix="3-Layer Pushforward Architecture (f_★F)")

# Render Direct 2-Layer Solve Animation
animate_comprehensive_escort(sim_data_dir, qstar_history_dir, target_history_dir, time_grid, r1, r2, r3, "direct_sheaf_comprehensive_animation.gif"; frame_step=4, fps=10, title_prefix="Direct 2-Layer Solve (F)")

println("Multilayer Escort simulation with generic architecture complete.")
