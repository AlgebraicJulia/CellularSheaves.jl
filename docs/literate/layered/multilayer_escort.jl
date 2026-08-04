# # Multilayer Hierarchical Control for Team-of-Teams Escort Formation

# This example implements a multi-ring formation coordination architecture using cellular sheaves.
#
# ## Scenario
# We have 15 quadrotor agents organized into three rings:
# - **Ring 1 (6 agents):** Escorting Target 1.
# - **Ring 2 (6 agents):** Escorting Target 2.
# - **Ring 3 (3 agents):** Providing communication support, tracking the midpoint between Ring 1 and Ring 2.

using CellularSheaves
using CellularSheaves.Formations
using CellularSheaves.AgentControllers
using CellularSheaves.DistributedLayeredControl
using CellularSheaves.Pushforwards
using CellularSheaves.GraphHomomorphisms
import CellularSheaves.NetworkSheaves.EuclideanSheaves: _harmonic_extension_restricted_laplacian
using LinearAlgebra
using Statistics
using Plots
using Printf
using Graphs

# --- Setup ---
default(framestyle=:box, grid=true, gridalpha=0.18, gridstyle=:dot,
    titlefontsize=10, guidefontsize=9, legendfontsize=8, tickfontsize=8,
    markerstrokewidth=0, size=(800, 400))

# --- Sheaf Topology ---
const NA1 = 6
const NA2 = 6
const NA3 = 3
const NA = NA1 + NA2 + NA3
const T1_NODE = 16
const T2_NODE = 17
const TOTAL_NODES = 17
const r_ring = 0.3

function build_multilayer_sheaf(na1, na2, na3, r1, r2, r3)
    F = EuclideanSheaf{Float64}(fill(4, TOTAL_NODES))

    # Ring 1 (agents 1..6) around Target 1 (node 16)
    ring1 = build_escort_ring(na1, na1 + 1, r1; observers=[1])
    for e in edges(ring1.underlying_graph)
        u_local, v_local = src(e), dst(e)
        u = u_local == na1 + 1 ? T1_NODE : u_local
        v = v_local == na1 + 1 ? T1_NODE : v_local
        add_sheaf_edge!(F, u, v, get_restriction_map(ring1, u_local, v_local), get_restriction_map(ring1, v_local, u_local))
    end

    # Ring 2 (agents 7..12) around Target 2 (node 17)
    ring2 = build_escort_ring(na2, na2 + 1, r2; observers=[1])
    for e in edges(ring2.underlying_graph)
        u_local, v_local = src(e), dst(e)
        u = u_local == na2 + 1 ? T2_NODE : u_local + na1
        v = v_local == na2 + 1 ? T2_NODE : v_local + na1
        add_sheaf_edge!(F, u, v, get_restriction_map(ring2, u_local, v_local), get_restriction_map(ring2, v_local, u_local))
    end

    # Ring 3 (agents 13..15): Consensus edges with identity restriction maps
    # Under the pushforward sheaf f_★F, Fiber 3 coarsens these 3 agents into a 4D consensus stalk.
    # Lifting via T⁺ (pinv(T)) projects references back into the 0-cochain subspace where all 3 agents
    # remain in exact consensus at the midpoint.
    r3_agents = (na1 + na2 + 1):(na1 + na2 + na3)
    for i in 1:na3
        u = r3_agents[i]
        v = r3_agents[i % na3 + 1]
        add_sheaf_edge!(F, u, v, Matrix{Float64}(I, 4, 4), Matrix{Float64}(I, 4, 4))
    end

    # Support tracking edges: Ring 3 agents observe Target 1 and Target 2 directly
    add_sheaf_edge!(F, 13, T1_NODE, Matrix{Float64}(I, 4, 4), Matrix{Float64}(I, 4, 4))
    add_sheaf_edge!(F, 14, T2_NODE, Matrix{Float64}(I, 4, 4), Matrix{Float64}(I, 4, 4))

    return F
end

F = build_multilayer_sheaf(NA1, NA2, NA3, r_ring, r_ring * 1.2, r_ring * 0.8)

# Target trajectories
target1_pos(t) = [1.0 + 0.05cos(0.5 * t), 0.05sin(0.5 * t), 1.5 + 0.01sin(1.0 * t), 1.0]
target2_pos(t) = [-1.0 - 0.05cos(0.5 * t), -0.05sin(0.5 * t), 1.5 + 0.01cos(1.0 * t), 1.0]

# --- 2. Define High-Level Graph H & Graph Homomorphism f ---
# Coarse graph H has 3 nodes: Ring 1 (1), Ring 2 (2), Support Ring 3 (3)
function build_homomorphism(na1, na2, na3)
    v_map = zeros(Int, TOTAL_NODES)
    v_map[1:na1] .= 1
    v_map[T1_NODE] = 1
    v_map[(na1+1):(na1+na2)] .= 2
    v_map[T2_NODE] = 2
    v_map[(na1+na2+1):(na1+na2+na3)] .= 3

    return GraphHomomorphism(v_map, 3)
end

f = build_homomorphism(NA1, NA2, NA3)

# --- 3. Construct Pushforward Sheaf PfF and Fiber Bases ---
PfF = pushforward_sheaf(f, F)
fiber_bases = all_fiber_bases(f, F)

# --- 4. 3-Layer Hierarchical Solvers ---
# Target 1 is at node 16 (Fiber 1, local rows 25:28 in fiber_bases[1])
# Target 2 is at node 17 (Fiber 2, local rows 25:28 in fiber_bases[2])
B1_T1 = fiber_bases[1][25:28, :]
B2_T2 = fiber_bases[2][25:28, :]

function solve_high_level_harmonic(PfF, p1, p2, B1_T1, B2_T2)
    # 1. Convert world target vectors p1, p2 into PfF stalk basis coordinates
    q_pf_1 = B1_T1 \ p1
    q_pf_2 = B2_T2 \ p2

    # 2. Solve harmonic extension on coarse pushforward sheaf PfF
    boundary = Dict(1 => q_pf_1, 2 => q_pf_2)
    _, _, H_mat, B_mat = _harmonic_extension_restricted_laplacian(PfF, boundary)
    q_pf_3 = vec(H_mat \ (-Matrix(B_mat) * [q_pf_1; q_pf_2]))
    
    return [q_pf_1, q_pf_2, q_pf_3]
end

function solve_mid_level_harmonic(q_H, fiber_bases)
    # Lift high-level stalk coordinates back to world coordinates using fiber_bases
    q_G_1 = fiber_bases[1] * q_H[1]  # Fiber 1 (Agents 1..6 + T1)
    q_G_2 = fiber_bases[2] * q_H[2]  # Fiber 2 (Agents 7..12 + T2)
    q_G_3 = fiber_bases[3] * q_H[3]  # Fiber 3 (Agents 13..15)

    # Extract 4D SE(3) positions for agents 1..15
    q_agents = zeros(15, 4)
    q_agents[1:6, :] .= reshape(q_G_1[1:24], 4, 6)'
    q_agents[7:12, :] .= reshape(q_G_2[1:24], 4, 6)'
    q_agents[13:15, :] .= reshape(q_G_3[1:12], 4, 3)'
    
    return q_agents
end

# --- Simulation Parameters ---
DT = 0.05
STEPS = 200
T_END = STEPS * DT
time_grid = 0:DT:T_END

dyns = [QuadrotorDynamics(m=0.5 + 0.01 * i, Ixx=0.01, Iyy=0.01) for i in 1:15]
init_states = [[0.0, 0.0, 1.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] for i in 1:15]

sim_data_direct = []
qstar_history_direct = []
target_history_direct = []

current_states_direct = copy(init_states)
for step in 1:STEPS
    t = time_grid[step]
    p1 = target1_pos(t)
    p2 = target2_pos(t)
    
    # 1. High-Level Planner: Solve harmonic extension on pushforward sheaf PfF
    q_H = solve_high_level_harmonic(PfF, p1, p2, B1_T1, B2_T2)
    
    # World reference targets for visualization (T1, T2, and Support center)
    p3_center = fiber_bases[3][1:4, :] * q_H[3]
    push!(target_history_direct, [p1, p2, p3_center])

    # 2. Mid-Level Planner: Lift high-level section to mid-level agent references q*
    qstar_agents = solve_mid_level_harmonic(q_H, fiber_bases)
    push!(qstar_history_direct, qstar_agents)

    step_states = []
    for i in 1:15
        x_actual = current_states_direct[i][1:3]
        x_ref = qstar_agents[i, 1:3]

        error = x_actual - x_ref
        current_states_direct[i][1:3] .-= 0.3 * error * DT
        push!(step_states, copy(current_states_direct[i]))
    end
    push!(sim_data_direct, step_states)
end

# --- Error Metrics ---
formation_errors_direct = zeros(STEPS, 3)
tracking_centroid_errors_direct = zeros(STEPS, 3)
teams = [1:6, 7:12, 13:15]
team_colors = [:steelblue :darkorange :darkseagreen]
expected_radii = [r_ring, r_ring * 1.2, r_ring * 0.8]

for step in 1:STEPS
    for t_idx in 1:3
        team = teams[t_idx]
        r_exp = expected_radii[t_idx]

        # Centroid of the team
        centroid = zeros(3)
        for i in team
            centroid .+= sim_data_direct[step][i][1:3]
        end
        centroid ./= length(team)

        # Formation error: Mean absolute deviation of agent distances from centroid relative to expected ring radius
        radii_devs = [abs(norm(sim_data_direct[step][i][1:3] - centroid) - r_exp) for i in team]
        formation_errors_direct[step, t_idx] = mean(radii_devs)

        # Tracking error: Centroid error relative to target
        target = target_history_direct[step][t_idx][1:3]
        tracking_centroid_errors_direct[step, t_idx] = norm(centroid - target)
    end
end

# --- Fixed Axis Limits Pre-computation ---
all_xs = vcat([[sim_data_direct[s][i][1] for i in 1:15] for s in 1:STEPS]...)
all_ys = vcat([[sim_data_direct[s][i][2] for i in 1:15] for s in 1:STEPS]...)
all_q_xs = vcat([qstar_history_direct[s][:, 1] for s in 1:STEPS]...)
all_q_ys = vcat([qstar_history_direct[s][:, 2] for s in 1:STEPS]...)
all_t1_xs = [target_history_direct[s][1][1] for s in 1:STEPS]
all_t1_ys = [target_history_direct[s][1][2] for s in 1:STEPS]
all_t2_xs = [target_history_direct[s][2][1] for s in 1:STEPS]
all_t2_ys = [target_history_direct[s][2][2] for s in 1:STEPS]

p1_min_x, p1_max_x = minimum(vcat(all_xs, all_q_xs, all_t1_xs, all_t2_xs)) - 0.5, maximum(vcat(all_xs, all_q_xs, all_t1_xs, all_t2_xs)) + 0.5
p1_min_y, p1_max_y = minimum(vcat(all_ys, all_q_ys, all_t1_ys, all_t2_ys)) - 0.5, maximum(vcat(all_ys, all_q_ys, all_t1_ys, all_t2_ys)) + 0.5
p1_xlims = (p1_min_x, p1_max_x)
p1_ylims = (p1_min_y, p1_max_y)

rel_lims = []
for (t_idx, team) in enumerate(teams)
    rel_xs_all = vcat([[sim_data_direct[s][i][1] - target_history_direct[s][t_idx][1] for i in team] for s in 1:STEPS]...)
    rel_ys_all = vcat([[sim_data_direct[s][i][2] - target_history_direct[s][t_idx][2] for i in team] for s in 1:STEPS]...)
    rel_q_xs_all = vcat([[qstar_history_direct[s][i, 1] - target_history_direct[s][t_idx][1] for i in team] for s in 1:STEPS]...)
    rel_q_ys_all = vcat([[qstar_history_direct[s][i, 2] - target_history_direct[s][t_idx][2] for i in team] for s in 1:STEPS]...)

    rx_min, rx_max = minimum(vcat(rel_xs_all, rel_q_xs_all, [0.0])) - 0.2, maximum(vcat(rel_xs_all, rel_q_xs_all, [0.0])) + 0.2
    ry_min, ry_max = minimum(vcat(rel_ys_all, rel_q_ys_all, [0.0])) - 0.2, maximum(vcat(rel_ys_all, rel_q_ys_all, [0.0])) + 0.2
    r_max = max(abs(rx_min), abs(rx_max), abs(ry_min), abs(ry_max))
    push!(rel_lims, (-r_max, r_max))
end

safe_form_err_full = max.(formation_errors_direct, 1e-6)
safe_track_err_full = max.(tracking_centroid_errors_direct, 1e-6)
p5_ylims = (10^floor(log10(minimum(safe_form_err_full))), 10^ceil(log10(maximum(safe_form_err_full))))
p6_ylims = (10^floor(log10(minimum(safe_track_err_full))), 10^ceil(log10(maximum(safe_track_err_full))))

# --- 6-Panel Visualization ---
proj2d(v) = v[1:2]

function make_comprehensive_frame(step)
    t = time_grid[step]
    states = sim_data_direct[step]
    qstar = qstar_history_direct[step]
    qH = target_history_direct[step]

    # Panel 1: Top-Down View with fixed limits
    p1 = plot(aspect_ratio=:equal, title="Multilayer Escort Top-Down View (t=$(round(t, digits=2))s)",
        xlabel="x [m]", ylabel="y [m]", legend=false, xlims=p1_xlims, ylims=p1_ylims)

    # Communication edges
    for edge in edges(F.underlying_graph)
        u, v = edge.src, edge.dst
        if u <= 15 && v <= 15
            x_u, x_v = proj2d(states[u]), proj2d(states[v])
            plot!(p1, [x_u[1], x_v[1]], [x_u[2], x_v[2]], color=:black, alpha=0.4, lw=1.5)
        end
    end

    # Targets 1 and 2
    scatter!(p1, [qH[1][1]], [qH[1][2]], marker=:star, color=:black, markersize=6)
    scatter!(p1, [qH[2][1]], [qH[2][2]], marker=:star, color=:black, markersize=6)

    # Harmonic extension positions
    for i in 1:15
        q_s = proj2d(qstar[i, :])
        scatter!(p1, [q_s[1]], [q_s[2]], marker=:square, color=:purple, alpha=0.3, markersize=4)
    end

    # Agents
    for (t_idx, team) in enumerate(teams)
        xs = [proj2d(states[i])[1] for i in team]
        ys = [proj2d(states[i])[2] for i in team]
        scatter!(p1, xs, ys, color=team_colors[t_idx], label="Team $(t_idx)", markersize=5)
    end

    # Relative Views with fixed limits
    p_rel = []
    for (t_idx, team) in enumerate(teams)
        target_pos = proj2d(qH[t_idx])
        p_target_view = plot(aspect_ratio=:equal, title="Team $(t_idx) Relative View",
            xlabel="rel x [m]", ylabel="rel y [m]", legend=false,
            xlims=rel_lims[t_idx], ylims=rel_lims[t_idx])
        scatter!(p_target_view, [0], [0], marker=:star, color=:black, markersize=6)

        rel_xs = [proj2d(states[i])[1] - target_pos[1] for i in team]
        rel_ys = [proj2d(states[i])[2] - target_pos[2] for i in team]
        scatter!(p_target_view, rel_xs, rel_ys, color=team_colors[t_idx], markersize=5)

        rel_q_xs = [proj2d(qstar[i, :])[1] - target_pos[1] for i in team]
        rel_q_ys = [proj2d(qstar[i, :])[2] - target_pos[2] for i in team]
        scatter!(p_target_view, rel_q_xs, rel_q_ys, marker=:square, color=:purple, alpha=0.3, markersize=4)

        # Communication edges
        for edge in edges(F.underlying_graph)
            u, v = edge.src, edge.dst
            if u in team && v in team
                x_u, x_v = proj2d(states[u]), proj2d(states[v])
                plot!(p_target_view,
                    [x_u[1] - target_pos[1], x_v[1] - target_pos[1]],
                    [x_u[2] - target_pos[2], x_v[2] - target_pos[2]],
                    color=:black, alpha=0.4, lw=1.5)
            end
        end

        push!(p_rel, p_target_view)
    end

    safe_form_err = max.(formation_errors_direct[1:step, :], 1e-6)
    safe_track_err = max.(tracking_centroid_errors_direct[1:step, :], 1e-6)

    ## Panel 5: Formation Error with fixed limits & xlims
    p5 = plot(time_grid[1:step], safe_form_err, yscale=:log10,
        title="Formation Shape Error", xlabel="Time [s]", ylabel="Error [m]",
        color=team_colors, lw=3, xlims=(0, time_grid[end]), ylims=p5_ylims,
        label=["Team 1" "Team 2" "Team 3"])

    ## Panel 6: Tracking Error with fixed limits & xlims
    p6 = plot(time_grid[1:step], safe_track_err, yscale=:log10,
        title="Target Tracking Error", xlabel="Time [s]", ylabel="Centroid Error [m]",
        color=team_colors, lw=3, xlims=(0, time_grid[end]), ylims=p6_ylims,
        label=["Team 1" "Team 2" "Team 3"])

    return plot(p1, p_rel[1], p_rel[2], p_rel[3], p5, p6, layout=(3, 2), size=(1200, 1500))
end

p = make_comprehensive_frame(STEPS)
savefig(p, "multilayer_comprehensive_snapshot.png")

anim = @animate for step in 1:4:STEPS
    make_comprehensive_frame(step)
end
gif(anim, "multilayer_comprehensive_animation.gif"; fps=10)

println("Multilayer Escort simulation complete.")

