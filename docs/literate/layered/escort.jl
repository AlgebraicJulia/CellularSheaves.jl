# # Layered Control Architecture for 12-Agent Escort Formation ($SE(3)$ Homogeneous Sheaf)
# 
# This example demonstrates scaling the layered control architecture to a 12-agent formation 
# escorting two moving targets in 3D space using an **$SE(3)$ Homogeneous Affine Cellular Sheaf**.
# 
# ## Non-Trivial Coordination Sheaf & Local Frame Rotations
# 
# Unlike identity sheaves ($I_D$) which decouple into independent scalar graph Laplacians ($L_{\mathcal{G}} \otimes I_D$), 
# this coordination sheaf operates on **4D homogeneous spatial cochains** $\tilde{\mathbf{x}} = [x, y, z, 1]^\top \in \mathbb{R}^4$. 
# 
# Each consensus edge enforces relative frame rotation $R_z(\theta_v) \in SO(3)$ ($\theta_v = \frac{2\pi(v-1)}{6}$) 
# and radial distance offset $\mathbf{d}_v = R_z(\theta_v) [r, 0, 0]^\top$ ($r = 1.2\,\text{m}$). 
# By pinning leader Agent 1 to Target 1 and leader Agent 7 to Target 2, consensus propagation across 
# the sheaf Laplacian produces **two perfect 3D spatial regular hexagonal escort rings** centered on the moving targets.
# 
# ## Theoretical Quadrotor Model & Local LQR Control
# 
# Low-level tracking is managed by 10D Discrete LQR controllers running on distributed Julia worker processes. 
# The physical quadrotor dynamics use the standard small-angle 10D linearization around hover 
# under fixed heading ($\psi = 0$), as formalized in seminal quadrotor control literature:
# 
# 1. **Mellinger, D., & Kumar, V. (2011)**. *Minimum snap trajectory generation and control for quadrotors*. *IEEE ICRA*, pp. 2520–2525.
# 2. **Bouabdallah, S., Noth, A., & Siegwart, R. (2004)**. *PID vs LQ control techniques applied to an indoor micro quadrotor*. *IEEE/RSJ IROS*, Vol. 3, pp. 2451–2456.


using CellularSheaves
import CellularSheaves.NetworkSheaves.EuclideanSheaves: _harmonic_extension_restricted_laplacian
using CellularSheaves.ControlSheaves.Tikhonov
using CellularSheaves.TrajectorySheaves: continuous_to_discrete_zoh
using CellularSheaves.NetworkSheaves.DistributedSolve
using CliqueTrees.Multifrontal
using LinearAlgebra
using SparseArrays
using Distributed
using Plots
using Printf

# A single house style for every figure below
default(framestyle = :box, grid = true, gridalpha = 0.18, gridstyle = :dot,
    titlefontsize = 10, guidefontsize = 9, legendfontsize = 8, tickfontsize = 8,
    markerstrokewidth = 0, size = (720, 380))

const RING_A, RING_B = :steelblue, :darkorange
const TARGET_A, TARGET_B = :gray35, :black

# ## Setup 10D Quadrotor Dynamics & DARE Solver

g = 9.81
m = 0.5
Ixx = 0.01
Iyy = 0.01

# Continuous-time state-space matrices (10D)
# State x = [x, y, z, phi, theta, x_dot, y_dot, z_dot, phi_dot, theta_dot]
Ac = zeros(10, 10)
Ac[1, 6] = 1.0
Ac[2, 7] = 1.0
Ac[3, 8] = 1.0
Ac[4, 9] = 1.0
Ac[5, 10] = 1.0
Ac[6, 5] = g
Ac[7, 4] = -g

Bc = zeros(10, 3)
Bc[8, 1] = 1.0 / m       # net thrust deviation
Bc[9, 2] = 1.0 / Ixx     # roll moment
Bc[10, 3] = 1.0 / Iyy    # pitch moment

DT = 0.05
Ad, Bd = continuous_to_discrete_zoh(Ac, Bc, DT)

nx = size(Ad, 1)
nu = size(Bd, 2)

# Compute Optimal LQR Gain via Discrete Algebraic Riccati Equation (DARE)
Q_diag = [150.0, 150.0, 150.0, 50.0, 50.0, 30.0, 30.0, 30.0, 1.0, 1.0]
Q_lqr = Matrix(Diagonal(Q_diag))
R_lqr = Matrix(Diagonal([0.01, 0.01, 0.01]))

function solve_dare(A, B, Q, R)
    P = Q
    for i in 1:200
        P_next = A' * P * A - (A' * P * B) * ((R + B' * P * B) \ (B' * P * A)) + Q
        if norm(P_next - P) < 1e-6
            break
        end
        P = P_next
    end
    return (R + B' * P * B) \ (B' * P * A)
end

K_lqr = solve_dare(Ad, Bd, Q_lqr, R_lqr)

A_cl = Ad - Bd * K_lqr
rho = maximum(abs.(eigvals(A_cl)))
@printf("Closed-loop spectral radius: %.4f\n", rho)

# ## SE(3) Homogeneous Coordination Sheaf Construction (D = 4)

# 12 agents in two escort rings of six, joined by a bridge edge (1, 7).
# Stalk dimension D = 4 for 3D spatial position + 1 homogeneous scale [x, y, z, 1].
const NA, NT = 12, 2
const TV1, TV2 = NA + 1, NA + 2
const D = 4
const I3 = Matrix{Float64}(I, 3, 3)

sheaf = EuclideanSheaf{Float64}(fill(D, NA + NT))

# SE(3) Rotation Matrix around Z-axis
Rz(theta) = [cos(theta) -sin(theta) 0.0;
             sin(theta)  cos(theta) 0.0;
             0.0         0.0        1.0]

r_ring = 1.2

consensus_edges = [(1,2),(2,3),(3,4),(4,5),(5,6),(6,1),
                   (7,8),(8,9),(9,10),(10,11),(11,12),(12,7)]

for (i, j) in consensus_edges
    angle_i = i <= 6 ? (i - 1) * 2π / 6 : (i - 7) * 2π / 6
    angle_j = j <= 6 ? (j - 1) * 2π / 6 : (j - 7) * 2π / 6
    Fi = [Rz(angle_i) zeros(3); 0 0 0 1]
    Fj = [Rz(angle_j) zeros(3); 0 0 0 1]
    add_sheaf_edge!(sheaf, i, j, Fi, Fj)
end

# Bridge edge (1, 7) connecting Ring A and Ring B
F1_bridge = [Rz(0.0) zeros(3); 0 0 0 1]
F7_bridge = [Rz(π) zeros(3); 0 0 0 1]
add_sheaf_edge!(sheaf, 1, 7, F1_bridge, F7_bridge)

# Pinning leader Agent 1 to Target 1 and leader Agent 7 to Target 2
# Consensus edges propagate the SE(3) ring geometry throughout the network
add_sheaf_edge!(sheaf, 1, TV1, [Rz(0.0) zeros(3); 0 0 0 1], [I3 -Rz(0.0)*[r_ring, 0, 0]; 0 0 0 1])
add_sheaf_edge!(sheaf, 7, TV2, [Rz(π) zeros(3); 0 0 0 1], [I3 -Rz(π)*[r_ring, 0, 0]; 0 0 0 1])

target1_pos(t) = [2cos(0.4t), 2sin(0.4t), 1.5 + 0.3sin(0.8t), 1.0]
target2_pos(t) = [2cos(-0.4t + pi), 2sin(-0.4t + pi), 1.5 + 0.3sin(-0.8t + pi), 1.0]

# Boundary conditions and factorisation
boundary0 = Dict(TV1 => target1_pos(0.0), TV2 => target2_pos(0.0))
_, _, Hraw, LIBraw = _harmonic_extension_restricted_laplacian(sheaf, boundary0)
H = Matrix(Hraw)
LIB = Matrix(LIBraw)

# Precompute chordal factorisation and tree partition across 12 agent processes
F = cholesky!(ChordalCholesky(sparse(H)), NoPivot())
Lfac = F.L
partition = partition_tree(Lfac, NA)

nchunk = length(partition.chunks)
@printf("Restricted Laplacian H size: %dx%d (rank %d)\n", size(H, 1), size(H, 2), rank(H))
@printf("Clique-tree nodes: %d supernodes, partitioned into %d chunks\n", length(partition.owner), nchunk)

# Provision worker processes for distributed parallel simulation (1 process per chunk)
workers_pids = addprocs(nchunk; exeflags = "--project=$(Base.active_project())")

# ## Load Distributed Agent Component
#
# Load the worker-side flight computer implementation (`_escort_agent_worker.jl`)
# on process 1 and all worker processes into `Main` using `Main.include`.

worker_file = joinpath(pkgdir(CellularSheaves), "docs", "literate", "layered", "_escort_agent_worker.jl")

Main.include(worker_file)
@everywhere workers_pids Main.include($worker_file)

# ## Simulation Framework (2 Full Orbital Revolutions)

function run_escort_simulation(mode=:distributed)
    STEPS = 630
    epsilon = 0.02
    init_states = [zeros(10) for _ in 1:NA]

    for i in 1:nchunk
        remotecall_fetch(Main.init_worker_agent!, workers_pids[i], init_states[i], K_lqr, Ad, Bd, epsilon)
    end

    sim_data = zeros(STEPS, NA, nx)
    qstar_history = zeros(STEPS, NA, D)

    H_pinv = pinv(H)

    for t_idx in 1:STEPS
        t = t_idx * DT
        b_t = [target1_pos(t); target2_pos(t)]

        rhs = Vector(-LIB * b_t)
        if mode == :centralised
            qstar_full = H_pinv * rhs
        else
            rhs_p = Vector(F.P' \ rhs)
            y_sol = distributed_tree_solve(Lfac, rhs_p, nchunk; pids = workers_pids)
            qstar_full = F.P \ y_sol
        end

        qstar = [qstar_full[(i-1)*D+1:i*D] for i in 1:NA]
        for i in 1:NA
            qstar_history[t_idx, i, :] = qstar[i]
        end

        step_futures = [remotecall(Main.step_worker_agent!, workers_pids[i], qstar[i], DT) for i in 1:nchunk]
        step_results = fetch.(step_futures)

        for i in 1:nchunk
            x_act, x_ref = step_results[i]
            sim_data[t_idx, i, :] = x_act
        end
    end

    return sim_data, qstar_history
end

# Run distributed and centralized simulations over 2 full revolutions
sim_d, q_d = run_escort_simulation(:distributed)
sim_c, q_c = run_escort_simulation(:centralised)

divergence = maximum(abs.(sim_d .- sim_c))
@printf("Max divergence between centralized and distributed 12-agent simulation over 2 revolutions: %.3e\n", divergence)

# Clean up worker processes
rmprocs(workers_pids)

# ## Multi-Projection Trajectory & Attitude Dynamics Visualization
#
# To demonstrate the physical quadrotor maneuvers across two full orbital revolutions (31.5 s), 
# we plot four complementary projections:
# 1. **Horizontal World Plane (x-y)**: Top-down view of the 12 quadrotors escorting targets T1 and T2 in SE(3) spatial rings.
# 2. **Target-Centered Relative Frame (x_rel - y_rel)**: Relative position of Ring A quadrotors centered on Target 1, demonstrating perfect 1.2 m regular hexagonal ring geometry.
# 3. **Roll Attitude Dynamics ϕ(t)**: Roll angle tilt (in degrees) driving lateral y-acceleration (y_ddot ≈ -g*ϕ).
# 4. **Pitch Attitude Dynamics θ(t)**: Pitch angle tilt (in degrees) driving forward x-acceleration (x_ddot ≈ g*θ).

STEPS = 630
orbit_t = range(0, STEPS * DT; length = 300)
lims_xy = (-3.5, 3.5)
lims_rel = (-1.8, 1.8)
trail = 40
ts = (1:STEPS) .* DT

anim = @animate for k in 1:6:STEPS
    t_curr = k * DT
    lo = max(1, k - trail)
    
    # 1. World Top-Down View
    p1 = plot(; aspect_ratio = 1, xlims = lims_xy, ylims = lims_xy,
              xlabel = "x position (m)", ylabel = "y position (m)",
              title = "World Top-Down View (x-y Plane)", legend = false)
    plot!(p1, [target1_pos(t)[1] for t in orbit_t], [target1_pos(t)[2] for t in orbit_t];
          color = :gray80, linewidth = 1, linestyle = :dot)
    plot!(p1, [target2_pos(t)[1] for t in orbit_t], [target2_pos(t)[2] for t in orbit_t];
          color = :gray80, linewidth = 1, linestyle = :dot)
    for i in 1:NA
        col = i <= 6 ? RING_A : RING_B
        sty = i <= 6 ? :solid : :dash
        plot!(p1, sim_d[lo:k, i, 1], sim_d[lo:k, i, 2];
              seriestype = :path, marker = :circle, markersize = 3, alpha = 0.6,
              linewidth = 1.4, color = col, linestyle = sty)
    end
    scatter!(p1, [target1_pos(t_curr)[1]], [target1_pos(t_curr)[2]]; marker = :star5, markersize = 8, color = TARGET_A)
    scatter!(p1, [target2_pos(t_curr)[1]], [target2_pos(t_curr)[2]]; marker = :star5, markersize = 8, color = TARGET_B)

    # 2. Target 1-Centered Relative Escort Ring Frame
    p2 = plot(; aspect_ratio = 1, xlims = lims_rel, ylims = lims_rel,
              xlabel = "rel x to T1 (m)", ylabel = "rel y to T1 (m)",
              title = "Target 1-Centered Escort Ring A (1.2m)", legend = false)
    # Draw reference 1.2m circle
    circ_ang = range(0, 2π; length = 100)
    plot!(p2, 1.2 .* cos.(circ_ang), 1.2 .* sin.(circ_ang); color = :gray80, linestyle = :dash, linewidth = 1)
    scatter!(p2, [0.0], [0.0]; marker = :star5, markersize = 10, color = TARGET_A)
    for i in 1:6
        rel_x = sim_d[k, i, 1] - target1_pos(t_curr)[1]
        rel_y = sim_d[k, i, 2] - target1_pos(t_curr)[2]
        scatter!(p2, [rel_x], [rel_y]; marker = :circle, markersize = 6, color = RING_A)
    end

    # 3. Roll Tilt Dynamics
    p3 = plot(; xlabel = "time (s)", ylabel = "roll angle ϕ (deg)",
              title = "Roll Tilt Dynamics ϕ(t) [y-accel]", legend = false, xlims = (0, 31.5), ylims = (-10.0, 10.0))
    for i in 1:NA
        col = i <= 6 ? RING_A : RING_B
        sty = i <= 6 ? :solid : :dash
        plot!(p3, ts[1:k], rad2deg.(sim_d[1:k, i, 4]); linewidth = 1.2, color = col, linestyle = sty)
    end

    # 4. Pitch Tilt Dynamics
    p4 = plot(; xlabel = "time (s)", ylabel = "pitch angle θ (deg)",
              title = "Pitch Tilt Dynamics θ(t) [x-accel]", legend = false, xlims = (0, 31.5), ylims = (-10.0, 10.0))
    for i in 1:NA
        col = i <= 6 ? RING_A : RING_B
        sty = i <= 6 ? :solid : :dash
        plot!(p4, ts[1:k], rad2deg.(sim_d[1:k, i, 5]); linewidth = 1.2, color = col, linestyle = sty)
    end

    plot(p1, p2, p3, p4; layout = (2, 2), size = (900, 700),
         plot_title = @sprintf("12-Agent SE(3) Escort Formation (2 Revolutions, t = %.2f s)", t_curr))
end
gif(anim, "layered_escort_tracking.gif"; fps = 15)
nothing # hide

# ![12-Agent Escort Tracking Projections](layered_escort_tracking.gif)
