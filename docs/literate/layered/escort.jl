# # Layered Control Architecture for 6-Agent Escort Formation ($SE(3)$ Homogeneous Sheaf)
# 
# This example demonstrates scaling the layered control architecture to a 6-agent formation 
# escorting a single stationary target in 3D space using an **$SE(3)$ Homogeneous Affine Cellular Sheaf**.
# 
# ## Non-Trivial Coordination Sheaf & Local Frame Rotations
# 
# Unlike identity sheaves ($I_D$) which decouple into independent scalar graph Laplacians ($L_{\mathcal{G}} \otimes I_D$), 
# this coordination sheaf operates on **4D homogeneous spatial cochains** $\tilde{\mathbf{x}} = [x, y, z, 1]^\top \in \mathbb{R}^4$. 
# 
# Each consensus edge enforces relative frame rotation $R_z(\theta_v) \in SO(3)$ ($\theta_v = \frac{2\pi(v-1)}{6}$) 
# and radial distance offset $\mathbf{d}_v = R_z(\theta_v) [r, 0, 0]^\top$ ($r = 1.2\,\text{m}$). 
# By pinning leader Agent 1 to the stationary target, consensus propagation across 
# the sheaf Laplacian produces a **perfect 3D spatial regular hexagonal escort ring** centered on the target.
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

const RING_COLOR = :steelblue
const TARGET_COLOR = :black

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

# 6 agents in a single escort ring.
# Stalk dimension D = 4 for 3D spatial position + 1 homogeneous scale [x, y, z, 1].
const NA = 6
const NT = 1
const TV1 = NA + 1
const D = 4
const I3 = Matrix{Float64}(I, 3, 3)

sheaf = EuclideanSheaf{Float64}(fill(D, NA + NT))

# SE(3) Rotation Matrix around Z-axis
Rz(theta) = [cos(theta) -sin(theta) 0.0;
             sin(theta)  cos(theta) 0.0;
             0.0         0.0        1.0]

r_ring = 1.2

consensus_edges = [(1,2),(2,3),(3,4),(4,5),(5,6),(6,1)]

for (i, j) in consensus_edges
    angle_i = (i - 1) * 2π / 6
    angle_j = (j - 1) * 2π / 6
    Fi = [Rz(angle_i) zeros(3); 0 0 0 1]
    Fj = [Rz(angle_j) zeros(3); 0 0 0 1]
    add_sheaf_edge!(sheaf, i, j, Fi, Fj)
end

# Pinning leader Agent 1 to Target 1
# Consensus edges propagate the SE(3) ring geometry throughout the network
add_sheaf_edge!(sheaf, 1, TV1, [Rz(0.0) zeros(3); 0 0 0 1], [I3 -Rz(0.0)*[r_ring, 0, 0]; 0 0 0 1])

target1_pos = [0.0, 0.0, 1.5, 1.0]

# Boundary conditions and factorisation
boundary0 = Dict(TV1 => target1_pos)
_, _, Hraw, LIBraw = _harmonic_extension_restricted_laplacian(sheaf, boundary0)
H = Matrix(Hraw)
LIB = Matrix(LIBraw)

# Precompute chordal factorisation and tree partition across agent processes
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

# ## Simulation Framework

function run_escort_simulation(mode=:distributed)
    STEPS = 200
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
        b_t = target1_pos

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

# Run distributed and centralized simulations
sim_d, q_d = run_escort_simulation(:distributed)
sim_c, q_c = run_escort_simulation(:centralised)

divergence = maximum(abs.(sim_d .- sim_c))
@printf("Max divergence between centralized and distributed 6-agent simulation: %.3e\n", divergence)

# Clean up worker processes
rmprocs(workers_pids)

# ## Multi-Projection Trajectory & Attitude Dynamics Visualization
#
# To demonstrate the physical quadrotor maneuvers around the stationary target, 
# we plot four complementary projections:
# 1. **Horizontal World Plane (x-y)**: Top-down view showing the 6 quadrotors forming the regular hexagonal escort ring.
# 2. **Altitude Profile (z vs t)**: Altitude convergence of all agents from the origin to the target z = 1.5 m.
# 3. **Roll Attitude Dynamics ϕ(t)**: Roll angle tilt (in degrees) driving lateral y-acceleration.
# 4. **Pitch Attitude Dynamics θ(t)**: Pitch angle tilt (in degrees) driving forward x-acceleration.

STEPS = 200
lims_xy = (-1.8, 1.8)
lims_z = (0.0, 1.8)
ts = (1:STEPS) .* DT

anim = @animate for k in 1:2:STEPS
    t_curr = k * DT
    
    p1 = plot(; aspect_ratio = 1, xlims = lims_xy, ylims = lims_xy,
              xlabel = "x position (m)", ylabel = "y position (m)",
              title = "World Top-Down View (x-y Plane)", legend = false)
    # Reference 1.2m circle
    circ_ang = range(0, 2π; length = 100)
    plot!(p1, r_ring .* cos.(circ_ang), r_ring .* sin.(circ_ang); color = :gray80, linestyle = :dash, linewidth = 1)
    scatter!(p1, [0.0], [0.0]; marker = :star5, markersize = 10, color = TARGET_COLOR)
    for i in 1:NA
        plot!(p1, sim_d[1:k, i, 1], sim_d[1:k, i, 2];
              seriestype = :path, marker = :circle, markersize = 3, alpha = 0.6,
              linewidth = 1.4, color = RING_COLOR)
    end

    p2 = plot(; xlims = (0, 10.0), ylims = lims_z,
              xlabel = "time (s)", ylabel = "altitude z (m)",
              title = "Altitude Profile z(t)", legend = false)
    for i in 1:NA
        plot!(p2, ts[1:k], sim_d[1:k, i, 3]; linewidth = 1.4, color = RING_COLOR)
    end
    plot!(p2, [0, 10.0], [1.5, 1.5]; color = :gray80, linestyle = :dot, linewidth = 1.2)

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
         plot_title = @sprintf("6-Agent SE(3) Stationary Escort Ring (t = %.2f s)", t_curr))
end
gif(anim, "layered_escort_tracking.gif"; fps = 15)
nothing # hide

# ![6-Agent Escort Tracking Projections](layered_escort_tracking.gif)
