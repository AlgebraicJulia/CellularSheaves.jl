# # Layered Control Architecture for Multi-Quadrotor Tracking (Scenario 5)
# 
# This example demonstrates a layered control architecture designed for multi-agent target tracking 
# where high-level coordination is handled by a static cellular sheaf and low-level tracking is 
# managed by independent LQR controllers running on distributed Julia worker processes.
# 
# ## Architecture Overview
# 
# The control pipeline is structured as follows:
# 
# 1. **Coordination Sheaf**: Resolves the conflict between target tracking and agent consensus 
#    by computing a harmonic extension $\mathbf{q}^*$ that minimizes a global energy functional.
# 2. **Distributed Tree Solve**: Solves for $\mathbf{q}^*$ across parallel worker processes using 
#    tree-parallel message passing on the clique tree of the restricted Laplacian.
# 3. **Tikhonov Filter**: Smooths the resulting harmonic reference on each worker process to ensure physical 
#    realizability and eliminate instantaneous jumps.
# 4. **LQR Control**: Drives each agent to the filtered reference on its dedicated worker process using
#    an optimal state-feedback gain computed via the Discrete Algebraic Riccati Equation (DARE).
# 
# This separation of concerns allows the system to handle complex coordination constraints 
# (like the "mirror $y$-consensus" in Scenario 5) independently of the low-level vehicle dynamics.


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

# ## Setup Dynamics & Constants

# Quadrotor physical parameters
g = 9.81
m_veh = 0.5
I_quad = 0.01
ell = 0.25

# Continuous-time dynamics: x_dot = Ac*x + Bc*u
# State x = [y, z, theta, y_dot, z_dot, theta_dot]
Ac = [0.0  0.0   0.0   1.0  0.0  0.0;
      0.0  0.0   0.0   0.0  1.0  0.0;
      0.0  0.0   0.0   0.0  0.0  1.0;
      0.0  0.0  -g     0.0  0.0  0.0;
      0.0  0.0   0.0   0.0  0.0  0.0;
      0.0  0.0   0.0   0.0  0.0  0.0]

Bc = [0.0               0.0;
      0.0               0.0;
      0.0               0.0;
      0.0               0.0;
      1.0 / m_veh       1.0 / m_veh;
      ell / (2I_quad)  -ell / (2I_quad)]

h = 0.05
Ad, Bd = continuous_to_discrete_zoh(Ac, Bc, h)

nx = size(Ad, 1)
nu = size(Bd, 2)

# Compute Optimal LQR Gain (Discrete-time)
# We use a non-uniform Q matrix to penalize position errors more than velocity errors,
# reducing overshoot while maintaining a fast response.
Q_diag = zeros(nx)
Q_diag[1:2] .= 10000.0
Q_diag[3] = 50.0
Q_diag[4] = 500.0
Q_diag[5:6] .= 10.0
Q_lqr = Matrix{Float64}(Diagonal(Q_diag))
R_lqr = Matrix{Float64}(I, nu, nu) * 0.0001

function solve_dare(A, B, Q, R)
    P = Q
    for i in 1:100
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

# ## Coordination Sheaf Construction

# Scenario 5 Configuration:
# - Agent 1 tracks Target 1 in z
# - Agent 2 tracks Target 2 in yz
# - Agents agree in y (Consensus)
D = nx
NA = 2
NT = 2
TotalV = NA + NT

sheaf = EuclideanSheaf{Float64}(fill(D, TotalV))

# Coordinate projection matrices
R_y  = zeros(1, D); R_y[1] = 1.0
R_z  = zeros(1, D); R_z[2] = 1.0
R_yz = zeros(2, D); R_yz[1, 1] = 1.0; R_yz[2, 2] = 1.0

# Setup coordination edges
add_sheaf_edge!(sheaf, 1, 2, R_y, R_y)
add_sheaf_edge!(sheaf, 1, 3, R_z, R_z)
add_sheaf_edge!(sheaf, 2, 4, R_yz, R_yz)

# Boundary conditions for restricted Laplacian
boundary0 = Dict{Int, Vector{Float64}}()
boundary0[3] = zeros(D)
boundary0[4] = zeros(D)

_, _, Hraw, LIBraw = _harmonic_extension_restricted_laplacian(sheaf, boundary0)
H = Matrix(Hraw)
LIB = Matrix(LIBraw)

# Precompute chordal factorisation for tree-parallel solve
H_reg = sparse(H) + 1e-8 * I
F = cholesky!(ChordalCholesky(H_reg), NoPivot())
Lfac = F.L

# Provision worker processes for distributed parallel simulation
nchunks = length(partition_tree(Lfac, NA).chunks)
nworkers = max(NA, nchunks)
workers_pids = addprocs(nworkers; exeflags = "--project=$(Base.active_project())")

# ## Load Distributed Agent Component
#
# Load the worker-side flight computer implementation (`_scenario5_agent_worker.jl`)
# on process 1 and all worker processes into `Main` using `Main.include`.

worker_file = joinpath(pkgdir(CellularSheaves), "docs", "literate", "layered", "_scenario5_agent_worker.jl")

Main.include(worker_file)
@everywhere workers_pids Main.include($worker_file)

# ## Simulation Framework

function run_layered_simulation(target_type=:bobbing, mode=:distributed)
    local t1_pos, t2_pos
    if target_type == :bobbing
        omega = 2π * 2 / (40 * h)
        t1_pos = t -> [0.0, 1.0 + 0.3sin(omega*t), 0.0, 0.0, 0.0, 0.0]
        t2_pos = t -> [1.5, 2.0 + 0.3sin(omega*t), 0.0, 0.0, 0.0, 0.0]
    else
        t1_pos = t -> [0.0, 1.0, 0.0, 0.0, 0.0, 0.0]
        t2_pos = t -> [1.5, 2.0, 0.0, 0.0, 0.0, 0.0]
    end

    T_end = 2.0
    steps = Int(T_end / h) + 1
    epsilon = 0.02
    
    init_states = [[-2.0, 0.5, 0.0, 0.0, 0.0, 0.0], 
                   [2.0, 1.0, 0.0, 0.0, 0.0, 0.0]]

    for i in 1:NA # Initialize agent flight computer state on each worker process
        remotecall_fetch(Main.init_worker_agent!, workers_pids[i], init_states[i], K_lqr, Ad, Bd, epsilon)
    end

    sim_data = zeros(steps, NA, nx)
    qstar_history = zeros(steps, NA, nx)
    filtered_ref_history = zeros(steps, NA, nx)

    H_pinv = pinv(H)

    for t_idx in 0:(steps-1)
        t = t_idx * h
        b_t = [t1_pos(t); t2_pos(t)]

        rhs = Vector(-LIB * b_t)
        if mode == :centralised
            qstar_full = H_pinv * rhs
        else
            rhs_p = Vector(F.P' \ rhs)
            y_sol = distributed_tree_solve(Lfac, rhs_p, NA; pids = workers_pids)
            qstar_full = F.P \ y_sol
        end
        
        qstar = [qstar_full[1:D], qstar_full[D+1:2D]]
        for i in 1:NA
            qstar_history[t_idx+1, i, :] = qstar[i]
        end

        step_futures = [remotecall(Main.step_worker_agent!, workers_pids[i], qstar[i], h) for i in 1:NA] # Dispatch LQR tracking step to workers
        step_results = fetch.(step_futures)

        for i in 1:NA
            x_act, x_ref = step_results[i]
            filtered_ref_history[t_idx+1, i, :] = x_ref
            sim_data[t_idx+1, i, :] = x_act
        end
    end
    
    return sim_data, qstar_history, filtered_ref_history, t1_pos, t2_pos
end

# ## Execution and Visualization

function plot_results(sim_fixed, q_fixed, ref_fixed, t1_fixed, t2_fixed, sim_bob, q_bob, ref_bob, t1_bob, t2_bob, filename)
    p = plot(layout=(2, 2), size=(1200, 800), plot_title="Layered Control: Scenario 5 Analysis (Distributed Workers)")
    t_axis = range(0, step=h, length=size(sim_fixed, 1))

    plot!(p[1], t_axis, sim_fixed[:, 1, 1], label="A1 Actual", color=:steelblue, lw=2)
    plot!(p[1], t_axis, sim_fixed[:, 2, 1], label="A2 Actual", color=:darkorange, lw=2)
    plot!(p[1], t_axis, q_fixed[:, 1, 1], label="A1 q*", color=:steelblue, alpha=0.5, ls=:dash, marker=:circle, ms=2)
    plot!(p[1], t_axis, q_fixed[:, 2, 1], label="A2 q*", color=:darkorange, alpha=0.5, ls=:dash, marker=:circle, ms=2)
    plot!(p[1], t_axis, ref_fixed[:, 1, 1], label="A1 Ref (Tikh)", color=:steelblue, ls=:dot, lw=1.5)
    plot!(p[1], t_axis, ref_fixed[:, 2, 1], label="A2 Ref (Tikh)", color=:darkorange, ls=:dot, lw=1.5)
    plot!(p[1], t_axis, [t1_fixed(t)[1] for t in t_axis], label="Target 1", color=:gray, ls=:dot)
    plot!(p[1], t_axis, [t2_fixed(t)[1] for t in t_axis], label="Target 2", color=:black, ls=:dot)
    title!(p[1], "Fixed: Lateral Position (y)")
    xlabel!(p[1], "Time (s)")
    ylabel!(p[1], "Position (m)")

    plot!(p[2], t_axis, sim_fixed[:, 1, 2], label="A1 Actual", color=:steelblue, lw=2)
    plot!(p[2], t_axis, sim_fixed[:, 2, 2], label="A2 Actual", color=:darkorange, lw=2)
    plot!(p[2], t_axis, [q_fixed[t, 1, 2] for t in 1:41], label="A1 q*", color=:steelblue, alpha=0.5, ls=:dash, marker=:circle, ms=2)
    plot!(p[2], t_axis, [q_fixed[t, 2, 2] for t in 1:41], label="A2 q*", color=:darkorange, alpha=0.5, ls=:dash, marker=:circle, ms=2)
    plot!(p[2], t_axis, [ref_fixed[t, 1, 2] for t in 1:41], label="A1 Ref (Tikh)", color=:steelblue, ls=:dot, lw=1.5)
    plot!(p[2], t_axis, [ref_fixed[t, 2, 2] for t in 1:41], label="A2 Ref (Tikh)", color=:darkorange, ls=:dot, lw=1.5)
    plot!(p[2], t_axis, [t1_fixed(t)[2] for t in t_axis], label="Target 1", color=:gray, ls=:dot)
    plot!(p[2], t_axis, [t2_fixed(t)[2] for t in t_axis], label="Target 2", color=:black, ls=:dot)
    title!(p[2], "Fixed: Altitude (z)")
    xlabel!(p[2], "Time (s)")
    ylabel!(p[2], "Altitude (m)")

    plot!(p[3], t_axis, sim_bob[:, 1, 1], label="A1 Actual", color=:steelblue, lw=2)
    plot!(p[3], t_axis, sim_bob[:, 2, 1], label="A2 Actual", color=:darkorange, lw=2)
    plot!(p[3], t_axis, [q_bob[t, 1, 1] for t in 1:41], label="A1 q*", color=:steelblue, alpha=0.5, ls=:dash, marker=:circle, ms=2)
    plot!(p[3], t_axis, [q_bob[t, 2, 1] for t in 1:41], label="A2 q*", color=:darkorange, alpha=0.5, ls=:dash, marker=:circle, ms=2)
    plot!(p[3], t_axis, [ref_bob[t, 1, 1] for t in 1:41], label="A1 Ref (Tikh)", color=:steelblue, ls=:dot, lw=1.5)
    plot!(p[3], t_axis, [ref_bob[t, 2, 1] for t in 1:41], label="A2 Ref (Tikh)", color=:darkorange, ls=:dot, lw=1.5)
    plot!(p[3], t_axis, [t1_bob(t)[1] for t in t_axis], label="Target 1", color=:gray, ls=:dot)
    plot!(p[3], t_axis, [t2_bob(t)[1] for t in t_axis], label="Target 2", color=:black, ls=:dot)
    title!(p[3], "Bobbing: Lateral Position (y)")
    xlabel!(p[3], "Time (s)")
    ylabel!(p[3], "Position (m)")

    plot!(p[4], t_axis, sim_bob[:, 1, 2], label="A1 Actual", color=:steelblue, lw=2)
    plot!(p[4], t_axis, sim_bob[:, 2, 2], label="A2 Actual", color=:darkorange, lw=2)
    plot!(p[4], t_axis, [q_bob[t, 1, 2] for t in 1:41], label="A1 q*", color=:steelblue, alpha=0.5, ls=:dash, marker=:circle, ms=2)
    plot!(p[4], t_axis, [q_bob[t, 2, 2] for t in 1:41], label="A2 q*", color=:darkorange, alpha=0.5, ls=:dash, marker=:circle, ms=2)
    plot!(p[4], t_axis, [ref_bob[t, 1, 2] for t in 1:41], label="A1 Ref (Tikh)", color=:steelblue, ls=:dot, lw=1.5)
    plot!(p[4], t_axis, [ref_bob[t, 2, 2] for t in 1:41], label="A2 Ref (Tikh)", color=:darkorange, ls=:dot, lw=1.5)
    plot!(p[4], t_axis, [t1_bob(t)[2] for t in t_axis], label="Target 1", color=:gray, ls=:dot)
    plot!(p[4], t_axis, [t2_bob(t)[2] for t in t_axis], label="Target 2", color=:black, ls=:dot)
    title!(p[4], "Bobbing: Altitude (z)")
    xlabel!(p[4], "Time (s)")
    ylabel!(p[4], "Altitude (m)")
    
    savefig(filename)
end

# Run distributed simulation scenarios
sim_fixed, q_fixed, ref_fixed, t1_fixed, t2_fixed = run_layered_simulation(:fixed, :distributed)
sim_bob, q_bob, ref_bob, t1_bob, t2_bob = run_layered_simulation(:bobbing, :distributed)

# Compare distributed against centralized simulation to verify precision
sim_fixed_c, _, _, _, _ = run_layered_simulation(:fixed, :centralised)
divergence = maximum(abs.(sim_fixed .- sim_fixed_c))
@printf("Max divergence between centralized and distributed simulation: %.3e\n", divergence)

# Clean up worker processes
rmprocs(workers_pids)

plot_results(sim_fixed, q_fixed, ref_fixed, t1_fixed, t2_fixed, sim_bob, q_bob, ref_bob, t1_bob, t2_bob, "layered_control_combined.png")

# ![Layered Control Analysis](layered_control_combined.png)
