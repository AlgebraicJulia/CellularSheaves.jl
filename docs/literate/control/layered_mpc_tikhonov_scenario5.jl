
using CellularSheaves
import CellularSheaves.NetworkSheaves.EuclideanSheaves: _harmonic_extension_restricted_laplacian
using CellularSheaves.ControlSheaves.Tikhonov
using CellularSheaves.TrajectorySheaves: continuous_to_discrete_zoh
using LinearAlgebra
using Plots
using Printf

# --- 1. Setup Dynamics & Constants ---
g = 9.81
m_veh = 0.5
I_quad = 0.01
ell = 0.25

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

IDX_Y = 1
IDX_Z = 2

# Compute Optimal LQR Gain (Discrete-time)
# Cost matrices: Q (state penalty) and R (control penalty)
# Penalize positions (1:3) heavily, velocities (4:6) moderately
Q_diag = zeros(nx)
Q_diag[1:2] .= 100000.0
Q_diag[3] = 10.0 # We don't need alignment in angle
Q_diag[4] = 50.0 # High damping for y-velocity to reduce overshoot
Q_diag[5:6] .= 10.0
Q_lqr = Matrix{Float64}(Diagonal(Q_diag))
R_lqr = Matrix{Float64}(I, nu, nu) * 0.0001

function solve_dare(A, B, Q, R)
    # Simple iterative solve for the Discrete Algebraic Riccati Equation (DARE)
    # P_{k+1} = A' P_k A - (A' P_k B)(R + B' P_k B)^-1 (B' P_k A) + Q
    P = Q
    for i in 1:100
        P_next = A' * P * A - (A' * P * B) * inv(R + B' * P * B) * (B' * P * A) + Q
        if norm(P_next - P) < 1e-6
            break
        end
        P = P_next
    end
    # Optimal Gain K = (R + B' P B)^-1 * B' P A
    return (R + B' * P * B) \ (B' * P * A)
end

K_lqr = solve_dare(Ad, Bd, Q_lqr, R_lqr)

# Stability Check
A_cl = Ad - Bd * K_lqr
rho = maximum(abs.(eigvals(A_cl)))
@printf("Closed-loop spectral radius: %.4f\n", rho)
if rho >= 1.0
    @warn "System is unstable! Spectral radius rho >= 1.0"
elseif rho > 0.9
    @info "System is stable but potentially slow or oscillatory (rho > 0.9)"
end

# --- 2. Static Coordination Sheaf (Scenario 5) ---
# Scenario 5:
# - Agent 1 tracks T1 in z
# - Agent 2 tracks T2 in yz
# - Agents agree in y (Consensus)
# - All agents have stalk R^nx (we only care about the state for the reference)

D = nx
NA = 2 # Agents
NT = 2 # Targets
TotalV = NA + NT

sheaf = EuclideanSheaf{Float64}(fill(D, TotalV))

# Projection matrices for coordination
R_y  = zeros(1, D)
R_y[1] = 1.0
R_z  = zeros(1, D)
R_z[2] = 1.0 # Projecting to z (index 2)
R_yz = zeros(2, D)
R_yz[1, 1] = 1.0; R_yz[2, 2] = 1.0

# Consensus Edge: A1 <-> A2 in y
add_sheaf_edge!(sheaf, 1, 2, R_y, R_y)

# Tracking Edges:
# A1 tracks T1 in z
add_sheaf_edge!(sheaf, 1, 3, R_z, R_z) 
# A2 tracks T2 in yz
add_sheaf_edge!(sheaf, 2, 4, R_yz, R_yz)

# Restricted Laplacian setup
boundary0 = Dict{Int, Vector{Float64}}()
boundary0[3] = zeros(D)
boundary0[4] = zeros(D)

_, _, Hraw, LIBraw = _harmonic_extension_restricted_laplacian(sheaf, boundary0)
H = Matrix(Hraw)
LIB = Matrix(LIBraw)

# --- 3. Simulation Loop ---

# Targets (Bobbing)
omega = 2π * 2 / (40 * h)
t1_pos(t) = [0.0, 1.0 + 0.3sin(omega*t), 0.0, 0.0, 0.0, 0.0]
t2_pos(t) = [1.5, 2.0 + 0.3sin(omega*t), 0.0, 0.0, 0.0, 0.0]

epsilon = 0.02
filters = [TikhonovFilter(zeros(D); epsilon) for _ in 1:2]
curr_states = [[-2.0, 0.5, 0.0, 0.0, 0.0, 0.0], 
               [2.0, 1.0, 0.0, 0.0, 0.0, 0.0]]

sim_a1 = zeros(41, nx)
sim_a2 = zeros(41, nx)
qstar_history = zeros(41, 2, nx)
filtered_ref_history = zeros(41, 2, nx)

for t_idx in 0:40
    t = t_idx * h
    
    # Update boundary data with current target positions
    b_t = [t1_pos(t); t2_pos(t)]
    
    # Use a robust solve for reference points qstar
    #Since H may be singular (due to the consensus/tracking compromise), 
    #we use the pseudo-inverse via pinv() or a regularized solve.
    #tikhonov_equilibrium currently uses \ which fails on singular H.
    qstar_full = pinv(H) * (-LIB * b_t)
    
    # Split qstar_full into agent references
    # qstar_full is a vector of size NA * D
    qstar = [qstar_full[1:D], qstar_full[D+1:2D]]
    for i in 1:2
        qstar_history[t_idx+1, i, :] = qstar[i]
    end
    
    # Filter references
    for i in 1:2
        tikhonov_step!(filters[i], qstar[i], qstar[i], h)
        filtered_ref_history[t_idx+1, i, :] = filters[i].x
    end
    
    # Local LQR Control: u = -K * (x - x_ref)
    for i in 1:2
        x = curr_states[i]
        x_ref = filters[i].x
        
        u = -K_lqr * (x - x_ref)
        
        # Apply dynamics
        curr_states[i] = Ad * x + Bd * u
    end
    
    sim_a1[t_idx+1, :] = curr_states[1]
    sim_a2[t_idx+1, :] = curr_states[2]
end

# Visualization
p = plot(layout=(1, 2), size=(1200, 400), plot_title="Layered Control: Scenario 5 Tracking")
plot!(p[1], 0:0.05:2, sim_a1[:, 1], label="A1 Actual", color=:steelblue, lw=2)
plot!(p[1], 0:0.05:2, sim_a2[:, 1], label="A2 Actual", color=:darkorange, lw=2)
plot!(p[1], 0:0.05:2, [qstar_history[t, 1, 1] for t in 1:41], label="A1 q*", color=:steelblue, alpha=0.5, ls=:dash, marker=:circle, ms=2)
plot!(p[1], 0:0.05:2, [qstar_history[t, 2, 1] for t in 1:41], label="A2 q*", color=:darkorange, alpha=0.5, ls=:dash, marker=:circle, ms=2)
plot!(p[1], 0:0.05:2, [filtered_ref_history[t, 1, 1] for t in 1:41], label="A1 Ref (Tikh)", color=:steelblue, ls=:dot, lw=1.5)
plot!(p[1], 0:0.05:2, [filtered_ref_history[t, 2, 1] for t in 1:41], label="A2 Ref (Tikh)", color=:darkorange, ls=:dot, lw=1.5)
plot!(p[1], 0:0.05:2, [t1_pos(t)[1] for t in 0:0.05:2], label="Target 1", color=:gray, ls=:dot)
plot!(p[1], 0:0.05:2, [t2_pos(t)[1] for t in 0:0.05:2], label="Target 2", color=:black, ls=:dot)
title!(p[1], "Lateral Position (y)")
xlabel!(p[1], "Time (s)")
ylabel!(p[1], "Position (m)")

plot!(p[2], 0:0.05:2, sim_a1[:, 2], label="A1 Actual", color=:steelblue, lw=2)
plot!(p[2], 0:0.05:2, sim_a2[:, 2], label="A2 Actual", color=:darkorange, lw=2)
plot!(p[2], 0:0.05:2, [qstar_history[t, 1, 2] for t in 1:41], label="A1 q*", color=:steelblue, alpha=0.5, ls=:dash, marker=:circle, ms=2)
plot!(p[2], 0:0.05:2, [qstar_history[t, 2, 2] for t in 1:41], label="A2 q*", color=:darkorange, alpha=0.5, ls=:dash, marker=:circle, ms=2)
plot!(p[2], 0:0.05:2, [filtered_ref_history[t, 1, 2] for t in 1:41], label="A1 Ref (Tikh)", color=:steelblue, ls=:dot, lw=1.5)
plot!(p[2], 0:0.05:2, [filtered_ref_history[t, 2, 2] for t in 1:41], label="A2 Ref (Tikh)", color=:darkorange, ls=:dot, lw=1.5)
plot!(p[2], 0:0.05:2, [t1_pos(t)[2] for t in 0:0.05:2], label="Target 1", color=:gray, ls=:dot)
plot!(p[2], 0:0.05:2, [t2_pos(t)[2] for t in 0:0.05:2], label="Target 2", color=:black, ls=:dot)
title!(p[2], "Altitude (z)")
xlabel!(p[2], "Time (s)")
ylabel!(p[2], "Altitude (m)")
savefig("layered_control_static_sheaf_fixed.png")
