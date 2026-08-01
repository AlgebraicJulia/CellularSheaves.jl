# # Layered Lagless Control: Harmonic Goals $\to$ Tikhonov $\to$ Control Affine
#
# This example demonstrates a layered control architecture designed to eliminate
# tracking lag. We move from a high-level harmonic goal to a filtered reference,
# and finally to a low-level control affine tracking law.
#
# ## The Layered Pipeline
#
# 1. **Harmonic Layer**: Solves for the ideal goal $q^\star(t)$ using the
#    sheaf-theoretic harmonic extension.
# 2. **Tikhonov Layer**: Smooths the goal through a normalized first-order filter.
#    We use **analytic feedforward** to cancel the $O(\epsilon)$ lag.
# 3. **Control Affine Layer**: A low-level controller tracks the filtered 
#    reference $x_f(t)$ for an agent with dynamics $\dot{x} = f(x) + B u$.
#
# ```math
# u = B^\dagger \left( \dot{x}_f - f(x) - K(x - x_f) \right)
# ```
#
# This ensures that neither the planner's regularization nor the agent's 
# inertia introduces significant phase lag.

using CellularSheaves
using LinearAlgebra
using Plots
using CellularSheaves.ControlSheaves
using CellularSheaves.ControlSheaves.Tikhonov

# ## Step 1: System Definition (Planar Quadrotor)
#
# We use the planar quadrotor model linearized around hover.
# Note: We define dynamics as functions to facilitate future nonlinear generalizations.

g  = 9.81    # (m/s²)
m  = 0.5     # (kg)
I_quad = 0.01 # (kg·m²)
ℓ  = 0.25    # (m)

# Linearized matrices for the quadrotor (from controlled_planar_quadrotor.jl)
Ac = [0.0  0.0   0.0   1.0  0.0  0.0;
      0.0  0.0   0.0   0.0  1.0  0.0;
      0.0  0.0   0.0   0.0  0.0  1.0;
      0.0  0.0  -g     0.0  0.0  0.0;
      0.0  0.0   0.0   0.0  0.0  0.0;
      0.0  0.0   0.0   0.0  0.0  0.0]

Bc = [0.0           0.0;
      0.0           0.0;
      0.0           0.0;
      0.0           0.0;
      1.0/m         1.0/m;
      ℓ/(2.0*I_quad)  -ℓ/(2.0*I_quad)]

# General dynamics wrapper. 
# FLAG: This currently uses a linear model (Ac*x), but is structured for f(x).
function agent_f(x)
    return Ac * x
end

# General control input map.
# FLAG: This currently uses a constant matrix B, but is structured for B(x).
function agent_B(x)
    return Bc
end

# ## Step 2: Reference Generation
#
# We create a time-varying goal $q^\star(t)$ and its derivative $\dot{q}^\star(t)$.
# In a full system, these would come from the Harmonic Extension.

t_span = 0.0:0.001:5.0
# A simple oscillatory goal in the y-axis (lateral position)
q_star = [ [0.5 * sin(t), 0.0, 0.0, 0.0, 0.0, 0.0] for t in t_span ]
# conceptual
# Real implementation below:
q_star_vals = [ [0.5 * sin(t), 0.0, 0.0, 0.0, 0.0, 0.0] for t in t_span ]
dq_star_vals = [ [0.5 * cos(t), 0.0, 0.0, 0.0, 0.0, 0.0] for t in t_span ]

# ## Step 3: Tikhonov Filtering (Reference Generation)
#
# We generate a lagless filtered reference $x_f(t)$ and its derivative $\dot{x}_f(t)$
# using the Tikhonov feedforward planner. This serves as the input to our
# low-level controllers.

epsilon = 0.1
x0 = zeros(6)

filter_ref = TikhonovFilter(x0; epsilon=epsilon)

xf_vals = Vector{Float64}[]
dx_f_vals = Vector{Float64}[]

dt = t_span[2] - t_span[1]

for i in 1:length(t_span)-1
    q0 = q_star_vals[i]
    q1 = q_star_vals[i+1]
    
    # Compute feedforward input: u_ff = q_star + epsilon * dq_star
    u0 = tikhonov_feedforward_reference(q0, dq_star_vals[i], epsilon)
    u1 = tikhonov_feedforward_reference(q1, dq_star_vals[i+1], epsilon)
    
    # Advance filter
    tikhonov_step!(filter_ref, u0, u1, dt)
    
    # Store filtered state and its instantaneous rate: dx_f = (-x_f + u_ff) / epsilon
    push!(xf_vals, copy(filter_ref.x))
    push!(dx_f_vals, (-filter_ref.x + u0) / epsilon)
end
push!(xf_vals, copy(filter_ref.x))
# Pad last derivative to match length
push!(dx_f_vals, dx_f_vals[end])

# ## Step 4: Inner Loop Simulation (LQR Tracking)
#
# We compare two low-level tracking laws:
# 1. **LQR Feedback Only**: Ignores the reference rate, resulting in lag.
# 2. **LQR Feedforward + Feedback**: Uses the reference rate to cancel tracking lag.
#
# For a linear system $\dot{x} = Ax + Bu$, the tracking law is:
# $u = -K(x - x_f) + B^\dagger(\dot{x}_f - Ax_f)$.

# Solve for K using the Algebraic Riccati Equation (ARE):
# A'P + PA - PBR⁻¹B'P + Q = 0
# K = R⁻¹B'P

function solve_lqr_gain(A, B, Q, R)
    nx = size(A, 1)
    nu = size(B, 2)
    
    R_inv = inv(R)
    B_invR = B * R_inv
    B_invR_Bt = B_invR * B'
    
    H = [A + B_invR_Bt * A   -B_invR_Bt;
         -Q                   A' - (B_invR_Bt * A)']
    
    vals, vecs = eigen(H)
    
    p = sortperm(real(vals))
    V_stable = vecs[:, p[1:nx]]
    
    X = V_stable[1:nx, :]
    Y = V_stable[nx+1:end, :]
    
    P = X * pinv(Y)
    
    # The resulting K might be complex due to eigen(H), take real part
    return real(R_inv * B' * P)
end

# Define LQR weights
Q_lqr = Diagonal([10.0, 10.0, 1.0, 1.0, 1.0, 1.0]) 
R_lqr = Matrix{Float64}(I, size(Bc, 2), size(Bc, 2)) * 1.0

# Calculate K dynamically
K_lqr = solve_lqr_gain(Ac, Bc, Q_lqr, R_lqr)

function compute_lqr_control(x, x_f, dx_f, use_ff=true)
    error = x - x_f
    u_fb = -K_lqr * error
    
    u_ff = zeros(size(Bc, 2))
    if use_ff
        res = dx_f - Ac * x_f
        u_ff = pinv(Bc) * res
    end
    return u_fb + u_ff
end

# Simulation loop for LQR tracking
function simulate_lqr_tracking(xf_traj, dx_f_traj, x_start, dt, use_ff)
    x = copy(x_start)
    traj = [copy(x)]
    
    for i in 1:length(xf_traj)-1
        xf = xf_traj[i]
        dx_f = dx_f_traj[i]
        
        u = compute_lqr_control(x, xf, dx_f, use_ff)
        
        dx = Ac * x + Bc * u
        x .+= dx * dt
        push!(traj, copy(x))
    end
    return traj
end

# Run both scenarios
traj_fb_only = simulate_lqr_tracking(xf_vals, dx_f_vals, x0, dt, false)
traj_ff_fb   = simulate_lqr_tracking(xf_vals, dx_f_vals, x0, dt, true)

# ## Step 5: Visualization
#
# We compare the trajectories against the filtered reference x_f.

function plot_tracking_comparison(times, xf_traj, traj_fb, traj_ff)
    y_ref = [s[1] for s in xf_traj]
    y_fb  = [s[1] for s in traj_fb]
    y_ff  = [s[1] for s in traj_ff]
    
    p = plot(times, y_ref, lw=3, color=:black, label="Filtered Ref x_f(t)", 
             title="LQR Low-Level Feedforward Comparison", xlabel="Time (s)", ylabel="Position y (m)")
    plot!(p, times, y_fb, lw=2, color=:red, label="LQR Feedback Only (Laggy)")
    plot!(p, times, y_ff, lw=2, color=:green, label="LQR FF + Feedback (Lagless)")
    
    return p
end

p_final = plot_tracking_comparison(t_span, xf_vals, traj_fb_only, traj_ff_fb)
display(p_final)

