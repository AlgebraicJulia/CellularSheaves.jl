# Compare optimal_control_trajectory vs IPM for vehicle platoon example

using CellularSheaves
using CellularSheaves.IPM
using LinearAlgebra
using SparseArrays
using BlockArrays
using BenchmarkTools

# ============================================================================
# Setup: Vehicle Platoon from literate example
# ============================================================================

Ac = [0.0  1.0  0.0  0.0;
      0.0  0.0  0.0  0.0;
      0.0  0.0  0.0  1.0;
      0.0  0.0  0.0  0.0]

Bc = [0.0  0.0;
      1.0  0.0;
      0.0  0.0;
      0.0  1.0]

h = 0.25   # sample period (seconds)
k = 500    # number of time steps (even larger)

F  = EuclideanSheaf{Float64}([2, 2])
ts = ControlledTrajectorySheaf(F, Ac, Bc, h, k)

n = ts.state_dim    # 4
m = ts.control_dim  # 2

x1  = [0.0, 0.0, 2.0, 0.0]
xk1 = [1.0, 0.0, 3.0, 0.0]

Q  = Matrix{Float64}(I, n, n)
Ru = Matrix{Float64}(I, m, m)
Qf = 10.0 * Matrix{Float64}(I, n, n)

H, f, _ = lqr_objective(ts, Q, Ru; Qf=Qf)

# ============================================================================
# Method 1: optimal_control_trajectory (current method)
# ============================================================================

println("="^60)
println("Method 1: optimal_control_trajectory")
println("="^60)

z_opt1, α_opt1, z_p1, N1 = optimal_control_trajectory(ts, x1, xk1, H, f)
z_opt1_vec = Array(z_opt1)

cost1 = 0.5 * dot(z_opt1_vec, H * z_opt1_vec) + dot(f, z_opt1_vec)
println("Optimal cost: ", cost1)
println("Null space dim: ", size(N1, 2))

# Verify constraints
println("Initial state error: ", norm(Array(z_opt1[Block(1)]) - x1))
println("Terminal state error: ", norm(Array(z_opt1[Block(k+1)]) - xk1))

# ============================================================================
# Method 2: IPM
# ============================================================================

println("\n" * "="^60)
println("Method 2: IPM")
println("="^60)

# Build constraint matrix for dynamics + boundary conditions
# Decision variable: z = [x_1, ..., x_{k+1}, u_1, ..., u_k]
#
# Dynamics constraints: x_{t+1} - A*x_t - B*u_t = 0  for t = 1,...,k
# Boundary constraints: x_1 = x1, x_{k+1} = xk1

p = (k + 1) * n + k * m  # total decision variable dimension

# Number of constraints: k dynamics + 2 boundary (initial and terminal)
num_dyn_constraints = k * n
num_bnd_constraints = 2 * n
num_constraints = num_dyn_constraints + num_bnd_constraints

# Build constraint matrix B and RHS g
# Using dense for simplicity first
B_constr = zeros(Float64, num_constraints, p)
g_constr = zeros(Float64, num_constraints)

Ad = ts.Ad
Bd = ts.Bd

# Dynamics: x_{t+1} - A*x_t - B*u_t = 0
for t in 1:k
    row = (t - 1) * n + 1
    state_col_t = (t - 1) * n + 1
    state_col_t1 = t * n + 1
    ctrl_col_t = (k + 1) * n + (t - 1) * m + 1

    # -A * x_t
    B_constr[row:row+n-1, state_col_t:state_col_t+n-1] = -Ad
    # I * x_{t+1}
    B_constr[row:row+n-1, state_col_t1:state_col_t1+n-1] = I(n)
    # -B * u_t
    B_constr[row:row+n-1, ctrl_col_t:ctrl_col_t+m-1] = -Bd
end

# Boundary: x_1 = x1
row = num_dyn_constraints + 1
B_constr[row:row+n-1, 1:n] = I(n)
g_constr[row:row+n-1] = x1

# Boundary: x_{k+1} = xk1
row = num_dyn_constraints + n + 1
B_constr[row:row+n-1, k*n+1:(k+1)*n] = I(n)
g_constr[row:row+n-1] = xk1

# Convert to BlockSparseMatrix format for IPM
# Q matrix (Hessian) - need to convert H to BlockSparseMatrix
# B matrix (constraints) - need to convert B_constr to BlockSparseMatrix

# For now, let's use a single block (one vertex with dimension p)
# This is the simplest case - no block structure

# Build Q as BlockSparseMatrix
Q_ipm = blocksparse(
    [1], [1],
    [Matrix(H)],
    [p], [p]
)

# Build B as BlockSparseMatrix
B_ipm = blocksparse(
    [1], [1],
    [B_constr],
    [num_constraints], [p]
)

# Vectors
c_ipm = collect(f)
g_ipm = g_constr

# Cone: CofreeCone (no cone constraint, just equality)
K_ipm = [CofreeCone()]

# Create IPM problem
prob = IPMProblem(Q_ipm, B_ipm, c_ipm, g_ipm, K_ipm)

# Solve
println("Solving with IPM...")
result = solve(prob; verbose=true, feas_tol=1e-10, gap_tol=1e-10)

println("\nIPM status: ", result.status)
println("IPM iterations: ", result.niter)

z_opt2 = result.p
cost2 = 0.5 * dot(z_opt2, H * z_opt2) + dot(f, z_opt2)
println("Optimal cost: ", cost2)

# Verify constraints
println("Initial state error: ", norm(z_opt2[1:n] - x1))
println("Terminal state error: ", norm(z_opt2[k*n+1:(k+1)*n] - xk1))

# ============================================================================
# Compare solutions
# ============================================================================

println("\n" * "="^60)
println("Comparison")
println("="^60)

println("Solution difference (norm): ", norm(z_opt1_vec - z_opt2))
println("Cost difference: ", abs(cost1 - cost2))

# ============================================================================
# Timing comparison
# ============================================================================

println("\n" * "="^60)
println("Timing (excluding setup)")
println("="^60)

println("\noptimal_control_trajectory:")
@btime optimal_control_trajectory($ts, $x1, $xk1, $H, $f)

println("\nIPM solve:")
@btime solve($prob; verbose=false, feas_tol=1e-10, gap_tol=1e-10)
