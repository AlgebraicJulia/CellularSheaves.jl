# Compare optimal_control_trajectory vs IPM with PROPER block structure
# Using natural sheaf structure: edge stalks → row blocks, vertex stalks → col blocks

using CellularSheaves
using CellularSheaves.NetworkSheaves: coboundary_map, vertex_stalks, edge_stalks
using CellularSheaves.IPM
using LinearAlgebra
using SparseArrays
using BenchmarkTools

# ============================================================================
# Setup: Vehicle Platoon (from docs)
# ============================================================================

Ac = [0.0  1.0  0.0  0.0;
      0.0  0.0  0.0  0.0;
      0.0  0.0  0.0  1.0;
      0.0  0.0  0.0  0.0]

Bc = [0.0  0.0;
      1.0  0.0;
      0.0  0.0;
      0.0  1.0]

h = 0.25
k = 100  # time steps

F  = EuclideanSheaf{Float64}([2, 2])
ts = ControlledTrajectorySheaf(F, Ac, Bc, h, k)

n = ts.state_dim    # 4
m = ts.control_dim  # 2

x1  = [0.0, 0.0, 2.0, 0.0]
xk1 = [1.0, 0.0, 3.0, 0.0]

Q_lqr  = Matrix{Float64}(I, n, n)
Ru = Matrix{Float64}(I, m, m)
Qf = 10.0 * Matrix{Float64}(I, n, n)

H, f, _ = lqr_objective(ts, Q_lqr, Ru; Qf=Qf)

# ============================================================================
# Method 1: optimal_control_trajectory
# ============================================================================

println("="^60)
println("Method 1: optimal_control_trajectory")
println("="^60)

z_opt1, α_opt1, z_p1, N1 = optimal_control_trajectory(ts, x1, xk1, H, f)
z_opt1_vec = Array(z_opt1)

cost1 = 0.5 * dot(z_opt1_vec, H * z_opt1_vec) + dot(f, z_opt1_vec)
println("Optimal cost: ", cost1)

# ============================================================================
# Method 2: IPM with PROPER block structure
# ============================================================================

println("\n" * "="^60)
println("Method 2: IPM (proper block structure from sheaf)")
println("="^60)

sheaf = ts.sheaf
vstks = vertex_stalks(sheaf)

println("Vertex stalks: ", vstks[1:min(5,end)], "... (", length(vstks), " total)")
println("Number of edges: ", length(edge_stalks(sheaf)))

# Interior vertices: 2, 3, ..., k+1 (total k vertices)
# These store (x_t, u_t) for t=1..k
interior_verts = 2:k+1
n_interior = length(interior_verts)

# Get the full coboundary map as a BlockSparseMatrix directly from sheaf
d_full = coboundary_map(sheaf)
println("Coboundary type: ", typeof(d_full))
println("Coboundary size: ", size(d_full))

# Extract interior columns while preserving block structure
# Column blocks: vertices → interior vertices only
# Row blocks: edges → all edges (they connect to interior)

# We need to select only interior vertex columns from d_full
# Use selectvtxs to get the proper block structure

using CellularSheaves.BlockSparseArrays: selectvtxs, nvtxs, nouts, narcs

d_int = selectvtxs(d_full, collect(interior_verts))
println("\nReduced coboundary (interior cols):")
println("  size: ", size(d_int))
println("  nvtxs (col blocks): ", nvtxs(d_int))
println("  nouts (row blocks): ", nouts(d_int))
println("  narcs (non-zero blocks): ", narcs(d_int))

# Compute boundary contribution: g = -d_bnd * x_bnd
# Boundary vertices are 1 and k+2
d_bnd = selectvtxs(d_full, [1, k+2])
x_bnd = vcat(x1, xk1)
g_reduced = -d_bnd * x_bnd

println("RHS g_reduced norm: ", norm(g_reduced))

# Build block-diagonal Q for interior vertices
# Each interior vertex has stalk size n+m = 6
Q_block = zeros(n+m, n+m)
Q_block[1:n, 1:n] = Q_lqr
Q_block[n+1:n+m, n+1:n+m] = Ru

I_Q = Int[]
J_Q = Int[]
V_Q = Matrix{Float64}[]
for v in 1:n_interior
    push!(I_Q, v)
    push!(J_Q, v)
    push!(V_Q, copy(Q_block))
end
col_dims = fill(n+m, n_interior)
row_dims = fill(n+m, n_interior)
Q_ipm = blocksparse(I_Q, J_Q, V_Q, row_dims, col_dims)

# Linear term (usually zero for standard LQR)
c_interior = zeros(n_interior * (n+m))

# Cones: one CofreeCone per interior vertex
K_ipm = [CofreeCone() for _ in 1:n_interior]

# Create IPM problem with proper block structure
prob = IPMProblem(Q_ipm, d_int, c_interior, vec(g_reduced), K_ipm)

println("\nIPM Problem structure:")
println("  Q: ", size(Q_ipm), " with ", narcs(Q_ipm), " blocks")
println("  B: ", size(d_int), " with ", narcs(d_int), " blocks")

println("\nSolving with IPM...")
result = solve(prob; verbose=true, feas_tol=1e-10, gap_tol=1e-10)

println("\nIPM status: ", result.status)
println("IPM iterations: ", result.niter)

# Extract solution and convert to public ordering
p_interior = result.p

# Reconstruct full solution
z_opt2_pub = zeros((k+1)*n + k*m)
z_opt2_pub[1:n] = x1
z_opt2_pub[k*n+1:(k+1)*n] = xk1

for t in 1:k
    int_start = (t-1)*(n+m) + 1
    x_t = p_interior[int_start : int_start+n-1]
    u_t = p_interior[int_start+n : int_start+n+m-1]

    z_opt2_pub[(t-1)*n+1 : t*n] = x_t
    z_opt2_pub[(k+1)*n + (t-1)*m+1 : (k+1)*n + t*m] = u_t
end

cost2 = 0.5 * dot(z_opt2_pub, H * z_opt2_pub) + dot(f, z_opt2_pub)
println("Optimal cost: ", cost2)

# ============================================================================
# Compare
# ============================================================================

println("\n" * "="^60)
println("Comparison")
println("="^60)

println("Solution difference (norm): ", norm(z_opt1_vec - z_opt2_pub))
println("Cost difference: ", abs(cost1 - cost2))

# ============================================================================
# Timing
# ============================================================================

println("\n" * "="^60)
println("Timing")
println("="^60)

println("\noptimal_control_trajectory:")
@btime optimal_control_trajectory($ts, $x1, $xk1, $H, $f)

println("\nIPM solve (excluding problem setup):")
@btime solve($prob; verbose=false, feas_tol=1e-10, gap_tol=1e-10)

println("\nIPM full (including problem setup):")
@btime begin
    prob_fresh = IPMProblem($Q_ipm, $d_int, $c_interior, vec($g_reduced), $K_ipm)
    solve(prob_fresh; verbose=false, feas_tol=1e-10, gap_tol=1e-10)
end
