# Compare optimal_control_trajectory vs IPM with proper block structure
# After reduction 1 (harmonic extension / boundary pinning) but before null space

using CellularSheaves
using CellularSheaves.NetworkSheaves: coboundary_map, vertex_stalks, harmonic_extension
using CellularSheaves.IPM
using LinearAlgebra
using SparseArrays
using BlockArrays
using BenchmarkTools

# ============================================================================
# Setup: Vehicle Platoon
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
println("Null space dim: ", size(N1, 2))

# ============================================================================
# Method 2: IPM with proper block structure (after boundary pinning)
# ============================================================================

println("\n" * "="^60)
println("Method 2: IPM (block-structured, reduced problem)")
println("="^60)

# The inner sheaf has k+2 vertices:
#   vertex 1: x₁ (dim n) - BOUNDARY (pinned)
#   vertices 2..k+1: (xₜ, uₜ) for t=1..k (dim n+m each) - INTERIOR
#   vertex k+2: x_{k+1} (dim n) - BOUNDARY (pinned)

sheaf = ts.sheaf
stalks = vertex_stalks(sheaf)
println("Sheaf stalks: ", stalks)

# Interior vertices: 2, 3, ..., k+1 (total k vertices)
interior_verts = 2:k+1
n_interior = length(interior_verts)
interior_stalks = stalks[interior_verts]
println("Interior vertices: ", n_interior, " with stalks ", interior_stalks[1:min(3,end)], "...")

# Get the full coboundary map
d_full = coboundary_map(sheaf)
println("Full coboundary size: ", size(d_full))

# We need to extract the reduced system after pinning boundary vertices
# The coboundary d maps 0-cochains to 1-cochains (edges)
# After pinning vertices 1 and k+2, we split: d * x = d_int * x_int + d_bnd * x_bnd
# The constraint becomes: d_int * x_int = -d_bnd * x_bnd = g_reduced

# Compute column ranges for each vertex in the 0-cochain
cumstalk = [0; cumsum(stalks)]
col_range(v) = cumstalk[v]+1 : cumstalk[v+1]

# Boundary columns
bnd_cols = vcat(collect(col_range(1)), collect(col_range(k+2)))
# Interior columns
int_cols = vcat([collect(col_range(v)) for v in interior_verts]...)

d_full_sparse = sparse(d_full)
d_int = d_full_sparse[:, int_cols]
d_bnd = d_full_sparse[:, bnd_cols]

# Boundary values
x_bnd = vcat(x1, xk1)
g_reduced = -d_bnd * x_bnd

println("Reduced coboundary size: ", size(d_int))
println("RHS g_reduced norm: ", norm(g_reduced))

# Now build the block-diagonal Q for interior vertices
# The LQR objective H is in "public" ordering: [x₁...x_{k+1}, u₁...uₖ]
# We need Q in "internal" ordering for interior: [(x₁,u₁), (x₂,u₂), ..., (xₖ,uₖ)]

# From lqr_objective, H is block-diagonal:
#   - Q_lqr on x₁, ..., xₖ (states 1 to k)
#   - Qf on x_{k+1} (terminal state)
#   - Ru on u₁, ..., uₖ (controls)

# For interior vertices 2..k+1 which store (xₜ, uₜ) for t=1..k:
# The cost for vertex t+1 (storing (xₜ, uₜ)) should be:
#   - Q_lqr on xₜ part
#   - Ru on uₜ part
# Exception: vertex k+1 stores (xₖ, uₖ), but x_{k+1} cost is Qf (but x_{k+1} is at boundary vertex k+2!)

# Actually looking more carefully at the internal structure:
# - vertex 2 stores (x₁, u₁)
# - vertex 3 stores (x₂, u₂)
# - ...
# - vertex k+1 stores (xₖ, uₖ)
# - vertex k+2 stores x_{k+1} (boundary, pinned)

# So for interior vertices, each has Q_block = [Q_lqr 0; 0 Ru]
Q_block = zeros(n+m, n+m)
Q_block[1:n, 1:n] = Q_lqr
Q_block[n+1:n+m, n+1:n+m] = Ru

# Build block-diagonal Q for IPM
# One block per interior vertex
Q_blocks = [copy(Q_block) for _ in 1:n_interior]

# Build BlockSparseMatrix for Q (block-diagonal)
I_Q = Int[]
J_Q = Int[]
V_Q = Matrix{Float64}[]

for v in 1:n_interior
    push!(I_Q, v)
    push!(J_Q, v)
    push!(V_Q, Q_blocks[v])
end

col_dims = fill(n+m, n_interior)
row_dims = fill(n+m, n_interior)
Q_ipm = blocksparse(I_Q, J_Q, V_Q, row_dims, col_dims)

# Build BlockSparseMatrix for B (reduced coboundary)
# d_int has rows = number of edges × edge_dim, cols = sum of interior stalks
# We need to convert to block structure

# The edges are:
#   edge 1: (1,2) - connects boundary to interior
#   edges 2..k: (t+1, t+2) for t=1..k-1 - connects interior to interior
#   edge k+1: (k+1, k+2) - connects interior to boundary

# Each edge has dimension n (the dynamics constraint dimension)
n_edges = k + 1
edge_dim = n

# For block structure, we need to identify which interior vertices each edge touches
# and build the appropriate blocks

# Actually, let's just convert d_int to a single-row-block BlockSparseMatrix
# where each column block corresponds to an interior vertex

# Simpler approach: use dense blocks for now
B_dense = Matrix(d_int)

# Build as BlockSparseMatrix with one row block and n_interior column blocks
# Actually the IPM expects B to have block structure matching the cones
# Let's check what structure is expected...

# From the IPM code: nvtxs(B) == nvtxs(Q) == length(K)
# So B should have the same column block structure as Q

# Let's build B with n_interior column blocks (one per interior vertex)
# and a single row block (all constraints together)

I_B = Int[]
J_B = Int[]
V_B = Matrix{Float64}[]

for v in 1:n_interior
    block = B_dense[:, (v-1)*(n+m)+1 : v*(n+m)]
    if norm(block) > 1e-15  # only add non-zero blocks
        push!(I_B, 1)  # single row block
        push!(J_B, v)
        push!(V_B, block)
    end
end

B_row_dims = [size(B_dense, 1)]  # single row block
B_col_dims = fill(n+m, n_interior)
B_ipm = blocksparse(I_B, J_B, V_B, B_row_dims, B_col_dims)

# Linear term c
# In public ordering: f = [f_x; f_u] where f_x has (k+1)*n and f_u has k*m
# We need to reorder to internal ordering for interior vertices

# f from lqr_objective is typically zero for standard LQR (no linear term)
# Let's check
println("Linear term f norm: ", norm(f))

# Reorder f to interior ordering
# Public: [x₁, x₂, ..., x_{k+1}, u₁, u₂, ..., uₖ]
# Interior internal: [(x₁,u₁), (x₂,u₂), ..., (xₖ,uₖ)]

c_interior = zeros(n_interior * (n+m))
for t in 1:k
    # xₜ is at public index t, goes to interior vertex t, first n components
    # uₜ is at public index (k+1)*n + t, goes to interior vertex t, last m components
    x_pub_start = (t-1)*n + 1
    u_pub_start = (k+1)*n + (t-1)*m + 1

    int_start = (t-1)*(n+m) + 1
    c_interior[int_start : int_start+n-1] = f[x_pub_start : x_pub_start+n-1]
    c_interior[int_start+n : int_start+n+m-1] = f[u_pub_start : u_pub_start+m-1]
end

# Cones: one CofreeCone per interior vertex
K_ipm = [CofreeCone() for _ in 1:n_interior]

# Create and solve IPM problem
prob = IPMProblem(Q_ipm, B_ipm, c_interior, vec(g_reduced), K_ipm)

println("\nSolving with IPM...")
result = solve(prob; verbose=true, feas_tol=1e-10, gap_tol=1e-10)

println("\nIPM status: ", result.status)
println("IPM iterations: ", result.niter)

# Extract solution and convert back to public ordering
p_interior = result.p

# Reconstruct full solution in public ordering
z_opt2_pub = zeros((k+1)*n + k*m)

# Boundary values
z_opt2_pub[1:n] = x1  # x₁
z_opt2_pub[k*n+1:(k+1)*n] = xk1  # x_{k+1}

# Interior values: extract from p_interior
for t in 1:k
    int_start = (t-1)*(n+m) + 1
    x_t = p_interior[int_start : int_start+n-1]
    u_t = p_interior[int_start+n : int_start+n+m-1]

    # xₜ goes to public index t (but x₁ is at boundary, so states 2..k go to indices 2..k)
    # Actually, public ordering is [x₁, x₂, ..., x_{k+1}, u₁, ..., uₖ]
    # Interior vertex t+1 stores (xₜ, uₜ)
    # So xₜ from interior should go to... wait, I need to be more careful

    # Interior vertex 2 stores (x₁, u₁) in internal ordering
    # But in public ordering x₁ is the first state
    # So x₁ from interior vertex 2 should go to public position 1

    # Hmm, but x₁ is also at boundary vertex 1... this is the "dummy" vertex
    # The dynamics are: x₁_dummy = x₁_traj (edge 1-2 enforces this)

    # So interior vertex 2 has the "real" x₁ and u₁
    # Interior vertex t+1 has xₜ and uₜ

    x_pub_idx = t  # xₜ at public position t
    u_pub_idx = t  # uₜ at public control position t

    z_opt2_pub[(x_pub_idx-1)*n+1 : x_pub_idx*n] = x_t
    z_opt2_pub[(k+1)*n + (u_pub_idx-1)*m+1 : (k+1)*n + u_pub_idx*m] = u_t
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
