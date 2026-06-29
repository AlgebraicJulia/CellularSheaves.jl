# Sweep raug parameter for IPM

using CellularSheaves
using CellularSheaves.NetworkSheaves: coboundary_map, vertex_stalks
using CellularSheaves.IPM
using LinearAlgebra
using SparseArrays

# Setup (same as before)
Ac = [0.0  1.0  0.0  0.0;
      0.0  0.0  0.0  0.0;
      0.0  0.0  0.0  1.0;
      0.0  0.0  0.0  0.0]

Bc = [0.0  0.0;
      1.0  0.0;
      0.0  0.0;
      0.0  1.0]

h = 0.25
k = 100

F  = EuclideanSheaf{Float64}([2, 2])
ts = ControlledTrajectorySheaf(F, Ac, Bc, h, k)

n = ts.state_dim
m = ts.control_dim

x1  = [0.0, 0.0, 2.0, 0.0]
xk1 = [1.0, 0.0, 3.0, 0.0]

Q_lqr  = Matrix{Float64}(I, n, n)
Ru = Matrix{Float64}(I, m, m)
Qf = 10.0 * Matrix{Float64}(I, n, n)

H, f, _ = lqr_objective(ts, Q_lqr, Ru; Qf=Qf)

# Build reduced problem (reusing code from v2)
sheaf = ts.sheaf
stalks = vertex_stalks(sheaf)

interior_verts = 2:k+1
n_interior = length(interior_verts)

d_full = coboundary_map(sheaf)
cumstalk = [0; cumsum(stalks)]
col_range(v) = cumstalk[v]+1 : cumstalk[v+1]

bnd_cols = vcat(collect(col_range(1)), collect(col_range(k+2)))
int_cols = vcat([collect(col_range(v)) for v in interior_verts]...)

d_full_sparse = sparse(d_full)
d_int = d_full_sparse[:, int_cols]
d_bnd = d_full_sparse[:, bnd_cols]

x_bnd = vcat(x1, xk1)
g_reduced = -d_bnd * x_bnd

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

B_dense = Matrix(d_int)
I_B = Int[]
J_B = Int[]
V_B = Matrix{Float64}[]
for v in 1:n_interior
    block = B_dense[:, (v-1)*(n+m)+1 : v*(n+m)]
    if norm(block) > 1e-15
        push!(I_B, 1)
        push!(J_B, v)
        push!(V_B, block)
    end
end
B_row_dims = [size(B_dense, 1)]
B_col_dims = fill(n+m, n_interior)
B_ipm = blocksparse(I_B, J_B, V_B, B_row_dims, B_col_dims)

c_interior = zeros(n_interior * (n+m))
K_ipm = [CofreeCone() for _ in 1:n_interior]

prob = IPMProblem(Q_ipm, B_ipm, c_interior, vec(g_reduced), K_ipm)

# Sweep raug
println("="^70)
println("Sweeping raug parameter")
println("="^70)

for raug in [1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8]
    println("\n--- raug = $raug ---")
    result = solve(prob; verbose=true, feas_tol=1e-10, gap_tol=1e-10,
                   kkt=UzawaSettings{Float64}(raug=raug))
    println("Status: ", result.status, ", Iterations: ", result.niter)
end
