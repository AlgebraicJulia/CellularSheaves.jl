# Profile make_kkt internals

using CellularSheaves
using CellularSheaves.NetworkSheaves: coboundary_map, vertex_stalks
using CellularSheaves.IPM
using CellularSheaves.IPM: weightedgraph, symbolic, selectvtxs, FChordalTriangular, UzawaWorkspace
using LinearAlgebra
using SparseArrays
using BenchmarkTools

# Setup (same as comparison test)
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

# Build reduced problem
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
    block_data = B_dense[:, (v-1)*(n+m)+1 : v*(n+m)]
    if norm(block_data) > 1e-15
        push!(I_B, 1)
        push!(J_B, v)
        push!(V_B, block_data)
    end
end
B_row_dims = [size(B_dense, 1)]
B_col_dims = fill(n+m, n_interior)
B_ipm = blocksparse(I_B, J_B, V_B, B_row_dims, B_col_dims)

println("="^60)
println("Profiling make_kkt internals")
println("="^60)

# Now profile each step of make_kkt
println("\nB_ipm size: ", size(B_ipm))
println("B_ipm structure: ", length(I_B), " non-zero blocks out of ", n_interior)

# Step 1: weightedgraph
println("\n1. weightedgraph(B):")
@btime weightedgraph($B_ipm)
weights, graph = weightedgraph(B_ipm)
println("   graph vertices: ", length(graph))

# Step 2: symbolic
println("\n2. symbolic(weights, graph):")
@btime symbolic($weights, $graph)
R, P, S = symbolic(weights, graph)
println("   R.perm length: ", length(R.perm))

# Step 3: selectvtxs
println("\n3. selectvtxs(B, R.perm):")
@btime selectvtxs($B_ipm, $(R.perm))
B_perm = selectvtxs(B_ipm, R.perm)

# Step 4: FChordalTriangular
println("\n4. FChordalTriangular{:N, :L, Float64, Int}(S):")
@btime FChordalTriangular{:N, :L, Float64, Int}($S)
F_ct = FChordalTriangular{:N, :L, Float64, Int}(S)

# Step 5: B' * B
println("\n5. B' * B:")
@btime $(B_perm)' * $B_perm
L = B_perm' * B_perm

# Step 6: UzawaWorkspace
println("\n6. UzawaWorkspace(F, L, B):")
@btime UzawaWorkspace($F_ct, $L, $B_perm)

# Total
println("\n" * "="^60)
println("Total make_kkt:")
@btime IPM.make_kkt(UzawaSettings{Float64}(), $B_ipm)
