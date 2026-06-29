# Profile B' * B internals

using CellularSheaves
using CellularSheaves.NetworkSheaves: coboundary_map, vertex_stalks
using CellularSheaves.IPM
using CellularSheaves.BlockSparseArrays: matprod_dest, gemm_impl!
using LinearAlgebra
using SparseArrays
using BenchmarkTools

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
println("Profiling B' * B internals")
println("="^60)
println("B size: ", size(B_ipm))

# Step 1: convert(BlockSparseMatrix, B') - transpose copy
println("\n1. convert(BlockSparseMatrix, B') [transpose copy]:")
@btime convert(BlockSparseMatrix, $(B_ipm)')
Bt = convert(BlockSparseMatrix, B_ipm')
println("   Bt size: ", size(Bt))

# Step 2: matprod_dest - allocate result with sparsity pattern
println("\n2. matprod_dest(Bt, B, T) [allocate result]:")
@btime matprod_dest($Bt, $B_ipm, Float64)
C = matprod_dest(Bt, B_ipm, Float64)
println("   C size: ", size(C))

# Step 3: gemm_impl! - actual multiplication
println("\n3. gemm_impl!(C, Val(:N), Bt, Val(:N), B, 1, 0) [multiply]:")
@btime gemm_impl!($C, Val(:N), $Bt, Val(:N), $B_ipm, 1.0, 0.0)

# Total
println("\n" * "="^60)
println("Total B' * B:")
@btime $(B_ipm)' * $B_ipm
