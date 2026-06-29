# Profile B' * B - understand structure

using CellularSheaves
using CellularSheaves.NetworkSheaves: coboundary_map, vertex_stalks
using CellularSheaves.IPM
using CellularSheaves.BlockSparseArrays: matprod_dest, gemm_impl!, nvtxs, narcs, nouts, vtxs, srcrange, ncols, nrows, block
using LinearAlgebra
using SparseArrays
using BenchmarkTools

# Setup
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

Bt = convert(BlockSparseMatrix, B_ipm')
C = matprod_dest(Bt, B_ipm, Float64)

println("="^60)
println("Block structure analysis")
println("="^60)

println("\nB (constraint matrix):")
println("  size: ", size(B_ipm))
println("  nvtxs (col blocks): ", nvtxs(B_ipm))
println("  nouts (row blocks): ", nouts(B_ipm))
println("  narcs (non-zero blocks): ", narcs(B_ipm))

println("\nBt (transposed):")
println("  size: ", size(Bt))
println("  nvtxs (col blocks): ", nvtxs(Bt))
println("  nouts (row blocks): ", nouts(Bt))
println("  narcs (non-zero blocks): ", narcs(Bt))

println("\nC = Bt * B (result):")
println("  size: ", size(C))
println("  nvtxs (col blocks): ", nvtxs(C))
println("  nouts (row blocks): ", nouts(C))
println("  narcs (non-zero blocks): ", narcs(C))

# Count mul! calls
mul_count = let cnt = 0
    for Cv in vtxs(C)
        for Be in srcrange(B_ipm, Cv)
            Bu = B_ipm.tgt[Be]
            for Ae in srcrange(Bt, Bu)
                cnt += 1
            end
        end
    end
    cnt
end
println("  Number of mul! calls: ", mul_count)

# Check block sizes
println("\nBlock sizes in C:")
let first_few = 0
    for Cv in vtxs(C)
        for Ce in srcrange(C, Cv)
            Cu = C.tgt[Ce]
            blk = block(C, Cu, Cv, Ce)
            if first_few < 5
                println("  Block ($Cu, $Cv): ", size(blk))
                first_few += 1
            end
        end
        first_few >= 5 && break
    end
end

# Time a single block multiply
println("\n" * "="^60)
println("Time for single block mul!")
println("="^60)

# Get one block from Bt and B
Bt_blk = block(Bt, 1, 1, 1)  # First block of Bt
B_blk = block(B_ipm, 1, 1, 1)  # First block of B
C_blk = similar(Bt_blk, size(Bt_blk, 1), size(B_blk, 2))
fill!(C_blk, 0.0)

println("Bt block size: ", size(Bt_blk))
println("B block size: ", size(B_blk))
println("C block size: ", size(C_blk))

println("\nSingle mul!(C_blk, Bt_blk, B_blk, 1.0, 1.0):")
@btime mul!($C_blk, $Bt_blk, $B_blk, 1.0, 1.0)

println("\nEstimated total from $(mul_count) mul! calls:")
t_single = @belapsed mul!($C_blk, $Bt_blk, $B_blk, 1.0, 1.0)
println("  $(mul_count) × $(round(t_single*1e6, digits=2)) μs = $(round(mul_count * t_single * 1000, digits=2)) ms")
