using CellularSheaves
using CliqueTrees.Multifrontal
using LinearAlgebra
using SparseArrays

# F is a cellular sheaf

n = 6; stalk_dim = 3

F = EuclideanSheaf{Float64}(Int[])

for i in 1:n
    add_vertex_stalk!(F, stalk_dim)
end

for i in 1:n
    v1 = i
    v2 = (i % n) + 1
    rm = Matrix{Float64}(I, stalk_dim, stalk_dim)
    add_sheaf_edge!(F, v1, v2, rm, rm)
end

# X is the Laplacian of F

C = sparse(coboundary_map(F))
X = C' * C

# M is an LDLt factorization of X
#
#   X = Pᵀ L D Lᵀ P
#
@time M = ldlt!(ChordalLDLt(X), RowMaximum())

L = M.L # lower triangular factor
D = M.D # diagonal factor
P = M.P # permutation

L

# V is a basis for the nullspace of X

ind = findall(i -> D[i, i] == 0, 1:18)
U = zeros(18, length(ind))

for j in eachindex(ind)
    U[ind[j], j] = 1
end

V = P \ (L' \ U)
