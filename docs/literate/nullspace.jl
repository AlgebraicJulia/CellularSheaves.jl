# # Computing the Nullspace of the Sheaf Laplacian

# ## Introduction
#
# This example demonstrates how to compute and visualize the nullspace of the Laplacian operator
# for a circular cellular sheaf. We construct a system of 6 agents, each with 3-dimensional stalks,
# arranged in a cycle topology. The nullspace of the Laplacian contains the global sections---the
# solutions that are preserved by the sheaf Laplacian operator. Understanding the structure of this
# nullspace is fundamental to understanding the harmonic analysis and signal processing on cellular sheaves.

using CellularSheaves
using LinearAlgebra
using Plots

# ## Building the Sheaf
#
# We start by creating a Euclidean sheaf with 6 vertices (agents), each equipped with a 3-dimensional
# vector space as its stalk. The stalks represent local data or state at each agent.

n = 6
stalk_dim = 3

F = EuclideanSheaf{Float64}(Int[])

for i in 1:n
    add_vertex_stalk!(F, stalk_dim)
end

# Next, we add edges connecting the agents in a cycle. Each edge is equipped with the identity
# restriction map, meaning that sections on adjacent agents must have the same value where they overlap.

for i in 1:n
    v1 = i
    v2 = (i % n) + 1
    rm = Matrix{Float64}(I, stalk_dim, stalk_dim)
    add_sheaf_edge!(F, v1, v2, rm, rm)
end

# ## Computing the Nullspace of the Laplacian
#
# We call [`nullspace_ldlt`](@ref) directly on the sheaf. Internally it forms the
# Laplacian ``X = d^\mathsf{T} d`` from the coboundary map ``d``, factorises it as
# ``X = P^\mathsf{T} L D L^\mathsf{T} P`` via a sparse Chordal LDLt decomposition,
# and returns a basis for the nullspace from the zero-diagonal entries of ``D``.
# The nullspace of the Laplacian is exactly the space of *global sections* of the sheaf.

V = nullspace_ldlt(F)

# ## Visualizing the Nullspace Basis
#
# Each column of V is a basis vector in the nullspace. We visualize all basis vectors simultaneously
# using three orthogonal projections (x vs y, y vs z, x vs z) arranged in a grid. Each basis vector
# is displayed in a distinct color.

num_basis_vectors = size(V, 2)
colors = palette(:tab10, num_basis_vectors)

# Extract coordinate components for each agent
# Coordinates are strided in V with period stalk_dim
x_coords = V[1:stalk_dim:end, :]  # first row in each stalk_dim-sized block
y_coords = V[2:stalk_dim:end, :]  # second row in each stalk_dim-sized block
z_coords = V[3:stalk_dim:end, :]  # third row in each stalk_dim-sized block

# Create three projection subplots
p1 = scatter(legend=false, xlabel="x", ylabel="y", title="x vs y projection",
    markersize=8, markerstrokewidth=1.5)
p2 = scatter(legend=false, xlabel="y", ylabel="z", title="y vs z projection",
    markersize=8, markerstrokewidth=1.5)
p3 = scatter(legend=false, xlabel="x", ylabel="z", title="x vs z projection",
    markersize=8, markerstrokewidth=1.5)

# Plot each basis vector with color coding
for b in 1:num_basis_vectors
    x_vals = x_coords[:, b]
    y_vals = y_coords[:, b]
    z_vals = z_coords[:, b]

    scatter!(p1, x_vals, y_vals, color=colors[b], label="basis $b",
        markersize=6, markerstrokewidth=0.5)

    scatter!(p2, y_vals, z_vals, color=colors[b], label="basis $b",
        markersize=6, markerstrokewidth=0.5)

    scatter!(p3, x_vals, z_vals, color=colors[b], label="basis $b",
        markersize=6, markerstrokewidth=0.5)
end

# Build a legend-only subplot for the bottom-right cell of the 2×2 grid
p4 = plot(; showaxis=false, grid=false, legend=:inside, border=:none, title="Legend")
for b in 1:num_basis_vectors
    scatter!(p4, [], []; color=colors[b], label="basis $b",
        markersize=6, markerstrokewidth=0.5)
end

# Arrange the three projections in a 2×2 grid layout with the legend in the bottom-right cell
plot_grid = @layout [a b; c d]
plot(p1, p2, p3, p4, layout=plot_grid, size=(900, 800))
