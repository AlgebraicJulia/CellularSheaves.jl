# # Sheaf Morphisms
#
# This example demonstrates constructing sheaf morphisms between Euclidean sheaves
# built over a cycle graph, including checking the naturality condition and composing morphisms.
#
# The types `SheafMorphism` and `ComplexMorphism`, together with
# `is_morphism` and `compose`, are part of the `CellularSheaves` package — they
# do not need to be defined here.

using CellularSheaves
using LinearAlgebra

# ## Build a cycle sheaf by adding stalks and edges
# Instead of first creating a Graph and calling `sheaf_from_graph`, we instantiate an
# empty Euclidean sheaf and populate it using `add_vertex_stalk!` and `add_sheaf_edge!`.

n = 6
stalk_dim = 3

# Start with an empty Euclidean sheaf (no vertex stalks)
F = EuclideanSheaf{Float64}(Int[])

# Add `n` vertex stalks of dimension `stalk_dim`.
for i in 1:n
    add_vertex_stalk!(F, stalk_dim)
end

# Add edges to form a cycle. Each edge has a 1-dimensional edge stalk; restriction
# maps are 1 x stalk_dim. Here we use identical maps on both endpoints for simplicity.
for i in 1:n
    v1 = i
    v2 = (i % n) + 1
    rm = Matrix{Float64}(I, stalk_dim, stalk_dim) # Identity map as restriction map
    add_sheaf_edge!(F, v1, v2, rm, rm)
end

println(F)

# ## Inspect the sheaf Laplacian and nearest global section

L = sheaf_laplacian_matrix(F)
println("Sheaf Laplacian is a ", size(L, 1), " x ", size(L, 2), " matrix.")

# Start from a random 0-cochain and project to the nearest global section
x0 = rand(sum(vertex_stalks(F)))
try
    gs = nearest_global_section(F, x0)
    println("Norm of coboundary on nearest global section: ", norm(coboundary_map(F) * gs))
    gs
catch e
    println("Could not compute nearest global section during docs build: ", e)
end

# The printed zero (or tiny numerical value) above indicates `gs` is close to a true global section.


# Start with an empty Euclidean sheaf (no vertex stalks)
G = EuclideanSheaf{Float64}(Int[])

# Add `n` vertex stalks of dimension `stalk_dim`.
for i in 1:n
    add_vertex_stalk!(G, stalk_dim)
end

# Add edges to form a cycle. Each edge has a 1-dimensional edge stalk; restriction
# maps are 1 x stalk_dim. Here we use identical maps on both endpoints for simplicity.
for i in 1:n
    v1 = i
    v2 = (i % n) + 1
    rm = Matrix{Float64}(I, stalk_dim - 1, stalk_dim)
    add_sheaf_edge!(G, v1, v2, rm, rm)
end

println(G)

Vmap = collect(1:n)
Emap = collect(1:n)
Vmaps = [Matrix{Float64}(I, stalk_dim, stalk_dim) for i in 1:n]
Emaps = [Matrix{Float64}(I, stalk_dim - 1, stalk_dim) for i in 1:n]

spec = SheafMorphism(Vmap, Emap, Vmaps, Emaps)
cm = ComplexMorphism(F, G, spec)

is_morphism(cm, F, G)

# The morphism condition holds, so the coboundary maps commute with the morphism matrices.

dF = coboundary_map(F)
dG = coboundary_map(G)

# Build a third sheaf H on the same graph topology but with 1-dimensional edge stalks.
# Vertex stalks match G (stalk_dim each). Edge restriction maps are 1 x stalk_dim
# projection matrices that pick the first coordinate.
H = EuclideanSheaf{Float64}(Int[])
for i in 1:n
  add_vertex_stalk!(H, stalk_dim)
end
for i in 1:n
  v1 = i
  v2 = (i % n) + 1
  rmH = zeros(1, stalk_dim)  # restriction: project vertex R^d -> R^1 by taking first coordinate
  rmH[1,1] = 1.0
  add_sheaf_edge!(H, v1, v2, rmH, rmH)
end

println(H)

# Build a SheafMorphism spec from G -> H. Vertex maps are identity (3x3),
# edge maps project G's edge stalk (size = stalk_dim - 1) to H's edge stalk (1)
Vmap_GH = collect(1:n)
Emap_GH = collect(1:n)
Vmaps_GH = [Matrix{Float64}(I, stalk_dim, stalk_dim) for i in 1:n]
edimsG = collect(values(edge_stalks(G)))
Emaps_GH = [zeros(1, edimsG[i]) for i in 1:length(edimsG)]
for i in 1:length(edimsG)
  Emaps_GH[i][1,1] = 1.0  # project to first coordinate of G's edge stalk
end

specGH = make_sheaf_morphism_spec(Vmap_GH, Emap_GH, Vmaps_GH, Emaps_GH)
cmGH = ComplexMorphism(G, H, specGH)

is_morphism(cmGH, G, H)

dH = coboundary_map(H)


compose(f::SheafMorphism, g::SheafMorphism) = SheafMorphism(
  [g.Vmap[i] for i in f.Vmap],
  [g.Emap[k] for k in f.Emap],
  [g.Vmaps[i] * f.Vmaps[g.Vmap[i]] for i in 1:length(f.Vmaps)],
  [g.Emaps[k] * f.Emaps[g.Emap[k]] for k in 1:length(f.Emaps)],
)

compose(f::ComplexMorphism, g::ComplexMorphism) = ComplexMorphism(
  g.V * f.V,
  g.E * f.E,
)

composite_morphism = ComplexMorphism(F, H, compose(spec, specGH)) 

# The V component is 

composite_morphism.V

# The E component is
composite_morphism.E

# The composite is supposed to be the same as the result of directly multiplying the ComplexMorphism matrices:

ϕ = compose(cm, cmGH)
norm(ϕ.V - composite_morphism.V)

# and the E component is

norm(ϕ.E - composite_morphism.E)