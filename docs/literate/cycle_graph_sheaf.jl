# # Cellular sheaf on a cycle graph
#
# This example builds a cycle graph, constructs a simple Euclidean sheaf on it using
# `sheaf_from_graph`, and demonstrates computing a nearest global section.

using CellularSheaves
using LinearAlgebra
using BlockArrays

# the V map should go C0(F) -> C0(G) and the E map should go C1(F) -> C1(G)
struct ComplexMorphism
  V::AbstractMatrix
  E::AbstractMatrix
end

# A small, structured specification type for describing a morphism from
# sheaf F -> sheaf G. This groups the mapping indices and per-stalk matrices
# into a single value that can be passed to constructors or stored.
struct SheafMorphism
  Vmap::Vector{Int}    # maps vertex indices of F -> vertex indices of G
  Emap::Vector{Int}    # maps edge indices of F -> edge indices of G
  Vmaps::Vector       # per-vertex matrices (one per source vertex)
  Emaps::Vector       # per-edge matrices (one per source edge)
end

# Pretty-print for the SheafMorphism object to aid interactive use and docs.
function Base.show(io::IO, spec::SheafMorphism)
  println(io, "SheafMorphism: $(length(spec.Vmap)) vertex mappings, $(length(spec.Emap)) edge mappings")
  println(io, "  Vmap (first 10): ", spec.Vmap[1:min(end,10)])
  println(io, "  Emap (first 10): ", spec.Emap[1:min(end,10)])
end

function is_morphism(cm::ComplexMorphism, dF::AbstractMatrix, dG::AbstractMatrix, tol=1e-8)
  norm(cm.E*dF - dG*cm.V) < tol
end

function  is_morphism(cm, F, G, tol=1e-8)
  is_morphism(cm, coboundary_map(F), coboundary_map(G), tol)
end

# We can construct a morphism of cellular sheaf complexes by building compatible vertex and edge maps. This is a morphism of chain complexes, so the edge map must satisfy the commutativity condition `E * dF = dG * V`. We can check this condition using `is_morphism`.  

# One way to specify a morphism is to build the vertex and edge maps as BlockArrays, where each block corresponds to a vertex or edge stalk. We can then check the morphism condition by multiplying these block matrices with the coboundary maps.

"""
    ComplexMorphism(F, G, Vmap::AbstractVector{<:Integer}, Emap::AbstractVector{<:Integer},
             Vmaps::AbstractVector, Emaps::AbstractVector)

Construct a ComplexMorphism between sheaf complexes F -> G.

Arguments:
- F, G : EuclideanSheaf (source and target)
- Vmap : Vector{Int} of length n_vertices(F). Vmap[i] == j means the i-th vertex stalk of F
     is mapped into the j-th vertex stalk of G.
- Emap : Vector{Int} of length n_edges(F). Emap[k] == ℓ means the k-th edge stalk of F
     is mapped into the ℓ-th edge stalk of G.
- Vmaps: Vector of matrices, one per source vertex stalk i, with size
     (vertex_stalks(G)[Vmap[i]], vertex_stalks(F)[i]).
- Emaps: Vector of matrices, one per source edge stalk k, with size
     (edge_stalks(G)[Emap[k]], edge_stalks(F)[k]).

Returns ComplexMorphism(V, E) where V and E are BlockArrays with the implied block partitions.
"""
function ComplexMorphism(F, G, Vmap::AbstractVector{<:Integer}, Emap::AbstractVector{<:Integer}, Vmaps::AbstractVector, Emaps::AbstractVector)
  vdimsF = vertex_stalks(F)  # vertex block dimensions (in the order used elsewhere in the repo)
  vdimsG = vertex_stalks(G)
  nvF = length(vdimsF)
  nvG = length(vdimsG)
  @assert length(Vmap) == nvF "Vmap must have length equal to number of vertex stalks in F"
  @assert length(Vmaps) == nvF "Vmaps must have one matrix per source vertex stalk"

  Vfull = zeros(sum(vdimsG), sum(vdimsF))  # build zero block matrix for vertex map: partitioned (vdimsG, vdimsF)
  V = BlockArray(Vfull, vdimsG, vdimsF)

  for i in 1:nvF
    tgt = Vmap[i]
    @assert 1 <= tgt <= nvG "Vmap[$i] = $tgt out of range for target vertices"
    A = Vmaps[i]
    expected = (vdimsG[tgt], vdimsF[i])
    @assert size(A) == expected "Vmaps[$i] must have size $(expected), got $(size(A))"
    blocks(V)[tgt, i] = A
  end

  edimsF = collect(values(edge_stalks(F)))  # edge block dimensions; following existing examples use collect(values(...))
  edimsG = collect(values(edge_stalks(G)))
  neF = length(edimsF)
  neG = length(edimsG)
  @assert length(Emap) == neF "Emap must have length equal to number of edge stalks in F"
  @assert length(Emaps) == neF "Emaps must have one matrix per source edge stalk"

  Efull = zeros(sum(edimsG), sum(edimsF))
  E = BlockArray(Efull, edimsG, edimsF)

  for k in 1:neF
    tgt = Emap[k]
    @assert 1 <= tgt <= neG "Emap[$k] = $tgt out of range for target edges"
    B = Emaps[k]
    expected = (edimsG[tgt], edimsF[k])
    @assert size(B) == expected "Emaps[$k] must have size $(expected), got $(size(B))"
    blocks(E)[tgt, k] = B
  end

  ComplexMorphism(V, E)
end

# This constructor implements thr functor from SheafMorphisms to ComplexMorphisms. 
# These are two categories on the same objects. A `SheafMorphism` is a more structured way to specify the mapping data, while a `ComplexMorphism` is the resulting pair of block matrices. By defining this constructor, we can build a `ComplexMorphism` directly from a `SheafMorphism` specification.
function ComplexMorphism(F, G, spec::SheafMorphism)
  return ComplexMorphism(F, G, spec.Vmap, spec.Emap, spec.Vmaps, spec.Emaps)
end

# Backwards-compatible helper: build a `SheafMorphism` from individual
# arguments and delegate to the main constructor above. This lets existing code
# continue to call the old form while preferring the new structured type.
function make_sheaf_morphism_spec(Vmap::AbstractVector{<:Integer}, Emap::AbstractVector{<:Integer},
                                 Vmaps::AbstractVector, Emaps::AbstractVector)
  return SheafMorphism(collect(Int, Vmap), collect(Int, Emap), collect(Vmaps), collect(Emaps))
end

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