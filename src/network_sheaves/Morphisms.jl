"""
Module for morphisms of cellular sheaf complexes.

A *sheaf morphism* `F -> G` consists of, for each vertex `v`, a linear map
`F(v) -> G(φ(v))` and, for each edge `e`, a linear map `F(e) -> G(ψ(e))`,
subject to the naturality condition `dG ∘ V = E ∘ dF` where `dF`, `dG` are the
coboundary maps of `F` and `G` respectively.

This module provides two representations:

- `SheafMorphism` — a *structured* specification: per-stalk matrices together with
  index maps that say which stalk of the target each source stalk maps into.  This
  is the natural form when building morphisms by hand.

- `ComplexMorphism` — the *matrix* representation: a pair `(V, E)` of block matrices
  `V : C⁰(F) → C⁰(G)` and `E : C¹(F) → C¹(G)`.  This is the form needed for
  linear-algebraic operations.

The functor `ComplexMorphism(F, G, spec::SheafMorphism)` assembles the block matrices
from the structured spec.
"""
module SheafMorphisms

export SheafMorphism, ComplexMorphism, is_morphism, compose, make_sheaf_morphism_spec, id

using LinearAlgebra
using BlockArrays
using Graphs

using ..SheafInterface: vertex_stalks, edge_stalks, edge_stalk_dimensions, coboundary_map,
                        underlying_graph, get_edge_stalk, get_restriction_map

# ---------------------------------------------------------------------------
# Types
# ---------------------------------------------------------------------------

"""
    ComplexMorphism

A morphism of cochain complexes `C*(F) -> C*(G)` represented as a pair of
matrices:

- `V :: AbstractMatrix` — the degree-0 component `C⁰(F) → C⁰(G)`.
- `E :: AbstractMatrix` — the degree-1 component `C¹(F) → C¹(G)`.

The naturality condition is `E * dF == dG * V` where `dF`, `dG` are the
coboundary maps.  Use [`is_morphism`](@ref) to check this numerically.
"""
struct ComplexMorphism
    V::AbstractMatrix
    E::AbstractMatrix
end

"""
    SheafMorphism

A structured specification of a cellular-sheaf morphism `F -> G`.

Fields
------
- `Vmap :: Vector{Int}` — `Vmap[i] = j` means the `i`-th vertex stalk of `F`
  maps *into* the `j`-th vertex stalk of `G`.
- `Emap :: Vector{Int}` — `Emap[k] = ℓ` means the `k`-th edge stalk of `F`
  maps into the `ℓ`-th edge stalk of `G`.
- `Vmaps :: Vector` — one matrix per source vertex stalk; the matrix at index
  `i` has size `(dim G_vertex[Vmap[i]], dim F_vertex[i])`.
- `Emaps :: Vector` — one matrix per source edge stalk; the matrix at index
  `k` has size `(dim G_edge[Emap[k]], dim F_edge[k])`.

Use [`ComplexMorphism(F, G, spec)`](@ref) to convert to block-matrix form.
"""
struct SheafMorphism
    Vmap::Vector{Int}
    Emap::Vector{Int}
    Vmaps::Vector
    Emaps::Vector
end

# ---------------------------------------------------------------------------
# Printing
# ---------------------------------------------------------------------------

function Base.show(io::IO, spec::SheafMorphism)
    println(io, "SheafMorphism: $(length(spec.Vmap)) vertex mappings, $(length(spec.Emap)) edge mappings")
    println(io, "  Vmap (first 10): ", spec.Vmap[1:min(length(spec.Vmap), 10)])
    println(io, "  Emap (first 10): ", spec.Emap[1:min(length(spec.Emap), 10)])
end

# ---------------------------------------------------------------------------
# Constructors
# ---------------------------------------------------------------------------

"""
    ComplexMorphism(F, G, Vmap, Emap, Vmaps, Emaps)

Assemble a [`ComplexMorphism`](@ref) `F -> G` from explicit index maps and
per-stalk matrices.

Arguments
---------
- `F`, `G` : source and target `AbstractNetworkSheaf`.
- `Vmap :: AbstractVector{<:Integer}` — vertex index map (length = `nv(F)`).
- `Emap :: AbstractVector{<:Integer}` — edge index map (length = `ne(F)`).
- `Vmaps :: AbstractVector` — per-source-vertex matrices.
- `Emaps :: AbstractVector` — per-source-edge matrices.

Returns a `ComplexMorphism(V, E)` where `V` and `E` are `BlockArray`s with the
block partitions induced by the stalk dimensions of `F` and `G`.
"""
function ComplexMorphism(F, G,
                         Vmap::AbstractVector{<:Integer},
                         Emap::AbstractVector{<:Integer},
                         Vmaps::AbstractVector,
                         Emaps::AbstractVector)
    vdimsF = vertex_stalks(F)
    vdimsG = vertex_stalks(G)
    nvF = length(vdimsF)
    nvG = length(vdimsG)
    @assert length(Vmap) == nvF "Vmap must have length equal to number of vertex stalks in F"
    @assert length(Vmaps) == nvF "Vmaps must have one matrix per source vertex stalk"

    T = if !isempty(Vmaps) || !isempty(Emaps)
        promote_type(
            isempty(Vmaps) ? Float64 : mapreduce(eltype, promote_type, Vmaps),
            isempty(Emaps) ? Float64 : mapreduce(eltype, promote_type, Emaps),
        )
    else
        Float64
    end
    Vfull = zeros(T, sum(vdimsG), sum(vdimsF))
    V = BlockArray(Vfull, vdimsG, vdimsF)
    for i in 1:nvF
        tgt = Vmap[i]
        @assert 1 <= tgt <= nvG "Vmap[$i] = $tgt out of range for target vertices"
        A = Vmaps[i]
        expected = (vdimsG[tgt], vdimsF[i])
        @assert size(A) == expected "Vmaps[$i] must have size $(expected), got $(size(A))"
        blocks(V)[tgt, i] = A
    end

    edimsF = edge_stalk_dimensions(F)
    edimsG = edge_stalk_dimensions(G)
    neF = length(edimsF)
    neG = length(edimsG)
    @assert length(Emap) == neF "Emap must have length equal to number of edge stalks in F"
    @assert length(Emaps) == neF "Emaps must have one matrix per source edge stalk"

    Efull = zeros(T, sum(edimsG), sum(edimsF))
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

"""
    ComplexMorphism(F, G, spec::SheafMorphism)

Convert a [`SheafMorphism`](@ref) spec into a [`ComplexMorphism`](@ref).

This is the functor from the category of structured sheaf-morphism
specifications to the category of cochain-complex morphisms on the same
objects.
"""
function ComplexMorphism(F, G, spec::SheafMorphism)
    return ComplexMorphism(F, G, spec.Vmap, spec.Emap, spec.Vmaps, spec.Emaps)
end

# ---------------------------------------------------------------------------
# Helper constructor
# ---------------------------------------------------------------------------

"""
    make_sheaf_morphism_spec(Vmap, Emap, Vmaps, Emaps) -> SheafMorphism

Convenience constructor: normalises argument types and returns a
[`SheafMorphism`](@ref).
"""
function make_sheaf_morphism_spec(Vmap::AbstractVector{<:Integer},
                                  Emap::AbstractVector{<:Integer},
                                  Vmaps::AbstractVector,
                                  Emaps::AbstractVector)
    return SheafMorphism(collect(Int, Vmap), collect(Int, Emap),
                         collect(Vmaps), collect(Emaps))
end

# ---------------------------------------------------------------------------
# Validity check
# ---------------------------------------------------------------------------

"""
    is_morphism(cm::ComplexMorphism, dF, dG; tol=1e-8) -> Bool

Return `true` if the naturality square `E * dF ≈ dG * V` holds within `tol`.
"""
function is_morphism(cm::ComplexMorphism, dF::AbstractMatrix, dG::AbstractMatrix;
                     tol=1e-8)
    return norm(cm.E * dF - dG * cm.V) < tol
end

"""
    is_morphism(cm, F, G; tol=1e-8) -> Bool

Convenience overload: computes `coboundary_map(F)` and `coboundary_map(G)` and
delegates to the matrix form.
"""
function is_morphism(cm, F, G; tol=1e-8)
    return is_morphism(cm, coboundary_map(F), coboundary_map(G); tol=tol)
end

# ---------------------------------------------------------------------------
# Composition
# ---------------------------------------------------------------------------

"""
    compose(f::SheafMorphism, g::SheafMorphism) -> SheafMorphism

Compose two compatible `SheafMorphism`s `f : F -> G` and `g : G -> H` to
produce `g ∘ f : F -> H`.
"""
function compose(f::SheafMorphism, g::SheafMorphism)
    return SheafMorphism(
        [g.Vmap[i] for i in f.Vmap],
        [g.Emap[k] for k in f.Emap],
        [g.Vmaps[f.Vmap[i]] * f.Vmaps[i] for i in 1:length(f.Vmaps)],
        [g.Emaps[f.Emap[k]] * f.Emaps[k] for k in 1:length(f.Emaps)],
    )
end

"""
    compose(f::ComplexMorphism, g::ComplexMorphism) -> ComplexMorphism

Compose two compatible `ComplexMorphism`s `f : C*(F) -> C*(G)` and
`g : C*(G) -> C*(H)` by matrix multiplication.
"""
function compose(f::ComplexMorphism, g::ComplexMorphism)
    return ComplexMorphism(g.V * f.V, g.E * f.E)
end

# ---------------------------------------------------------------------------
# Identity morphisms
# ---------------------------------------------------------------------------

# Helper: infer element type from the first available restriction map.
function _sheaf_scalar_type(X)
    g = underlying_graph(X)
    for e in edges(g)
        return eltype(get_restriction_map(X, src(e), dst(e)))
    end
    return Float64
end

"""
    id(::Type{SheafMorphism}, X) -> SheafMorphism

Return the identity [`SheafMorphism`](@ref) on sheaf `X`.

Each vertex stalk maps to itself via `Vmap[i] = i` with a square identity
matrix, and likewise for every edge stalk.
"""
function id(::Type{SheafMorphism}, X)
    g = underlying_graph(X)
    vdims = vertex_stalks(X)
    edims = edge_stalk_dimensions(X)
    nv = length(vdims)
    ne = length(edims)
    T = _sheaf_scalar_type(X)
    return SheafMorphism(
        collect(1:nv),
        collect(1:ne),
        [Matrix{T}(I, vdims[i], vdims[i]) for i in 1:nv],
        [Matrix{T}(I, edims[k], edims[k]) for k in 1:ne],
    )
end

"""
    id(::Type{ComplexMorphism}, X) -> ComplexMorphism

Return the identity [`ComplexMorphism`](@ref) on sheaf `X`.

`V` is the identity matrix on `C⁰(X)` (total vertex-stalk dimension) and
`E` is the identity matrix on `C¹(X)` (total edge-stalk dimension).
"""
function id(::Type{ComplexMorphism}, X)
    g = underlying_graph(X)
    vdims = vertex_stalks(X)
    edims = edge_stalk_dimensions(X)
    T = _sheaf_scalar_type(X)
    V = Matrix{T}(I, sum(vdims), sum(vdims))
    E = Matrix{T}(I, sum(edims), sum(edims))
    return ComplexMorphism(V, E)
end

end # module SheafMorphisms
