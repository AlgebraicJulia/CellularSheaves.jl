# Overlapping Schwarz domain decomposition phrased on a cellular sheaf.
#
# A cover of the unknowns of a discretized PDE by overlapping subdomains
# Ω₁, …, Ω_N determines a network sheaf (the "overlap sheaf" of the cover):
# one vertex per subdomain with stalk ℝ^{Ω_i}, one edge per nonempty pairwise
# overlap with stalk ℝ^{Ω_i ∩ Ω_j}, and coordinate projections as restriction
# maps. A Schwarz iterate is a 0-cochain of this sheaf — every subdomain keeps
# its own copy of the solution — and the method has converged exactly when
# that cochain is a global section that glues to the solution of A u = f.
module SchwarzMethods

export OverlapCover, ghost_layer_cover, Ownership, LocalProblem, RobinFace, SchwarzDecomposition,
    overlap_sheaf, overlapping_subdomains, localize, glue, overlap_disagreement,
    TransmissionCondition, DirichletTransmission, RobinTransmission, optimized_robin_parameter,
    SchwarzSweep, MultiplicativeSweep, MulticolorSweep, ParallelSweep, schwarz_step!,
    AbstractCoarseSpace, TruncatedPushforwardCoarseSpace, ExactPushforwardCoarseSpace,
    coarse_dimension, coarse_correct!,
    SchwarzProblem, SchwarzIteration, SchwarzCG, SchwarzResult, SchwarzPreconditioner, solve,
    SheafADMM, LocalObjective, local_objectives, SchwarzGMRES, SchwarzSweepPreconditioner

using ArgCheck: @argcheck
using BlockArrays: BlockVector, mortar, blocks
using Graphs: SimpleGraph, add_edge!, has_edge, neighbors, ne, nv, vertices, degree
using LinearAlgebra
using LinearAlgebra: ldlt!, RowMaximum
using SparseArrays
using CliqueTrees.Multifrontal: ChordalLDLt
using Krylov: cg, gmres
import CommonSolve: solve

using ..SheafInterface: add_sheaf_edge!
using ..EuclideanSheaves: EuclideanSheaf
using ..RestrictionMaps: SelectionRestriction
using ..GraphHomomorphisms: GraphHomomorphism, fiber_vertices

# ===== Sparse-matrix helpers =====

# Structural neighbours of dof k in the (symmetric) matrix S.
_adjacent(S::SparseMatrixCSC, k::Int) =
    (rowvals(S)[p] for p in nzrange(S, k) if rowvals(S)[p] != k && !iszero(nonzeros(S)[p]))

# Γ: dofs outside `dofs` that S couples to `dofs` (the discrete Dirichlet boundary).
function _boundary(S::SparseMatrixCSC, dofs::Vector{Int})
    inside = falses(size(S, 1))
    inside[dofs] .= true
    Γ = Int[]
    for c in dofs, r in _adjacent(S, c)
        inside[r] || push!(Γ, r)
    end
    return sort!(unique!(Γ))
end

# The sparsity pattern of S + Sᵀ. The ghost layers, overlaps and ownership of a
# cover only depend on which dofs are coupled, so for a nonsymmetric matrix they
# are computed from this symmetrized pattern.
_structure(S::SparseMatrixCSC) = issymmetric(S) ? S : spones(S) + spones(sparse(transpose(S)))

# Local factorization: ChordalLDLt for symmetric matrices (with a positive
# definiteness check), sparse LU otherwise.
function _local_factor(Ai::SparseMatrixCSC, symmetric::Bool, i::Int)
    symmetric || return lu(Ai)
    factor = ldlt!(ChordalLDLt(Ai), RowMaximum(); check=false)
    @argcheck all(>(0), factor.D.diag) "the local matrix of subdomain $i is not positive definite; increase the Robin parameter"
    return factor
end

_factor_solve(M::ChordalLDLt, b::AbstractVector) = _ldlt_solve(M, b)
_factor_solve(M, b::AbstractVector) = M \ b

# Solve M v = b for a ChordalLDLt factor with X = P' L D L' P.
function _ldlt_solve(M, b::AbstractVector)
    c = M.P' \ b
    z = M.L \ c
    w = M.D \ z
    y = M.L' \ w
    return M.P \ y
end

# ===== The cover and its overlap sheaf =====

"""
    OverlapCover(subdomains, ndofs=maximum dof index; ghosts=nothing)
    ghost_layer_cover(A, subdomains) -> OverlapCover

The combinatorics of a cover of the dofs `1:ndofs` by subdomains
``\\Omega_i`` with optional *ghost layers* ``\\Gamma_i``, bundled together. The
vertex stalk of subdomain ``i`` is the *closed* subdomain
``\\overline\\Omega_i = \\Omega_i \\cup \\Gamma_i``. Its interior values are what
the subdomain computes, and its ghost values are copies received from
neighbours.

- `subdomains[i]`: the sorted dofs of ``\\Omega_i`` (the interior);
- `stalks[i]`: the sorted dofs of ``\\overline\\Omega_i``, the coordinates of
  the vertex stalk;
- `interior[i]`: a mask over `stalks[i]`, `true` on ``\\Omega_i``;
- `members[k]`: the pairs `(i, ℓ)` with dof `k` at position `ℓ` of `stalks[i]`;
- `graph`: the 1-skeleton of the nerve of the closed cover (an edge for every
  nonempty ``\\overline\\Omega_i \\cap \\overline\\Omega_j``);
- `overlaps[i => j]`: the positions in `stalks[i]` of
  ``\\overline\\Omega_i \\cap \\overline\\Omega_j``, ordered so that
  `overlaps[i => j]` and `overlaps[j => i]` address the same dofs entry by
  entry. These are the restriction maps of [`overlap_sheaf`](@ref) in index form.

[`ghost_layer_cover`](@ref) takes ``\\Gamma_i`` to be the discrete boundary of
``\\Omega_i`` in the sparsity graph of `A`, the dofs a local solve on
``\\Omega_i`` needs as boundary data. Every such value then lies in an edge stalk
shared with the subdomain that computes it. As a result, each message of a
Schwarz method is a restriction map applied to a vertex stalk.
"""
struct OverlapCover
    subdomains::Vector{Vector{Int}}
    stalks::Vector{Vector{Int}}
    interior::Vector{BitVector}
    members::Vector{Vector{Tuple{Int,Int}}}
    graph::SimpleGraph{Int}
    overlaps::Dict{Pair{Int,Int},Vector{Int}}
end

function OverlapCover(subdomains::AbstractVector{<:AbstractVector{<:Integer}},
                      ndofs::Integer=maximum(d -> isempty(d) ? 0 : maximum(d), subdomains; init=0);
                      ghosts::Union{Nothing,AbstractVector{<:AbstractVector{<:Integer}}}=nothing)
    @argcheck !isempty(subdomains) "need at least one subdomain"
    doms = [sort!(unique(Vector{Int}(d))) for d in subdomains]
    layers = ghosts === nothing ? [Int[] for _ in doms] : [setdiff(Vector{Int}(g), d) for (g, d) in zip(ghosts, doms)]
    @argcheck length(layers) == length(doms) "need one ghost layer per subdomain"
    stalks = [sort!(union(d, g)) for (d, g) in zip(doms, layers)]
    for (i, s) in enumerate(stalks)
        @argcheck isempty(s) || (1 <= first(s) && last(s) <= ndofs) "subdomain $i has dof indices outside 1:$ndofs"
    end
    interior = [BitVector(insorted(k, d) for k in s) for (s, d) in zip(stalks, doms)]
    members = [Tuple{Int,Int}[] for _ in 1:ndofs]
    for (i, s) in enumerate(stalks), (ℓ, k) in enumerate(s)
        push!(members[k], (i, ℓ))
    end
    graph = SimpleGraph(length(stalks))
    overlaps = Dict{Pair{Int,Int},Vector{Int}}()
    for m in members, a in eachindex(m), b in (a + 1):lastindex(m)
        (i, ℓi), (j, ℓj) = m[a], m[b]
        add_edge!(graph, i, j)
        push!(get!(overlaps, i => j, Int[]), ℓi)
        push!(get!(overlaps, j => i, Int[]), ℓj)
    end
    return OverlapCover(doms, stalks, interior, members, graph, overlaps)
end

"""
    ghost_layer_cover(A, subdomains) -> OverlapCover

The cover by closed subdomains ``\\overline\\Omega_i = \\Omega_i \\cup \\Gamma_i``,
where the ghost layer ``\\Gamma_i`` is the set of dofs outside ``\\Omega_i``
that `A` couples to ``\\Omega_i``. See [`OverlapCover`](@ref).
"""
function ghost_layer_cover(A::AbstractMatrix, subdomains::AbstractVector{<:AbstractVector{<:Integer}})
    S = _structure(dropzeros(sparse(A)))
    doms = [sort!(unique(Vector{Int}(d))) for d in subdomains]
    return OverlapCover(doms, size(S, 1); ghosts=[_boundary(S, d) for d in doms])
end

# Position of dof k in the interior of subdomain i, or nothing.
function _interior_position(cover::OverlapCover, i::Int, k::Int)
    for (j, ℓ) in cover.members[k]
        j == i && cover.interior[i][ℓ] && return ℓ
    end
    return nothing
end

_interior_members(cover::OverlapCover, k::Int) =
    ((i, ℓ) for (i, ℓ) in cover.members[k] if cover.interior[i][ℓ])

_contains(cover::OverlapCover, i::Int, k::Int) = _interior_position(cover, i, k) !== nothing

"""
    overlap_sheaf(cover::OverlapCover, T=Float64) -> EuclideanSheaf{T,SelectionRestriction{T}}
    overlap_sheaf(subdomains, T=Float64) -> EuclideanSheaf{T,SelectionRestriction{T}}
    overlap_sheaf(dd::SchwarzDecomposition) -> EuclideanSheaf

The cellular sheaf of a cover of a set of degrees of freedom (dofs).

- **Vertices** are subdomains, with stalk ``F(i) = \\mathbb{R}^{\\overline\\Omega_i}``
  (one coordinate per dof of the closed subdomain, in increasing global order).
- **Edges** are the nonempty pairwise overlaps
  ``\\overline\\Omega_i \\cap \\overline\\Omega_j``, with stalk
  ``F(ij) = \\mathbb{R}^{\\overline\\Omega_i \\cap \\overline\\Omega_j}``.
- **Restriction maps** ``F(i) \\to F(ij)`` are coordinate selections
  ([`SelectionRestriction`](@ref)): they read off the values a subdomain holds
  on the overlap. They are stored as index lists, never as matrices.

For a plain list of subdomains there are no ghost layers and
``\\overline\\Omega_i = \\Omega_i``. The sheaf of a [`SchwarzDecomposition`](@ref)
uses its ghost-layer cover, so its edge stalks also carry the boundary data
the local solves exchange.

The underlying graph is the 1-skeleton of the nerve of the cover. A 0-cochain
``x = (x_i)`` assigns each subdomain its own local function; the coboundary
``(\\delta x)_{ij} = x_i|_{\\overline\\Omega_i\\cap\\overline\\Omega_j} - x_j|_{\\overline\\Omega_i\\cap\\overline\\Omega_j}``
measures how much neighbouring subdomains disagree. Because every nonempty
overlap is an edge, the global sections ``H^0`` are exactly the cochains coming
from a single function on ``\\bigcup_i \\overline\\Omega_i``.

This is the sheaf on which the Schwarz methods of this module run: iterates are
0-cochains, and convergence means they become a global section.
"""
function overlap_sheaf(cover::OverlapCover, ::Type{T}=Float64) where {T}
    stalks = cover.stalks
    s = EuclideanSheaf{T,SelectionRestriction{T}}(length.(stalks))
    for i in eachindex(stalks), j in neighbors(cover.graph, i)
        i < j || continue
        add_sheaf_edge!(s, i, j,
            SelectionRestriction{T}(cover.overlaps[i => j], length(stalks[i])),
            SelectionRestriction{T}(cover.overlaps[j => i], length(stalks[j])))
    end
    return s
end

overlap_sheaf(subdomains::AbstractVector{<:AbstractVector{<:Integer}}, ::Type{T}=Float64) where {T} =
    overlap_sheaf(OverlapCover(subdomains), T)

"""
    overlapping_subdomains(A, parts; overlap=1) -> Vector{Vector{Int}}

Grow a non-overlapping partition of the dofs into an overlapping cover.

`parts[k] ∈ 1:N` labels the part containing dof `k`. Each part is enlarged by
`overlap` layers of neighbours in the adjacency graph of the sparse matrix `A`
(dofs `k`, `m` are adjacent when `A[k, m] ≠ 0`), the usual algebraic way to
build the subdomains ``\\Omega_i`` of an overlapping Schwarz method
(Smith–Bjørstad–Gropp 1996, §1.3). `overlap = 0` returns the partition itself.

Pass `owner = parts` to [`SchwarzDecomposition`](@ref) to make the original
partition decide which subdomain supplies each dof's value.
"""
function overlapping_subdomains(A::AbstractMatrix, parts::AbstractVector{<:Integer}; overlap::Integer=1)
    n = size(A, 1)
    @argcheck size(A, 2) == n "A must be square"
    @argcheck length(parts) == n "need one part label per dof"
    @argcheck overlap >= 0
    @argcheck all(>=(1), parts) "part labels must be positive"
    S = _structure(dropzeros(sparse(A)))
    subdomains = Vector{Vector{Int}}(undef, maximum(parts))
    for i in eachindex(subdomains)
        mark = parts .== i
        @argcheck any(mark) "part $i is empty"
        for _ in 1:overlap
            grown = copy(mark)
            for c in findall(mark), r in _adjacent(S, c)
                grown[r] = true
            end
            mark = grown
        end
        subdomains[i] = findall(mark)
    end
    return subdomains
end

# ===== Ownership =====

"""
    Ownership(cover, owner)
    Ownership(cover, S::SparseMatrixCSC)

Which subdomain's copy supplies each dof: `owner[k]` is a subdomain whose interior contains
dof `k`, and `local_index[k]` is the position of `k` in its stalk. The owner's copy is
read whenever another subdomain needs `k` as boundary data, and when a cochain
is glued into a global vector ([`glue`](@ref)).

Given the matrix instead of an owner vector, each dof is owned by a subdomain
that contains it together with all of its matrix neighbours, when one exists.
In a ghost-layer cover a boundary dof of one subdomain lies in the ghost layer
of that subdomain and in the interior of its owner, so it is always in a shared
edge stalk.
"""
struct Ownership
    owner::Vector{Int}
    local_index::Vector{Int}
end

function Ownership(cover::OverlapCover, owner::AbstractVector{<:Integer})
    n = length(cover.members)
    @argcheck length(owner) == n "need one owner per dof"
    local_index = Vector{Int}(undef, n)
    for k in 1:n
        ℓ = _interior_position(cover, Int(owner[k]), k)
        @argcheck ℓ !== nothing "owner[$k] = $(owner[k]), but dof $k is not in the interior of subdomain $(owner[k])"
        local_index[k] = ℓ
    end
    return Ownership(Vector{Int}(owner), local_index)
end

function Ownership(cover::OverlapCover, S::SparseMatrixCSC)
    owner = map(eachindex(cover.members)) do k
        candidates = [i for (i, _) in _interior_members(cover, k)]
        @argcheck !isempty(candidates) "dof $k lies in no subdomain interior"
        deep = findfirst(i -> all(r -> _contains(cover, i, r), _adjacent(S, k)), candidates)
        candidates[something(deep, 1)]
    end
    return Ownership(cover, owner)
end

# ===== Transmission conditions =====

"""
    TransmissionCondition

How neighbouring subdomains pass data to each other in a local solve:
[`DirichletTransmission`](@ref) (classical Schwarz) or
[`RobinTransmission`](@ref) (optimized Schwarz).
"""
abstract type TransmissionCondition end

"""
    DirichletTransmission()

Each subdomain takes the values of its neighbours on its boundary
``\\Gamma_i`` as Dirichlet data. The local matrix is
``\\tilde A_i = A_{\\Omega_i\\Omega_i}``. This is classical Schwarz.
"""
struct DirichletTransmission <: TransmissionCondition end

"""
    RobinTransmission(p)

Subdomains exchange Robin data ``(\\partial_n + p) u`` instead of values, which
gives optimized Schwarz methods (Gander 2006). `p` is a positive number, or a
function `(i, j) -> p_ij` giving one parameter per edge of the overlap graph.

The algebraic form follows St-Cyr–Gander–Thomas (2007). The *interface* dofs of
``\\Omega_i`` are those coupled to ``\\Gamma_i``, and the local matrix is

```math
\\tilde A_i = A_{\\Omega_i\\Omega_i} - N_i + P_i .
```

The modification is made face by face. A *face* is a coupling ``A_{mk} \\neq 0``
between an interface dof ``m \\in \\Omega_i`` and a boundary dof
``k \\in \\Gamma_i`` owned by the neighbour ``j`` (a [`RobinFace`](@ref)).
``N_i = \\mathrm{diag}\\big(\\sum_{k} |A_{mk}|\\big)`` removes the Dirichlet
coupling of every face and leaves an algebraic Neumann condition. ``P_i`` adds
``p_{ij}`` once per face. The face's Robin data, the discrete
``(\\partial_n + p_{ij})\\, u_j``, is read entirely from neighbour ``j``. The
modification only touches interface rows, so the exact solution is still the
fixed point.

# Cross points

At a subdomain corner an interface dof has several faces, often towards
different neighbours. Discrete optimized Schwarz methods are known to be
fragile there and can diverge (Gander–Kwok 2012/2013). Three choices make
the corners consistent, and together they removed the divergence on box
decompositions in our tests:

1. a Robin term ``p_{ij}`` for *each* face whose Dirichlet coupling is removed,
   so that a corner dof is not under-penalized;
2. each face's Robin data taken from the neighbour across that face, never a
   mix of copies from different subdomains;
3. no pushes in the alternating sweeps: each subdomain keeps its own copy on
   the overlaps, as optimized Schwarz requires (see [`MultiplicativeSweep`](@ref)).

A large `p` approaches Dirichlet transmission. `p` is in the units of `A`: for a
stencil scaled by ``h^{-2}``, a continuous Robin parameter ``p`` corresponds to
``p/h`` (see [`optimized_robin_parameter`](@ref)). Every ``\\tilde A_i`` must be
positive definite, which holds for diagonally dominant `A` (e.g. M-matrices)
and `p > 0`. [`SchwarzDecomposition`](@ref) checks this.
"""
struct RobinTransmission{P} <: TransmissionCondition
    parameter::P

    function RobinTransmission(p::P) where {P<:Union{Real,Function}}
        p isa Real && @argcheck p > 0 "the Robin parameter must be positive"
        return new{P}(p)
    end
end

_robin_parameter(t::RobinTransmission{<:Real}, i, j) = t.parameter
_robin_parameter(t::RobinTransmission, i, j) = t.parameter(i, j)

"""
    optimized_robin_parameter(L; kmin=π, η=0) -> Real

The optimized Robin parameter of zeroth order for overlapping Schwarz on
``(\\eta - \\Delta) u = f``, to leading order in the overlap width ``L``:

```math
p^* = \\frac{(k_{\\min}^2 + \\eta)^{1/3}}{(2L)^{1/3}}
```

(Gander 2006, overlapping OO0). ``k_{\\min}`` is the lowest frequency along the
interface, ``\\pi / \\ell`` for an interface of length ``\\ell`` with Dirichlet
ends. With Robin transmission the contraction factor is ``1 - O(L^{1/3})``
instead of ``1 - O(L)`` for Dirichlet transmission.

This is a parameter of the continuous problem. For a finite-difference matrix
scaled by ``h^{-2}``, use `RobinTransmission(p / h)`. The formula is asymptotic
and for model problems, so treat it as a starting point. Where subdomain
boundaries meet at cross points, larger values are safer (Gander–Kwok 2013).
"""
function optimized_robin_parameter(L::Real; kmin::Real=π, η::Real=0)
    @argcheck L > 0 "the overlap width must be positive"
    @argcheck kmin > 0 && η >= 0
    return cbrt(kmin^2 + η) / cbrt(2L)
end

# ===== Local problems =====

"""
    RobinFace

One face of a subdomain's interface under [`RobinTransmission`](@ref): the
coupling between an interface dof ``m \\in \\Omega_i`` and a boundary dof
``k \\in \\Gamma_i`` across it.

- `dof`: the local index of ``m`` in ``\\Omega_i``;
- `weight`: ``p_{ij} - |A_{mk}|``, the Robin term added for this face minus
  the Dirichlet coupling it replaces;
- `source`: the neighbour ``j`` across the face (the owner of ``k``);
- `source_dof`: the position of ``m`` in the stalk of ``j``, where ``m`` is an
  interior dof of ``j``.

The face's Robin data is read entirely from the neighbour ``j``: its values
at ``k`` (through the ghost layer) and at ``m`` (through `weight`). Together
they form the discrete ``(\\partial_n + p_{ij})\\, u_j`` on that face. Both
``m`` and ``k`` lie in ``\\overline\\Omega_i \\cap \\overline\\Omega_j``, so the
data is in the edge stalk between ``i`` and ``j``. An interface dof at a
subdomain corner has one face per outside neighbour, each with its own
neighbour, parameter and data.
"""
struct RobinFace{T}
    dof::Int
    weight::T
    source::Int
    source_dof::Int
end

"""
    LocalProblem

Everything one subdomain needs for its local solve:

- `dofs`: the global dofs of ``\\Omega_i``;
- `boundary`: its ghost layer, the discrete boundary
  ``\\Gamma_i = \\{k \\notin \\Omega_i : A_{km} \\neq 0 \\text{ for some } m \\in \\Omega_i\\}``;
- `interior`, `ghosts`: the positions of ``\\Omega_i`` and ``\\Gamma_i`` in the
  vertex stalk ``\\mathbb R^{\\overline\\Omega_i}``;
- `coupling`: the block ``A_{\\Omega_i \\Gamma_i}``;
- `faces`: the [`RobinFace`](@ref)s of the interface (empty for
  [`DirichletTransmission`](@ref));
- `factor`: a sparse `ChordalLDLt` factorization of the local matrix
  ``\\tilde A_i``, computed once and reused by every local solve.
"""
struct LocalProblem{T,F}
    dofs::Vector{Int}
    boundary::Vector{Int}
    interior::Vector{Int}
    ghosts::Vector{Int}
    coupling::SparseMatrixCSC{T,Int}
    faces::Vector{RobinFace{T}}
    factor::F
end

# The global data every local problem is assembled from.
struct _Assembly{T,C<:TransmissionCondition}
    A::SparseMatrixCSC{T,Int}
    cover::OverlapCover
    ownership::Ownership
    transmission::C
end

function LocalProblem(asm::_Assembly, i::Int)
    S = asm.A
    dofs = asm.cover.subdomains[i]
    interior = findall(asm.cover.interior[i])
    ghosts = findall(!, asm.cover.interior[i])
    boundary = asm.cover.stalks[i][ghosts]
    faces = _robin_faces(asm, i, boundary)
    Ai = S[dofs, dofs]
    if !isempty(faces)
        local_dofs = [face.dof for face in faces]
        Ai = Ai + sparse(local_dofs, local_dofs, [face.weight for face in faces], length(dofs), length(dofs))
    end
    factor = _local_factor(Ai, issymmetric(S), i)
    return LocalProblem(dofs, boundary, interior, ghosts, S[dofs, boundary], faces, factor)
end

_robin_faces(::_Assembly{T,DirichletTransmission}, i, boundary) where {T} = RobinFace{T}[]

function _robin_faces(asm::_Assembly{T,<:RobinTransmission}, i, boundary) where {T}
    S = asm.A
    on_boundary = falses(size(S, 1))
    on_boundary[boundary] .= true
    faces = RobinFace{T}[]
    for (ℓ, m) in enumerate(asm.cover.subdomains[i]), idx in nzrange(S, m)
        k = rowvals(S)[idx]
        (on_boundary[k] && !iszero(nonzeros(S)[idx])) || continue
        j = asm.ownership.owner[k]
        p = T(_robin_parameter(asm.transmission, i, j))
        @argcheck p > 0 "Robin parameters must be positive (got $p on edge $(i)–$(j))"
        position = _interior_position(asm.cover, j, m)
        src, src_dof = position === nothing ?
            (asm.ownership.owner[m], asm.ownership.local_index[m]) : (j, position)
        @argcheck src != i "Robin transmission needs overlapping subdomains: the interface dof $m of subdomain $i is not in the interior of any neighbour (use overlap >= 1)"
        push!(faces, RobinFace(ℓ, p - abs(nonzeros(S)[idx]), src, src_dof))
    end
    return faces
end

# ===== Decomposition =====

"""
    SchwarzDecomposition(A, subdomains; owner=nothing, transmission=DirichletTransmission())

A domain decomposition of the sparse system ``A u = f``, prepared for Schwarz
iteration. `A` is typically symmetric positive definite (a finite-difference or
finite-element discretization of an elliptic PDE). With
[`DirichletTransmission`](@ref), `A` may also be nonsymmetric, for example the
M-matrix of an upwind discretization of a transport or Hamilton–Jacobi–Bellman
equation: the local problems are then factored with a sparse LU, and the ghost
layers come from the symmetrized sparsity pattern of `A`. Multiplicative and
additive Schwarz converge for nonsingular M-matrices (Frommer and Szyld,
*Weighted max norm estimates for additive Schwarz methods*, 1999). Robin
transmission, the coarse spaces, [`SchwarzCG`](@ref) and [`SheafADMM`](@ref)
need a symmetric `A`. `subdomains[i]` lists the
global dofs of ``\\Omega_i``; together they must cover `1:size(A, 1)`. They may
overlap, or not: with `overlap = 0` in [`overlapping_subdomains`](@ref) the
sweeps are block Gauss–Seidel and block Jacobi. The decomposition bundles

- `A`, as a sparse matrix;
- `cover`: the ghost-layer [`OverlapCover`](@ref)
  (see [`ghost_layer_cover`](@ref)) with its overlap graph and restriction
  maps (see [`overlap_sheaf`](@ref));
- `ownership`: the [`Ownership`](@ref) of each dof, from the vector `owner`
  if given and from the sparsity of `A` otherwise;
- `locals`: one [`LocalProblem`](@ref) per subdomain, built with the given
  [`TransmissionCondition`](@ref);
- `colors`: a coloring of the subdomains for [`MulticolorSweep`](@ref);
- `transmission`: the [`TransmissionCondition`](@ref) the local problems were built with.

Each vertex stalk is the closed subdomain ``\\overline\\Omega_i``. A local
solve receives its ghost values ``\\Gamma_i`` from their owners. Each of those
values lies in the edge stalk ``\\overline\\Omega_i \\cap \\overline\\Omega_j``
shared with its owner ``j``, so all communication is along edges of the overlap
graph, through the restriction maps.

Two subdomains get different colors when they *conflict*: they share a dof that
is interior to at least one of them, i.e. their interiors overlap or one's
ghost layer meets the other's interior. Subdomains of one color neither read
nor write each other's data, so they can be solved concurrently. The coloring
is greedy, largest conflict degree first.
"""
struct SchwarzDecomposition{T,F}
    A::SparseMatrixCSC{T,Int}
    cover::OverlapCover
    ownership::Ownership
    locals::Vector{LocalProblem{T,F}}
    colors::Vector{Vector{Int}}
    transmission::TransmissionCondition
end

function SchwarzDecomposition(A::AbstractMatrix, subdomains::AbstractVector{<:AbstractVector{<:Integer}};
                              owner::Union{Nothing,AbstractVector{<:Integer}}=nothing,
                              transmission::TransmissionCondition=DirichletTransmission())
    n = size(A, 1)
    @argcheck size(A, 2) == n "A must be square"
    S = dropzeros(sparse(float.(A)))
    @argcheck issymmetric(S) || transmission isa DirichletTransmission "Robin transmission needs a symmetric A"
    @argcheck !isempty(subdomains) "need at least one subdomain"
    for (i, d) in enumerate(subdomains)
        @argcheck !isempty(d) "subdomain $i is empty"
        @argcheck all(k -> 1 <= k <= n, d) "subdomain $i has dof indices outside 1:$n"
    end

    cover = ghost_layer_cover(S, subdomains)
    uncovered = findfirst(k -> isempty(_interior_members(cover, k)), 1:n)
    @argcheck uncovered === nothing "dof $uncovered lies in no subdomain; the subdomains must cover 1:$n"

    ownership = owner === nothing ? Ownership(cover, _structure(S)) : Ownership(cover, owner)
    asm = _Assembly(S, cover, ownership, transmission)
    locals = [LocalProblem(asm, i) for i in eachindex(cover.subdomains)]
    colors = _greedy_coloring(_conflict_graph(cover))
    return SchwarzDecomposition(S, cover, ownership, locals, colors, transmission)
end

# Subdomains conflict when they share a dof that is interior to at least one of
# them: a concurrent update of one could change data the other reads or writes.
function _conflict_graph(cover::OverlapCover)
    conflicts = SimpleGraph(length(cover.stalks))
    for m in cover.members, a in eachindex(m), b in (a + 1):lastindex(m)
        (i, ℓi), (j, ℓj) = m[a], m[b]
        (cover.interior[i][ℓi] || cover.interior[j][ℓj]) && add_edge!(conflicts, i, j)
    end
    return conflicts
end

# Greedy coloring, largest degree first (ties by vertex id). Returns the color
# classes, each in increasing vertex order.
function _greedy_coloring(g::SimpleGraph{Int})
    color = zeros(Int, nv(g))
    for v in sort(collect(vertices(g)); by=v -> (-degree(g, v), v))
        used = Set(color[u] for u in neighbors(g, v))
        c = 1
        while c in used
            c += 1
        end
        color[v] = c
    end
    return [findall(==(c), color) for c in 1:maximum(color; init=0)]
end

overlap_sheaf(dd::SchwarzDecomposition{T}) where {T} = overlap_sheaf(dd.cover, T)

function Base.show(io::IO, dd::SchwarzDecomposition)
    print(io, "SchwarzDecomposition(", length(dd.locals), " subdomains, ",
        size(dd.A, 1), " dofs, ", ne(dd.cover.graph), " overlaps)")
end

# ===== Cochains =====

function _cochain_blocks(dd::SchwarzDecomposition, x::BlockVector)
    xs = blocks(x)
    @argcheck length(xs) == length(dd.locals) && all(length.(xs) .== length.(dd.cover.stalks)) "cochain blocks must match the stalk sizes"
    return xs
end

function _cochain_blocks(dd::SchwarzDecomposition, x::AbstractVector)
    sizes = length.(dd.cover.stalks)
    @argcheck length(x) == sum(sizes) "cochain length must equal the total size of the stalks"
    offsets = [0; cumsum(sizes)]
    return [x[offsets[i]+1:offsets[i+1]] for i in eachindex(sizes)]
end

_owned_value(dd::SchwarzDecomposition, xs, k::Int) =
    xs[dd.ownership.owner[k]][dd.ownership.local_index[k]]

"""
    localize(dd::SchwarzDecomposition, u) -> BlockVector

Restrict a global vector `u` to every closed subdomain, giving the 0-cochain
``(u|_{\\overline\\Omega_1}, \\dots, u|_{\\overline\\Omega_N})`` of the overlap sheaf (interior and ghost
values). This is the
isomorphism from functions on the domain onto the global sections ``H^0``.
"""
function localize(dd::SchwarzDecomposition, u::AbstractVector)
    @argcheck length(u) == size(dd.A, 1)
    return mortar([u[d] for d in dd.cover.stalks])
end

"""
    glue(dd::SchwarzDecomposition, x) -> Vector

Glue a 0-cochain of the overlap sheaf into a single global vector, taking the
value of each dof from the copy held by its owner subdomain. On a global
section all copies agree, so this inverts [`localize`](@ref); on a general
cochain it is the restricted (owner-weighted) partition of unity used by
restricted additive Schwarz (Cai–Sarkis 1999).
"""
function glue(dd::SchwarzDecomposition, x::AbstractVector)
    xs = _cochain_blocks(dd, x)
    return [_owned_value(dd, xs, k) for k in eachindex(dd.ownership.owner)]
end

"""
    overlap_disagreement(dd::SchwarzDecomposition, x) -> Real

The norm ``\\|\\delta x\\|_2`` of the coboundary of the 0-cochain `x` in the
overlap sheaf: the root-sum-square of the differences between the values
neighbouring subdomains hold on their overlaps. It vanishes exactly on global
sections, so it measures how far a Schwarz iterate is from gluing.
Equal to `norm(coboundary_map(overlap_sheaf(dd)) * x)`, without assembling the
sheaf.
"""
function overlap_disagreement(dd::SchwarzDecomposition, x::AbstractVector)
    xs = _cochain_blocks(dd, x)
    acc = zero(real(eltype(eltype(xs))))
    for ((i, j), idx) in dd.cover.overlaps
        i < j || continue
        acc += sum(abs2, view(xs[i], idx) .- view(xs[j], dd.cover.overlaps[j => i]))
    end
    return sqrt(acc)
end

# ===== Problems =====

"""
    SchwarzProblem(dd::SchwarzDecomposition, f; u0=zeros(n))

The linear system ``A u = f`` together with the decomposition `dd` of its
dofs and an initial guess `u0`. Solve it with
[`solve`](@ref)`(problem, algorithm)`, where the algorithm is a
[`SchwarzIteration`](@ref) or a [`SchwarzCG`](@ref).
"""
struct SchwarzProblem{T,D<:SchwarzDecomposition{T}}
    decomposition::D
    rhs::Vector{T}
    u0::Vector{T}
end

function SchwarzProblem(dd::SchwarzDecomposition{T}, f::AbstractVector;
                        u0::AbstractVector=zeros(T, size(dd.A, 1))) where {T}
    @argcheck length(f) == size(dd.A, 1) "the right-hand side must have one entry per dof"
    @argcheck length(u0) == size(dd.A, 1) "the initial guess must have one entry per dof"
    return SchwarzProblem(dd, Vector{T}(f), Vector{T}(u0))
end

_relative_residual(prob::SchwarzProblem, u) =
    norm(prob.rhs - prob.decomposition.A * u) / _rhs_scale(prob)

function _rhs_scale(prob::SchwarzProblem{T}) where {T}
    scale = norm(prob.rhs)
    return iszero(scale) ? one(T) : scale
end

# ===== Sweeps =====

"""
    SchwarzSweep

One sweep of local solves over all subdomains, applied by
[`schwarz_step!`](@ref): [`MultiplicativeSweep`](@ref),
[`MulticolorSweep`](@ref) or [`ParallelSweep`](@ref).
"""
abstract type SchwarzSweep end

"""
    MultiplicativeSweep(order=nothing)

Schwarz's *alternating* method (Schwarz 1870). It visits the subdomains in
`order` (default `1:N`), each solve reading the newest data of its neighbours.

- With [`DirichletTransmission`](@ref), each local solve is followed by a push
  of the new values through the restriction maps onto every overlapping
  neighbour, ``x_j|_{\\Omega_i\\cap\\Omega_j} \\leftarrow x_i|_{\\Omega_i\\cap\\Omega_j}``.
  Starting from a global section the iterate stays a section, and its glued
  vector is exactly the classical multiplicative Schwarz iterate
  ``u \\leftarrow u + R_i^\\mathsf{T} A_i^{-1} R_i (f - A u)``. It converges for
  every SPD `A`, since it is block Gauss–Seidel over overlapping blocks
  (Toselli–Widlund 2005, ch. 2).
- With [`RobinTransmission`](@ref) nothing is pushed: each subdomain keeps its
  own copy on the overlaps, which is the alternating optimized Schwarz method.
  Pushing would overwrite the neighbours' Robin solutions with values that do
  not satisfy their transmission conditions; with cross points that made the
  iteration diverge for small `p`.
"""
struct MultiplicativeSweep <: SchwarzSweep
    order::Union{Nothing,Vector{Int}}
end

MultiplicativeSweep() = MultiplicativeSweep(nothing)

"""
    MulticolorSweep()

The alternating method with the subdomains visited one color class of
`dd.colors` at a time. The subdomains in a class do not conflict, so they are
solved (and, for Dirichlet transmission, pushed) concurrently
(`Threads.@threads`; start Julia with several threads to benefit). The result equals `MultiplicativeSweep(reduce(vcat,
dd.colors))`. It keeps that method's convergence guarantee while the number of
sequential steps per sweep drops from the number of subdomains to the number
of colors (Smith–Bjørstad–Gropp 1996, §1.4).
"""
struct MulticolorSweep <: SchwarzSweep end

"""
    ParallelSweep()

Lions' *parallel* Schwarz method (Lions 1988). All subdomains solve
simultaneously from the previous cochain and no values are pushed. Copies on
overlaps disagree during the iteration (the cochain is not a section) and agree
in the limit. With Dirichlet transmission the glued iterate coincides with
restricted additive Schwarz (RAS) for the owner partition
(Efstathiou–Gander 2003), and it converges when `A` is an M-matrix, such as
standard discretizations of ``-\\Delta`` (Frommer–Szyld 2001), but not for
every SPD matrix. With [`RobinTransmission`](@ref) it is the discrete parallel
optimized Schwarz method. A step from a global section equals one step of
optimized RAS (St-Cyr–Gander–Thomas 2007). Later steps differ, because each
face reads the neighbour across it rather than the owner of each dof.
"""
struct ParallelSweep <: SchwarzSweep end

# Receive: fill the ghost layer of subdomain i from the owners of its ghost
# dofs. Each value lies in the edge stalk between i and its owner, so this is
# the restriction map of the owner followed by the adjoint restriction of i.
function _receive!(xs, dd::SchwarzDecomposition, i::Int)
    lp = dd.locals[i]
    for (g, k) in zip(lp.ghosts, lp.boundary)
        xs[i][g] = _owned_value(dd, xs, k)
    end
    return nothing
end

# Subproblem on Ω_i: Ã_i x_i = f|Ω_i − A_{Ω_i Γ_i} g_i + Σ_faces w x_j(m), with
# g_i the ghost values just received and each Robin face reading the interface
# value from the neighbour across it. The face sum vanishes for Dirichlet
# transmission. Returns the new interior values.
function _local_solve(prob::SchwarzProblem, xs, i::Int)
    lp = prob.decomposition.locals[i]
    b = prob.rhs[lp.dofs] - lp.coupling * xs[i][lp.ghosts]
    for face in lp.faces
        b[face.dof] += face.weight * xs[face.source][face.source_dof]
    end
    return _factor_solve(lp.factor, b)
end

# Receive, solve on Ω_i, then publish the result to the neighbours.
function _solve_and_publish!(xs, prob::SchwarzProblem, i::Int)
    dd = prob.decomposition
    _receive!(xs, dd, i)
    xs[i][dd.locals[i].interior] .= _local_solve(prob, xs, i)
    _publish!(xs, dd.transmission, dd.cover, i)
    return nothing
end

# Dirichlet transmission: push the new interior values along the restriction
# maps onto every neighbour's copy, so a section stays a section.
function _publish!(xs, ::DirichletTransmission, cover::OverlapCover, i::Int)
    for j in neighbors(cover.graph, i)
        mine, theirs = cover.overlaps[i => j], cover.overlaps[j => i]
        for (a, b) in zip(mine, theirs)
            cover.interior[i][a] && (xs[j][b] = xs[i][a])
        end
    end
    return nothing
end

# Robin transmission: each subdomain keeps its own copy on the overlaps, as
# optimized Schwarz requires; neighbours read the new values on their next solve.
_publish!(xs, ::RobinTransmission, cover::OverlapCover, i::Int) = nothing

"""
    schwarz_step!(x::BlockVector, prob::SchwarzProblem, sweep::SchwarzSweep) -> x

Apply one `sweep` to the 0-cochain `x` of the overlap sheaf, in place. On
subdomain ``i`` the local solve is

```math
\\tilde A_i\\, x_i = f|_{\\Omega_i} - A_{\\Omega_i \\Gamma_i}\\, g_i
    + \\sum_{\\text{faces } (m, k)} (p_{ij} - |A_{mk}|)\\, x_j(m),
```

where ``\\tilde A_i`` is the local matrix of the
[`TransmissionCondition`](@ref). The boundary data ``g_i`` on ``\\Gamma_i`` is
read from the copies of the owners, which are neighbours in the overlap graph.
The face sum is the Robin data (see [`RobinFace`](@ref)) and is empty for
Dirichlet transmission.

With Robin transmission ``\\tilde A_i`` no longer dominates
``A_{\\Omega_i\\Omega_i}``, so convergence of the stationary iteration is not
proven for every parameter, although with the corner treatment of
[`RobinTransmission`](@ref) it converged for every `p` we tried. Use
[`SchwarzCG`](@ref) when a guarantee is needed.
"""
function schwarz_step!(x::BlockVector, prob::SchwarzProblem, sweep::MultiplicativeSweep)
    xs = _cochain_blocks(prob.decomposition, x)
    order = something(sweep.order, eachindex(xs))
    for i in order
        _solve_and_publish!(xs, prob, i)
    end
    return x
end

function schwarz_step!(x::BlockVector, prob::SchwarzProblem, ::MulticolorSweep)
    xs = _cochain_blocks(prob.decomposition, x)
    for class in prob.decomposition.colors
        Threads.@threads for i in class
            _solve_and_publish!(xs, prob, i)
        end
    end
    return x
end

function schwarz_step!(x::BlockVector, prob::SchwarzProblem, ::ParallelSweep)
    xs = _cochain_blocks(prob.decomposition, x)
    dd = prob.decomposition
    foreach(i -> _receive!(xs, dd, i), eachindex(xs))
    updates = [_local_solve(prob, xs, i) for i in eachindex(xs)]
    for (xi, lp, update) in zip(xs, dd.locals, updates)
        xi[lp.interior] .= update
    end
    return x
end

# ===== Coarse spaces =====

"""
    AbstractCoarseSpace

A second level for Schwarz methods. In a [`SchwarzIteration`](@ref) every fine
sweep is followed by a coarse correction ([`coarse_correct!`](@ref)); in
[`SchwarzCG`](@ref) the coarse term is added to the preconditioner.

Both implementations start from a graph homomorphism ``\\varphi : G \\to H``
that groups the subdomains (vertices of the overlap graph ``G``) into
aggregates, and from the pushforward ``\\varphi_* F`` of the overlap sheaf
``F``. The stalk ``(\\varphi_* F)(h)`` is the space of global sections of ``F``
over the fiber ``\\varphi^{-1}(h)``, i.e. functions on the aggregate
``\\widehat\\Omega_h = \\bigcup_{i \\in \\varphi^{-1}(h)} \\Omega_i``.

- [`TruncatedPushforwardCoarseSpace`](@ref) keeps a few modes of each stalk
  and solves a small Galerkin problem. It is cheap and scalable but inexact.
- [`ExactPushforwardCoarseSpace`](@ref) keeps the whole stalk, i.e. it solves
  on the aggregates themselves. It is more accurate per sweep while there are
  few aggregates, but its cost grows with the aggregate size.
"""
abstract type AbstractCoarseSpace end

function _check_hom(dd::SchwarzDecomposition, hom::GraphHomomorphism)
    @argcheck length(hom.vertex_map) == length(dd.locals) "the graph homomorphism must have one source vertex per subdomain"
    for h in 1:hom.n_target
        @argcheck !isempty(fiber_vertices(hom, h)) "aggregate $h has an empty fiber; every target vertex must receive a subdomain"
    end
end

"""
    TruncatedPushforwardCoarseSpace(dd, hom=identity; modes=ones(n, 1))

A low-dimensional coarse space built from the pushforward of the overlap sheaf
along ``\\varphi`` = `hom` by keeping only a few modes of every stalk.

Let ``\\mu_k`` be the number of subdomains containing dof ``k`` and
``D_h = \\mathrm{diag}\\big(\\#\\{i \\in \\varphi^{-1}(h) : k \\in \\Omega_i\\} / \\mu_k\\big)``,
so ``\\sum_h D_h = I`` is a partition of unity subordinate to the aggregates.
For every aggregate ``h`` and every column ``z`` of `modes` (a near-nullspace of
``A``, e.g. constants for the Laplacian, rigid body modes for elasticity), the
coarse basis contains ``D_h z``. This vector lives in the stalk
``(\\varphi_* F)(h)``, so the coarse space is a sub-cochain space
``V_0 \\subset C^0(\\varphi_* F)`` with `size(modes, 2)` dimensions per vertex of
``H``. Basis vectors that vanish are dropped.

With ``\\Phi`` the matrix of basis vectors, the restriction is
``R_0 = \\Phi^\\mathsf{T}``, the prolongation is ``\\Phi``, and the coarse matrix
``A_0 = \\Phi^\\mathsf{T} A \\Phi`` is sparse on the graph ``H`` and factored
once with `ChordalLDLt`. The correction is the ``A``-orthogonal projection
``e = \\Phi A_0^{-1} \\Phi^\\mathsf{T}(f - A u)``, so adding it after
multiplicative sweeps still converges for every SPD `A`.

With the identity homomorphism and constant modes this is the Nicolaides coarse
space (Nicolaides 1987), which makes the iteration count roughly independent of
the number of subdomains (Toselli–Widlund 2005, §3.10).
"""
struct TruncatedPushforwardCoarseSpace{T,F} <: AbstractCoarseSpace
    hom::GraphHomomorphism
    basis::SparseMatrixCSC{T,Int}
    stalks::Vector{Int}
    matrix::SparseMatrixCSC{T,Int}
    factor::F
end

function TruncatedPushforwardCoarseSpace(dd::SchwarzDecomposition{T},
                                         hom::GraphHomomorphism=GraphHomomorphism(collect(eachindex(dd.locals)));
                                         modes::AbstractVecOrMat=ones(T, size(dd.A, 1), 1)) where {T}
    _check_hom(dd, hom)
    @argcheck issymmetric(dd.A) "coarse spaces need a symmetric A"
    n = size(dd.A, 1)
    Z = reshape(modes, size(modes, 1), :)
    @argcheck size(Z, 1) == n "modes must have one row per dof"
    multiplicity = [count(Returns(true), _interior_members(dd.cover, k)) for k in 1:n]

    I, J, V = Int[], Int[], T[]
    stalks = zeros(Int, hom.n_target)
    for h in 1:hom.n_target
        weight = zeros(T, n)
        for i in fiber_vertices(hom, h)
            weight[dd.cover.subdomains[i]] .+= 1
        end
        support = findall(!iszero, weight)
        weight[support] ./= multiplicity[support]
        for c in axes(Z, 2)
            vals = weight[support] .* Z[support, c]
            nz = findall(!iszero, vals)
            isempty(nz) && continue
            stalks[h] += 1
            append!(I, support[nz])
            append!(J, fill(sum(stalks), length(nz)))
            append!(V, vals[nz])
        end
    end
    Φ = sparse(I, J, V, n, sum(stalks))
    A0 = Φ' * dd.A * Φ
    A0 = (A0 + A0') / 2
    factor = ldlt!(ChordalLDLt(A0), RowMaximum())
    return TruncatedPushforwardCoarseSpace(hom, Φ, stalks, A0, factor)
end

"""
    ExactPushforwardCoarseSpace(dd, hom; sweep=MulticolorSweep(), transmission=DirichletTransmission())

The coarse level given by the full pushforward ``\\varphi_* F`` of the overlap
sheaf along ``\\varphi`` = `hom`. Its stalk at an aggregate ``h`` is all
functions on ``\\widehat\\Omega_h = \\bigcup_{i \\in \\varphi^{-1}(h)} \\Omega_i``,
so its vertex stalks have the same dimensions as those of
`pushforward_sheaf(hom, overlap_sheaf(dd))`. Since no modes are discarded,
``H^0(\\varphi_* F) \\cong H^0(F)``.

Concretely this is a second [`SchwarzDecomposition`](@ref) whose subdomains
are the aggregates ``\\widehat\\Omega_h``, with the given `transmission`
between aggregates. A dof is owned by the aggregate containing its fine owner.
The correction is one `sweep` on that decomposition, started from the current
glued iterate. Each local solve is exact on a larger region, so while there are
few aggregates a sweep reduces the error more than a truncated coarse solve.
With ``\\varphi`` to a single vertex the coarse level is a direct solve of the
whole problem (one iteration, no scalability). The level carries no global
information beyond its aggregates, so with a fixed aggregation ratio it is a
one-level method on larger subdomains: the iteration count grows with the
number of aggregates and eventually exceeds that of
[`TruncatedPushforwardCoarseSpace`](@ref), while its factorizations cost as
much as the whole problem plus overlaps.
"""
struct ExactPushforwardCoarseSpace{D<:SchwarzDecomposition,S<:SchwarzSweep} <: AbstractCoarseSpace
    hom::GraphHomomorphism
    decomposition::D
    sweep::S
end

function ExactPushforwardCoarseSpace(dd::SchwarzDecomposition, hom::GraphHomomorphism;
                                     sweep::SchwarzSweep=MulticolorSweep(),
                                     transmission::TransmissionCondition=DirichletTransmission())
    _check_hom(dd, hom)
    @argcheck issymmetric(dd.A) "coarse spaces need a symmetric A"
    aggregates =[reduce(vcat, (dd.cover.subdomains[i] for i in fiber_vertices(hom, h)))
                  for h in 1:hom.n_target]
    owner = hom.vertex_map[dd.ownership.owner]
    return ExactPushforwardCoarseSpace(hom, SchwarzDecomposition(dd.A, aggregates; owner, transmission), sweep)
end

"""
    coarse_dimension(coarse::AbstractCoarseSpace) -> Int

Number of unknowns the coarse level solves for in each correction: the
dimension of the truncated coarse space, or the total size of the aggregates
(the sum of the vertex stalk dimensions of the pushforward sheaf) for the
exact pushforward.
"""
coarse_dimension(c::TruncatedPushforwardCoarseSpace) = size(c.basis, 2)
coarse_dimension(c::ExactPushforwardCoarseSpace) = sum(length, c.decomposition.cover.subdomains)

function Base.show(io::IO, c::Union{TruncatedPushforwardCoarseSpace,ExactPushforwardCoarseSpace})
    print(io, nameof(typeof(c)), "(", c.hom.n_target, " aggregates, dimension ",
        coarse_dimension(c), ")")
end

_coarse_correction(c::TruncatedPushforwardCoarseSpace, prob::SchwarzProblem, u) =
    c.basis * _ldlt_solve(c.factor, c.basis' * (prob.rhs - prob.decomposition.A * u))

function _coarse_correction(c::ExactPushforwardCoarseSpace, prob::SchwarzProblem, u)
    coarse_prob = SchwarzProblem(c.decomposition, prob.rhs, u)
    x = localize(c.decomposition, u)
    schwarz_step!(x, coarse_prob, c.sweep)
    return glue(c.decomposition, x) - u
end

"""
    coarse_correct!(x::BlockVector, prob::SchwarzProblem, coarse::AbstractCoarseSpace) -> x

Apply one coarse-level correction to the 0-cochain `x`, in place. The glued
iterate ``u`` = `glue(dd, x)` gives the residual ``f - A u``. The coarse level
turns it into a global correction ``e``, which is added to every subdomain's
copy: ``x \\leftarrow x + \\mathrm{localize}(e)``. The glued iterate becomes
``u + e``, and the disagreement ``\\|\\delta x\\|`` is unchanged.
"""
function coarse_correct!(x::BlockVector, prob::SchwarzProblem, coarse::AbstractCoarseSpace)
    dd = prob.decomposition
    xs = _cochain_blocks(dd, x)
    e = _coarse_correction(coarse, prob, glue(dd, x))
    for (xi, d) in zip(xs, dd.cover.stalks)
        xi .+= view(e, d)
    end
    return x
end

# ===== Algorithms =====

"""
    SchwarzResult

Output of [`solve`](@ref) for a [`SchwarzProblem`](@ref).

- `u`: the glued global solution.
- `x`: the final 0-cochain of the overlap sheaf (one local solution per
  subdomain).
- `residuals`: relative residual ``\\|f - A u_k\\| / \\|f\\|`` of each iterate
  (entry 1 is the initial guess).
- `disagreements`: [`overlap_disagreement`](@ref) of each iterate. Krylov
  iterates are global vectors, i.e. sections, so for [`SchwarzCG`](@ref)
  these are all zero.
- `iterations`, `converged`.
"""
struct SchwarzResult{T}
    u::Vector{T}
    x::BlockVector{T}
    residuals::Vector{T}
    disagreements::Vector{T}
    iterations::Int
    converged::Bool
end

"""
    SchwarzIteration(; sweep=MultiplicativeSweep(), coarse=nothing, tol=1e-8, maxiter=1000)

The stationary Schwarz iteration on the overlap sheaf. The iterate is a
0-cochain ``x = (x_1, \\dots, x_N)``: each subdomain holds its own copy of the
solution on ``\\Omega_i``. It starts from the global section `localize(dd,
u0)`. Each iteration applies one `sweep` ([`schwarz_step!`](@ref)) followed,
when a `coarse` space is given, by [`coarse_correct!`](@ref), which gives a
two-level method. Iteration stops once the glued iterate satisfies
``\\|f - A u\\| \\le \\mathrm{tol}\\,\\|f\\|``, or after `maxiter` iterations. At
the fixed point every local solve is consistent with its neighbours, so ``x``
is a global section and its gluing solves ``A u = f``.
"""
Base.@kwdef struct SchwarzIteration{S<:SchwarzSweep,C<:Union{Nothing,AbstractCoarseSpace}}
    sweep::S = MultiplicativeSweep()
    coarse::C = nothing
    tol::Float64 = 1e-8
    maxiter::Int = 1000
end

"""
    SchwarzCG(; coarse=nothing, tol=1e-8, maxiter=1000)

Conjugate gradients (Krylov.jl) preconditioned with the additive Schwarz
operator [`SchwarzPreconditioner`](@ref), with an optional `coarse` space.
Unlike [`SchwarzIteration`](@ref), this converges for every SPD `A`, with every
coarse space and with Robin transmission. The two-level version needs a number
of iterations bounded independently of the number of subdomains
(Toselli–Widlund 2005, Thm. 3.13). It stops once the relative residual of the
iterate drops below `tol`, the same criterion as [`SchwarzIteration`](@ref).
"""
Base.@kwdef struct SchwarzCG{C<:Union{Nothing,AbstractCoarseSpace}}
    coarse::C = nothing
    tol::Float64 = 1e-8
    maxiter::Int = 1000
end

"""
    solve(prob::SchwarzProblem, alg::SchwarzIteration) -> SchwarzResult
    solve(prob::SchwarzProblem, alg::SchwarzCG) -> SchwarzResult

Solve the [`SchwarzProblem`](@ref) with a stationary [`SchwarzIteration`](@ref)
or with Schwarz-preconditioned conjugate gradients ([`SchwarzCG`](@ref)).
"""
function solve(prob::SchwarzProblem{T}, alg::SchwarzIteration) where {T}
    @argcheck alg.maxiter >= 0
    dd = prob.decomposition
    x = localize(dd, prob.u0)
    u = glue(dd, x)
    residuals = T[_relative_residual(prob, u)]
    disagreements = T[overlap_disagreement(dd, x)]
    iterations = 0
    while last(residuals) > alg.tol && iterations < alg.maxiter
        schwarz_step!(x, prob, alg.sweep)
        alg.coarse === nothing || coarse_correct!(x, prob, alg.coarse)
        iterations += 1
        u = glue(dd, x)
        push!(residuals, _relative_residual(prob, u))
        push!(disagreements, overlap_disagreement(dd, x))
    end
    return SchwarzResult(u, x, residuals, disagreements, iterations, last(residuals) <= alg.tol)
end

function solve(prob::SchwarzProblem{T}, alg::SchwarzCG) where {T}
    @argcheck alg.maxiter >= 1
    @argcheck issymmetric(prob.decomposition.A) "SchwarzCG needs a symmetric A; use SchwarzIteration or SchwarzGMRES"
    dd = prob.decomposition
    residuals = T[_relative_residual(prob, prob.u0)]
    # CG solves for the correction e in A e = f − A u0 from e = 0, so the
    # iterate is u0 + e.
    function monitor(workspace)
        push!(residuals, _relative_residual(prob, prob.u0 + workspace.x))
        return last(residuals) <= alg.tol
    end
    u, iterations = if last(residuals) <= alg.tol
        copy(prob.u0), 0
    else
        e, stats = cg(dd.A, prob.rhs - dd.A * prob.u0; M=SchwarzPreconditioner(dd, alg.coarse),
            atol=zero(T), rtol=zero(T), itmax=alg.maxiter, callback=monitor)
        prob.u0 + e, stats.niter
    end
    return SchwarzResult(u, localize(dd, u), residuals, zeros(T, length(residuals)),
        iterations, last(residuals) <= alg.tol)
end

# ===== Krylov preconditioner =====

"""
    SchwarzPreconditioner(dd::SchwarzDecomposition, coarse=nothing)

The additive Schwarz preconditioner on the overlap sheaf of `dd`,

```math
M^{-1} r = \\sum_i R_i^\\mathsf{T} \\tilde A_i^{-1} R_i\\, r + (\\text{coarse term}),
```

where ``R_i`` restricts to ``\\Omega_i`` and ``\\tilde A_i`` is the local matrix
of `dd` (from its [`TransmissionCondition`](@ref)). All local solves are
independent and run concurrently (`Threads.@threads`). The coarse term is

- ``\\Phi A_0^{-1} \\Phi^\\mathsf{T} r`` for a
  [`TruncatedPushforwardCoarseSpace`](@ref), the two-level additive Schwarz
  method (Toselli–Widlund 2005, ch. 3);
- ``\\sum_h \\widehat R_h^\\mathsf{T} \\widehat A_h^{-1} \\widehat R_h r`` over the
  aggregates of an [`ExactPushforwardCoarseSpace`](@ref).

The full restrictions ``R_i^\\mathsf{T}`` (not the owner-restricted ones of
RAS) keep ``M^{-1}`` symmetric. It is positive definite because each
``\\tilde A_i`` is, so it can precondition conjugate gradients
([`SchwarzCG`](@ref)). Apply it with `P * r` or `mul!(z, P, r)`.
"""
struct SchwarzPreconditioner{T,D<:SchwarzDecomposition{T},C<:Union{Nothing,AbstractCoarseSpace}}
    decomposition::D
    coarse::C
end

SchwarzPreconditioner(dd::SchwarzDecomposition) = SchwarzPreconditioner(dd, nothing)

Base.size(P::SchwarzPreconditioner) = size(P.decomposition.A)
Base.size(P::SchwarzPreconditioner, d::Integer) = size(P.decomposition.A, d)
Base.eltype(::SchwarzPreconditioner{T}) where {T} = T

function _additive!(y::AbstractVector, dd::SchwarzDecomposition{T}, r::AbstractVector) where {T}
    corrections = Vector{Vector{T}}(undef, length(dd.locals))
    Threads.@threads for i in eachindex(dd.locals)
        lp = dd.locals[i]
        corrections[i] = _factor_solve(lp.factor, r[lp.dofs])
    end
    fill!(y, zero(eltype(y)))
    for (lp, z) in zip(dd.locals, corrections)
        view(y, lp.dofs) .+= z
    end
    return y
end

_add_coarse!(y, ::Nothing, r) = y
_add_coarse!(y, c::TruncatedPushforwardCoarseSpace, r) =
    (y .+= c.basis * _ldlt_solve(c.factor, c.basis' * r); y)
_add_coarse!(y, c::ExactPushforwardCoarseSpace, r) =
    (y .+= _additive!(similar(y), c.decomposition, r); y)

function LinearAlgebra.mul!(y::AbstractVector, P::SchwarzPreconditioner, r::AbstractVector)
    @argcheck length(y) == length(r) == size(P, 1)
    _additive!(y, P.decomposition, r)
    return _add_coarse!(y, P.coarse, r)
end

Base.:*(P::SchwarzPreconditioner{T}, r::AbstractVector) where {T} =
    mul!(similar(r, promote_type(T, eltype(r))), P, r)

# ===== GMRES acceleration =====

"""
    SchwarzSweepPreconditioner(dd, sweep, coarse=nothing)

One step of a stationary Schwarz iteration as a linear map: `P * r` starts from
the zero cochain, applies one `sweep` of [`schwarz_step!`](@ref) to
``A e = r``, then the `coarse` correction ([`coarse_correct!`](@ref)) if any,
and returns the glued result ``e = M^{-1} r``. The map is linear in ``r`` for
every sweep, transmission condition and coarse space:

- with [`ParallelSweep`](@ref) and Dirichlet transmission ``M^{-1}`` is
  restricted additive Schwarz (RAS);
- with [`ParallelSweep`](@ref) and [`RobinTransmission`](@ref) it is optimized
  RAS (ORAS, St-Cyr–Gander–Thomas 2007);
- with a coarse space it is the corresponding hybrid two-level method,
  ``M^{-1} = M_1^{-1} + \\Phi A_0^{-1}\\Phi^\\mathsf{T}(I - A M_1^{-1})`` for a
  [`TruncatedPushforwardCoarseSpace`](@ref).

It is not symmetric in general, so it preconditions GMRES ([`SchwarzGMRES`](@ref)).
"""
struct SchwarzSweepPreconditioner{T,D<:SchwarzDecomposition{T},S<:SchwarzSweep,C<:Union{Nothing,AbstractCoarseSpace}}
    decomposition::D
    sweep::S
    coarse::C
end

SchwarzSweepPreconditioner(dd::SchwarzDecomposition{T}, sweep::SchwarzSweep, coarse=nothing) where {T} =
    SchwarzSweepPreconditioner{T,typeof(dd),typeof(sweep),typeof(coarse)}(dd, sweep, coarse)

Base.size(P::SchwarzSweepPreconditioner) = size(P.decomposition.A)
Base.size(P::SchwarzSweepPreconditioner, d::Integer) = size(P.decomposition.A, d)
Base.eltype(::SchwarzSweepPreconditioner{T}) where {T} = T

function LinearAlgebra.mul!(y::AbstractVector, P::SchwarzSweepPreconditioner{T}, r::AbstractVector) where {T}
    dd = P.decomposition
    @argcheck length(y) == length(r) == size(dd.A, 1)
    n = size(dd.A, 1)
    prob = SchwarzProblem(dd, Vector{T}(r), zeros(T, n))
    x = localize(dd, prob.u0)
    schwarz_step!(x, prob, P.sweep)
    P.coarse === nothing || coarse_correct!(x, prob, P.coarse)
    y .= glue(dd, x)
    return y
end

Base.:*(P::SchwarzSweepPreconditioner{T}, r::AbstractVector) where {T} =
    mul!(similar(r, promote_type(T, eltype(r))), P, r)

"""
    SchwarzGMRES(; sweep=ParallelSweep(), coarse=nothing, tol=1e-8, maxiter=1000, memory=100)

GMRES (Krylov.jl, restarted every `memory` iterations) right-preconditioned
with one step of the stationary Schwarz iteration
([`SchwarzSweepPreconditioner`](@ref)). It accelerates any
[`SchwarzIteration`](@ref) with the same sweep and coarse space, including
combinations whose stationary iteration diverges. The standard example is
optimized (Robin) RAS with a coarse level, the two-level ORAS method
(Dolean–Jolivet–Nataf 2015, ch. 5).

Right preconditioning makes GMRES minimize the true residual. Iteration stops
when ``\\|f - A u\\| \\le \\mathrm{tol}\\,\\|f\\|``, the same criterion as
[`SchwarzIteration`](@ref) and [`SchwarzCG`](@ref). The relative residual of
the final iterate is checked and recorded.
"""
Base.@kwdef struct SchwarzGMRES{S<:SchwarzSweep,C<:Union{Nothing,AbstractCoarseSpace}}
    sweep::S = ParallelSweep()
    coarse::C = nothing
    tol::Float64 = 1e-8
    maxiter::Int = 1000
    memory::Int = 100
end

function solve(prob::SchwarzProblem{T}, alg::SchwarzGMRES) where {T}
    @argcheck alg.maxiter >= 1 && alg.memory >= 1
    dd = prob.decomposition
    r0 = prob.rhs - dd.A * prob.u0
    initial = _relative_residual(prob, prob.u0)
    if initial <= alg.tol || iszero(norm(r0))
        return SchwarzResult(copy(prob.u0), localize(dd, prob.u0), T[initial], T[0], 0, true)
    end
    # GMRES solves A e = f − A u0 for the correction e from e = 0.
    P = SchwarzSweepPreconditioner(dd, alg.sweep, alg.coarse)
    e, stats = gmres(dd.A, r0; N=P, memory=alg.memory, restart=true, atol=zero(T),
        rtol=T(alg.tol) * _rhs_scale(prob) / norm(r0), itmax=alg.maxiter, history=true)
    u = prob.u0 + e
    residuals = T.(stats.residuals) ./ _rhs_scale(prob)
    residuals[end] = _relative_residual(prob, u)
    return SchwarzResult(u, localize(dd, u), residuals, zeros(T, length(residuals)),
        stats.niter, last(residuals) <= alg.tol)
end

# ===== Sheaf ADMM =====

"""
    LocalObjective

The quadratic local objective ``f_i(x_i) = \\tfrac12 x_i^\\mathsf{T} K_i x_i - b_i^\\mathsf{T} x_i``
of one subdomain, on its vertex stalk ``\\mathbb R^{\\overline\\Omega_i}``, with
`matrix` ``K_i`` and `rhs` ``b_i``. [`local_objectives`](@ref) builds them so
that on every global section ``x = \\mathrm{localize}(u)``

```math
\\sum_i f_i(x_i) = \\tfrac12\\, u^\\mathsf{T} A u - f^\\mathsf{T} u .
```
"""
struct LocalObjective{T}
    matrix::SparseMatrixCSC{T,Int}
    rhs::Vector{T}
end

"""
    local_objectives(prob::SchwarzProblem) -> Vector{LocalObjective}

Split the energy ``\\tfrac12 u^\\mathsf{T} A u - f^\\mathsf{T} u`` of `prob` into one
[`LocalObjective`](@ref) per vertex stalk of the ghost-layer cover. The energy is
a sum of edge and vertex terms,

```math
u^\\mathsf{T} A u = \\sum_{k<l} |a_{kl}|\\,(u_k + \\operatorname{sign}(a_{kl})\\, u_l)^2
    + \\sum_k s_k\\, u_k^2, \\qquad s_k = a_{kk} - \\sum_{l \\neq k} |a_{kl}| .
```

Each edge term is shared equally among the stalks containing both ``k`` and
``l``, and each vertex term and ``f_k`` among the stalks containing ``k``. Every
edge of `A` lies in some closed subdomain, so nothing is lost. When `A` is
diagonally dominant (``s_k \\ge 0``, e.g. an M-matrix) every term is a square,
so every ``K_i`` is positive semidefinite.
"""
function local_objectives(prob::SchwarzProblem{T}) where {T}
    dd = prob.decomposition
    S, cover = dd.A, dd.cover
    entries = [(Int[], Int[], T[]) for _ in cover.stalks]
    rhs = [zeros(T, length(s)) for s in cover.stalks]
    function add!(i, a, b, v)
        push!(entries[i][1], a); push!(entries[i][2], b); push!(entries[i][3], v)
    end
    for k in axes(S, 2)
        holders = cover.members[k]
        excess = S[k, k] - sum(abs(nonzeros(S)[p]) for p in nzrange(S, k) if rowvals(S)[p] != k; init=zero(T))
        for (i, ℓ) in holders
            add!(i, ℓ, ℓ, excess / length(holders))
            rhs[i][ℓ] += prob.rhs[k] / length(holders)
        end
        for p in nzrange(S, k)
            l, a = rowvals(S)[p], nonzeros(S)[p]
            (l > k && !iszero(a)) || continue
            common = [(i, ℓk, ℓl) for (i, ℓk) in holders for (j, ℓl) in cover.members[l] if i == j]
            w, σ = abs(a) / length(common), sign(a)
            for (i, ℓk, ℓl) in common
                add!(i, ℓk, ℓk, w); add!(i, ℓl, ℓl, w)
                add!(i, ℓk, ℓl, σ * w); add!(i, ℓl, ℓk, σ * w)
            end
        end
    end
    return [LocalObjective(sparse(I, J, V, length(s), length(s)), b)
            for ((I, J, V), b, s) in zip(entries, rhs, cover.stalks)]
end

"""
    SheafADMM(; rho, penalty=:stalk, projection_steps=nothing, tol=1e-8, maxiter=1000)

Sheaf ADMM (Hanks, Riess, Cohen, Gross, Hale, Fairbanks 2025, Algorithm 1) on
the ghost-layer overlap sheaf. It solves the homological program

```math
\\min_{x \\in C^0} \\sum_i f_i(x_i) \\quad \\text{subject to} \\quad x \\in H^0 ,
```

with the [`local_objectives`](@ref) of the PDE, so its solution glues to the
solution of ``A u = f``. With copies ``z`` and scaled multipliers ``y`` on the
vertex stalks, each iteration is

```math
\\begin{aligned}
x_i &\\leftarrow \\operatorname{argmin}_{x_i} f_i(x_i) + \\tfrac{\\rho}{2}\\lVert x_i - z_i + y_i\\rVert_{P_i}^2
     = (K_i + \\rho P_i)^{-1}\\big(b_i + \\rho P_i (z_i - y_i)\\big),\\\\
z &\\leftarrow \\Pi_{H^0}(x + y),\\\\
y_i &\\leftarrow y_i + x_i - z_i .
\\end{aligned}
```

- `penalty = :stalk` is the method of the paper: ``P_i = I``, the penalty acts on
  the whole stalk. `penalty = :shared` penalizes only the dofs held by more than
  one stalk (the overlaps and ghost layers), which is where the consensus
  constraint actually binds. It is closer to a Robin transmission condition,
  whose penalty acts only on the interface (see the comparison in the feature
  guide).
- ``\\Pi_{H^0}`` is the orthogonal projection onto global sections. For this
  sheaf it averages the copies of every dof, one exchange with the neighbours
  (`projection_steps = nothing`). With `projection_steps = k` it is replaced by
  ``k`` explicit sheaf-diffusion steps ``z \\leftarrow z - \\alpha L_F z`` with
  ``\\alpha = 1/c_{\\max}`` (``c_{\\max}`` the largest number of copies of a dof),
  as in Algorithm 1, at one neighbour exchange per step.

The local matrices ``K_i + \\rho P_i`` are factored once with `ChordalLDLt`.
Convergence (Boyd et al. 2011, §3.2) needs convex ``f_i`` and an exact
projection, so `A` should be diagonally dominant; the local matrices are
checked for positive definiteness.
"""
Base.@kwdef struct SheafADMM
    rho::Float64
    penalty::Symbol = :stalk
    projection_steps::Union{Nothing,Int} = nothing
    tol::Float64 = 1e-8
    maxiter::Int = 1000
end

# Π_{H^0}: average the copies of every dof.
function _project_sections!(zs, dd::SchwarzDecomposition{T}, vs) where {T}
    total = zeros(T, size(dd.A, 1))
    copies = zeros(Int, size(dd.A, 1))
    for (v, s) in zip(vs, dd.cover.stalks)
        total[s] .+= v
        copies[s] .+= 1
    end
    total ./= copies
    for (z, s) in zip(zs, dd.cover.stalks)
        z .= view(total, s)
    end
    return zs
end

# k explicit sheaf-diffusion steps z ← z − α L z on the overlap sheaf.
function _diffuse_sections!(zs, dd::SchwarzDecomposition, vs, steps::Int)
    cover = dd.cover
    α = 1 / maximum(length, cover.members)
    foreach(copyto!, zs, vs)
    for _ in 1:steps
        Lz = [zero(z) for z in zs]
        for ((i, j), mine) in cover.overlaps
            theirs = cover.overlaps[j => i]
            Lz[i][mine] .+= zs[i][mine] .- zs[j][theirs]
        end
        foreach((z, l) -> (z .-= α .* l), zs, Lz)
    end
    return zs
end

function solve(prob::SchwarzProblem{T}, alg::SheafADMM) where {T}
    @argcheck issymmetric(prob.decomposition.A) "sheaf ADMM needs a symmetric A"
    @argcheck alg.rho > 0 "the ADMM penalty must be positive"
    @argcheck alg.penalty in (:stalk, :shared) "penalty must be :stalk or :shared"
    @argcheck alg.projection_steps === nothing || alg.projection_steps >= 1
    dd = prob.decomposition
    cover = dd.cover
    objectives = local_objectives(prob)
    masks = [alg.penalty === :stalk ? trues(length(s)) : BitVector(length(cover.members[k]) > 1 for k in s)
             for s in cover.stalks]
    factors = map(enumerate(objectives)) do (i, obj)
        Ki = obj.matrix + alg.rho * spdiagm(T.(masks[i]))
        F = ldlt!(ChordalLDLt(Ki), RowMaximum(); check=false)
        @argcheck all(>(0), F.D.diag) "the ADMM local matrix of subdomain $i is not positive definite; increase rho or use penalty=:stalk"
        F
    end

    z = localize(dd, prob.u0)
    zs = blocks(z)
    xs = [copy(zi) for zi in zs]
    ys = [zero(zi) for zi in zs]
    u = glue(dd, z)
    residuals = T[_relative_residual(prob, u)]
    disagreements = T[0]
    iterations = 0
    while last(residuals) > alg.tol && iterations < alg.maxiter
        Threads.@threads for i in eachindex(xs)
            b = objectives[i].rhs .+ alg.rho .* masks[i] .* (zs[i] .- ys[i])
            xs[i] .= _ldlt_solve(factors[i], b)
        end
        vs = [x .+ y for (x, y) in zip(xs, ys)]
        alg.projection_steps === nothing ? _project_sections!(zs, dd, vs) :
            _diffuse_sections!(zs, dd, vs, alg.projection_steps)
        for (y, x, zi) in zip(ys, xs, zs)
            y .+= x .- zi
        end
        iterations += 1
        u = glue(dd, z)
        push!(residuals, _relative_residual(prob, u))
        push!(disagreements, overlap_disagreement(dd, mortar(xs)))
    end
    return SchwarzResult(u, mortar(xs), residuals, disagreements, iterations, last(residuals) <= alg.tol)
end

end
