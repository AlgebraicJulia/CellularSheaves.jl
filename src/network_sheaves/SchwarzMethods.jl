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

export OverlapCover, Ownership, LocalProblem, RobinFace, SchwarzDecomposition,
    overlap_sheaf, overlapping_subdomains, localize, glue, overlap_disagreement,
    TransmissionCondition, DirichletTransmission, RobinTransmission, optimized_robin_parameter,
    SchwarzSweep, MultiplicativeSweep, MulticolorSweep, ParallelSweep, schwarz_step!,
    AbstractCoarseSpace, TruncatedPushforwardCoarseSpace, ExactPushforwardCoarseSpace,
    coarse_dimension, coarse_correct!,
    SchwarzProblem, SchwarzIteration, SchwarzCG, SchwarzResult, SchwarzPreconditioner, solve

using ArgCheck: @argcheck
using BlockArrays: BlockVector, mortar, blocks
using Graphs: SimpleGraph, add_edge!, has_edge, neighbors, ne, nv, vertices, degree
using LinearAlgebra
using LinearAlgebra: ldlt!, RowMaximum
using SparseArrays
using CliqueTrees.Multifrontal: ChordalLDLt
using Krylov: cg
import CommonSolve: solve

using ..SheafInterface: add_sheaf_edge!
using ..EuclideanSheaves: EuclideanSheaf
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
    OverlapCover(subdomains, ndofs=maximum dof index)

The combinatorics of an overlapping cover ``\\{\\Omega_i\\}`` of the dofs
`1:ndofs`, bundled together:

- `subdomains[i]`: the sorted global dof indices of ``\\Omega_i``;
- `members[k]`: the pairs `(i, ℓ)` with dof `k` the `ℓ`-th dof of ``\\Omega_i``;
- `graph`: the overlap graph, i.e. the 1-skeleton of the nerve of the cover
  (an edge for every nonempty ``\\Omega_i \\cap \\Omega_j``);
- `overlaps[i => j]`: the local indices in ``\\Omega_i`` of
  ``\\Omega_i \\cap \\Omega_j``, ordered so that `overlaps[i => j]` and
  `overlaps[j => i]` address the same dofs entry by entry. These are the
  restriction maps of [`overlap_sheaf`](@ref) in index form.
"""
struct OverlapCover
    subdomains::Vector{Vector{Int}}
    members::Vector{Vector{Tuple{Int,Int}}}
    graph::SimpleGraph{Int}
    overlaps::Dict{Pair{Int,Int},Vector{Int}}
end

function OverlapCover(subdomains::AbstractVector{<:AbstractVector{<:Integer}},
                      ndofs::Integer=maximum(d -> isempty(d) ? 0 : maximum(d), subdomains; init=0))
    @argcheck !isempty(subdomains) "need at least one subdomain"
    doms = [sort!(unique(Vector{Int}(d))) for d in subdomains]
    for (i, d) in enumerate(doms)
        @argcheck isempty(d) || (1 <= first(d) && last(d) <= ndofs) "subdomain $i has dof indices outside 1:$ndofs"
    end
    members = [Tuple{Int,Int}[] for _ in 1:ndofs]
    for (i, d) in enumerate(doms), (ℓ, k) in enumerate(d)
        push!(members[k], (i, ℓ))
    end
    graph = SimpleGraph(length(doms))
    overlaps = Dict{Pair{Int,Int},Vector{Int}}()
    for m in members, a in eachindex(m), b in (a + 1):lastindex(m)
        (i, ℓi), (j, ℓj) = m[a], m[b]
        add_edge!(graph, i, j)
        push!(get!(overlaps, i => j, Int[]), ℓi)
        push!(get!(overlaps, j => i, Int[]), ℓj)
    end
    return OverlapCover(doms, members, graph, overlaps)
end

_contains(cover::OverlapCover, i::Int, k::Int) = any(p -> first(p) == i, cover.members[k])

"""
    overlap_sheaf(cover::OverlapCover, T=Float64) -> EuclideanSheaf{T}
    overlap_sheaf(subdomains, T=Float64) -> EuclideanSheaf{T}
    overlap_sheaf(dd::SchwarzDecomposition) -> EuclideanSheaf

The cellular sheaf of an overlapping cover ``\\{\\Omega_i\\}`` of a set of
degrees of freedom (dofs).

- **Vertices** are subdomains, with stalk ``F(i) = \\mathbb{R}^{\\Omega_i}``
  (one coordinate per dof of ``\\Omega_i``, in increasing global order).
- **Edges** are the nonempty pairwise overlaps ``\\Omega_i \\cap \\Omega_j``,
  with stalk ``F(ij) = \\mathbb{R}^{\\Omega_i \\cap \\Omega_j}``.
- **Restriction maps** ``F(i) \\to F(ij)`` are coordinate projections: they
  read off the values a subdomain holds on the overlap.

The underlying graph is the 1-skeleton of the nerve of the cover. A 0-cochain
``x = (x_i)`` assigns each subdomain its own local function; the coboundary
``(\\delta x)_{ij} = x_i|_{\\Omega_i\\cap\\Omega_j} - x_j|_{\\Omega_i\\cap\\Omega_j}``
measures how much neighbouring subdomains disagree. Because every nonempty
overlap is an edge, the global sections ``H^0`` are exactly the cochains coming
from a single function on ``\\bigcup_i \\Omega_i``, so
``\\dim H^0 = |\\bigcup_i \\Omega_i|``.

This is the sheaf on which the Schwarz methods of this module run: iterates are
0-cochains, and convergence means they become a global section.

The restriction maps are stored densely (as `EuclideanSheaf` requires), so build
this for inspection and small problems; the solvers use the index form held by
[`OverlapCover`](@ref).
"""
function overlap_sheaf(cover::OverlapCover, ::Type{T}=Float64) where {T}
    doms = cover.subdomains
    s = EuclideanSheaf{T}(length.(doms))
    for i in eachindex(doms), j in neighbors(cover.graph, i)
        i < j || continue
        add_sheaf_edge!(s, i, j,
            _selection(T, cover.overlaps[i => j], length(doms[i])),
            _selection(T, cover.overlaps[j => i], length(doms[j])))
    end
    return s
end

overlap_sheaf(subdomains::AbstractVector{<:AbstractVector{<:Integer}}, ::Type{T}=Float64) where {T} =
    overlap_sheaf(OverlapCover(subdomains), T)

function _selection(::Type{T}, rows::Vector{Int}, ncols::Int) where {T}
    P = zeros(T, length(rows), ncols)
    for (r, c) in enumerate(rows)
        P[r, c] = one(T)
    end
    return P
end

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
    S = dropzeros(sparse(A))
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

Which subdomain's copy supplies each dof: `owner[k]` is a subdomain containing
dof `k`, and `local_index[k]` is the position of `k` in it. The owner's copy is
read whenever another subdomain needs `k` as boundary data, and when a cochain
is glued into a global vector ([`glue`](@ref)).

Given the matrix instead of an owner vector, each dof is owned by a subdomain
that contains it together with all of its matrix neighbours, when one exists.
Any subdomain reading `k` as boundary data then shares a neighbour of `k` with
the owner, so the data travels along an edge of the overlap graph.
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
        idx = findfirst(p -> first(p) == owner[k], cover.members[k])
        @argcheck idx !== nothing "owner[$k] = $(owner[k]), but subdomain $(owner[k]) does not contain dof $k"
        local_index[k] = last(cover.members[k][idx])
    end
    return Ownership(Vector{Int}(owner), local_index)
end

function Ownership(cover::OverlapCover, S::SparseMatrixCSC)
    owner = map(eachindex(cover.members)) do k
        candidates = first.(cover.members[k])
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
- `source_dof`: the local index of ``m`` in ``\\Omega_j``.

The face's Robin data is read entirely from the neighbour ``j``: its values
at ``k`` (through the coupling block) and at ``m`` (through `weight`). Together
they form the discrete ``(\\partial_n + p_{ij})\\, u_j`` on that face. An
interface dof at a subdomain corner has one face per outside neighbour, each
with its own neighbour, parameter and data.
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
- `boundary`: its discrete boundary
  ``\\Gamma_i = \\{k \\notin \\Omega_i : A_{km} \\neq 0 \\text{ for some } m \\in \\Omega_i\\}``;
- `coupling`: the block ``A_{\\Omega_i \\Gamma_i}``;
- `faces`: the [`RobinFace`](@ref)s of the interface (empty for
  [`DirichletTransmission`](@ref));
- `factor`: a sparse `ChordalLDLt` factorization of the local matrix
  ``\\tilde A_i``, computed once and reused by every local solve.
"""
struct LocalProblem{T,F}
    dofs::Vector{Int}
    boundary::Vector{Int}
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
    boundary = _boundary(S, dofs)
    faces = _robin_faces(asm, i, boundary)
    Ai = S[dofs, dofs]
    if !isempty(faces)
        local_dofs = [face.dof for face in faces]
        Ai = Ai + sparse(local_dofs, local_dofs, [face.weight for face in faces], length(dofs), length(dofs))
    end
    factor = ldlt!(ChordalLDLt(Ai), RowMaximum(); check=false)
    @argcheck all(>(0), factor.D.diag) "the local matrix of subdomain $i is not positive definite; increase the Robin parameter"
    return LocalProblem(dofs, boundary, S[dofs, boundary], faces, factor)
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
        source = findfirst(q -> first(q) == j, asm.cover.members[m])
        src, src_dof = source === nothing ?
            (asm.ownership.owner[m], asm.ownership.local_index[m]) : asm.cover.members[m][source]
        push!(faces, RobinFace(ℓ, p - abs(nonzeros(S)[idx]), src, src_dof))
    end
    return faces
end

# ===== Decomposition =====

"""
    SchwarzDecomposition(A, subdomains; owner=nothing, transmission=DirichletTransmission())

An overlapping domain decomposition of the sparse symmetric positive-definite
system ``A u = f`` (typically a finite-difference or finite-element
discretization of an elliptic PDE), prepared for Schwarz iteration.
`subdomains[i]` lists the global dofs of ``\\Omega_i``; together they must
cover `1:size(A, 1)`. The decomposition bundles

- `A`, as a sparse matrix;
- `cover`: the [`OverlapCover`](@ref), i.e. the overlap graph and its
  restriction maps (see [`overlap_sheaf`](@ref));
- `ownership`: the [`Ownership`](@ref) of each dof, from the vector `owner`
  if given and from the sparsity of `A` otherwise;
- `locals`: one [`LocalProblem`](@ref) per subdomain, built with the given
  [`TransmissionCondition`](@ref);
- `colors`: a coloring of the subdomains for [`MulticolorSweep`](@ref);
- `transmission`: the [`TransmissionCondition`](@ref) the local problems were built with.

Every boundary dof ``k \\in \\Gamma_i`` must be owned by a subdomain that
overlaps ``\\Omega_i``, so that all communication runs along edges of the
overlap graph. A decomposition built with `overlap >= 1` by
[`overlapping_subdomains`](@ref) always satisfies this.

Two subdomains get different colors when they *conflict*: they overlap, or
one's boundary ``\\Gamma_i`` meets the other subdomain (they are coupled by
`A`). Subdomains of one color neither read nor write each other's data, so
they can be solved concurrently. The coloring is greedy, largest conflict
degree first.
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
    @argcheck issymmetric(S) "A must be symmetric"

    cover = OverlapCover(subdomains, n)
    for (i, d) in enumerate(cover.subdomains)
        @argcheck !isempty(d) "subdomain $i is empty"
    end
    uncovered = findfirst(isempty, cover.members)
    @argcheck uncovered === nothing "dof $uncovered lies in no subdomain; the subdomains must cover 1:$n"

    ownership = owner === nothing ? Ownership(cover, S) : Ownership(cover, owner)
    asm = _Assembly(S, cover, ownership, transmission)
    locals = [LocalProblem(asm, i) for i in eachindex(cover.subdomains)]
    for (i, lp) in enumerate(locals), k in lp.boundary
        j = ownership.owner[k]
        @argcheck has_edge(cover.graph, i, j) "subdomain $i needs boundary dof $k from its owner, subdomain $j, but the two do not overlap; increase the overlap or choose a different owner"
    end
    colors = _greedy_coloring(_conflict_graph(cover, locals))
    return SchwarzDecomposition(S, cover, ownership, locals, colors, transmission)
end

# Subdomains conflict when they overlap or when one's Dirichlet boundary lies
# in the other: either way a concurrent update could change data the other
# reads or writes.
function _conflict_graph(cover::OverlapCover, locals::Vector{<:LocalProblem})
    conflicts = copy(cover.graph)
    for (i, lp) in enumerate(locals), k in lp.boundary, (j, _) in cover.members[k]
        i == j || add_edge!(conflicts, i, j)
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
    @argcheck length(xs) == length(dd.locals) && all(length.(xs) .== length.(dd.cover.subdomains)) "cochain blocks must match the subdomain sizes"
    return xs
end

function _cochain_blocks(dd::SchwarzDecomposition, x::AbstractVector)
    sizes = length.(dd.cover.subdomains)
    @argcheck length(x) == sum(sizes) "cochain length must equal the total size of the subdomains"
    offsets = [0; cumsum(sizes)]
    return [x[offsets[i]+1:offsets[i+1]] for i in eachindex(sizes)]
end

_owned_value(dd::SchwarzDecomposition, xs, k::Int) =
    xs[dd.ownership.owner[k]][dd.ownership.local_index[k]]

"""
    localize(dd::SchwarzDecomposition, u) -> BlockVector

Restrict a global vector `u` to every subdomain, giving the 0-cochain
``(u|_{\\Omega_1}, \\dots, u|_{\\Omega_N})`` of the overlap sheaf. This is the
isomorphism from functions on the domain onto the global sections ``H^0``.
"""
function localize(dd::SchwarzDecomposition, u::AbstractVector)
    @argcheck length(u) == size(dd.A, 1)
    return mortar([u[d] for d in dd.cover.subdomains])
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

# Subproblem on Ω_i: Ã_i x_i = f|Ω_i − A_{Ω_i Γ_i} g_i + Σ_faces w x_j(m), where
# g_i is read from the owners' copies and each Robin face reads the interface
# value from the neighbour across it. The face sum vanishes for Dirichlet
# transmission.
function _local_solve(prob::SchwarzProblem, xs, i::Int)
    dd = prob.decomposition
    lp = dd.locals[i]
    g = [_owned_value(dd, xs, k) for k in lp.boundary]
    b = prob.rhs[lp.dofs] - lp.coupling * g
    for face in lp.faces
        b[face.dof] += face.weight * xs[face.source][face.source_dof]
    end
    return _ldlt_solve(lp.factor, b)
end

# Local solve on Ω_i, then publish the result to the neighbours.
function _solve_and_publish!(xs, prob::SchwarzProblem, i::Int)
    xs[i] .= _local_solve(prob, xs, i)
    _publish!(xs, prob.decomposition.transmission, prob.decomposition.cover, i)
    return nothing
end

# Dirichlet transmission: push the new values along the restriction maps onto
# every overlapping neighbour, so a section stays a section.
function _publish!(xs, ::DirichletTransmission, cover::OverlapCover, i::Int)
    for j in neighbors(cover.graph, i)
        xs[j][cover.overlaps[j => i]] .= view(xs[i], cover.overlaps[i => j])
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
    updates = [_local_solve(prob, xs, i) for i in eachindex(xs)]
    foreach(copyto!, xs, updates)
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
    n = size(dd.A, 1)
    Z = reshape(modes, size(modes, 1), :)
    @argcheck size(Z, 1) == n "modes must have one row per dof"
    multiplicity = length.(dd.cover.members)

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
    aggregates = [reduce(vcat, (dd.cover.subdomains[i] for i in fiber_vertices(hom, h)))
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
    for (xi, d) in zip(xs, dd.cover.subdomains)
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
        corrections[i] = _ldlt_solve(lp.factor, r[lp.dofs])
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

end
