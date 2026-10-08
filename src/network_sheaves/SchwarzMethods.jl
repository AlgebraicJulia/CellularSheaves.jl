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

export SchwarzDecomposition, SchwarzResult, overlap_sheaf, overlapping_subdomains,
    AbstractCoarseSpace, TruncatedPushforwardCoarseSpace, ExactPushforwardCoarseSpace,
    coarse_dimension, coarse_correct!,
    localize, glue, overlap_disagreement, schwarz_step!, schwarz_solve

using ArgCheck: @argcheck
using BlockArrays: BlockVector, mortar, blocks
using Graphs: SimpleGraph, add_edge!, has_edge, neighbors, ne, nv, vertices, degree
using LinearAlgebra
using LinearAlgebra: ldlt!, RowMaximum
using SparseArrays
using CliqueTrees.Multifrontal: ChordalLDLt

using ..SheafInterface: add_sheaf_edge!
using ..EuclideanSheaves: EuclideanSheaf
using ..GraphHomomorphisms: GraphHomomorphism, fiber_vertices

# ===== Cover combinatorics =====

# members[k] lists (i, ℓ) for every subdomain i containing dof k, where ℓ is
# the local index of k in Ω_i. Subdomain ids within each list are increasing.
function _memberships(doms::Vector{Vector{Int}}, n::Int)
    members = [Tuple{Int,Int}[] for _ in 1:n]
    for (i, d) in enumerate(doms), (ℓ, k) in enumerate(d)
        push!(members[k], (i, ℓ))
    end
    return members
end

# The 1-skeleton of the nerve of the cover, together with the restriction maps
# in index form: overlaps[i => j] holds the local indices (in Ω_i) of
# Ω_i ∩ Ω_j, ordered by global index so that overlaps[i => j] and
# overlaps[j => i] address the same dofs entry by entry.
function _overlap_structure(members::Vector{Vector{Tuple{Int,Int}}}, N::Int)
    graph = SimpleGraph(N)
    overlaps = Dict{Pair{Int,Int},Vector{Int}}()
    for m in members, a in eachindex(m), b in (a + 1):lastindex(m)
        (i, ℓi), (j, ℓj) = m[a], m[b]
        add_edge!(graph, i, j)
        push!(get!(overlaps, i => j, Int[]), ℓi)
        push!(get!(overlaps, j => i, Int[]), ℓj)
    end
    return graph, overlaps
end

_normalize_subdomains(subdomains) = [sort!(unique(Vector{Int}(d))) for d in subdomains]

# Structural neighbours of dof k in the (symmetric) matrix S.
_adjacent(S::SparseMatrixCSC, k::Int) =
    (rowvals(S)[p] for p in nzrange(S, k) if rowvals(S)[p] != k && !iszero(nonzeros(S)[p]))

# Γ_i: dofs outside Ω_i that A couples to Ω_i (the discrete Dirichlet boundary).
function _boundary(S::SparseMatrixCSC, dom::Vector{Int})
    inside = falses(size(S, 1))
    inside[dom] .= true
    Γ = Int[]
    for c in dom, r in _adjacent(S, c)
        inside[r] || push!(Γ, r)
    end
    return sort!(unique!(Γ))
end

# Prefer, for each dof, a subdomain that contains it *together with* all of its
# A-neighbours. Any subdomain reading k as Dirichlet data then shares an
# A-neighbour of k with the owner, so the data travels along an edge of the
# overlap graph.
function _default_owner(S::SparseMatrixCSC, members)
    n = length(members)
    owner = Vector{Int}(undef, n)
    for k in 1:n
        owner[k] = first(first(members[k]))
        for (i, _) in members[k]
            if all(r -> any(p -> first(p) == i, members[r]), _adjacent(S, k))
                owner[k] = i
                break
            end
        end
    end
    return owner
end

# Subdomains conflict when they overlap or when one's Dirichlet boundary lies
# in the other: either way a concurrent update could change data the other
# reads or writes.
function _conflict_graph(graph::SimpleGraph{Int}, boundary, members)
    conflicts = copy(graph)
    for (i, Γ) in enumerate(boundary), k in Γ, (j, _) in members[k]
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

# Solve M v = b for a ChordalLDLt factor with X = P' L D L' P.
function _ldlt_solve(M, b::AbstractVector)
    c = M.P' \ b
    z = M.L \ c
    w = M.D \ z
    y = M.L' \ w
    return M.P \ y
end

# ===== The overlap sheaf =====

"""
    overlap_sheaf(subdomains, T=Float64) -> EuclideanSheaf{T}
    overlap_sheaf(dd::SchwarzDecomposition) -> EuclideanSheaf

The cellular sheaf of an overlapping cover ``\\{\\Omega_i\\}`` of a set of
degrees of freedom (dofs). Each `subdomains[i]` is a collection of global dof
indices.

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

This is the sheaf on which [`schwarz_solve`](@ref) runs: Schwarz iterates are
0-cochains, and convergence means they become a global section.

The restriction maps are stored densely (as `EuclideanSheaf` requires), so build
this for inspection and small problems; the solver itself uses an index
representation of the same maps.
"""
function overlap_sheaf(subdomains::AbstractVector{<:AbstractVector{<:Integer}}, ::Type{T}=Float64) where {T}
    @argcheck !isempty(subdomains) "need at least one subdomain"
    doms = _normalize_subdomains(subdomains)
    @argcheck all(d -> isempty(d) || first(d) >= 1, doms) "dof indices must be positive"
    n = maximum(d -> isempty(d) ? 0 : last(d), doms)
    graph, overlaps = _overlap_structure(_memberships(doms, n), length(doms))
    return _overlap_sheaf(doms, graph, overlaps, T)
end

function _overlap_sheaf(doms, graph, overlaps, ::Type{T}) where {T}
    s = EuclideanSheaf{T}(length.(doms))
    for i in eachindex(doms), j in neighbors(graph, i)
        i < j || continue
        add_sheaf_edge!(s, i, j,
            _selection(T, overlaps[i => j], length(doms[i])),
            _selection(T, overlaps[j => i], length(doms[j])))
    end
    return s
end

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
(Smith–Bjørstad–Gropp 1996, §1.3). `overlap = 0` returns the partition itself
(block Jacobi / block Gauss–Seidel when used with [`schwarz_solve`](@ref)).

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
    N = maximum(parts)
    subdomains = Vector{Vector{Int}}(undef, N)
    for i in 1:N
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

# ===== Decomposition =====

"""
    SchwarzDecomposition(A, subdomains; owner=nothing)

An overlapping domain decomposition of the sparse symmetric positive-definite
system ``A u = f`` (typically a finite-difference or finite-element
discretization of an elliptic PDE), prepared for Schwarz iteration.

`subdomains[i]` lists the global dof indices of ``\\Omega_i``; together they
must cover `1:size(A, 1)`. The decomposition stores

- the overlap graph and its restriction maps (see [`overlap_sheaf`](@ref)):
  vertices are subdomains, edges are nonempty overlaps;
- for each subdomain its *discrete boundary*
  ``\\Gamma_i = \\{k \\notin \\Omega_i : A_{km} \\neq 0 \\text{ for some } m \\in \\Omega_i\\}``
  and the coupling block ``A_{\\Omega_i \\Gamma_i}``;
- a sparse `ChordalLDLt` factorization of the local stiffness matrix
  ``A_i = A_{\\Omega_i \\Omega_i}`` (SPD, as a principal submatrix of an SPD
  matrix), computed once and reused by every local solve.

`owner[k]` names the subdomain whose copy supplies dof `k` when another
subdomain needs it as Dirichlet data, and when a cochain is glued back to a
global vector ([`glue`](@ref)). By default each dof is owned by a subdomain
containing it together with all of its matrix neighbours, when one exists.

Every boundary dof ``k \\in \\Gamma_i`` must be owned by a subdomain that
overlaps ``\\Omega_i``, so that all communication runs along edges of the
overlap graph; a decomposition built with `overlap >= 1` by
[`overlapping_subdomains`](@ref) always satisfies this.

`colors` partitions the subdomains into color classes for the `:multicolor`
sweep of [`schwarz_step!`](@ref). Two subdomains get different colors when they
*conflict*: they overlap (an edge of the overlap graph), or one's boundary
``\\Gamma_i`` meets the other subdomain (they are coupled by `A`). Subdomains
of one color neither read nor write each other's data, so they can be solved
concurrently. The coloring is greedy, largest conflict degree first.
"""
struct SchwarzDecomposition{T,F}
    A::SparseMatrixCSC{T,Int}
    subdomains::Vector{Vector{Int}}
    owner::Vector{Int}
    owner_local::Vector{Int}
    graph::SimpleGraph{Int}
    overlaps::Dict{Pair{Int,Int},Vector{Int}}
    boundary::Vector{Vector{Int}}
    coupling::Vector{SparseMatrixCSC{T,Int}}
    factors::Vector{F}
    colors::Vector{Vector{Int}}
end

function SchwarzDecomposition(A::AbstractMatrix, subdomains::AbstractVector{<:AbstractVector{<:Integer}};
                              owner::Union{Nothing,AbstractVector{<:Integer}}=nothing)
    n = size(A, 1)
    @argcheck size(A, 2) == n "A must be square"
    S = dropzeros(sparse(float.(A)))
    @argcheck issymmetric(S) "A must be symmetric"
    @argcheck !isempty(subdomains) "need at least one subdomain"

    doms = _normalize_subdomains(subdomains)
    for (i, d) in enumerate(doms)
        @argcheck !isempty(d) "subdomain $i is empty"
        @argcheck 1 <= first(d) && last(d) <= n "subdomain $i has dof indices outside 1:$n"
    end
    members = _memberships(doms, n)
    uncovered = findfirst(isempty, members)
    @argcheck uncovered === nothing "dof $uncovered lies in no subdomain; the subdomains must cover 1:$n"
    graph, overlaps = _overlap_structure(members, length(doms))

    own = owner === nothing ? _default_owner(S, members) : Vector{Int}(owner)
    @argcheck length(own) == n "need one owner per dof"
    owner_local = Vector{Int}(undef, n)
    for k in 1:n
        idx = findfirst(p -> first(p) == own[k], members[k])
        @argcheck idx !== nothing "owner[$k] = $(own[k]), but subdomain $(own[k]) does not contain dof $k"
        owner_local[k] = last(members[k][idx])
    end

    boundary = [_boundary(S, d) for d in doms]
    for (i, Γ) in enumerate(boundary), k in Γ
        @argcheck has_edge(graph, i, own[k]) "subdomain $i needs boundary dof $k from its owner, subdomain $(own[k]), but the two do not overlap; increase the overlap or choose a different owner"
    end

    coupling = [S[d, Γ] for (d, Γ) in zip(doms, boundary)]
    factors = [ldlt!(ChordalLDLt(S[d, d]), RowMaximum()) for d in doms]
    colors = _greedy_coloring(_conflict_graph(graph, boundary, members))
    return SchwarzDecomposition{eltype(S),eltype(factors)}(
        S, doms, own, owner_local, graph, overlaps, boundary, coupling, factors, colors)
end

overlap_sheaf(dd::SchwarzDecomposition{T}) where {T} =
    _overlap_sheaf(dd.subdomains, dd.graph, dd.overlaps, T)

function Base.show(io::IO, dd::SchwarzDecomposition)
    print(io, "SchwarzDecomposition(", length(dd.subdomains), " subdomains, ",
        size(dd.A, 1), " dofs, ", ne(dd.graph), " overlaps)")
end

# ===== Cochains =====

function _cochain_blocks(dd::SchwarzDecomposition, x::BlockVector)
    xs = blocks(x)
    @argcheck length(xs) == length(dd.subdomains) && all(length.(xs) .== length.(dd.subdomains)) "cochain blocks must match the subdomain sizes"
    return xs
end

function _cochain_blocks(dd::SchwarzDecomposition, x::AbstractVector)
    sizes = length.(dd.subdomains)
    @argcheck length(x) == sum(sizes) "cochain length must equal the total size of the subdomains"
    offsets = [0; cumsum(sizes)]
    return [x[offsets[i]+1:offsets[i+1]] for i in eachindex(sizes)]
end

"""
    localize(dd::SchwarzDecomposition, u) -> BlockVector

Restrict a global vector `u` to every subdomain, giving the 0-cochain
``(u|_{\\Omega_1}, \\dots, u|_{\\Omega_N})`` of the overlap sheaf. This is the
isomorphism from functions on the domain onto the global sections ``H^0``.
"""
function localize(dd::SchwarzDecomposition, u::AbstractVector)
    @argcheck length(u) == size(dd.A, 1)
    return mortar([u[d] for d in dd.subdomains])
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
    return [xs[dd.owner[k]][dd.owner_local[k]] for k in eachindex(dd.owner)]
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
    for ((i, j), idx) in dd.overlaps
        i < j || continue
        acc += sum(abs2, view(xs[i], idx) .- view(xs[j], dd.overlaps[j => i]))
    end
    return sqrt(acc)
end

# ===== Schwarz iteration =====

# Dirichlet subproblem on Ω_i: A_i x_i = f|Ω_i − A_{Ω_i Γ_i} g_i, where g_i
# reads each boundary dof from its owner's copy.
function _local_solve(dd::SchwarzDecomposition, i::Int, xs, f::AbstractVector)
    g = [xs[dd.owner[k]][dd.owner_local[k]] for k in dd.boundary[i]]
    return _ldlt_solve(dd.factors[i], f[dd.subdomains[i]] - dd.coupling[i] * g)
end

# Local solve on Ω_i, then push the result along the restriction maps onto
# every overlapping neighbour.
function _solve_and_push!(xs, dd::SchwarzDecomposition, i::Int, f::AbstractVector)
    xs[i] .= _local_solve(dd, i, xs, f)
    for j in neighbors(dd.graph, i)
        xs[j][dd.overlaps[j => i]] .= view(xs[i], dd.overlaps[i => j])
    end
    return nothing
end

const SCHWARZ_METHODS = (:multiplicative, :multicolor, :parallel)

"""
    schwarz_step!(x::BlockVector, dd::SchwarzDecomposition, f;
                  method=:multiplicative, order=1:N) -> x

One sweep of an overlapping Schwarz method on the 0-cochain `x` of the overlap
sheaf, in place. On subdomain ``i`` the *local solve* is the discrete Dirichlet
problem

```math
A_{\\Omega_i\\Omega_i}\\, x_i = f|_{\\Omega_i} - A_{\\Omega_i \\Gamma_i}\\, g_i,
```

whose boundary data ``g_i`` on ``\\Gamma_i`` is read from neighbouring
subdomains' copies (each boundary dof from its owner, a neighbour in the
overlap graph).

- `method = :multiplicative` is Schwarz's *alternating* method (Schwarz 1870).
  It visits subdomains in `order`. After each local solve it pushes the new
  values through the restriction maps onto every overlapping neighbour,
  ``x_j|_{\\Omega_i\\cap\\Omega_j} \\leftarrow x_i|_{\\Omega_i\\cap\\Omega_j}``.
  Starting from a global section the iterate stays a section, and its glued
  vector is exactly the classical multiplicative Schwarz iterate
  ``u \\leftarrow u + R_i^\\mathsf{T} A_i^{-1} R_i (f - A u)``. It converges for
  every SPD `A`, since it is block Gauss–Seidel over overlapping blocks
  (Toselli–Widlund 2005, ch. 2).
- `method = :multicolor` is the same alternating method, with the subdomains
  visited one color class of `dd.colors` at a time. The subdomains in a class
  do not conflict, so they are solved and pushed concurrently
  (`Threads.@threads`; start Julia with several threads to benefit). The result
  equals `:multiplicative` with `order = reduce(vcat, dd.colors)`, so it keeps
  that method's convergence guarantee while the number of sequential steps per
  sweep drops from the number of subdomains to the number of colors
  (Smith–Bjørstad–Gropp 1996, §1.4).
- `method = :parallel` is Lions' *parallel* Schwarz method (Lions 1988). All
  subdomains solve simultaneously from the previous cochain and no values are
  pushed. Copies on overlaps disagree during the iteration (the cochain is not
  a section) and agree in the limit. The glued iterate coincides with
  restricted additive Schwarz (RAS) for the owner partition
  (Efstathiou–Gander 2003). It converges when `A` is an M-matrix, such as
  standard discretizations of ``-\\Delta`` (Frommer–Szyld 2001), but not for
  every SPD matrix.
"""
function schwarz_step!(x::BlockVector, dd::SchwarzDecomposition, f::AbstractVector;
                       method::Symbol=:multiplicative, order=eachindex(dd.subdomains))
    @argcheck method in SCHWARZ_METHODS "method must be one of $SCHWARZ_METHODS"
    @argcheck length(f) == size(dd.A, 1)
    xs = _cochain_blocks(dd, x)
    if method === :parallel
        updates = [_local_solve(dd, i, xs, f) for i in eachindex(xs)]
        foreach(copyto!, xs, updates)
    elseif method === :multicolor
        for class in dd.colors
            Threads.@threads for i in class
                _solve_and_push!(xs, dd, i, f)
            end
        end
    else
        for i in order
            _solve_and_push!(xs, dd, i, f)
        end
    end
    return x
end

# ===== Coarse spaces =====

"""
    AbstractCoarseSpace

A second level for [`schwarz_solve`](@ref). After every fine sweep the glued
iterate ``u`` receives a correction ``e`` computed from the residual
``f - A u``, and ``e`` is added to every subdomain's copy (the cochain
`localize(dd, e)`, which leaves the overlap disagreement unchanged).

Both implementations start from a graph homomorphism ``\\varphi : G \\to H``
that groups the subdomains (vertices of the overlap graph ``G``) into
aggregates, and from the pushforward ``\\varphi_* F`` of the overlap sheaf
``F``. The stalk ``(\\varphi_* F)(h)`` is the space of global sections of ``F``
over the fiber ``\\varphi^{-1}(h)``, i.e. functions on the aggregate
``\\widehat\\Omega_h = \\bigcup_{i \\in \\varphi^{-1}(h)} \\Omega_i``.

- [`TruncatedPushforwardCoarseSpace`](@ref) keeps a few modes of each stalk
  and solves a small Galerkin problem. It is cheap and scalable but inexact.
- [`ExactPushforwardCoarseSpace`](@ref) keeps the whole stalk, i.e. it solves
  on the aggregates themselves. It is more accurate per sweep but its cost
  grows with the aggregate size.
"""
abstract type AbstractCoarseSpace end

function _check_hom(dd::SchwarzDecomposition, hom::GraphHomomorphism)
    @argcheck length(hom.vertex_map) == length(dd.subdomains) "the graph homomorphism must have one source vertex per subdomain"
    for h in 1:hom.n_target
        @argcheck !isempty(fiber_vertices(hom, h)) "aggregate $h has an empty fiber; every target vertex must receive a subdomain"
    end
end

_identity_hom(dd::SchwarzDecomposition) = GraphHomomorphism(collect(eachindex(dd.subdomains)))

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

function TruncatedPushforwardCoarseSpace(dd::SchwarzDecomposition{T}, hom::GraphHomomorphism=_identity_hom(dd);
                                         modes::AbstractVecOrMat=ones(T, size(dd.A, 1), 1)) where {T}
    _check_hom(dd, hom)
    n = size(dd.A, 1)
    Z = reshape(modes, size(modes, 1), :)
    @argcheck size(Z, 1) == n "modes must have one row per dof"
    multiplicity = zeros(Int, n)
    for d in dd.subdomains
        multiplicity[d] .+= 1
    end

    I, J, V = Int[], Int[], T[]
    stalks = zeros(Int, hom.n_target)
    col = 0
    for h in 1:hom.n_target
        weight = zeros(T, n)
        for i in fiber_vertices(hom, h)
            weight[dd.subdomains[i]] .+= 1
        end
        support = findall(!iszero, weight)
        weight[support] ./= multiplicity[support]
        for c in axes(Z, 2)
            vals = weight[support] .* Z[support, c]
            nz = findall(!iszero, vals)
            isempty(nz) && continue
            col += 1
            stalks[h] += 1
            append!(I, support[nz])
            append!(J, fill(col, length(nz)))
            append!(V, vals[nz])
        end
    end
    Φ = sparse(I, J, V, n, col)
    A0 = Φ' * dd.A * Φ
    A0 = (A0 + A0') / 2
    factor = ldlt!(ChordalLDLt(A0), RowMaximum())
    return TruncatedPushforwardCoarseSpace{T,typeof(factor)}(hom, Φ, stalks, A0, factor)
end

"""
    ExactPushforwardCoarseSpace(dd, hom; method=:multicolor)

The coarse level given by the full pushforward ``\\varphi_* F`` of the overlap
sheaf along ``\\varphi`` = `hom`. Its stalk at an aggregate ``h`` is all
functions on ``\\widehat\\Omega_h = \\bigcup_{i \\in \\varphi^{-1}(h)} \\Omega_i``,
so its vertex stalks have the same dimensions as those of
`pushforward_sheaf(hom, overlap_sheaf(dd))`. Since no modes are discarded,
``H^0(\\varphi_* F) \\cong H^0(F)``.

Concretely this is a second [`SchwarzDecomposition`](@ref) whose subdomains
are the aggregates ``\\widehat\\Omega_h``. A dof is owned by the aggregate
containing its fine owner. The correction is one sweep of `method` on that
decomposition, started from the current glued iterate. Each local solve is
exact on a larger region, so a sweep reduces the error far more than a
truncated coarse solve. The cost is in the factorizations: with ``\\varphi``
to a single vertex the coarse level is a direct solve of the whole problem
(one iteration, no scalability). With a fixed aggregation ratio, the iteration
count still grows with the number of aggregates, as for any one-level method.
"""
struct ExactPushforwardCoarseSpace{D<:SchwarzDecomposition} <: AbstractCoarseSpace
    hom::GraphHomomorphism
    decomposition::D
    method::Symbol
end

function ExactPushforwardCoarseSpace(dd::SchwarzDecomposition, hom::GraphHomomorphism; method::Symbol=:multicolor)
    _check_hom(dd, hom)
    @argcheck method in SCHWARZ_METHODS "method must be one of $SCHWARZ_METHODS"
    aggregates = [reduce(vcat, (dd.subdomains[i] for i in fiber_vertices(hom, h)))
                  for h in 1:hom.n_target]
    owner = hom.vertex_map[dd.owner]
    return ExactPushforwardCoarseSpace(hom, SchwarzDecomposition(dd.A, aggregates; owner), method)
end

"""
    coarse_dimension(coarse::AbstractCoarseSpace) -> Int

Number of unknowns the coarse level solves for in each correction: the
dimension of the truncated coarse space, or the total size of the aggregates
(the sum of the vertex stalk dimensions of the pushforward sheaf) for the
exact pushforward.
"""
coarse_dimension(c::TruncatedPushforwardCoarseSpace) = size(c.basis, 2)
coarse_dimension(c::ExactPushforwardCoarseSpace) = sum(length, c.decomposition.subdomains)

function Base.show(io::IO, c::TruncatedPushforwardCoarseSpace)
    print(io, "TruncatedPushforwardCoarseSpace(", c.hom.n_target, " aggregates, dimension ",
        coarse_dimension(c), ")")
end

function Base.show(io::IO, c::ExactPushforwardCoarseSpace)
    print(io, "ExactPushforwardCoarseSpace(", c.hom.n_target, " aggregates, dimension ",
        coarse_dimension(c), ")")
end

_coarse_correction(c::TruncatedPushforwardCoarseSpace, A, f, u) =
    c.basis * _ldlt_solve(c.factor, c.basis' * (f - A * u))

function _coarse_correction(c::ExactPushforwardCoarseSpace, A, f, u)
    dd = c.decomposition
    x = localize(dd, u)
    schwarz_step!(x, dd, f; method=c.method)
    return glue(dd, x) - u
end

"""
    coarse_correct!(x::BlockVector, dd::SchwarzDecomposition, coarse, f) -> x

Apply one coarse-level correction to the 0-cochain `x`, in place. The glued
iterate ``u`` = `glue(dd, x)` gives the residual ``f - A u``. The coarse level
turns it into a global correction ``e``, which is added to every subdomain's
copy: ``x \\leftarrow x + \\mathrm{localize}(e)``. The glued iterate becomes
``u + e``, and the disagreement ``\\|\\delta x\\|`` is unchanged.
"""
function coarse_correct!(x::BlockVector, dd::SchwarzDecomposition, coarse::AbstractCoarseSpace, f::AbstractVector)
    @argcheck length(f) == size(dd.A, 1)
    xs = _cochain_blocks(dd, x)
    e = _coarse_correction(coarse, dd.A, f, glue(dd, x))
    for (xi, d) in zip(xs, dd.subdomains)
        xi .+= view(e, d)
    end
    return x
end

"""
    SchwarzResult

Output of [`schwarz_solve`](@ref).

- `u` — the glued global solution.
- `x` — the final 0-cochain of the overlap sheaf (one local solution per
  subdomain).
- `residuals` — relative residual ``\\|f - A u_k\\| / \\|f\\|`` of the glued
  iterate after each sweep (entry 1 is the initial guess).
- `disagreements` — [`overlap_disagreement`](@ref) of each iterate.
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
    schwarz_solve(dd::SchwarzDecomposition, f; method=:multiplicative, coarse=nothing,
                  u0=zeros(n), tol=1e-8, maxiter=1000, order=1:N) -> SchwarzResult

Solve ``A u = f`` by overlapping Schwarz iteration on the overlap sheaf of `dd`.

The iterate is a 0-cochain ``x = (x_1, \\dots, x_N)``: each subdomain holds its
own copy of the solution on ``\\Omega_i``. It starts from the global section
`localize(dd, u0)`. Each iteration is one [`schwarz_step!`](@ref) sweep with
the given `method` (`:multiplicative`, `:multicolor` or `:parallel`). When a
`coarse` space is given (an [`AbstractCoarseSpace`](@ref)), the sweep is
followed by [`coarse_correct!`](@ref), giving a two-level method. Iteration
stops once the glued iterate satisfies
``\\|f - A u\\| \\le \\mathrm{tol}\\,\\|f\\|``. At the fixed point every local
solve is consistent with its neighbours, so ``x`` is a global section and its
gluing solves ``A u = f``.
"""
function schwarz_solve(dd::SchwarzDecomposition{T}, f::AbstractVector;
                       method::Symbol=:multiplicative,
                       coarse::Union{Nothing,AbstractCoarseSpace}=nothing,
                       u0::AbstractVector=zeros(T, size(dd.A, 1)),
                       tol::Real=1e-8, maxiter::Integer=1000,
                       order=eachindex(dd.subdomains)) where {T}
    @argcheck method in SCHWARZ_METHODS "method must be one of $SCHWARZ_METHODS"
    @argcheck length(f) == size(dd.A, 1)
    @argcheck maxiter >= 0
    x = localize(dd, Vector{T}(u0))
    scale = norm(f)
    scale = iszero(scale) ? one(T) : scale
    relres(u) = T(norm(f - dd.A * u) / scale)

    u = glue(dd, x)
    residuals = [relres(u)]
    disagreements = [T(overlap_disagreement(dd, x))]
    iterations = 0
    while last(residuals) > tol && iterations < maxiter
        schwarz_step!(x, dd, f; method, order)
        coarse === nothing || coarse_correct!(x, dd, coarse, f)
        iterations += 1
        u = glue(dd, x)
        push!(residuals, relres(u))
        push!(disagreements, T(overlap_disagreement(dd, x)))
    end
    return SchwarzResult{T}(u, x, residuals, disagreements, iterations, last(residuals) <= tol)
end

end
