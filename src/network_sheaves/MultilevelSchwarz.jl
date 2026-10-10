"""
    MultilevelSchwarz

When do Schwarz domain decomposition and the pushforward grid hierarchy of
[`GridMultigrid`](@ref CellularSheaves.NetworkSheaves.GridMultigrid) commute?
This module is a reference implementation, with explicit sparse matrices for small
grids, of both orders of composition, so that the statements below can be checked.

**Setting.** A grid of `points` (vertices numbered column-major) carries an operator
``A``. A *cover* by boxes ``U_i`` gives the inclusions ``ι_i : U_i → G``, the
restriction ``R_i = ι_i^*`` to a box and its extension by zero ``R_i^T``, and the
local (zero-Dirichlet) problems ``A_i = R_i A R_i^T``. One step of coarsening is the
aggregation homomorphism ``ψ : G → G_c``
([`aggregation_homomorphism`](@ref CellularSheaves.NetworkSheaves.GridMultigrid.aggregation_homomorphism)),
with the pullback ``E = ψ^*`` ([`prolongation_matrix`](@ref)), the pushforward
transfer ``T`` (the block average, [`transfer_matrix`](@ref)) and the Galerkin coarse
operator ``A_c = T A E``.

**Base change.** Call a box *aligned* ([`is_aligned`](@ref)) when it is a union of
whole blocks, ``U_i = ψ^{-1}(V_i)`` with ``V_i = ψ(U_i)`` ([`coarse_box`](@ref)). Then
the square formed by ``ι_i``, ``ψ``, the coarse inclusion ``j_i : V_i → G_c`` and the
restriction ``ψ_i : U_i → V_i`` is cartesian, and base change says that restriction to
the cover commutes with the pushforward and the pullback, and extension by zero with
both:

```math
R_{V_i} T = T_i R_{U_i}, \\quad R_{U_i} E = E_i R_{V_i}, \\quad
T R_{U_i}^T = R_{V_i}^T T_i, \\quad E R_{V_i}^T = R_{U_i}^T E_i,
```

where ``E_i, T_i`` are the transfers of the box's own grid. Consequently

```math
T_i (R_{U_i} A R_{U_i}^T) E_i = R_{V_i} (T A E) R_{V_i}^T :
```

the Galerkin coarse operator of a local problem is the local problem of the Galerkin
coarse operator on the pushed-forward box. By induction this holds on every level of
a hierarchy when the cover is aligned at every level (box edges on multiples of
``2^L`` for ``L`` halvings). The coarse cover ``\\{ψ(U_i)\\}`` has half the overlap.
For a box that cuts a block the square is not cartesian and the identity fails at
the box boundary.

**The two solvers.** Both are subspace corrections (Xu 1992) over the same
subspaces ``W_{i,ℓ}``, the images of box ``i`` on level ``ℓ``, with the same local
operators by the identity above:

- [`schwarz_of_multigrid`](@ref): a multigrid of each box on its own grid, summed
  over the cover: ``\\sum_i R_i^T \\mathrm{MG}_i R_i``;
- [`multigrid_of_schwarz`](@ref): a multigrid of ``A`` on the pushforward hierarchy
  whose correction on each level is additive Schwarz over the pushed-forward cover.

With additive composition over the levels (BPX, Bramble–Pasciak–Xu 1990; the
multilevel additive Schwarz of Dryja–Widlund) both are the double sum
``\\sum_{i,ℓ}`` of the same corrections, so for an aligned cover they are **the same
operator**: exchanging the sums over boxes and levels is exactly base change. With
multiplicative composition over the levels (V-cycles) they differ. Write
``C_{i,ℓ}`` for box ``i``'s correction on level ``ℓ`` lifted to the fine grid
(the same in both orders). A V-cycle is a noncommutative polynomial
``p(C_0, C_1, …; A)`` of its corrections, built by ``B ← B + C_k (I - A B)``, and
since ``R_i^T X A_i Y R_i = (R_i^T X R_i) A (R_i^T Y R_i)`` for ``A_i = R_i A R_i^T``,

```math
\\texttt{schwarz\\_of\\_multigrid} = \\sum_i p(C_{i,0}, C_{i,1}, …; A), \\qquad
\\texttt{multigrid\\_of\\_schwarz} = p\\Big(\\sum_i C_{i,0}, \\sum_i C_{i,1}, …; A\\Big).
```

A polynomial of degree one (additive composition) commutes with the sum; the
products ``C_{i,ℓ} A C_{j,ℓ'}`` of a V-cycle leave the cross terms ``i ≠ j``, each
through the coupling ``R_i A R_j^T`` between two boxes. So the box-wise V-cycles
equal the global one exactly when no two boxes are coupled (disjoint boxes and no
stencil across their faces, or a single box), and the difference is linear in the
couplings. Already with one level the second of two exact Schwarz sweeps solves
nothing new inside a box (``C_i A C_i = C_i``) but on the whole grid picks up the
neighbours' corrections. See `docs/scripts/multilevel_schwarz_orders.jl`.

Neither construction has a global coarse problem: the coarsest level is still
covered by the pushed-forward boxes, so each box's coarsest problem is local.
"""
module MultilevelSchwarz

export prolongation_matrix, transfer_matrix, box_restriction, coarse_box, is_aligned,
    schwarz_of_multigrid, multigrid_of_schwarz

using ArgCheck: @argcheck
using LinearAlgebra
using SparseArrays
using ..GridMultigrid: aggregation_homomorphism, coarse_points, _block

"""
    prolongation_matrix(points, factors) -> SparseMatrixCSC

The pullback ``E = ψ^*`` along the aggregation homomorphism of a grid of `points`
coarsened by `factors` (see
[`coarsening_factors`](@ref CellularSheaves.NetworkSheaves.GridMultigrid.coarsening_factors)):
``E_{pb} = 1`` when fine point ``p`` lies in block ``b``. It copies each coarse value
onto its block.
"""
function prolongation_matrix(points::NTuple{D,Integer}, factors::NTuple{D,Integer}) where {D}
    ψ = aggregation_homomorphism(points, factors)
    N = length(ψ.vertex_map)
    return sparse(1:N, ψ.vertex_map, ones(N), N, ψ.n_target)
end

"""
    transfer_matrix(points, factors) -> SparseMatrixCSC

The pushforward transfer ``T`` of the constant sheaf along the aggregation
homomorphism: the average over each block, ``T = (E^T E)^{-1} E^T`` with
``E`` = [`prolongation_matrix`](@ref). ``T E = I``.
"""
function transfer_matrix(points::NTuple{D,Integer}, factors::NTuple{D,Integer}) where {D}
    E = prolongation_matrix(points, factors)
    return sparse(Diagonal(1 ./ vec(sum(E; dims=1))) * sparse(E'))
end

"""
    box_restriction(points, box) -> SparseMatrixCSC

The restriction ``R = ι^*`` from a grid of `points` to the box `box` (a range of
points per dimension), rows numbered column-major in the box. ``R^T`` is the
extension by zero.
"""
function box_restriction(points::NTuple{D,Integer}, box::NTuple{D,AbstractUnitRange{<:Integer}}) where {D}
    @argcheck all(k -> !isempty(box[k]) && first(box[k]) >= 1 && last(box[k]) <= points[k], 1:D) "a box is a nonempty range of grid points in every dimension"
    L = LinearIndices(Tuple(points))
    cols = vec([L[I] for I in CartesianIndices(box)])
    m = length(cols)
    return sparse(1:m, cols, ones(m), m, prod(points))
end

"""
    coarse_box(box, factors) -> NTuple

The image ``ψ(U)`` of the box `box` under the aggregation homomorphism with
`factors`: the blocks it meets.
"""
coarse_box(box::NTuple{D,AbstractUnitRange{<:Integer}}, factors::NTuple{D,Integer}) where {D} =
    ntuple(k -> _block(first(box[k]), factors[k]):_block(last(box[k]), factors[k]), D)

"""
    is_aligned(points, box, factors) -> Bool

Whether the box `box` of a grid of `points` is a union of whole blocks of the
aggregation with `factors`, ``U = ψ^{-1}(ψ(U))``: in every coarsened dimension it
starts at the first point of a block (odd) and ends at the last (even, or the last
grid point, a block of one when the count is odd). Then restriction to the box
commutes with the transfers (see the module documentation).
"""
is_aligned(points::NTuple{D,Integer}, box::NTuple{D,AbstractUnitRange{<:Integer}}, factors::NTuple{D,Integer}) where {D} =
    all(k -> factors[k] == 1 || (isodd(first(box[k])) && (iseven(last(box[k])) || last(box[k]) == points[k])), 1:D)

_jacobi(M) = Diagonal(1 ./ diag(M))
_exact(M) = inv(Matrix(M))

# The preconditioner of the corrections Cₖ applied in order from x = 0,
# x ← x + Cₖ (b - A x): the sum for additive composition, B ← B + Cₖ (I - A B)
# for multiplicative.
function _compose(A, corrections, composition)
    composition === :additive && return sum(corrections)
    B = zeros(size(A))
    for C in corrections
        B = B + C * (I - A * B)
    end
    return B
end

# The multilevel preconditioner of A on a grid of `points` coarsened by each of
# `factors` in turn, with Galerkin coarse operators; level_solve(ℓ, Aℓ) is the
# correction on level ℓ (0 the finest). Multiplicative composition is a V-cycle:
# levels 0, 1, …, L, …, 1, 0.
function _multilevel(A, points, factors, level_solve, composition)
    N = prod(points)
    sizes = [Tuple(Int.(points))]
    P, Q, As = [sparse(1.0I, N, N)], [sparse(1.0I, N, N)], Any[A]
    for r in factors
        E, T = prolongation_matrix(sizes[end], r), transfer_matrix(sizes[end], r)
        push!(P, P[end] * E)
        push!(Q, T * Q[end])
        push!(As, T * As[end] * E)
        push!(sizes, coarse_points(sizes[end], r))
    end
    L = length(factors)
    C = [Matrix(P[ℓ + 1] * level_solve(ℓ, As[ℓ + 1]) * Q[ℓ + 1]) for ℓ in 0:L]
    order = composition === :additive ? collect(0:L) : [0:L; (L - 1):-1:0]
    return _compose(A, [C[ℓ + 1] for ℓ in order], composition)
end

function _check(A, points, boxes, composition)
    @argcheck size(A) == (prod(points), prod(points)) "A must act on the grid of $(points) points"
    @argcheck !isempty(boxes) "the cover needs at least one box"
    @argcheck composition in (:additive, :multiplicative) "composition must be :additive or :multiplicative"
end

"""
    schwarz_of_multigrid(A, points, boxes, factors; composition = :additive,
        smoother = Jacobi, coarsest_solve = exact) -> Matrix

The preconditioner ``\\sum_i R_i^T \\mathrm{MG}_i R_i``: on every box ``U_i`` of
`boxes`, a multigrid of the local problem ``A_i = R_i A R_i^T`` on the box's own
grid, coarsened by each of `factors` in turn with Galerkin coarse operators, then
summed over the cover. Each level's correction is `smoother(Aℓ)` (an approximate
inverse, Jacobi by default), the coarsest `coarsest_solve(Aℓ)` (exact by default);
`composition` combines the levels additively or as a V-cycle (`:multiplicative`).

For a cover aligned at every level ([`is_aligned`](@ref)) and additive composition
this equals [`multigrid_of_schwarz`](@ref) (see the module documentation). A dense
reference implementation for small grids.
"""
function schwarz_of_multigrid(A::AbstractMatrix, points::NTuple{D,Integer}, boxes,
        factors::AbstractVector{<:NTuple{D,Integer}}; composition::Symbol = :additive,
        smoother = _jacobi, coarsest_solve = _exact) where {D}
    _check(A, points, boxes, composition)
    L = length(factors)
    solve(ℓ, M) = ℓ < L ? smoother(M) : coarsest_solve(M)
    B = zeros(size(A))
    for box in boxes
        R = box_restriction(points, box)
        B += R' * _multilevel(R * A * R', length.(box), factors, solve, composition) * R
    end
    return B
end

"""
    multigrid_of_schwarz(A, points, boxes, factors; composition = :additive,
        smoother = Jacobi, coarsest_solve = exact) -> Matrix

The multigrid of `A` on the pushforward hierarchy (coarsened by each of `factors` in
turn, Galerkin coarse operators ``A_ℓ``) whose correction on level ``ℓ`` is additive
Schwarz over the pushed-forward cover ``\\{U_{i,ℓ}\\}`` (``U_{i,ℓ+1} = ψ_ℓ(U_{i,ℓ})``,
[`coarse_box`](@ref)): ``\\sum_i R_{i,ℓ}^T S(R_{i,ℓ} A_ℓ R_{i,ℓ}^T) R_{i,ℓ}``, with
``S`` = `smoother` below the coarsest level and `coarsest_solve` on it.
`composition` combines the levels additively or as a V-cycle (`:multiplicative`).

For a cover aligned at every level ([`is_aligned`](@ref)) and additive composition
this equals [`schwarz_of_multigrid`](@ref) (see the module documentation). A dense
reference implementation for small grids.
"""
function multigrid_of_schwarz(A::AbstractMatrix, points::NTuple{D,Integer}, boxes,
        factors::AbstractVector{<:NTuple{D,Integer}}; composition::Symbol = :additive,
        smoother = _jacobi, coarsest_solve = _exact) where {D}
    _check(A, points, boxes, composition)
    L = length(factors)
    sizes = [Tuple(Int.(points))]
    covers = [collect(boxes)]
    for r in factors
        push!(covers, [coarse_box(U, r) for U in covers[end]])
        push!(sizes, coarse_points(sizes[end], r))
    end
    function solve(ℓ, M)
        S = ℓ < L ? smoother : coarsest_solve
        return sum(covers[ℓ + 1]) do U
            R = box_restriction(sizes[ℓ + 1], U)
            Matrix(R' * S(R * M * R') * R)
        end
    end
    return _multilevel(A, points, factors, solve, composition)
end

end # module
