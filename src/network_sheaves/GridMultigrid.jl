"""
    GridMultigrid

Hierarchies of structured grids for multigrid solvers of the implicit stencil
operators of [`GridSchwarz`](@ref CellularSheaves.NetworkSheaves.GridSchwarz),
with the maps between levels written as a pushforward–pullback pair.

**Coarsening.** Each dimension ``k`` of a grid is coarsened by a factor
``r_k \\in \\{1, 2\\}`` ([`coarsening_factors`](@ref)): with ``r_k = 2`` the
points ``2b - 1, 2b`` form block ``b`` (an odd count leaves a last block of one
point), so the coarse grid has ``\\lceil n_k / 2 \\rceil`` points; with
``r_k = 1`` the dimension is not coarsened. Every coarse dimension keeps at
least three points. A non-periodic dimension is
coarsened while it has at least five points: the coarse points are the block
centres (spacing ``2h``, the first at ``\\mathrm{lower} + h/2``) for an even
count, and every other point (same ends, spacing ``2h``) for an odd one, so
every level is a uniform grid. A periodic dimension is coarsened only while
its count is a multiple of 4 (and at least 8), which keeps every coarse
level's count even, as
red–black Gauss–Seidel around the circle needs (otherwise it is left alone:
semicoarsening). Halving the grid at every level gives ``O(\\log n)`` levels
on grids of ``n`` points per dimension, whatever the counts. The map sending
each fine point to its block is a graph homomorphism
``ψ : G_\\text{fine} → G_\\text{coarse}`` of the grid graphs
([`aggregation_homomorphism`](@ref)): an edge inside a block is collapsed, an
edge between blocks goes to the edge between their points.

**Transfers as pushforward and pullback.** Grid functions are 0-cochains of
the constant sheaf ``\\underline{ℝ}``. Its pushforward ``ψ_* \\underline{ℝ}``
([`pushforward_sheaf`](@ref CellularSheaves.NetworkSheaves.Pushforwards.pushforward_sheaf))
has as stalk over a coarse point the global sections of ``\\underline{ℝ}`` over
the fibre, the constants on the block: it is the constant sheaf of the coarse
grid. The two maps between the cochain spaces are

- the *pullback* (prolongation) ``E = ψ^*``: the inclusion of fibre sections,
  copying each coarse value onto its block ([`prolong_add!`](@ref));
- the *pushforward transfer* (restriction) ``T``
  ([`pushforward_transfer_map`](@ref CellularSheaves.NetworkSheaves.Pushforwards.pushforward_transfer_map)),
  the fibre-wise pseudoinverse of ``E``: the average over each block
  ([`restrict_average!`](@ref)). ``T E = I``.

The **Galerkin coarse operator** of an operator ``A`` is the composite
``A_c = T A E``, pushforward ∘ operator ∘ pullback ([`galerkin_coefficients!`](@ref)).
Of an axis stencil it is again an axis stencil on the coarse grid, and of an
M-matrix with positive row sums again such an M-matrix (Nicolaides 1987's
aggregation; Braess 1995).

**Rediscretization.** For the first-order upwind operators of the HJB solvers,
``A_c = T A E`` is again a first-order upwind discretization on the coarse
grid. Its coupling across a block face in dimension ``j`` is the upwind rate
``f_j^\\pm / h`` summed over the children on that (outflow) face and averaged
over all the children of the block: ``\\bar f_j^\\pm / (2h)`` for a block of two points
per dimension,
with ``\\bar f_j`` the drift averaged over the outflow face; its diagonal is
``ρ`` plus the outflow rates. Rediscretizing the scheme on the coarse grid
instead takes the drift at the block centre. The two coincide exactly when
each drift component ``f_j`` is constant along its own axis ``x_j`` within the
block and keeps its sign there (constant drift; the position rows ``\\dot q = v``
of a mechanical system written in velocities). Otherwise they differ at first
order in ``h``: by the variation of ``f_j`` along ``x_j`` across half a block, and
where ``f_j`` changes sign inside the block, by averaging the upwind rates
rather than upwinding the averaged drift. Both are consistent coarse
discretizations of the same transport operator.
"""
module GridMultigrid

export coarsening_factors, coarse_points, coarse_layout, aggregation_homomorphism, restrict_average!,
    prolong_add!, stencil_coefficients!, galerkin_coefficients!

using ArgCheck: @argcheck
using KernelAbstractions
using KernelAbstractions: @kernel, @index, @Const, get_backend, synchronize
using ..GridSchwarz: GridOperator, BoxLayout, _coefficients, _shift
using ..GraphHomomorphisms: GraphHomomorphism

"""
    coarsening_factors(points, periodic) -> NTuple

The factor ``r_k \\in \\{1, 2\\}`` by which each dimension of a grid of `points`
is coarsened (see the module documentation): 2 for a non-periodic dimension
of at least five points and for a periodic dimension whose count is a
multiple of 4 and at least 8 (every coarse dimension keeps at least three
points), 1 otherwise. All ones means the grid cannot be coarsened.
"""
coarsening_factors(points::NTuple{D,Integer}, periodic::NTuple{D,Bool}) where {D} =
    ntuple(k -> (periodic[k] ? points[k] % 4 == 0 && points[k] >= 8 : points[k] >= 5) ? 2 : 1, D)

"""
    coarse_points(points, factors) -> NTuple

The number of points per dimension of the coarse grid: ``\\lceil n_k / 2 \\rceil``
where the factor is 2, ``n_k`` where it is 1.
"""
function coarse_points(points::NTuple{D,Integer}, factors::NTuple{D,Integer}) where {D}
    @argcheck all(in((1, 2)), factors) "coarsening factors must be 1 or 2"
    return ntuple(k -> factors[k] == 2 ? cld(Int(points[k]), 2) : Int(points[k]), D)
end

"""
    coarse_layout(layout::BoxLayout) -> BoxLayout

The single-box layout of the coarse grid of the single-box `layout`
(factors from [`coarsening_factors`](@ref)), with the same periodic dimensions
and overlap.
"""
function coarse_layout(layout::BoxLayout{D}) where {D}
    @argcheck all(==(1), layout.ranks) "grid hierarchies are single-box"
    factors = coarsening_factors(layout.points, layout.periodic)
    return BoxLayout(coarse_points(layout.points, factors), ntuple(_ -> 1, D), 0; overlap=layout.overlap,
        periodic=layout.periodic)
end

# The factors relating a fine grid of nf points to a coarse grid of nc points.
function _factors(nf::NTuple{D,Int}, nc::NTuple{D,Int}) where {D}
    factors = ntuple(k -> nc[k] == nf[k] ? 1 : 2, D)
    @argcheck nc == coarse_points(nf, factors) "a grid of $(nc) points is not a coarse grid of one of $(nf) points"
    return factors
end

"""
    aggregation_homomorphism(points, factors) -> GraphHomomorphism

The graph homomorphism ``ψ`` from the grid graph of `points` (vertices
numbered column-major) to the grid graph of its coarse grid, sending each
point to its block.
"""
function aggregation_homomorphism(points::NTuple{D,Integer}, factors::NTuple{D,Integer}) where {D}
    nc = coarse_points(points, factors)
    L = LinearIndices(nc)
    vmap = [L[CartesianIndex(ntuple(k -> _block(I[k], factors[k]), D))] for I in CartesianIndices(Tuple(points))]
    return GraphHomomorphism(vec(vmap), prod(nc))
end

# The block of fine index i, and the fine indices of block b.
@inline _block(i, r) = r == 2 ? (i + 1) ÷ 2 : i
@inline _children(b, r, n) = r == 2 ? ((2b - 1):min(2b, n)) : (b:b)
@inline _children(Ic, r::NTuple{D}, n::NTuple{D}, ::Val{D}) where {D} =
    CartesianIndices(ntuple(d -> _children(Ic[d], r[d], n[d]), Val(D)))

@kernel function _restrict_average_kernel!(rc, @Const(rf), gc, gf, r, nf, ::Val{D}) where {D}
    Ic = @index(Global, Cartesian)
    acc = zero(eltype(rc))
    count = 0
    for I in _children(Ic, r, nf, Val(D))
        acc += rf[I + _shift(gf, Val(D))]
        count += 1
    end
    rc[Ic + _shift(gc, Val(D))] = acc / count
end

@kernel function _prolong_add_kernel!(xf, @Const(xc), gf, gc, r, ::Val{D}) where {D}
    I = @index(Global, Cartesian)
    Ic = CartesianIndex(ntuple(d -> _block(I[d], r[d]), Val(D)))
    xf[I + _shift(gf, Val(D))] += xc[Ic + _shift(gc, Val(D))]
end

"""
    restrict_average!(rc, opc, rf, opf) -> rc

The pushforward transfer ``T``: set the interior of the coarse vector `rc`
(operator `opc`) to the averages of the fine vector `rf` (operator `opf`) over
the blocks.
"""
function restrict_average!(rc::AbstractArray, opc::GridOperator{T,D}, rf::AbstractArray, opf::GridOperator) where {T,D}
    r = _factors(size(opf), size(opc))
    backend = get_backend(opc)
    _restrict_average_kernel!(backend)(rc, rf, opc.origin, opf.origin, r, size(opf), Val(D); ndrange=size(opc))
    synchronize(backend)
    return rc
end

"""
    prolong_add!(xf, opf, xc, opc) -> xf

The pullback ``E = ψ^*``: add to every point of the interior of the fine
vector `xf` the value of the coarse vector `xc` at its block.
"""
function prolong_add!(xf::AbstractArray, opf::GridOperator{T,D}, xc::AbstractArray, opc::GridOperator) where {T,D}
    r = _factors(size(opf), size(opc))
    backend = get_backend(opf)
    _prolong_add_kernel!(backend)(xf, xc, opf.origin, opc.origin, r, Val(D); ndrange=size(opf))
    synchronize(backend)
    return xf
end

@kernel function _coefficients_kernel!(coef, stencil, ::Val{D}) where {D}
    I = @index(Global, Cartesian)
    c0, cm, cp = _coefficients(stencil, I)
    coef[I, 1] = c0
    for j in 1:D
        coef[I, 1 + j] = cm[j]
        coef[I, 1 + D + j] = cp[j]
    end
end

"""
    stencil_coefficients!(coef, op) -> coef

The coefficients of the stencil of `op` at every interior point, in the layout
of a [`CoefficientStencil`](@ref CellularSheaves.NetworkSheaves.GridSchwarz.CoefficientStencil):
`coef[I, 1]` ``= c_0``, `coef[I, 1 + j]` ``= c^-_j``, `coef[I, 1 + D + j]` ``= c^+_j``.
"""
function stencil_coefficients!(coef::AbstractArray, op::GridOperator{T,D}) where {T,D}
    @argcheck size(coef) == (size(op)..., 2D + 1)
    backend = get_backend(op)
    _coefficients_kernel!(backend)(coef, op.stencil, Val(D); ndrange=size(op))
    synchronize(backend)
    return coef
end

# T A E at coarse point Ic: the average over the block's children of (c₀ minus
# the couplings to children of the same block), and per face the couplings of
# the children on that face to the neighbouring block.
@kernel function _galerkin_kernel!(coefc, stencil, r, nf, ::Val{D}) where {D}
    Ic = @index(Global, Cartesian)
    children = _children(Ic, r, nf, Val(D))
    lo = first(children)
    hi = last(children)
    diag = 0.0
    low = ntuple(_ -> 0.0, Val(D))
    high = ntuple(_ -> 0.0, Val(D))
    count = 0
    for I in children
        count += 1
        c0, cm, cp = _coefficients(stencil, I)
        diag += c0
        for j in 1:D
            if I[j] == lo[j]
                low = Base.setindex(low, low[j] + cm[j], j)     # to the block below
            else
                diag -= cm[j]                                    # to a child of this block
            end
            if I[j] == hi[j]
                high = Base.setindex(high, high[j] + cp[j], j)
            else
                diag -= cp[j]
            end
        end
    end
    coefc[Ic, 1] = diag / count
    for j in 1:D
        coefc[Ic, 1 + j] = low[j] / count
        coefc[Ic, 1 + D + j] = high[j] / count
    end
end

"""
    galerkin_coefficients!(coefc, opf) -> coefc

The Galerkin coarse operator ``T A E`` (pushforward ∘ operator ∘ pullback; see
the module documentation) of the stencil of `opf`, as stored coefficients of
the coarse grid whose size `coefc` gives (the layout of
[`stencil_coefficients!`](@ref)). A coupling of a fine point to a neighbour
outside a non-periodic grid (a zero ghost value) becomes a coupling of its
block to the coarse ghost, also zero.
"""
function galerkin_coefficients!(coefc::AbstractArray, opf::GridOperator{T,D}) where {T,D}
    @argcheck ndims(coefc) == D + 1 && size(coefc, D + 1) == 2D + 1
    r = _factors(size(opf), ntuple(k -> size(coefc, k), D))
    backend = get_backend(opf)
    _galerkin_kernel!(backend)(coefc, opf.stencil, r, size(opf), Val(D); ndrange=size(coefc)[1:D])
    synchronize(backend)
    return coefc
end

end # module
