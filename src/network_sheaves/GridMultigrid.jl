"""
    GridMultigrid

Hierarchies of structured grids for multigrid solvers of the implicit stencil
operators of [`GridSchwarz`](@ref CellularSheaves.NetworkSheaves.GridSchwarz),
with the maps between levels written as a pushforward–pullback pair.

**Coarsening.** A grid of ``n_1 \\times \\dots \\times n_D`` points (every
``n_k`` even) is cut into blocks of ``2^D`` points; the coarse grid has one
point per block, at its centre (cell-centred coarsening: spacing ``2h``,
``n_k / 2`` points, the first at ``\\mathrm{lower} + h/2``). Periodic
dimensions stay periodic. The map sending each fine point to its block is a
graph homomorphism ``ψ : G_\\text{fine} → G_\\text{coarse}`` of the grid
graphs ([`aggregation_homomorphism`](@ref)): an edge inside a block is
collapsed, an edge between blocks goes to the edge between their centres.

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
``A_c = T A E`` is the same scheme discretized again on the coarse grid
wherever the drift does not change sign inside a block: a coupling ``f_j / h``
of each of the ``2^{D-1}`` points on a block face, averaged over the ``2^D``
points of the block, is ``f_j / (2h)``, and the diagonal is
``ρ + \\sum_j |f_j| / (2h)``; a drift varying linearly across the block averages
to its value at the centre. The two coarse operators differ only on blocks
that straddle a switching surface of the drift, where Galerkin averages the
upwind directions and rediscretization upwinds the averaged drift.
"""
module GridMultigrid

export coarse_points, coarse_layout, aggregation_homomorphism, restrict_average!, prolong_add!,
    stencil_coefficients!, galerkin_coefficients!

using ArgCheck: @argcheck
using KernelAbstractions
using KernelAbstractions: @kernel, @index, @Const, get_backend, synchronize
using ..GridSchwarz: GridOperator, BoxLayout, _coefficients, _shift
using ..GraphHomomorphisms: GraphHomomorphism

"""
    coarse_points(points) -> NTuple

The number of points per dimension of the cell-centred coarse grid of a grid
of `points` (every entry even): `points .÷ 2`.
"""
function coarse_points(points::NTuple{D,Integer}) where {D}
    @argcheck all(iseven, points) "cell-centred coarsening needs an even number of points in every dimension"
    return Int.(points) .÷ 2
end

"""
    coarse_layout(layout::BoxLayout) -> BoxLayout

The single-box layout of the coarse grid of the single-box `layout`, with the
same periodic dimensions and overlap.
"""
function coarse_layout(layout::BoxLayout{D}) where {D}
    @argcheck all(==(1), layout.ranks) "grid hierarchies are single-box"
    return BoxLayout(coarse_points(layout.points), ntuple(_ -> 1, D), 0; overlap=layout.overlap,
        periodic=layout.periodic)
end

"""
    aggregation_homomorphism(points) -> GraphHomomorphism

The graph homomorphism ``ψ`` from the grid graph of `points` (vertices
numbered column-major) to the grid graph of its cell-centred coarse grid,
sending each point to its block of ``2^D`` points.
"""
function aggregation_homomorphism(points::NTuple{D,Integer}) where {D}
    nc = coarse_points(points)
    L = LinearIndices(nc)
    vmap = [L[CartesianIndex(ntuple(k -> (I[k] + 1) ÷ 2, D))] for I in CartesianIndices(Tuple(points))]
    return GraphHomomorphism(vec(vmap), prod(nc))
end

@inline _child(Ic, k, ::Val{D}) where {D} = CartesianIndex(ntuple(d -> 2Ic[d] - 1 + ((k >> (d - 1)) & 1), Val(D)))

@kernel function _restrict_average_kernel!(rc, @Const(rf), gc, gf, ::Val{D}) where {D}
    Ic = @index(Global, Cartesian)
    acc = zero(eltype(rc))
    for k in 0:(2^D - 1)
        acc += rf[_child(Ic, k, Val(D)) + _shift(gf, Val(D))]
    end
    rc[Ic + _shift(gc, Val(D))] = acc / 2^D
end

@kernel function _prolong_add_kernel!(xf, @Const(xc), gf, gc, ::Val{D}) where {D}
    I = @index(Global, Cartesian)
    Ic = CartesianIndex(ntuple(d -> (I[d] + 1) ÷ 2, Val(D)))
    xf[I + _shift(gf, Val(D))] += xc[Ic + _shift(gc, Val(D))]
end

function _check_levels(opf::GridOperator{T,D}, opc::GridOperator{S,D}) where {T,S,D}
    @argcheck size(opc) == coarse_points(size(opf)) "the coarse operator must be on the coarse grid of the fine one"
end

"""
    restrict_average!(rc, opc, rf, opf) -> rc

The pushforward transfer ``T``: set the interior of the coarse vector `rc`
(operator `opc`) to the averages of the fine vector `rf` (operator `opf`) over
the blocks of ``2^D`` points.
"""
function restrict_average!(rc::AbstractArray, opc::GridOperator{T,D}, rf::AbstractArray, opf::GridOperator) where {T,D}
    _check_levels(opf, opc)
    backend = get_backend(opc)
    _restrict_average_kernel!(backend)(rc, rf, opc.origin, opf.origin, Val(D); ndrange=size(opc))
    synchronize(backend)
    return rc
end

"""
    prolong_add!(xf, opf, xc, opc) -> xf

The pullback ``E = ψ^*``: add to every point of the interior of the fine
vector `xf` the value of the coarse vector `xc` at its block.
"""
function prolong_add!(xf::AbstractArray, opf::GridOperator{T,D}, xc::AbstractArray, opc::GridOperator) where {T,D}
    _check_levels(opf, opc)
    backend = get_backend(opf)
    _prolong_add_kernel!(backend)(xf, xc, opf.origin, opc.origin, Val(D); ndrange=size(opf))
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

# T A E at coarse point Ic: average over the children of (c₀ minus the
# couplings to children of the same block), and per face the couplings of the
# children on that face to the neighbouring block.
@kernel function _galerkin_kernel!(coefc, stencil, ::Val{D}) where {D}
    Ic = @index(Global, Cartesian)
    diag = 0.0
    low = ntuple(_ -> 0.0, Val(D))
    high = ntuple(_ -> 0.0, Val(D))
    for k in 0:(2^D - 1)
        c0, cm, cp = _coefficients(stencil, _child(Ic, k, Val(D)))
        diag += c0
        for j in 1:D
            if (k >> (j - 1)) & 1 == 1          # upper child in dimension j
                diag -= cm[j]                    # its lower neighbour is in the block
                high = Base.setindex(high, high[j] + cp[j], j)
            else
                diag -= cp[j]
                low = Base.setindex(low, low[j] + cm[j], j)
            end
        end
    end
    coefc[Ic, 1] = diag / 2^D
    for j in 1:D
        coefc[Ic, 1 + j] = low[j] / 2^D
        coefc[Ic, 1 + D + j] = high[j] / 2^D
    end
end

"""
    galerkin_coefficients!(coefc, opf) -> coefc

The Galerkin coarse operator ``T A E`` (pushforward ∘ operator ∘ pullback; see
the module documentation) of the stencil of `opf`, as stored coefficients of
the coarse grid (the layout of [`stencil_coefficients!`](@ref)). A coupling of a
fine point to a neighbour outside a non-periodic grid (a zero ghost value)
becomes a coupling of its block to the coarse ghost, also zero.
"""
function galerkin_coefficients!(coefc::AbstractArray, opf::GridOperator{T,D}) where {T,D}
    nc = coarse_points(size(opf))
    @argcheck size(coefc) == (nc..., 2D + 1)
    backend = get_backend(opf)
    _galerkin_kernel!(backend)(coefc, opf.stencil, Val(D); ndrange=nc)
    synchronize(backend)
    return coefc
end

end # module
