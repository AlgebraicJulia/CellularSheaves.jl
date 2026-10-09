# Implicit (matrix-free) operators on structured grids, for multithreaded,
# distributed and GPU Schwarz solvers.
#
# A box of a D-dimensional grid carries an operator whose stencil couples each
# point only to its neighbours along the axes (first-order upwind schemes,
# the 2D+1-point Laplacian, …):
#
#   (A x)[I] = c₀[I] x[I] − Σ_j ( c⁻_j[I] x[I − e_j] + c⁺_j[I] x[I + e_j] ).
#
# No matrix is ever formed: every application evaluates the stencil inside a
# kernel, from a stencil object that either stores the coefficients or
# computes them from the problem data on the fly. Vectors live in arrays
# padded with ghost cells around the box, where a box reads its neighbours'
# values (or zero Dirichlet data). Every kernel is a KernelAbstractions kernel,
# so the same code runs on CPU threads and GPUs. Coupling only axis neighbours
# makes the grid graph bipartite: the red–black coloring is exact, and
# Gauss–Seidel within a color is fully parallel.
module GridSchwarz

export AbstractStencil, CoefficientStencil, GridOperator, grid_zeros, interior, apply!, red_black_sgs!,
    grid_dot, grid_reduce, GridWorkspace, grid_bicgstab!, BoxLayout, box_coordinates, balanced_ranks, box_operator,
    BoxCommunicator, SerialBoxes, mpi_boxes, box_count, box_rank, box_allreduce, exchange!, box_allgather, gather_boxes

using Adapt: Adapt
using ArgCheck: @argcheck
using KernelAbstractions
using KernelAbstractions: @kernel, @index, @Const, get_backend, synchronize

"""
    AbstractStencil{D}

The coefficients of an axis stencil on a box of a ``D``-dimensional grid,
evaluated inside kernels. A subtype implements
`GridSchwarz._coefficients(s, I) -> (c₀, c⁻, c⁺)` for the interior point `I`
(a `CartesianIndex`, 1-based in the box), with `c⁻`, `c⁺` `NTuple{D}`s, so that

```math
(A x)_I = c_0\\, x_I - \\sum_{j=1}^{D} \\big( c^-_j\\, x_{I - e_j} + c^+_j\\, x_{I + e_j} \\big),
```

and `KernelAbstractions.get_backend(s)`. Stencils holding arrays implement
`Adapt.adapt_structure` so they can be passed to GPU kernels. Computing the
coefficients from the problem data (as the HJB upwind stencil does from the
policy) keeps the kernels free of stored matrices and reads less memory.
"""
abstract type AbstractStencil{D} end

"""
    CoefficientStencil(coef)

A stencil with stored coefficients: `coef` has size ``(n_1, \\dots, n_D, 2D + 1)``
with `coef[I, 1]` ``= c_0``, `coef[I, 1 + j]` ``= c^-_j`` and
`coef[I, 1 + D + j]` ``= c^+_j``.
"""
struct CoefficientStencil{D,C<:AbstractArray} <: AbstractStencil{D}
    coef::C
end

function CoefficientStencil(coef::AbstractArray)
    D = ndims(coef) - 1
    @argcheck D >= 1 && size(coef, D + 1) == 2D + 1 "coef must have size (n₁, …, n_D, 2D + 1)"
    return CoefficientStencil{D,typeof(coef)}(coef)
end

function Adapt.adapt_structure(to, s::CoefficientStencil{D}) where {D}
    coef = Adapt.adapt(to, s.coef)
    return CoefficientStencil{D,typeof(coef)}(coef)
end
KernelAbstractions.get_backend(s::CoefficientStencil) = get_backend(s.coef)

Base.@propagate_inbounds _coefficients(s::CoefficientStencil{D}, I) where {D} =
    (s.coef[I, 1], ntuple(j -> s.coef[I, 1 + j], Val(D)), ntuple(j -> s.coef[I, 1 + D + j], Val(D)))

"""
    GridOperator(stencil::AbstractStencil, n; ghost = 1, origin = ghost, padded = n .+ 2ghost,
                 parity = 0, eltype = Float64)
    GridOperator(coef::AbstractArray; ghost = 1, parity = 0)

The implicit operator of `stencil` on a box of `n = (n₁, …, n_D)` interior
points (the second form wraps stored coefficients in a
[`CoefficientStencil`](@ref)). Vectors are arrays of size `padded`
([`grid_zeros`](@ref)), by default the box plus `ghost` layers on every side;
interior point `I` sits at `I + origin` (`origin` a number or an `NTuple`),
so several operators can act on overlapping regions of the same storage
(the owned box of a subdomain and its overlap-extended box). The operator acts
on the [`interior`](@ref) and reads the surrounding cells for the neighbours of
boundary points. A coefficient pointing out of the global domain should be
zero, with the boundary data moved into the right-hand side.

`parity` is the parity of the global index of the first interior point, so
boxes of one global grid agree on the red–black coloring.
"""
struct GridOperator{T,D,S<:AbstractStencil{D}}
    stencil::S
    size::NTuple{D,Int}
    origin::NTuple{D,Int}
    padded::NTuple{D,Int}
    parity::Int
end

function GridOperator(stencil::AbstractStencil{D}, n::NTuple{D,Integer}; ghost::Integer=1,
                      origin::Union{Integer,NTuple{D,Integer}}=ghost, padded::NTuple{D,Integer}=n .+ 2ghost,
                      parity::Integer=0, eltype::Type=Float64) where {D}
    o = origin isa Integer ? ntuple(_ -> Int(origin), D) : Int.(origin)
    @argcheck all(>=(1), n) "the box needs at least one point per dimension"
    @argcheck all(>=(1), o) && all(o .+ n .< padded) "need at least one cell around the box in the padded array"
    return GridOperator{eltype,D,typeof(stencil)}(stencil, Int.(n), o, Int.(padded), Int(parity) & 1)
end

function GridOperator(coef::AbstractArray{T}; ghost::Integer=1, parity::Integer=0) where {T}
    s = CoefficientStencil(coef)
    D = ndims(coef) - 1
    return GridOperator(s, ntuple(k -> size(coef, k), D); ghost, parity, eltype=T)
end

Base.size(op::GridOperator) = op.size
Base.size(op::GridOperator, k::Integer) = op.size[k]
Base.ndims(::GridOperator{T,D}) where {T,D} = D
Base.eltype(::GridOperator{T}) where {T} = T
KernelAbstractions.get_backend(op::GridOperator) = get_backend(op.stencil)

"""
    grid_zeros(op::GridOperator) -> array

A zero vector of the padded size of `op`, on the backend of `op`.
"""
grid_zeros(op::GridOperator{T}) where {T} = KernelAbstractions.zeros(get_backend(op), T, op.padded...)

"""
    interior(x, op::GridOperator)

The view of the interior of the padded vector `x`.
"""
interior(x::AbstractArray, op::GridOperator) =
    view(x, ntuple(k -> (op.origin[k] + 1):(op.origin[k] + size(op, k)), ndims(op))...)

@inline _unit(j, ::Val{D}) where {D} = CartesianIndex(ntuple(k -> k == j ? 1 : 0, Val(D)))
@inline _shift(g::Integer, ::Val{D}) where {D} = CartesianIndex(ntuple(_ -> g, Val(D)))
@inline _shift(o::NTuple{D,Int}, ::Val{D}) where {D} = CartesianIndex(o)

# Each kernel has one body per point, launched in one of two shapes:
# one point per work item on GPUs, and one grid row (along the first, contiguous
# dimension) per work item on the CPU backend, where decoding a Cartesian index
# per point would cost more than the stencil itself.
_rowwise(backend) = backend isa KernelAbstractions.CPU
_rows(op::GridOperator{T,D}) where {T,D} = CartesianIndices(D == 1 ? (1,) : Base.tail(size(op)))
@inline _point(i, K, ::Val{D}) where {D} = CartesianIndex(i, ntuple(k -> K[k], Val(D - 1))...)

function _launch(kernel_points, kernel_rows, op::GridOperator, args...)
    backend = get_backend(op)
    if _rowwise(backend)
        rows = _rows(op)
        groups = cld(length(rows), 8 * Threads.nthreads())
        kernel_rows(backend, groups)(args..., rows, size(op, 1); ndrange=length(rows))
    else
        kernel_points(backend)(args...; ndrange=size(op))
    end
    synchronize(backend)
end

@inline function _apply_point!(y, stencil, x, g, I, ::Val{D}) where {D}
    @inbounds begin
        J = I + _shift(g, Val(D))
        c0, cm, cp = _coefficients(stencil, I)
        acc = c0 * x[J]
        for j in 1:D
            e = _unit(j, Val(D))
            acc -= cm[j] * x[J - e] + cp[j] * x[J + e]
        end
        y[J] = acc
    end
end

@kernel function _apply_kernel!(y, stencil, @Const(x), g, ::Val{D}) where {D}
    _apply_point!(y, stencil, x, g, @index(Global, Cartesian), Val(D))
end

@kernel function _apply_rows_kernel!(y, stencil, @Const(x), g, ::Val{D}, rows, n1) where {D}
    row = @index(Global, Linear)
    K = rows[row]
    for i in 1:n1
        _apply_point!(y, stencil, x, g, _point(i, K, Val(D)), Val(D))
    end
end

"""
    apply!(y, op::GridOperator, x) -> y

Set the interior of `y` to ``A x`` by evaluating the stencil at every point
(no matrix is formed), reading the ghost layer of `x` for the neighbours of
boundary points.
"""
function apply!(y::AbstractArray, op::GridOperator{T,D}, x::AbstractArray) where {T,D}
    _launch(_apply_kernel!, _apply_rows_kernel!, op, y, op.stencil, x, op.origin, Val(D))
    return y
end

# One Gauss–Seidel pass over the points of one color: each point of the color
# is updated from its row with the current values of its neighbours, which all
# have the other color, so the points of a color are independent.
@inline function _color_point!(y, stencil, r, g, I, ::Val{D}) where {D}
    @inbounds begin
        J = I + _shift(g, Val(D))
        c0, cm, cp = _coefficients(stencil, I)
        acc = r[J]
        for j in 1:D
            e = _unit(j, Val(D))
            acc += cm[j] * y[J - e] + cp[j] * y[J + e]
        end
        y[J] = acc / c0
    end
end

@kernel function _color_kernel!(y, stencil, @Const(r), g, shift::Int, color::Int, ::Val{D}) where {D}
    I = @index(Global, Cartesian)
    if (sum(Tuple(I)) + shift) & 1 == color
        _color_point!(y, stencil, r, g, I, Val(D))
    end
end

# Along a row the colors alternate: visit only the points of `color`.
@kernel function _color_rows_kernel!(y, stencil, @Const(r), g, shift::Int, color::Int, ::Val{D}, rows, n1) where {D}
    row = @index(Global, Linear)
    K = rows[row]
    first_i = 1 + ((color - (1 + sum(Tuple(K)) + shift)) & 1)
    for i in first_i:2:n1
        _color_point!(y, stencil, r, g, _point(i, K, Val(D)), Val(D))
    end
end

function _color_pass!(y, op::GridOperator{T,D}, r, color::Int) where {T,D}
    _launch(_color_kernel!, _color_rows_kernel!, op, y, op.stencil, r, op.origin, op.parity - D, color, Val(D))
    return y
end

"""
    red_black_sgs!(y, op::GridOperator, r; exchange = identity) -> y

One symmetric Gauss–Seidel sweep for ``A y = r`` from ``y = 0`` in red–black
order, i.e. ``y = (D + U)^{-1} D (D + L)^{-1} r`` with the red points (even
global index sum) ordered first: a forward pass red then black, and a backward
pass black then red. The backward black pass would recompute the values the
forward pass just wrote (black points only see red neighbours), so the sweep is
red, black, red.
Equal to
`SymmetricGaussSeidel` with the red–black coloring on the assembled matrix.

When `op` is one box of a distributed grid, `exchange(y)` refreshes the ghost
layer from the neighbouring boxes after each color: the sweep is then the
global red–black sweep, identical to a single-box sweep of the whole grid, at
the price of two halo exchanges (the first and last passes need none).
With the default `identity` the ghost layer stays zero: the local problem has
zero Dirichlet data, as in restricted additive Schwarz.
"""
function red_black_sgs!(y::AbstractArray, op::GridOperator, r::AbstractArray; exchange=identity)
    fill!(y, zero(eltype(y)))
    _color_pass!(y, op, r, 0)
    exchange(y)
    _color_pass!(y, op, r, 1)
    exchange(y)
    _color_pass!(y, op, r, 0)
    return y
end

# ===== Reductions and vector updates on the interior =====

# Partial reductions of f(x[J], y[J]) with ⊕ along the first dimension; the
# host or device then reduces the partial array (deterministic for a fixed
# backend). The reduction starts from zero, so ⊕ is + or max of nonnegatives.
@kernel function _rowreduce_kernel!(partial, f, op, @Const(x), @Const(y), o, n1::Int, ::Val{D}) where {D}
    K = @index(Global, Cartesian)
    acc = zero(eltype(partial))
    for i in 1:n1
        J = CartesianIndex(i + o[1], ntuple(k -> K[k] + o[k + 1], Val(D - 1))...)
        acc = op(acc, f(x[J], y[J]))
    end
    partial[K] = acc
end

_reduction_shape(op::GridOperator{T,D}) where {T,D} = D == 1 ? (1,) : Base.tail(size(op))

"""
    grid_dot(x, y, op, partial) -> number

``\\sum_I x_I y_I`` over the interior of `op`, with `partial` an array of
size `size(op)[2:end]` on the backend of `op` for the partial sums.
"""
grid_dot(x, y, op::GridOperator, partial) = grid_reduce(*, +, x, y, op, partial)

"""
    grid_reduce(f, ⊕, x, y, op, partial) -> number

``\\bigoplus_I f(x_I, y_I)`` over the interior of `op`, starting from zero
(so `⊕` is `+`, or `max` of nonnegative values), with `partial` as in
[`grid_dot`](@ref).
"""
function grid_reduce(f, ⊕, x, y, op::GridOperator{T,D}, partial) where {T,D}
    backend = get_backend(op)
    _rowreduce_kernel!(backend)(partial, f, ⊕, x, y, op.origin, size(op, 1), Val(D); ndrange=_reduction_shape(op))
    synchronize(backend)
    return reduce(⊕, partial)
end

@kernel function _lincomb_kernel!(out, a, x, b, y, c, z, g, ::Val{D}) where {D}
    J = @index(Global, Cartesian) + _shift(g, Val(D))
    out[J] = a * x[J] + b * y[J] + c * z[J]
end

@kernel function _lincomb_rows_kernel!(out, a, x, b, y, c, z, g, ::Val{D}, rows, n1) where {D}
    row = @index(Global, Linear)
    K = rows[row]
    o = _shift(g, Val(D))
    @inbounds @simd for i in 1:n1
        J = _point(i, K, Val(D)) + o
        out[J] = a * x[J] + b * y[J] + c * z[J]
    end
end

# out = a x + b y + c z on the interior (out may alias any of x, y, z).
function _lincomb!(out, op::GridOperator{T,D}, a, x, b, y, c, z) where {T,D}
    _launch(_lincomb_kernel!, _lincomb_rows_kernel!, op, out, T(a), x, T(b), y, T(c), z, op.origin, Val(D))
    return out
end

"""
    GridWorkspace(op::GridOperator)

The vectors [`grid_bicgstab!`](@ref) needs, allocated once on the backend of `op`.
"""
struct GridWorkspace{A,P}
    r::A
    rhat::A
    p::A
    v::A
    s::A
    t::A
    phat::A
    shat::A
    partial::P
end

function GridWorkspace(op::GridOperator{T}) where {T}
    vectors = ntuple(_ -> grid_zeros(op), 8)
    partial = KernelAbstractions.zeros(get_backend(op), T, _reduction_shape(op)...)
    return GridWorkspace(vectors..., partial)
end

"""
    grid_bicgstab!(x, op, M!, b, ws; tol, maxiter, reduce = identity,
                   exchange = identity, operator = (y, v) -> apply!(y, op, v))
        -> (iterations, converged, residual norms)

Right-preconditioned BiCGStab (van der Vorst 1992) for ``A x = b`` on the
interior of `op`, from the initial `x`, overwritten with the solution.
`M!(y, v)` sets ``y = M^{-1} v``. Ghost layers are the caller's:
`exchange(v)` is called before every application of `operator` to fill them
(a no-op on a single box), and `reduce` combines the inner products of all
boxes (`MPI.Allreduce` in a distributed solve). Every vector operation is a
kernel on the backend of `op`. Stops when ``\\|b - A x\\| \\le \\mathrm{tol}\\,\\|b\\|``;
a breakdown restarts the shadow residual.
"""
function grid_bicgstab!(x, op::GridOperator{T}, M!, b, ws::GridWorkspace; tol::Real, maxiter::Integer,
                        reduce=identity, exchange=identity, operator=((y, v) -> apply!(y, op, v))) where {T}
    dot(u, w) = reduce(grid_dot(u, w, op, ws.partial))
    (; r, rhat, p, v, s, t, phat, shat) = ws
    exchange(x)
    operator(r, x)
    _lincomb!(r, op, 1, b, -1, r, 0, r)
    copyto!(rhat, r)
    fill!(p, zero(T))
    fill!(v, zero(T))
    target = tol * max(sqrt(dot(b, b)), eps(T))
    residuals = T[sqrt(dot(r, r))]
    last(residuals) <= target && return 0, true, residuals
    ρ = α = ω = one(T)
    for it in 1:maxiter
        ρnew = dot(rhat, r)
        if iszero(ρnew)
            copyto!(rhat, r); fill!(p, zero(T)); fill!(v, zero(T))
            ρ = α = ω = one(T)
            ρnew = dot(r, r)
        end
        β = (ρnew / ρ) * (α / ω)
        _lincomb!(p, op, 1, r, β, p, -β * ω, v)
        M!(phat, p)
        exchange(phat)
        operator(v, phat)
        σ = dot(rhat, v)
        if iszero(σ)
            copyto!(rhat, r); fill!(p, zero(T)); fill!(v, zero(T))
            ρ = α = ω = one(T)
            push!(residuals, last(residuals))
            continue
        end
        α = ρnew / σ
        _lincomb!(s, op, 1, r, -α, v, 0, v)
        snorm = sqrt(dot(s, s))
        if snorm <= target
            _lincomb!(x, op, 1, x, α, phat, 0, phat)
            push!(residuals, snorm)
            return it, true, residuals
        end
        M!(shat, s)
        exchange(shat)
        operator(t, shat)
        ω = dot(t, s) / dot(t, t)
        _lincomb!(x, op, 1, x, α, phat, ω, shat)
        _lincomb!(r, op, 1, s, -ω, t, 0, t)
        ρ = ρnew
        push!(residuals, sqrt(dot(r, r)))
        last(residuals) <= target && return it, true, residuals
    end
    return maxiter, false, residuals
end


# ===== Boxes of a distributed grid =====

"""
    BoxLayout(points, ranks, coords; overlap = 1)

Box `coords` (0-based, one per dimension) of the partition of a grid of
`points` into `ranks[1] × … × ranks[D]` boxes of nearly equal index ranges,
one per rank: the vertices of the overlap sheaf of the cover, the boxes that
share a face being its edges. Fields:

- `owned`: the global index ranges the box owns;
- `extended`: `owned` grown by `overlap` points into each neighbour (clipped
  at the grid boundary), the subdomain of a restricted additive Schwarz solve;
- `width = overlap + 1`: the ghost cells stored around `owned`, enough for the
  extended box and its stencil;
- `neighbors[d] = (lower, upper)`: the ranks of the boxes sharing a face
  across dimension `d`, or `-1`.

Local vectors are padded arrays of size `length.(owned) .+ 2width`
([`box_operator`](@ref)). Ranks are numbered column-major over `coords`.
"""
struct BoxLayout{D}
    points::NTuple{D,Int}
    ranks::NTuple{D,Int}
    coords::NTuple{D,Int}
    owned::NTuple{D,UnitRange{Int}}
    extended::NTuple{D,UnitRange{Int}}
    overlap::Int
    width::Int
    neighbors::NTuple{D,NTuple{2,Int}}
end

function BoxLayout(points::NTuple{D,Integer}, ranks::NTuple{D,Integer}, coords::NTuple{D,Integer};
                   overlap::Integer=1) where {D}
    @argcheck overlap >= 0 "the overlap must be nonnegative"
    @argcheck all(1 .<= ranks .<= points) "need between 1 and points[d] boxes in dimension d"
    @argcheck all(0 .<= coords .< ranks) "box coordinates must lie in 0:ranks[d]-1"
    owned = ntuple(d -> _split(points[d], ranks[d], coords[d]), D)
    w = overlap + 1
    for d in 1:D, c in 0:(ranks[d] - 1)
        @argcheck length(_split(points[d], ranks[d], c)) >= w "boxes of fewer than overlap + 1 = $w points in dimension $d"
    end
    extended = ntuple(d -> max(1, first(owned[d]) - overlap):min(points[d], last(owned[d]) + overlap), D)
    rank(c) = sum(c[d] * prod(ranks[1:(d - 1)]; init=1) for d in 1:D)
    neighbors = ntuple(D) do d
        lower = coords[d] > 0 ? rank(ntuple(k -> k == d ? coords[k] - 1 : coords[k], D)) : -1
        upper = coords[d] < ranks[d] - 1 ? rank(ntuple(k -> k == d ? coords[k] + 1 : coords[k], D)) : -1
        (lower, upper)
    end
    return BoxLayout{D}(Int.(points), Int.(ranks), Int.(coords), owned, extended, Int(overlap), w, neighbors)
end

BoxLayout(points::NTuple{D,Integer}, ranks::NTuple{D,Integer}, rank::Integer; overlap::Integer=1) where {D} =
    BoxLayout(points, ranks, box_coordinates(ranks, rank); overlap)

_split(n, k, c) = (c * n ÷ k + 1):((c + 1) * n ÷ k)

"""
    box_coordinates(ranks, rank) -> NTuple

The 0-based coordinates of box `rank` (0-based, column-major) in a
`ranks[1] × … × ranks[D]` partition.
"""
box_coordinates(ranks::NTuple{D,Integer}, rank::Integer) where {D} =
    Tuple(CartesianIndices(ranks)[rank + 1]) .- 1

"""
    balanced_ranks(nboxes, points) -> NTuple

A factorization of `nboxes` into boxes per dimension that keeps the boxes of a
grid of `points` close to cubes (fewest ghost cells per owned cell): each prime
factor, largest first, splits the dimension whose boxes are currently longest.
"""
function balanced_ranks(nboxes::Integer, points::NTuple{D,Integer}) where {D}
    @argcheck nboxes >= 1
    factors = Int[]
    m, p = Int(nboxes), 2
    while m > 1
        while m % p == 0
            push!(factors, p)
            m ÷= p
        end
        p += 1
    end
    ranks = ones(Int, D)
    for f in sort!(factors; rev=true)
        d = argmax(ntuple(k -> points[k] / ranks[k], D))
        ranks[d] *= f
    end
    return Tuple(ranks)
end

"""
    box_operator(layout, stencil, region = :owned; eltype = Float64) -> GridOperator

The [`GridOperator`](@ref) of `stencil` on the `owned` or `extended` box of
`layout`, both acting on the same padded local arrays, with the red–black
parity of the global grid.
"""
function box_operator(layout::BoxLayout{D}, stencil::AbstractStencil{D}, region::Symbol=:owned;
                      eltype::Type=Float64) where {D}
    @argcheck region in (:owned, :extended)
    box = region === :owned ? layout.owned : layout.extended
    origin = ntuple(d -> layout.width - (first(layout.owned[d]) - first(box[d])), D)
    padded = length.(layout.owned) .+ 2layout.width
    parity = sum(first.(box)) - D
    return GridOperator(stencil, length.(box); origin, padded, parity, eltype)
end

"""
    BoxCommunicator

How the boxes of a distributed grid talk to each other. Implementations:
[`SerialBoxes`](@ref) (one box, nothing to exchange) and, with MPI.jl loaded,
[`mpi_boxes`](@ref). Each implements `box_count(c)`, `box_rank(c)`,
`box_allreduce(c, x, op)` (`op` is `+` or `max`), `exchange!(c, layout, x)`
(fill the `layout.width` ghost cells of `x` from the neighbouring boxes) and
`box_allgather(c, v)` (the vectors `v` of all ranks, in rank order).
"""
abstract type BoxCommunicator end

"""
    SerialBoxes()

The communicator of a grid held as a single box by one process.
"""
struct SerialBoxes <: BoxCommunicator end

box_count(::SerialBoxes) = 1
box_rank(::SerialBoxes) = 0
box_allreduce(::SerialBoxes, x, op) = x
exchange!(::SerialBoxes, layout::BoxLayout, x) = x
box_allgather(::SerialBoxes, v::AbstractVector) = [Vector(v)]

"""
    mpi_boxes(comm) -> BoxCommunicator

The [`BoxCommunicator`](@ref) of the MPI communicator `comm`: one box per
rank, halo exchanges with `MPI.Sendrecv!` between face neighbours, one
dimension at a time (so edge and corner ghost cells arrive in ``D`` rounds),
inner products with `MPI.Allreduce`. Needs `using MPI` (a package extension).
"""
function mpi_boxes end

# The slabs of the padded array exchanged across dimension d: the first and
# last `width` owned layers are sent, the ghost layers below and above them are
# received; the other dimensions span the whole padded range, so ghost cells
# filled across earlier dimensions are passed on (edges and corners).
function _slabs(layout::BoxLayout{D}, d::Int) where {D}
    w = layout.width
    n = length(layout.owned[d])
    padded = length.(layout.owned) .+ 2w
    slab(r) = ntuple(k -> k == d ? r : (1:padded[k]), D)
    return (send_lower=slab((w + 1):(2w)), recv_lower=slab(1:w),
            send_upper=slab((n + 1):(n + w)), recv_upper=slab((n + w + 1):(n + 2w)))
end

"""
    gather_boxes(c::BoxCommunicator, layout, x, op) -> Array
    gather_boxes(c::BoxCommunicator, layout, block) -> Array

The global array assembled, on every rank, from the owned interiors of `x` (a
padded local array, `op` its owned operator), or from the owned `block`s (of
size `length.(layout.owned)`), of every box.
"""
gather_boxes(c::BoxCommunicator, layout::BoxLayout, x::AbstractArray, op::GridOperator) =
    gather_boxes(c, layout, interior(x, op))

function gather_boxes(c::BoxCommunicator, layout::BoxLayout{D}, block::AbstractArray) where {D}
    @argcheck size(block) == length.(layout.owned) "the block must have the size of the owned box"
    blocks = box_allgather(c, vec(Array(block)))
    out = zeros(eltype(block), layout.points...)
    for (r, block) in enumerate(blocks)
        other = BoxLayout(layout.points, layout.ranks, r - 1; overlap=layout.overlap)
        out[other.owned...] .= reshape(block, length.(other.owned))
    end
    return out
end

end # module
