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
# padded with `ghost` layers on every side, where a box reads its neighbours'
# values (or zero Dirichlet data). Every kernel is a KernelAbstractions kernel,
# so the same code runs on CPU threads and GPUs. Coupling only axis neighbours
# makes the grid graph bipartite: the red–black coloring is exact, and
# Gauss–Seidel within a color is fully parallel.
module GridSchwarz

export AbstractStencil, CoefficientStencil, GridOperator, grid_zeros, interior, apply!, red_black_sgs!,
    grid_dot, grid_reduce, GridWorkspace, grid_bicgstab!

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

@inline _coefficients(s::CoefficientStencil{D}, I) where {D} =
    (s.coef[I, 1], ntuple(j -> s.coef[I, 1 + j], Val(D)), ntuple(j -> s.coef[I, 1 + D + j], Val(D)))

"""
    GridOperator(stencil::AbstractStencil, n; ghost = 1, parity = 0, eltype = Float64)
    GridOperator(coef::AbstractArray; ghost = 1, parity = 0)

The implicit operator of `stencil` on a box of `n = (n₁, …, n_D)` interior
points (the second form wraps stored coefficients in a
[`CoefficientStencil`](@ref)). Vectors are arrays of size ``n_k + 2\\,``
`ghost` per dimension ([`grid_zeros`](@ref)); the operator acts on the
[`interior`](@ref) and reads the ghost layer for the neighbours of boundary
points. A coefficient pointing out of the global domain should be zero, with
the boundary data moved into the right-hand side.

`parity` is the parity of the global index of the first interior point, so
boxes of one global grid agree on the red–black coloring.
"""
struct GridOperator{T,D,S<:AbstractStencil{D}}
    stencil::S
    size::NTuple{D,Int}
    ghost::Int
    parity::Int
end

function GridOperator(stencil::AbstractStencil{D}, n::NTuple{D,Integer}; ghost::Integer=1, parity::Integer=0,
                      eltype::Type=Float64) where {D}
    @argcheck ghost >= 1 "need at least one ghost layer"
    @argcheck all(>=(1), n) "the box needs at least one point per dimension"
    return GridOperator{eltype,D,typeof(stencil)}(stencil, Int.(n), Int(ghost), Int(parity) & 1)
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

A zero vector for `op`: its interior plus `op.ghost` layers on every side, on
the backend of `op`.
"""
grid_zeros(op::GridOperator{T}) where {T} = KernelAbstractions.zeros(get_backend(op), T, (size(op) .+ 2op.ghost)...)

"""
    interior(x, op::GridOperator)

The view of the interior of the padded vector `x`.
"""
interior(x::AbstractArray, op::GridOperator) =
    view(x, ntuple(k -> (op.ghost + 1):(op.ghost + size(op, k)), ndims(op))...)

@inline _unit(j, ::Val{D}) where {D} = CartesianIndex(ntuple(k -> k == j ? 1 : 0, Val(D)))
@inline _shift(g, ::Val{D}) where {D} = CartesianIndex(ntuple(_ -> g, Val(D)))

@kernel function _apply_kernel!(y, stencil, @Const(x), g::Int, ::Val{D}) where {D}
    I = @index(Global, Cartesian)
    J = I + _shift(g, Val(D))
    c0, cm, cp = _coefficients(stencil, I)
    acc = c0 * x[J]
    for j in 1:D
        e = _unit(j, Val(D))
        acc -= cm[j] * x[J - e] + cp[j] * x[J + e]
    end
    y[J] = acc
end

"""
    apply!(y, op::GridOperator, x) -> y

Set the interior of `y` to ``A x`` by evaluating the stencil at every point
(no matrix is formed), reading the ghost layer of `x` for the neighbours of
boundary points.
"""
function apply!(y::AbstractArray, op::GridOperator{T,D}, x::AbstractArray) where {T,D}
    backend = get_backend(op)
    _apply_kernel!(backend)(y, op.stencil, x, op.ghost, Val(D); ndrange=size(op))
    synchronize(backend)
    return y
end

# One Gauss–Seidel pass over the points of one color: each point of the color
# is updated from its row with the current values of its neighbours, which all
# have the other color, so the points of a color are independent.
@kernel function _color_kernel!(y, stencil, @Const(r), g::Int, shift::Int, color::Int, ::Val{D}) where {D}
    I = @index(Global, Cartesian)
    if (sum(Tuple(I)) + shift) & 1 == color
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

function _color_pass!(y, op::GridOperator{T,D}, r, color::Int) where {T,D}
    backend = get_backend(op)
    _color_kernel!(backend)(y, op.stencil, r, op.ghost, op.parity - D, color, Val(D); ndrange=size(op))
    synchronize(backend)
    return y
end

"""
    red_black_sgs!(y, op::GridOperator, r) -> y

One symmetric Gauss–Seidel sweep for ``A y = r`` from ``y = 0`` in red–black
order, i.e. ``y = (D + U)^{-1} D (D + L)^{-1} r`` with the red points (even
global index sum) ordered first: a forward pass red then black, and a backward
pass black then red. The backward black pass would recompute the values the
forward pass just wrote (black points only see red neighbours), so the sweep is
red, black, red. The ghost layer of `y` is zero throughout: the local problem
has zero Dirichlet data, as in restricted additive Schwarz. Equal to
`SymmetricGaussSeidel` with the red–black coloring on the assembled matrix.
"""
function red_black_sgs!(y::AbstractArray, op::GridOperator, r::AbstractArray)
    fill!(y, zero(eltype(y)))
    _color_pass!(y, op, r, 0)
    _color_pass!(y, op, r, 1)
    _color_pass!(y, op, r, 0)
    return y
end

# ===== Reductions and vector updates on the interior =====

# Partial reductions of f(x[J], y[J]) with ⊕ along the first dimension; the
# host or device then reduces the partial array (deterministic for a fixed
# backend). The reduction starts from zero, so ⊕ is + or max of nonnegatives.
@kernel function _rowreduce_kernel!(partial, f, op, @Const(x), @Const(y), g::Int, n1::Int, ::Val{D}) where {D}
    K = @index(Global, Cartesian)
    acc = zero(eltype(partial))
    for i in 1:n1
        J = CartesianIndex(i + g, ntuple(k -> K[k] + g, Val(D - 1))...)
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
    _rowreduce_kernel!(backend)(partial, f, ⊕, x, y, op.ghost, size(op, 1), Val(D); ndrange=_reduction_shape(op))
    synchronize(backend)
    return reduce(⊕, partial)
end

@kernel function _lincomb_kernel!(out, a, x, b, y, c, z, g::Int, ::Val{D}) where {D}
    I = @index(Global, Cartesian)
    J = I + _shift(g, Val(D))
    out[J] = a * x[J] + b * y[J] + c * z[J]
end

# out = a x + b y + c z on the interior (out may alias any of x, y, z).
function _lincomb!(out, op::GridOperator{T,D}, a, x, b, y, c, z) where {T,D}
    backend = get_backend(op)
    _lincomb_kernel!(backend)(out, T(a), x, T(b), y, T(c), z, op.ghost, Val(D); ndrange=size(op))
    synchronize(backend)
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

end # module
