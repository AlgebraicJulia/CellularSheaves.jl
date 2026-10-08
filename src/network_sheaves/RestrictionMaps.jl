# Concrete restriction maps: linear maps from a vertex stalk to an edge stalk
# that need not be stored as dense matrices. Each map implements `size`,
# `mul!(y, R, x)` and `mul!(y, R', x)`; everything else is derived.
module RestrictionMaps

export DenseRestriction, SparseRestriction, SelectionRestriction, FunctionRestriction,
    restriction_map

using ArgCheck: @argcheck
using LinearAlgebra
using LinearAlgebra: Adjoint
using SparseArrays
import LinearOperators: LinearOperator

using ..SheafInterface: AbstractRestrictionMap

# ===== Concrete maps =====

"""
    DenseRestriction(matrix)

A restriction map stored as a dense matrix.
"""
struct DenseRestriction{T} <: AbstractRestrictionMap{T}
    matrix::Matrix{T}
end

"""
    SparseRestriction(matrix)

A restriction map stored as a sparse matrix.
"""
struct SparseRestriction{T} <: AbstractRestrictionMap{T}
    matrix::SparseMatrixCSC{T,Int}
end

"""
    SelectionRestriction{T}(rows, ncols)

The coordinate selection ``(R x)_r = x_{\\mathrm{rows}[r]}`` from
``\\mathbb R^{\\mathrm{ncols}}`` onto the coordinates listed in `rows`: the
restriction map of a cover onto an overlap. It stores only the index list.
"""
struct SelectionRestriction{T} <: AbstractRestrictionMap{T}
    rows::Vector{Int}
    ncols::Int

    function SelectionRestriction{T}(rows::AbstractVector{<:Integer}, ncols::Integer) where {T}
        @argcheck all(r -> 1 <= r <= ncols, rows) "selected coordinates must lie in 1:$ncols"
        return new{T}(Vector{Int}(rows), Int(ncols))
    end
end

"""
    FunctionRestriction{T}(apply!, apply_adjoint!, nrows, ncols)

A matrix-free restriction map known only through its action:
`apply!(y, x)` sets ``y = R x`` and `apply_adjoint!(y, x)` sets
``y = R^\\mathsf{T} x``. Use it for maps that are cheap to apply but expensive or
impossible to store, for example a trace operator, an interpolation, or the
output of a simulation.
"""
struct FunctionRestriction{T,F,G} <: AbstractRestrictionMap{T}
    apply!::F
    apply_adjoint!::G
    nrows::Int
    ncols::Int
end

FunctionRestriction{T}(apply!::F, apply_adjoint!::G, nrows::Integer, ncols::Integer) where {T,F,G} =
    FunctionRestriction{T,F,G}(apply!, apply_adjoint!, Int(nrows), Int(ncols))

"""
    restriction_map(A) -> AbstractRestrictionMap
    restriction_map(T, A) -> AbstractRestrictionMap{T}

Wrap a matrix as a restriction map: sparse matrices become
[`SparseRestriction`](@ref)s, other matrices [`DenseRestriction`](@ref)s. A
restriction map is returned unchanged (converted to element type `T` when
given).
"""
restriction_map(A::SparseMatrixCSC{T}) where {T} = SparseRestriction{T}(A)
restriction_map(A::AbstractMatrix{T}) where {T} = DenseRestriction{T}(Matrix(A))
restriction_map(R::AbstractRestrictionMap) = R
restriction_map(::Type{T}, A::AbstractMatrix) where {T} = restriction_map(T.(A))
restriction_map(::Type{T}, R::AbstractRestrictionMap{T}) where {T} = R
restriction_map(::Type{T}, R::SelectionRestriction) where {T} = SelectionRestriction{T}(R.rows, R.ncols)
restriction_map(::Type{T}, R::AbstractRestrictionMap) where {T} = SparseRestriction{T}(SparseMatrixCSC{T,Int}(sparse(R)))

# ===== Interface =====

Base.size(R::DenseRestriction) = size(R.matrix)
Base.size(R::SparseRestriction) = size(R.matrix)
Base.size(R::SelectionRestriction) = (length(R.rows), R.ncols)
Base.size(R::FunctionRestriction) = (R.nrows, R.ncols)
Base.size(R::AbstractRestrictionMap, d::Integer) = d <= 2 ? size(R)[d] : 1
Base.eltype(::AbstractRestrictionMap{T}) where {T} = T
Base.eltype(::Type{<:AbstractRestrictionMap{T}}) where {T} = T

LinearAlgebra.mul!(y::AbstractVector, R::Union{DenseRestriction,SparseRestriction}, x::AbstractVector) =
    mul!(y, R.matrix, x)
LinearAlgebra.mul!(y::AbstractVector, R::Adjoint{<:Any,<:Union{DenseRestriction,SparseRestriction}}, x::AbstractVector) =
    mul!(y, parent(R).matrix', x)

function LinearAlgebra.mul!(y::AbstractVector, R::SelectionRestriction, x::AbstractVector)
    @argcheck length(y) == length(R.rows) && length(x) == R.ncols
    @inbounds for (r, c) in enumerate(R.rows)
        y[r] = x[c]
    end
    return y
end

function LinearAlgebra.mul!(y::AbstractVector, Rt::Adjoint{<:Any,<:SelectionRestriction}, x::AbstractVector)
    R = parent(Rt)
    @argcheck length(y) == R.ncols && length(x) == length(R.rows)
    fill!(y, zero(eltype(y)))
    @inbounds for (r, c) in enumerate(R.rows)
        y[c] += x[r]
    end
    return y
end

LinearAlgebra.mul!(y::AbstractVector, R::FunctionRestriction, x::AbstractVector) = (R.apply!(y, x); y)
LinearAlgebra.mul!(y::AbstractVector, R::Adjoint{<:Any,<:FunctionRestriction}, x::AbstractVector) =
    (parent(R).apply_adjoint!(y, x); y)

# ===== Derived operations =====

Base.adjoint(R::AbstractRestrictionMap) = Adjoint(R)
Base.size(R::Adjoint{<:Any,<:AbstractRestrictionMap}) = reverse(size(parent(R)))
Base.size(R::Adjoint{<:Any,<:AbstractRestrictionMap}, d::Integer) = d <= 2 ? size(R)[d] : 1

const _MaybeAdjointRestriction = Union{AbstractRestrictionMap,Adjoint{<:Any,<:AbstractRestrictionMap}}

_eltype(R::AbstractRestrictionMap) = eltype(R)
_eltype(R::Adjoint{<:Any,<:AbstractRestrictionMap}) = eltype(parent(R))

function Base.:*(R::_MaybeAdjointRestriction, x::AbstractVector)
    return mul!(zeros(promote_type(_eltype(R), eltype(x)), size(R, 1)), R, x)
end

function Base.:*(R::_MaybeAdjointRestriction, X::AbstractMatrix)
    @argcheck size(X, 1) == size(R, 2) "dimension mismatch"
    Y = zeros(promote_type(_eltype(R), eltype(X)), size(R, 1), size(X, 2))
    for j in axes(X, 2)
        mul!(view(Y, :, j), R, X[:, j])
    end
    return Y
end

Base.:*(X::AbstractMatrix, R::AbstractRestrictionMap) = Matrix((R' * Matrix(X'))')
Base.:*(X::AbstractMatrix, R::Adjoint{<:Any,<:AbstractRestrictionMap}) = Matrix((parent(R) * Matrix(X'))')

Base.Matrix(R::AbstractRestrictionMap{T}) where {T} = R * Matrix{T}(I, size(R, 2), size(R, 2))
Base.Matrix(R::DenseRestriction) = copy(R.matrix)
Base.Matrix(R::SparseRestriction) = Matrix(R.matrix)
Base.Matrix{T}(R::AbstractRestrictionMap) where {T} = Matrix{T}(Matrix(R))

SparseArrays.sparse(R::AbstractRestrictionMap) = sparse(Matrix(R))
SparseArrays.sparse(R::SparseRestriction) = copy(R.matrix)
SparseArrays.sparse(R::SelectionRestriction{T}) where {T} =
    sparse(eachindex(R.rows), R.rows, ones(T, length(R.rows)), length(R.rows), R.ncols)

"""
    LinearOperator(R::AbstractRestrictionMap)

View a restriction map as a `LinearOperators.LinearOperator`, for use with
Krylov.jl and other matrix-free solvers.
"""
function LinearOperator(R::AbstractRestrictionMap{T}) where {T}
    prod!(y, x, α, β) = _axpby_apply!(y, R, x, α, β)
    tprod!(y, x, α, β) = _axpby_apply!(y, R', x, α, β)
    m, n = size(R)
    return LinearOperator(T, m, n, false, false, prod!, tprod!, tprod!)
end

function _axpby_apply!(y, R, x, α, β)
    if iszero(β)
        mul!(y, R, x)
        isone(α) || (y .*= α)
    else
        y .= α .* (R * x) .+ β .* y
    end
    return y
end

Base.:(==)(A::DenseRestriction, B::DenseRestriction) = A.matrix == B.matrix
Base.:(==)(A::SparseRestriction, B::SparseRestriction) = A.matrix == B.matrix
Base.:(==)(A::SelectionRestriction, B::SelectionRestriction) = A.rows == B.rows && A.ncols == B.ncols
Base.hash(R::DenseRestriction, h::UInt) = hash(R.matrix, hash(:DenseRestriction, h))
Base.hash(R::SparseRestriction, h::UInt) = hash(R.matrix, hash(:SparseRestriction, h))
Base.hash(R::SelectionRestriction, h::UInt) = hash(R.rows, hash(R.ncols, hash(:SelectionRestriction, h)))

function Base.show(io::IO, R::AbstractRestrictionMap)
    m, n = size(R)
    print(io, nameof(typeof(R)), "{", eltype(R), "}(", m, "×", n, ")")
end

end
