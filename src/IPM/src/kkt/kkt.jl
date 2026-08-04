abstract type KKTSolver{T} end

#
# compute the KKT matrix-vector product
#
#   [ u ] = [ A  -Bᵀ ] [ x ]
#   [ v ]   [ B   0  ] [ y ]
#
function mulkkt!(
        u::AbstractVector{T},
        v::AbstractVector{T},
        A::AbstractMatrix{T},
        B::BlockSparseMatrix{T},
        x::AbstractVector{T},
        y::AbstractVector{T},
    ) where {T}
    mul!(u, Symmetric(A, :L), x)
    mul!(u, B', y, -one(T), one(T))
    mul!(v, B, x)
    return
end

include("uzawa.jl")
include("bordered.jl")
