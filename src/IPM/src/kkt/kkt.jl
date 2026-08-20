abstract type KKTSolver{T} end

@enum KKTStatus KKT_CONTINUE KKT_SOLVED KKT_STAGNATED KKT_ITMAX KKT_NUMERICAL_FAILURE KKT_RFMAX

const INTERIOR_THETA = 0.1

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
    mul!(u, A, x)
    mul!(u, B', y, -one(T), one(T))
    mul!(v, B, x)
    return
end

include("cg.jl")
include("uzawa.jl")
include("bordered.jl")
