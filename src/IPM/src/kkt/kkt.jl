@enum KKTStatus KKT_CONTINUE KKT_SOLVED KKT_STAGNATED KKT_ITMAX KKT_NUMERICAL_FAILURE KKT_IRMAX

#
# compute the KKT matrix-vector product
#
#   [ u ] = [ A  -Bᵀ ] [ x ]
#   [ v ]   [ B   0  ] [ y ]
#
function mulkkt!(
        u::AbstractVector,
        v::AbstractVector,
        A::AbstractMatrix,
        B::AbstractMatrix,
        x::AbstractVector,
        y::AbstractVector,
    )
    mul!(u, A, x)
    mul!(u, B', y, -1, 1)
    mul!(v, B, x)
    return
end

function symbkkt(
        Q::AbstractMatrix{T},
        B::AbstractMatrix{T},
        f::AbstractVector{T},
        g::AbstractVector{T},
        K::AbstractVector,
        P1::AbstractMatrix,
        P2::AbstractMatrix,
        alg::EliminationAlgorithm,
    ) where {T}
    m, n = size(B)

    weights, graph = weightedgraph(B, Q)
    D, P, S = symbolic(weights, graph; alg)

    sQ = halfselectvtxs(halfselectvtxs(Q, D.perm), D.perm)
    sB = selectvtxs(B, D.perm)

    sf = FVector{T}(undef, n)
    sg = FVector{T}(undef, m)

    mul!(sf, P, f)
    copyto!(sg, g)

    sK = tounion(K, D.perm)

    sP1 =     P1
    sP2 = P * P2

    return S, sQ, sB, sf, sg, sK, sP1, sP2
end

include("cg.jl")
include("uzawa.jl")
include("skew.jl")
