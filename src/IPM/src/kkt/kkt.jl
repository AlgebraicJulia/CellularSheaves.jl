abstract type KKTSolver{T} end

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

#
# symbolic KKT setup: pick a fill-reducing elimination ordering for the
# KKT system and permute the problem data into that ordering.
#
# returns the symbolic factorization S, the permuted problem data
# (B, Q, f, g, cones), and the column/row permutations (C, R) composed
# against the caller's own.
#
function symbkkt(
        B::AbstractMatrix,
        Q::AbstractMatrix,
        f::AbstractVector,
        g::AbstractVector,
        K::AbstractVector,
        C::AbstractMatrix,
        R::AbstractMatrix,
        alg::EliminationAlgorithm,
    )
    weights, graph = weightedgraph(B, Q)
    D, P, S = symbolic(weights, graph; alg)

    B = selectvtxs(B, D.perm)
    Q = halfselectvtxs(halfselectvtxs(Q, D.perm), D.perm)
    f = P * f
    g = copy(g)
    K = tounion(K, D.perm)
    C = P * C

    return S, B, Q, f, g, K, C, R
end

include("cg.jl")
include("uzawa.jl")
include("bordered.jl")
