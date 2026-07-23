@kwdef struct UzawaSettings{T, E <: EliminationAlgorithm} <: KKTSettings{T}
    aaug::T = zero(T)
    raug::T = 1e6
    itmax::Int = 1000
    elim::E = DEFAULT_ELIMINATION_ALGORITHM
end

function UzawaSettings{T}(; elim::E = DEFAULT_ELIMINATION_ALGORITHM, kwargs...) where {T, E <: EliminationAlgorithm}
    return UzawaSettings{T, E}(; elim=elim, kwargs...)
end

function showsettings(io::IO, set::UzawaSettings; indent::Integer=0)
    pad = " "^indent
    @printf(io, "%saaug:  %8.2e  raug:  %8.2e\n", pad, set.aaug, set.raug)
    @printf(io, "%sitmax: %8d\n", pad, set.itmax)
    return
end

function Base.show(io::IO, ::MIME"text/plain", set::T) where {T <: UzawaSettings}
    println(io, T, ":")
    return showsettings(io, set; indent=2)
end

struct UzawaWorkspace{UPLO, T, I <: Integer} <: KKTWorkspace{T}
    F::FChordalTriangular{:N, UPLO, T, I}
    L::BlockSparseMatrix{T, I}
    facwrk::FactorizationWorkspace{T, I}
    divwrk::DivisionWorkspace{T, I}
    itrwrk::CraigWorkspace{T, T, Vector{T}}
    r::Vector{T}
    α::Scalar{T}
end

function UzawaWorkspace(F::FChordalTriangular{:N, UPLO, T, I}, L::BlockSparseMatrix{T, I}, B::BlockSparseMatrix{T, I}) where {UPLO, T, I <: Integer}
    m, n = size(B)
    facwrk = FactorizationWorkspace(F)
    divwrk = DivisionWorkspace(F, 1)
    itrwrk = CraigWorkspace(m, n, Vector{T})
    r = zeros(T, m)
    α = ones(T)
    return UzawaWorkspace(F, L, facwrk, divwrk, itrwrk, r, α)
end

function make_kkt(settings::UzawaSettings{T}, B::BlockSparseMatrix{T, I}) where {T, I}
    weights, graph = weightedgraph(B)

    R, P, S = symbolic(weights, graph; alg=settings.elim)

    B = selectvtxs(B, R.perm)

    F = FChordalTriangular{:N, :L, T, I}(S)
    L = B' * B

    wrk = UzawaWorkspace(F, L, B)

    return R, P, B, wrk
end

function initkkt!(wrk::UzawaWorkspace{UPLO, T}, set::UzawaSettings{T}, A::BlockSparseMatrix, nH::T, nB::T, rgmin::T, rgmax::T) where {UPLO, T}
    wrk.α[] = α = set.aaug + set.raug * nH / nB^2
    return init_uzw!(wrk.facwrk, wrk.F, wrk.L, A, α, rgmin, rgmax)
end

function init_uzw!(
        facwrk::FactorizationWorkspace{T},
        F::ChordalTriangular{:N, UPLO, T},
        L::BlockSparseMatrix{T},
        A::BlockSparseMatrix{T},
        α::T,
        rgmin::T,
        rgmax::T
    ) where {UPLO, T}
    @assert size(F, 1) == size(L, 1) == size(A, 1)

    β = inv(α)
    ρ = rgmin
    #
    # assemble and factor the augmented matrix
    #
    #   F = β A + Bᵀ B
    #
    copyto!(F, L)
    axpy!(β, A, F)
    info = cholesky!(facwrk, F; check=false)
    #
    # on failure, add a diagonal shift ρ I and retry:
    #
    #   F ← β A + Bᵀ B + ρ I
    #
    if !iszero(info)
        copyto!(F, L)
        axpy!(β, A, F)
        axpy!(ρ, I, F)
        info = cholesky!(facwrk, F; check=false)
    end

    while !iszero(info) && 8ρ ≤ rgmax
        ρ *= 8
        copyto!(F, L)
        axpy!(β, A, F)
        axpy!(ρ, I, F)
        info = cholesky!(facwrk, F; check=false)
    end

    return iszero(info), ρ
end

function solve_kkt!(
    wrk::UzawaWorkspace{UPLO, T},
    set::UzawaSettings{T},
    x::AbstractVector{T},
    y::AbstractVector{T},
    A::BlockSparseMatrix{T},
    B::BlockSparseMatrix{T},
    f::AbstractVector{T},
    g::AbstractVector{T},
    y0 = nothing;
    atol::T,
) where {UPLO, T}
    return solve_uzw!(wrk.divwrk, wrk.itrwrk, x, y, wrk.r, wrk.F, A, B,
                      f, g, wrk.α[], atol, set.itmax, y0)
end

#
# Solve the KKT system
#
#   [ A -Bᵀ ] [ x ] = [ f ]
#   [ B  0  ] [ y ]   [ g ]
#
function solve_uzw!(
        divwrk::DivisionWorkspace{T},
        itrwrk::CraigWorkspace{T},
        x::AbstractVector{T},
        y::AbstractVector{T},
        r::AbstractVector{T},
        F::ChordalTriangular{:N, UPLO, T},
        A::BlockSparseMatrix{T},
        B::BlockSparseMatrix{T},
        f::AbstractVector{T},
        g::AbstractVector{T},
        α::T,
        atol::T,
        itmax::Int,
        y0 = nothing,
    ) where {UPLO, T}
    m, n = size(B)

    @assert length(x) == n
    @assert length(y) == m
    @assert length(r) == m
    @assert length(f) == n
    @assert length(g) == m
    @assert size(F, 1) == n
    @assert size(A, 1) == n
    niter = 1; β = inv(α)
    #
    # solve for x:
    #
    #   (β A + Bᵀ B) x =  β f + Bᵀ (g + β y₀)
    #
    copyto!(r, g)

    if !isnothing(y0)
        axpy!(β, y0, r)
    end

    copyto!(x, f)
    mul!(x, B', r, one(T), β)
    ldiv!(divwrk, F,  x)
    ldiv!(divwrk, F', x)
    #
    # compute the residual
    #
    #   r = g - B x
    #
    copyto!(r, g)
    mul!(r, B, x, -one(T), one(T))
    #
    # solve for δx and δy:
    #
    #   [ β A + Bᵀ B  -Bᵀ ] [ δx ] = [ 0 ]
    #   [ B            0  ] [ δy ]   [ r ]
    #
    function prec!(u, v)
        #
        # solve for u:
        #
        #   (β A + Bᵀ B) u = v
        #
        copyto!(u, v)
        ldiv!(divwrk, F,  u)
        ldiv!(divwrk, F', u)
    end

    N = LinearOperator(T, n, n, true, true, prec!)
    craig!(itrwrk, B, r; ldiv = false, btol = zero(T), N, atol, itmax)
    niter += itrwrk.stats.niter
    #
    # update x:
    #
    #   x = x + δx
    #
    axpy!(one(T), itrwrk.x, x)
    #
    # recover y:
    #
    #   y = y₀ + 1/β (δy + (g - B x))
    #
    copyto!(r, g)
    mul!(r, B, x, -one(T), one(T))
    copyto!(y, itrwrk.y)
    axpy!(one(T), r, y)
    lmul!(α, y)

    if !isnothing(y0)
        axpy!(one(T), y0, y)
    end

    return niter
end
