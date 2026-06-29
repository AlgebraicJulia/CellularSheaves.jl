@kwdef struct UzawaSettings{T} <: KKTSettings{T}
    aaug::T = zero(T)
    raug::T = 1e6
    atol::T = √eps(T)
    rtol::T = √eps(T)
    itmax::Int = 1000
    rgmin::T = 1e-9
    rgmax::T = 1e-6
end

function showsettings(io::IO, set::UzawaSettings; indent::Integer=0)
    pad = " "^indent
    @printf(io, "%saaug:  %8.2e  raug:  %8.2e\n", pad, set.aaug, set.raug)
    @printf(io, "%satol:  %8.2e  rtol:  %8.2e\n", pad, set.atol, set.rtol)
    @printf(io, "%srgmin: %8.2e  rgmax: %8.2e\n", pad, set.rgmin, set.rgmax)
    @printf(io, "%sitmax: %8d\n", pad, set.itmax)
    return
end

function Base.show(io::IO, ::MIME"text/plain", set::T) where {T <: UzawaSettings}
    println(io, T, ":")
    return showsettings(io, set; indent=2)
end

struct UzawaWorkspace{UPLO, T, I <: Integer, ItrWrk <: IterationWorkspace{T}} <: KKTWorkspace{T}
    F::FChordalTriangular{:N, UPLO, T, I}
    L::BlockSparseMatrix{T, I}
    facwrk::FactorizationWorkspace{T, I}
    divwrk::DivisionWorkspace{T, I}
    itrwrk::ItrWrk
    r::Vector{T}
    α::Scalar{T}
    nrm::T
end

function UzawaWorkspace(F::FChordalTriangular{:N, UPLO, T, I}, L::BlockSparseMatrix{T, I}, B::BlockSparseMatrix{T, I}) where {UPLO, T, I <: Integer}
    m = size(B, 1)
    facwrk = FactorizationWorkspace(F)
    divwrk = DivisionWorkspace(F, 1)
    itrwrk = CgWorkspace(m, m, Vector{T})
    r = zeros(T, m)
    α = ones(T)
    nrm = norm(B)^2
    return UzawaWorkspace(F, L, facwrk, divwrk, itrwrk, r, α, nrm)
end

function make_kkt(::UzawaSettings{T}, B::BlockSparseMatrix{T, I}) where {T, I}
    weights, graph = weightedgraph(B)

    R, P, S = symbolic(weights, graph)

    B = selectvtxs(B, R.perm)

    F = FChordalTriangular{:N, :L, T, I}(S)
    L = B' * B

    wrk = UzawaWorkspace(F, L, B)

    return R, P, B, wrk
end

function initkkt!(wrk::UzawaWorkspace{UPLO, T}, set::UzawaSettings{T}, A::BlockSparseMatrix) where {UPLO, T}
    wrk.α[] = α = set.aaug + set.raug * norm(Symmetric(A, :L)) / wrk.nrm
    return init_uzw!(wrk.facwrk, wrk.F, wrk.L, A, α, set.rgmin, set.rgmax)
end

# form the augmented block
#
#   F = A + α Bᵀ B
#
# and factorize it. If the
# factorization fails, increasing
# diagonal perturbations
#
#   rgmin ≤ ρ ≤ rgmin
#
# are applied until it succeeds.
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

    ρ = rgmin

    copyto!(F, A)
    axpy!(α, L, F)
    info = cholesky!(facwrk, F; check=false)

    while !iszero(info) && ρ ≤ rgmax
        copyto!(F, A)
        axpy!(α, L, F)
        axpy!(ρ, I, F)
        info = cholesky!(facwrk, F; check=false)
        ρ *= 2
    end

    return iszero(info)
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
    rtol::T = set.rtol,
) where {UPLO, T}
    rtol = max(set.rtol, rtol)
    return solve_uzw!(wrk.divwrk, wrk.itrwrk, x, y, wrk.r, wrk.F, B, f, g, wrk.α[], set.atol, rtol, set.itmax, y0)
end

#
# Solve the KKT system
#
#   [ A -Bᵀ ] [ x ] = [ f ]
#   [ B  0  ] [ y ]   [ g ]
#
function solve_uzw!(
        divwrk::DivisionWorkspace{T},
        itrwrk::IterationWorkspace{T},
        x::AbstractVector{T},
        y::AbstractVector{T},
        r::AbstractVector{T},
        F::ChordalTriangular{:N, UPLO, T},
        B::BlockSparseMatrix{T},
        f::AbstractVector{T},
        g::AbstractVector{T},
        α::T,
        atol::T,
        rtol::T,
        itmax::Int,
        y0 = nothing
    ) where {UPLO, T}
    m, n = size(B)

    @assert length(x) == n
    @assert length(y) == m
    @assert length(r) == m
    @assert length(f) == n
    @assert length(g) == m
    @assert size(F, 1) == n
    #
    # solve for x:
    #
    #   F Fᵀ x = f + α Bᵀ g
    #
    copyto!(x, f)
    mul!(x, B', g, α, 1)
    ldiv!(divwrk, F, x)
    ldiv!(divwrk, F', x)
    #
    # compute the residual
    #
    #   r = g - B x
    #
    copyto!(r, g)
    mul!(r, B, x, -1, 1)

    function schur!(u, b)
        #
        # compute
        #
        #   u = B x
        #
        # where x solves
        #
        #   F F' x = Bᵀ b
        #
        mul!(x, B', b)
        ldiv!(divwrk, F, x)
        ldiv!(divwrk, F', x)
        mul!(u, B, x)
    end
    #
    # S is the augmented Schur complement:
    #
    #   S = B (F Fᵀ)⁻¹ Bᵀ
    #
    S = LinearOperator(T, m, m, true, true, schur!)
    #
    # solve for y:
    #
    #   S y = r
    #
    it!(itrwrk, S, r, y0; α, atol, rtol, itmax)

    copyto!(y, solution(itrwrk))
    #
    # compute the dual correction
    #
    #   r = y + α g
    #
    copyto!(r, y)
    axpy!(α, g, r)
    #
    # solve for x:
    #
    #   F Fᵀ x = f + Bᵀ r
    #
    copyto!(x, f)
    mul!(x, B', r, 1, 1)
    ldiv!(divwrk, F, x)
    ldiv!(divwrk, F', x)
    #
    # update y:
    #
    #   y = y + α (g - B x)
    #
    copyto!(r, g)
    mul!(r, B, x, -1, 1)
    axpy!(α, r, y)

    return niter(itrwrk)
end
