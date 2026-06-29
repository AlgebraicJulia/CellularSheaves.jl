"""
    SemidefiniteCone <: AbstractCone

A cone of n × n positive-semidefinite
matrices.
"""
struct SemidefiniteCone <: AbstractCone end

# SDP workspace layout in data:
#   M1: n×n at offset 0
#   M2: n×n at offset n²
#   M3: n×n at offset 2n²
#   M4: n×n at offset 3n²
# Total: 4n²
function workspacesize(::SemidefiniteCone, n::Integer)
    d = triroot(n)
    return 4 * d^2
end

struct SemidefiniteConeCache{T} <: AbstractCache{SemidefiniteCone}
    cone::SemidefiniteCone
    #
    # the lower triangular Cholesky factor of
    # the primal variable P:
    #
    #   P = LP LPᵀ
    #
    LP::FMatrixView{T}
    #
    # the lower triangular Cholesky factor of
    # the dual variable D:
    #
    #   D = LD LDᵀ
    #
    LD::FMatrixView{T}
    #
    # the orthogonal factor U in the singular
    # value decomposition
    #
    #   LPᵀ LD = U Σ Vᵀ
    #
    U::FMatrixView{T}
    #
    # the diagonal factor Σ in the singular
    # value decomposition
    #
    #   LPᵀ LD = U Σ Vᵀ
    #
    s::FVectorView{T}
end

function triroot(n::Integer)
    return (isqrt(1 + 8n) - 1) ÷ 2
end

function roottwo(::Type{T}) where {T}
    return sqrt(two(T))
end

function symmetrize!(M::AbstractMatrix)
    for j in axes(M, 1)
        for i in 1:j - 1
            M[i, j] = M[j, i]
        end
    end

    return M
end

function svec!(v::AbstractVector{T}, M::AbstractMatrix{T}) where {T}
    n = size(M, 1); k = 0
    α = roottwo(T)

    for j in 1:n
        k += 1; v[k] = M[j, j]

        for i in j + 1:n
            k += 1; v[k] = α * M[i, j]
        end
    end

    return v
end

function smat!(M::AbstractMatrix{T}, v::AbstractVector{T}) where {T}
    n = size(M, 1); k = 0
    α = roottwo(T)

    for j in 1:n
        k += 1; M[j, j] = v[k]

        for i in j + 1:n
            k += 1; M[i, j] = v[k] / α
        end
    end

    return M
end

# compute the symmetric Kronecker product
#
#   H = A ⊗ A
#
function skron!(H::AbstractMatrix{T}, A::AbstractMatrix{T}) where {T}
    n = size(A, 1)
    α = roottwo(T)
    tll = 1

    @inbounds for l in 1:n
        tij = 0

        for j in 1:n
            Ajl = A[j, l]

            tij += 1; H[tij, tll] = Ajl^2

            for i in j + 1:n
                tij += 1; H[tij, tll] = α * A[i, l] * Ajl
            end
        end

        tkl = tll

        for k in l + 1:n
            tkl += 1; tij = 0

            for j in 1:n
                Ajk = A[j, k]
                Ajl = A[j, l]

                tij += 1; H[tij, tkl] = α * Ajk * Ajl

                for i in j + 1:n
                    tij += 1; H[tij, tkl] = A[i, k] * Ajl + A[i, l] * Ajk
                end
            end
        end

        tll += n - l + 1
    end

    return H
end

# compute the symmetric Kronecker product
#
#   H = W⁻¹ ⊗ W⁻¹
#
# where W is the Nesterov-Todd scaling point
#
#   W = √P √(√P D √P)⁻¹ √P
#     = √D √(√D P √D)⁻¹ √D
#
function sdpscale!(
        H::AbstractMatrix{T},
        LP::AbstractMatrix{T},
        LD::AbstractMatrix{T},
        U::AbstractMatrix{T},
        s::AbstractVector{T},
        p::AbstractVector,
        d::AbstractVector,
        wrk::ConeWorkspace{T},
    ) where {T}
    n = size(LP, 1); m = n * n

    V = reshape(view(wrk.data, 0m + 1:1m), n, n)
    W = reshape(view(wrk.data, 1m + 1:2m), n, n)

    work  = wrk.work
    iwork = wrk.iwork

    smat!(LP, p)
    smat!(LD, d)
    #
    # factorize P:
    #
    #   P = LP LPᵀ
    #
    cholesky!(Symmetric(LP, :L))
    #
    # factorize D:
    #
    #   D = LD LDᵀ
    #
    cholesky!(Symmetric(LD, :L))
    #
    # factorize the product LPᵀ LD
    #
    #   LPᵀ LD = U Σ Vᵀ
    #
    tril!(LD)
    mul!(U, LowerTriangular(LP)', LD)
    svd!(s, U, V, work, iwork)
    #
    # compute the inverse W⁻¹ of the Nesterov-Todd
    # scaling point
    #
    #   W = LP U Σ⁻¹ Uᵀ LPᵀ
    #
    mul!(V, U, Diagonal(s))
    mul!(W, V, U')
    ldiv!(LowerTriangular(LP)', W)
    rdiv!(W, LowerTriangular(LP))
    #
    # compute the symmetric Kronecker product
    #
    #   H = W⁻¹ ⊗ W⁻¹
    #
    skron!(H, W)

    return
end

# Compute the corrector term
#
#   σμ Σ⁻¹ - Σ - 𝓛⁻¹(X)
#
# where
#
#   R = L U √Σ⁻¹
#
# is a factor of the Nesterov-Todd
# scaling point W = R Rᵀ, X is the sum
#
#   X = R⁻¹ ΔP ΔD R + Rᵀ ΔD ΔP R⁻ᵀ,
#
# and 𝓛 is the Lyapunov operator
#
#   𝓛(Y) = ΣY + YΣ.
#
function sdpcorr!(
        r::AbstractVector{T},
        L::AbstractMatrix{T},
        U::AbstractMatrix{T},
        s::AbstractVector{T},
        Δp::AbstractVector{T},
        Δd::AbstractVector{T},
        σμ::Real,
        wrk::ConeWorkspace{T},
    ) where {T}
    n = size(L, 1); m = n * n

    ΔP = reshape(view(wrk.data, 0m + 1:1m), n, n)
    ΔD = reshape(view(wrk.data, 1m + 1:2m), n, n)
    W  = reshape(view(wrk.data, 2m + 1:3m), n, n)
    X  = reshape(view(wrk.data, 3m + 1:4m), n, n)

    smat!(ΔP, Δp)
    smat!(ΔD, Δd)
    symmetrize!(ΔD)
    #
    # compute the product
    #
    #   X = Uᵀ L⁻¹ ΔP ΔD L U
    #
    mul!(W, ΔD, LowerTriangular(L))
    mul!(X, Symmetric(ΔP, :L), W)

    ldiv!(LowerTriangular(L), X)

    mul!(W, U', X)
    mul!(X, W,  U)
    #
    # compute
    #
    #   W = # TODO
    #
    for j in 1:n
        sj = s[j]

        for i in 1:j - 1
            W[i, j] = W[j, i] = -weightedmean(s[i], sj, X[i, j], X[j, i])
        end

        W[j, j] = σμ - sj^2 - X[j, j]
    end
    #
    # compute the product
    #
    #   W = L⁻ᵀ U W Uᵀ L⁻¹
    #
    # and write it to r.
    #
    mul!(X, W, U')
    mul!(W, U, X)

    ldiv!(LowerTriangular(L)', W)
    rdiv!(W, LowerTriangular(L))

    svec!(r, W)
    return r
end

# Find the largest number 0 < τ ≤ 1 such that
#
#   L Lᵀ + τ ΔX = L (I + τ M) Lᵀ
#
# is positive definite, where
#
#   M = L⁻¹ ΔX L⁻ᵀ.
#
# This matrix is positive definite if and
# only if M is, so the solution is given by
#
#   τ⁻¹ = max {1, -λ},
#
# where λ is the smallest eigenvalue of L⁻¹ ΔX L⁻ᵀ.
# construct the identity matrix
#
#   I
#
function sdpid!(x::AbstractVector{T}) where {T}
    d = triroot(length(x))
    k = 1

    fill!(x, zero(T))

    for j in 1:d
        x[k] = one(T); k += d - j + 1
    end

    return x
end

function sdpmaxstep(L::LowerTriangular{T}, Δx::AbstractVector{T}, wrk::ConeWorkspace{T}) where {T}
    n = size(L, 1); m = n * n

    M = reshape(view(wrk.data,     1:m    ), n, n)
    W =         view(wrk.data, m + 1:m + n)

    work  = wrk.work
    iwork = wrk.iwork
    #
    # compute the product
    #
    #   M = L⁻¹ ΔX L⁻ᵀ
    #
    smat!(M, Δx)
    symmetrize!(M)
    ldiv!(L, M)
    rdiv!(M, L')
    #
    # λ is the smallest eigenvalue of M
    #
    λ = eigmin!(M, W, work, iwork)

    return inv(max(one(T), -λ))
end

#
# AbstractCone Interface
#

function degree(::SemidefiniteCone, n::Integer)
    return triroot(n)
end

function cachesize(::SemidefiniteCone, n::Integer)
    d = triroot(n)
    return 3d^2 + d
end

function cache(c::Caches, i::Integer, cone::SemidefiniteCone)
    n = c.xcol[i + 1] - c.xcol[i]
    d = triroot(n)

    data = cachedata(c, i)

    LP = reshape(view(data, 0d^2 + 1:1d^2    ), d, d)
    LD = reshape(view(data, 1d^2 + 1:2d^2    ), d, d)
    U  = reshape(view(data, 2d^2 + 1:3d^2    ), d, d)
    s  =         view(data, 3d^2 + 1:3d^2 + d)

    SemidefiniteConeCache(cone, LP, LD, U, s)
end

function identity!(x::AbstractVector, ::SemidefiniteCone)
    return sdpid!(x)
end

function scale!(H::AbstractMatrix{T}, p::AbstractVector{T}, d::AbstractVector{T}, cache::SemidefiniteConeCache{T}, wrk::ConeWorkspace{T}) where {T}
    return sdpscale!(H, cache.LP, cache.LD, cache.U, cache.s, p, d, wrk)
end

function corr!(
        r::AbstractVector{T},
        ::AbstractVector{T},
        ::AbstractVector{T},
        Δp::AbstractVector{T},
        Δd::AbstractVector{T},
        σμ::Real,
        cache::SemidefiniteConeCache{T},
        wrk::ConeWorkspace{T},
    ) where {T}
    return sdpcorr!(r, cache.LP, cache.U, cache.s, Δp, Δd, σμ, wrk)
end

function maxsteps(::AbstractVector{T}, Δp::AbstractVector{T}, ::AbstractVector{T}, Δd::AbstractVector{T}, cache::SemidefiniteConeCache{T}, wrk::ConeWorkspace{T}) where {T}
    τp = sdpmaxstep(LowerTriangular(cache.LP), Δp, wrk)
    τd = sdpmaxstep(LowerTriangular(cache.LD), Δd, wrk)
    return τp, τd
end
