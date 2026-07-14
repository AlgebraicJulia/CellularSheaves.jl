"""
    SemidefiniteCone <: AbstractCone

A cone of n × n positive-semidefinite
matrices.
"""
struct SemidefiniteCone <: AbstractCone end

# SDP workspace layout in data:
#   Correction path (sdpcorr!): 4 × n² matrices
#   Lanczos path:
#     M: n² (the matrix L⁻¹ΔXL⁻ᵀ)
#     C: n² (Cholesky workspace, avoids rebuilding M on recovery)
#     Q: n × (k + 1)
#     α: k
#     β: k
#     ritzwork: 6k
#   Total: 2n² + n(k+1) + 8k
function workspacesize(::SemidefiniteCone, n::Integer)
    d = triroot(n)
    k = min(LANCZOS_ITMAX, d)
    return max(4d^2, 2d^2 + d * (k + 1) + 8k)
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
    #
    # warm start vectors for Lanczos iteration
    #
    qp::FVectorView{T}
    qd::FVectorView{T}
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

@propagate_inbounds function svec!(v::AbstractVector{T}, M::AbstractMatrix{T}) where {T}
    n = size(M, 1)
    @boundscheck checkbounds(M, n, n)
    @boundscheck checkbounds(v, n * (n + 1) ÷ 2)

    k = 0
    α = roottwo(T)

    @inbounds for j in 1:n
        k += 1; v[k] = M[j, j]

        for i in j + 1:n
            k += 1; v[k] = α * M[i, j]
        end
    end

    return v
end

@propagate_inbounds function smat!(M::AbstractMatrix{T}, v::AbstractVector{T}) where {T}
    n = size(M, 1)
    @boundscheck checkbounds(M, n, n)
    @boundscheck checkbounds(v, n * (n + 1) ÷ 2)

    k = 0
    α = roottwo(T)

    @inbounds for j in 1:n
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

# n = 1: scalar case (SVD is trivial)
function sdpscale1!(
        H::AbstractMatrix{T},
        LP::AbstractMatrix{T},
        LD::AbstractMatrix{T},
        U::AbstractMatrix{T},
        s::AbstractVector{T},
        p::AbstractVector,
        d::AbstractVector,
    ) where {T}
    @inbounds smat!(LP, p)
    @inbounds smat!(LD, d)

    @inbounds FP = cholesky!(Symmetric(LP, :L); check=false)
    issuccess(FP) || return false

    @inbounds FD = cholesky!(Symmetric(LD, :L); check=false)
    issuccess(FD) || return false

    @inbounds m = LP[1,1] * LD[1,1]
    @inbounds s[1] = abs(m)
    @inbounds U[1,1] = m ≥ zero(T) ? one(T) : -one(T)

    @inbounds w = s[1] / LP[1,1]^2
    @inbounds H[1,1] = w * w

    return true
end

# n = 2: Jacobi SVD
function sdpscale2!(
        H::AbstractMatrix{T},
        LP::AbstractMatrix{T},
        LD::AbstractMatrix{T},
        U::AbstractMatrix{T},
        s::AbstractVector{T},
        p::AbstractVector,
        d::AbstractVector,
        wrk::ConeWorkspace{T},
    ) where {T}
    A = reshape(view(wrk.data, 1:4), 2, 2)
    B = reshape(view(wrk.data, 5:8), 2, 2)
    W = reshape(view(wrk.data, 9:12), 2, 2)

    @inbounds smat!(LP, p)
    @inbounds smat!(LD, d)

    @inbounds chol2!(Symmetric(LP, :L)) || return false
    @inbounds chol2!(Symmetric(LD, :L)) || return false

    @inbounds mul2!(A, LowerTriangular(LP)', LowerTriangular(LD))
    @inbounds svd2!(U, s, A)

    @inbounds mul2!(B, U, Diagonal(s))
    @inbounds mul2!(W, B, U')
    @inbounds ldiv2!(LowerTriangular(LP)', Symmetric(W, :L))
    @inbounds rdiv2!(W, LowerTriangular(LP))

    @inbounds skron!(H, W)

    return true
end

# n = 3: Jacobi SVD
function sdpscale3!(
        H::AbstractMatrix{T},
        LP::AbstractMatrix{T},
        LD::AbstractMatrix{T},
        U::AbstractMatrix{T},
        s::AbstractVector{T},
        p::AbstractVector,
        d::AbstractVector,
        wrk::ConeWorkspace{T},
    ) where {T}
    A = reshape(view(wrk.data, 1:9), 3, 3)
    B = reshape(view(wrk.data, 10:18), 3, 3)
    W = reshape(view(wrk.data, 19:27), 3, 3)

    @inbounds smat!(LP, p)
    @inbounds smat!(LD, d)

    @inbounds chol3!(Symmetric(LP, :L)) || return false
    @inbounds chol3!(Symmetric(LD, :L)) || return false

    @inbounds mul3!(A, LowerTriangular(LP)', LowerTriangular(LD))
    @inbounds svd3!(U, s, A)

    @inbounds mul3!(B, U, Diagonal(s))
    @inbounds mul3!(W, B, U')
    @inbounds ldiv3!(LowerTriangular(LP)', Symmetric(W, :L))
    @inbounds rdiv3!(W, LowerTriangular(LP))

    @inbounds skron!(H, W)

    return true
end

# n = 4: Jacobi SVD
function sdpscale4!(
        H::AbstractMatrix{T},
        LP::AbstractMatrix{T},
        LD::AbstractMatrix{T},
        U::AbstractMatrix{T},
        s::AbstractVector{T},
        p::AbstractVector,
        d::AbstractVector,
        wrk::ConeWorkspace{T},
    ) where {T}
    A = reshape(view(wrk.data, 1:16), 4, 4)
    B = reshape(view(wrk.data, 17:32), 4, 4)
    W = reshape(view(wrk.data, 33:48), 4, 4)

    @inbounds smat!(LP, p)
    @inbounds smat!(LD, d)

    @inbounds chol4!(Symmetric(LP, :L)) || return false
    @inbounds chol4!(Symmetric(LD, :L)) || return false

    @inbounds mul4!(A, LowerTriangular(LP)', LowerTriangular(LD))
    @inbounds svd4!(U, s, A)

    @inbounds mul4!(B, U, Diagonal(s))
    @inbounds mul4!(W, B, U')
    @inbounds ldiv4!(LowerTriangular(LP)', Symmetric(W, :L))
    @inbounds rdiv4!(W, LowerTriangular(LP))

    @inbounds skron!(H, W)

    return true
end

# n ≥ 5: LAPACK SVD
function sdpscalelapack!(
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

    FP = cholesky!(Symmetric(LP, :L); check=false)
    issuccess(FP) || return false

    FD = cholesky!(Symmetric(LD, :L); check=false)
    issuccess(FD) || return false

    tril!(LD)
    mul!(U, LowerTriangular(LP)', LD)
    svd!(s, U, V, work, iwork)

    mul!(V, U, Diagonal(s))
    mul!(W, V, U')
    ldiv!(LowerTriangular(LP)', W)
    rdiv!(W, LowerTriangular(LP))

    skron!(H, W)

    return true
end

# Dispatcher
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
    n = size(LP, 1)
    if n == 1
        sdpscale1!(H, LP, LD, U, s, p, d)
    elseif n == 2
        sdpscale2!(H, LP, LD, U, s, p, d, wrk)
    elseif n == 3
        sdpscale3!(H, LP, LD, U, s, p, d, wrk)
    elseif n == 4
        sdpscale4!(H, LP, LD, U, s, p, d, wrk)
    else
        sdpscalelapack!(H, LP, LD, U, s, p, d, wrk)
    end
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

#
# Tiered step-length computation: dispatch by matrix size
#
# n = 1,2,3,4: closed-form eigmin (inline trsm + analytic formula)
# n ≥ 5: matrix-free Lanczos
#

const LANCZOS_ITMAX = 32

# n = 1: scalar case
function sdpmaxstep1(L::LowerTriangular{T}, Δx::AbstractVector{T}, wrk::ConeWorkspace{T}) where T
    M = reshape(view(wrk.data, 1:1), 1, 1)
    @inbounds smat!(M, Δx)
    @inbounds ldiv1!(L, Symmetric(M, :L))
    @inbounds rdiv1!(M, L')
    @inbounds λ = eigmin1(M)
    return inv(max(one(T), -λ))
end

# n = 2: quadratic eigmin
function sdpmaxstep2(L::LowerTriangular{T}, Δx::AbstractVector{T}, wrk::ConeWorkspace{T}) where T
    M = reshape(view(wrk.data, 1:4), 2, 2)
    @inbounds smat!(M, Δx)
    @inbounds ldiv2!(L, Symmetric(M, :L))
    @inbounds rdiv2!(M, L')
    @inbounds λ = eigmin2(M)
    return inv(max(one(T), -λ))
end

# n = 3: trigonometric eigmin
function sdpmaxstep3(L::LowerTriangular{T}, Δx::AbstractVector{T}, wrk::ConeWorkspace{T}) where T
    M = reshape(view(wrk.data, 1:9), 3, 3)
    @inbounds smat!(M, Δx)
    @inbounds ldiv3!(L, Symmetric(M, :L))
    @inbounds rdiv3!(M, L')
    @inbounds λ = eigmin3(M)
    return inv(max(one(T), -λ))
end

# n = 4: quartic eigmin
function sdpmaxstep4(L::LowerTriangular{T}, Δx::AbstractVector{T}, wrk::ConeWorkspace{T}) where T
    M = reshape(view(wrk.data, 1:16), 4, 4)
    @inbounds smat!(M, Δx)
    @inbounds ldiv4!(L, Symmetric(M, :L))
    @inbounds rdiv4!(M, L')
    @inbounds λ = eigmin4(M)
    return inv(max(one(T), -λ))
end

# n ≥ 5: Lanczos with Cholesky gate and deflated recovery
#
# blend: noise scale relative to ‖Q0‖. Injects an O(η/√n) component on every
# eigenvector, defeating structured deficiency where the warm vector is
# orthogonal to a well-separated minimum.
#
function sdpmaxsteplanczos(
        L::LowerTriangular{T},
        Δx::AbstractVector{T},
        wrk::ConeWorkspace{T},
        Q0::AbstractVector{T};
        tol::T = T(1e-3),
        blend::T = T(0.05),
        itmax::Int = 32
    ) where {T}
    n = size(L, 1)
    jmax = min(itmax, n)

    o = 0
    M        = reshape(view(wrk.data, o + 1:o +  n * n         ), n, n       ); o += n * n
    C        = reshape(view(wrk.data, o + 1:o +  n * n         ), n, n       ); o += n * n
    Q        = reshape(view(wrk.data, o + 1:o +  n * (jmax + 1)), n, jmax + 1); o += n * (jmax + 1)
    α        =         view(wrk.data, o + 1:o +  jmax);                         o += jmax
    β        =         view(wrk.data, o + 1:o +  jmax);                         o += jmax
    ritzwork =         view(wrk.data, o + 1:o + 6jmax)
    #
    # form M = L⁻¹ ΔX L⁻ᵀ
    #
    smat!(M, Δx)
    symmetrize!(M)
    rdiv!(M, L')
    ldiv!(L, M)
    #
    # initialize Q[:,1] with blended warm start
    #
    nQ0sq = normsq(Q0)

    if nQ0sq > zero(T) && isfinite(nQ0sq)
        η = blend * sqrt(nQ0sq)
    else
        η = one(T)
    end

    Q1 = view(Q, :, 1)

    @inbounds for i in 1:n
        Q1[i] = Q0[i] + η * randn(wrk.rng, T)
    end

    ldiv!(sqrt(normsq(Q1)), Q1)

    return eigmin!(Symmetric(M, :L), C, Q, α, β, ritzwork, Q0, wrk.work, wrk.iwork; tol=tol)
end

# Dispatcher
function sdpmaxstep(L::LowerTriangular{T}, Δx::AbstractVector{T}, wrk::ConeWorkspace{T}, Q0::AbstractVector{T}) where T
    n = size(L, 1)

    if n == 1
        sdpmaxstep1(L, Δx, wrk)
    elseif n == 2
        sdpmaxstep2(L, Δx, wrk)
    elseif n == 3
        sdpmaxstep3(L, Δx, wrk)
    elseif n == 4
        sdpmaxstep4(L, Δx, wrk)
    else
        sdpmaxsteplanczos(L, Δx, wrk, Q0; tol = T(1e-3), blend = T(0.05), itmax = LANCZOS_ITMAX)
    end
end

#
# AbstractCone Interface
#

function degree(::SemidefiniteCone, n::Integer)
    return triroot(n)
end

function cachesize(::SemidefiniteCone, n::Integer)
    d = triroot(n)
    return 3d^2 + 3d  # LP, LD, U (d² each) + s, qp, qd (d each)
end

function cache(c::Caches, i::Integer, cone::SemidefiniteCone)
    n = c.xcol[i + 1] - c.xcol[i]
    d = triroot(n)

    data = cachedata(c, i)

    LP = reshape(view(data, 0d^2 + 1:1d^2      ), d, d)
    LD = reshape(view(data, 1d^2 + 1:2d^2      ), d, d)
    U  = reshape(view(data, 2d^2 + 1:3d^2      ), d, d)
    s  =         view(data, 3d^2 + 1:3d^2 +  d)
    qp =         view(data, 3d^2 + d + 1:3d^2 + 2d)
    qd =         view(data, 3d^2 + 2d + 1:3d^2 + 3d)

    SemidefiniteConeCache(cone, LP, LD, U, s, qp, qd)
end

function initcache!(c::SemidefiniteConeCache)
    fill!(c.qp, false)
    fill!(c.qd, false)
    return c
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
    τp = sdpmaxstep(LowerTriangular(cache.LP), Δp, wrk, cache.qp)
    τd = sdpmaxstep(LowerTriangular(cache.LD), Δd, wrk, cache.qd)
    return τp, τd
end
