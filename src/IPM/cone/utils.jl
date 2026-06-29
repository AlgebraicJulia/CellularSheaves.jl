#
# 3×3 linear algebra helpers
#

# Cross product of two 3-vectors: z = x × y
function cross3!(z::AbstractVector, x::AbstractVector, y::AbstractVector)
    z[1] = x[2] * y[3] - x[3] * y[2]
    z[2] = x[3] * y[1] - x[1] * y[3]
    z[3] = x[1] * y[2] - x[2] * y[1]
    return z
end

# 3×3 gemm: C = α*A*B + β*C (works for matrix-vector and matrix-matrix)
function mul3!(C::AbstractArray, A::AbstractArray, B::AbstractArray)
    return mul3!(C, A, B, true, false)
end

function mul3!(C::AbstractVector, A::AbstractMatrix, B::AbstractVector, α::Number, β::Number)
    if iszero(β)
        C[1] = α * (A[1,1] * B[1] + A[1,2] * B[2] + A[1,3] * B[3])
        C[2] = α * (A[2,1] * B[1] + A[2,2] * B[2] + A[2,3] * B[3])
        C[3] = α * (A[3,1] * B[1] + A[3,2] * B[2] + A[3,3] * B[3])
    else
        C[1] = α * (A[1,1] * B[1] + A[1,2] * B[2] + A[1,3] * B[3]) + β * C[1]
        C[2] = α * (A[2,1] * B[1] + A[2,2] * B[2] + A[2,3] * B[3]) + β * C[2]
        C[3] = α * (A[3,1] * B[1] + A[3,2] * B[2] + A[3,3] * B[3]) + β * C[3]
    end

    return C
end

function mul3!(C::AbstractMatrix, A::AbstractMatrix, B::AbstractMatrix, α::Number, β::Number)
    if iszero(β)
        C[1,1] = α * (A[1,1] * B[1,1] + A[1,2] * B[2,1] + A[1,3] * B[3,1])
        C[2,1] = α * (A[2,1] * B[1,1] + A[2,2] * B[2,1] + A[2,3] * B[3,1])
        C[3,1] = α * (A[3,1] * B[1,1] + A[3,2] * B[2,1] + A[3,3] * B[3,1])
        C[1,2] = α * (A[1,1] * B[1,2] + A[1,2] * B[2,2] + A[1,3] * B[3,2])
        C[2,2] = α * (A[2,1] * B[1,2] + A[2,2] * B[2,2] + A[2,3] * B[3,2])
        C[3,2] = α * (A[3,1] * B[1,2] + A[3,2] * B[2,2] + A[3,3] * B[3,2])
        C[1,3] = α * (A[1,1] * B[1,3] + A[1,2] * B[2,3] + A[1,3] * B[3,3])
        C[2,3] = α * (A[2,1] * B[1,3] + A[2,2] * B[2,3] + A[2,3] * B[3,3])
        C[3,3] = α * (A[3,1] * B[1,3] + A[3,2] * B[2,3] + A[3,3] * B[3,3])
    else
        C[1,1] = α * (A[1,1] * B[1,1] + A[1,2] * B[2,1] + A[1,3] * B[3,1]) + β * C[1,1]
        C[2,1] = α * (A[2,1] * B[1,1] + A[2,2] * B[2,1] + A[2,3] * B[3,1]) + β * C[2,1]
        C[3,1] = α * (A[3,1] * B[1,1] + A[3,2] * B[2,1] + A[3,3] * B[3,1]) + β * C[3,1]
        C[1,2] = α * (A[1,1] * B[1,2] + A[1,2] * B[2,2] + A[1,3] * B[3,2]) + β * C[1,2]
        C[2,2] = α * (A[2,1] * B[1,2] + A[2,2] * B[2,2] + A[2,3] * B[3,2]) + β * C[2,2]
        C[3,2] = α * (A[3,1] * B[1,2] + A[3,2] * B[2,2] + A[3,3] * B[3,2]) + β * C[3,2]
        C[1,3] = α * (A[1,1] * B[1,3] + A[1,2] * B[2,3] + A[1,3] * B[3,3]) + β * C[1,3]
        C[2,3] = α * (A[2,1] * B[1,3] + A[2,2] * B[2,3] + A[2,3] * B[3,3]) + β * C[2,3]
        C[3,3] = α * (A[3,1] * B[1,3] + A[3,2] * B[2,3] + A[3,3] * B[3,3]) + β * C[3,3]
    end

    return C
end

# M ← α x yᵀ + β M (3×3 rank-1 update)
function ger3!(M::AbstractMatrix, x::AbstractVector, y::AbstractVector, α::Number, β::Number)
    if iszero(β)
        M[1,1] = α * x[1] * y[1]
        M[2,1] = α * x[2] * y[1]
        M[3,1] = α * x[3] * y[1]
        M[1,2] = α * x[1] * y[2]
        M[2,2] = α * x[2] * y[2]
        M[3,2] = α * x[3] * y[2]
        M[1,3] = α * x[1] * y[3]
        M[2,3] = α * x[2] * y[3]
        M[3,3] = α * x[3] * y[3]
    else
        M[1,1] = α * x[1] * y[1] + β * M[1,1]
        M[2,1] = α * x[2] * y[1] + β * M[2,1]
        M[3,1] = α * x[3] * y[1] + β * M[3,1]
        M[1,2] = α * x[1] * y[2] + β * M[1,2]
        M[2,2] = α * x[2] * y[2] + β * M[2,2]
        M[3,2] = α * x[3] * y[2] + β * M[3,2]
        M[1,3] = α * x[1] * y[3] + β * M[1,3]
        M[2,3] = α * x[2] * y[3] + β * M[2,3]
        M[3,3] = α * x[3] * y[3] + β * M[3,3]
    end

    return M
end

# 3-element BLAS-style operations
function copy3!(y::AbstractVector, x::AbstractVector)
    y[1] = x[1]
    y[2] = x[2]
    y[3] = x[3]
    return y
end

function axpy3!(a, x::AbstractVector, y::AbstractVector)
    y[1] += a * x[1]
    y[2] += a * x[2]
    y[3] += a * x[3]
    return y
end

function axpby3!(a, x::AbstractVector, b, y::AbstractVector)
    y[1] = a * x[1] + b * y[1]
    y[2] = a * x[2] + b * y[2]
    y[3] = a * x[3] + b * y[3]
    return y
end

function lmul3!(a::Number, x::AbstractVector)
    x[1] *= a
    x[2] *= a
    x[3] *= a
    return x
end

function lmul3!(a::Number, M::AbstractMatrix)
    M[1,1] *= a; M[2,1] *= a; M[3,1] *= a
    M[1,2] *= a; M[2,2] *= a; M[3,2] *= a
    M[1,3] *= a; M[2,3] *= a; M[3,3] *= a
    return M
end

function ldiv3!(a::Number, x::AbstractVector)
    x[1] /= a
    x[2] /= a
    x[3] /= a
    return x
end

function fill3!(x::AbstractVector, a)
    x[1] = a
    x[2] = a
    x[3] = a
    return x
end

function fill3!(M::AbstractMatrix, a)
    M[1,1] = a; M[2,1] = a; M[3,1] = a
    M[1,2] = a; M[2,2] = a; M[3,2] = a
    M[1,3] = a; M[2,3] = a; M[3,3] = a
    return M
end

function dot3(x::AbstractVector, y::AbstractVector)
    return x[1] * y[1] + x[2] * y[2] + x[3] * y[3]
end

function norm3(x::AbstractVector)
    return sqrt(x[1]^2 + x[2]^2 + x[3]^2)
end

# Solve 3×3 system A x = b using Cramer's rule, store result in x
function ldiv3!(x::AbstractVector{T}, A::AbstractMatrix{T}, b::AbstractVector{T}) where {T}
    a11, a12, a13 = A[1,1], A[1,2], A[1,3]
    a21, a22, a23 = A[2,1], A[2,2], A[2,3]
    a31, a32, a33 = A[3,1], A[3,2], A[3,3]
    b1, b2, b3 = b[1], b[2], b[3]

    # 2×2 minors (cofactors of row 1)
    m11 = a22 * a33 - a23 * a32
    m12 = a21 * a33 - a23 * a31
    m13 = a21 * a32 - a22 * a31

    det = a11 * m11 - a12 * m12 + a13 * m13

    # x1: det of A with col 1 replaced by b
    n1 = b1 * m11 - a12 * (b2 * a33 - a23 * b3) + a13 * (b2 * a32 - a22 * b3)

    # x2: det of A with col 2 replaced by b
    n2 = a11 * (b2 * a33 - a23 * b3) - b1 * m12 + a13 * (a21 * b3 - b2 * a31)

    # x3: det of A with col 3 replaced by b
    n3 = a11 * (a22 * b3 - b2 * a32) - a12 * (a21 * b3 - b2 * a31) + b1 * m13

    x[1], x[2], x[3] = n1 / det, n2 / det, n3 / det
    return x
end

# 2-arg version: solve in-place (b overwritten with solution)
ldiv3!(A::AbstractMatrix{T}, b::AbstractVector{T}) where {T} = ldiv3!(b, A, b)

# In-place Cholesky factorization of 3×3 symmetric positive definite matrix
# Overwrites lower triangle of A with L such that A = L L'
function chol3!(A::AbstractMatrix{T}) where {T}
    A[1,1] = sqrt(A[1,1])
    A[2,1] = A[2,1] / A[1,1]
    A[3,1] = A[3,1] / A[1,1]
    A[2,2] = sqrt(A[2,2] - A[2,1]^2)
    A[3,2] = (A[3,2] - A[3,1] * A[2,1]) / A[2,2]
    A[3,3] = sqrt(A[3,3] - A[3,1]^2 - A[3,2]^2)
    return A
end

# Forward substitution: solve L b = b in-place (L lower triangular, stored in lower part of A)
function ldiv3!(L::LowerTriangular, b::AbstractVector)
    A = parent(L)
    b[1] = b[1] / A[1,1]
    b[2] = (b[2] - A[2,1] * b[1]) / A[2,2]
    b[3] = (b[3] - A[3,1] * b[1] - A[3,2] * b[2]) / A[3,3]
    return b
end

# Back substitution: solve L' b = b in-place
function ldiv3!(Lt::Adjoint{<:Any, <:LowerTriangular}, b::AbstractVector)
    A = parent(parent(Lt))
    b[3] = b[3] / A[3,3]
    b[2] = (b[2] - A[3,2] * b[3]) / A[2,2]
    b[1] = (b[1] - A[2,1] * b[2] - A[3,1] * b[3]) / A[1,1]
    return b
end

# Binary search for last t in [lo, hi] where f(t) is true
function binarysearchlast(f, lo::T, hi::T, tol::T, itmax::Integer) where {T}
    for _ in 1:itmax
        mid = (lo + hi) / 2

        if f(mid)
            lo = mid
        else
            hi = mid
        end

        if hi - lo < tol
            break
        end
    end

    return lo
end

#
# Safeguarded Newton on a scalar root of f in [lo, hi],
# where f(lo) < 0 < f(hi) (increasing across bracket).
# Takes the Newton step r - f/f' when it lands inside
# the bracket, otherwise a bisection step. Converges on
# the RELATIVE step size.
#
# For decreasing functions, negate f and f' at the call site.
#
function rtsafe(f, fp, lo::T, hi::T, r0::T;
                tol::T = T(1e-12), maxit::Int = 60) where {T}
    r = clamp(r0, lo, hi)

    for _ in 1:maxit
        rp = r
        fr = f(r)

        # tighten the bracket toward the root
        if fr >= 0
            hi = r
        else
            lo = r
        end

        # Newton step, safeguarded into a bisection
        # step if it would leave the bracket
        r = r - fr / fp(r)

        if !(lo < r < hi)
            r = (lo + hi) / 2
        end

        # converge on the relative step size
        if abs(r - rp) < tol * (one(T) + abs(r))
            break
        end
    end

    return r
end