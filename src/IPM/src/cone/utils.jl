#
# Small linear algebra helpers for IPM cone operations
#
# Organization:
#   1. 3×3 vector/matrix helpers (cross, dot, norm, copy, axpy, lmul, etc.)
#   2. mul - matrix multiplication (dense loops, triangular inlined)
#   3. chol - Cholesky factorization (inlined)
#   4. ldiv - left triangular solves (inlined)
#   5. rdiv - right triangular solves (inlined)
#   6. Newton solvers (rtsafe, nflast)
#

# ============================================================================
# 3×3 vector/matrix helpers
# ============================================================================

function cross3!(z::AbstractVector, x::AbstractVector, y::AbstractVector)
    z[1] = x[2] * y[3] - x[3] * y[2]
    z[2] = x[3] * y[1] - x[1] * y[3]
    z[3] = x[1] * y[2] - x[2] * y[1]
    return z
end

function dot3(x::AbstractVector, y::AbstractVector)
    return x[1] * y[1] + x[2] * y[2] + x[3] * y[3]
end

function norm3(x::AbstractVector)
    return sqrt(x[1]^2 + x[2]^2 + x[3]^2)
end

function normsq(x::AbstractVector{T}) where {T}
    out = zero(T)

    @inbounds for i in eachindex(x)
        out += x[i]^2
    end

    return out
end

@propagate_inbounds function copy3!(y::AbstractVector, x::AbstractVector)
    @boundscheck checkbounds(y, 3)
    @boundscheck checkbounds(x, 3)

    @inbounds for i in 1:3
        y[i] = x[i]
    end

    return y
end

@propagate_inbounds function axpy3!(a, x::AbstractVector, y::AbstractVector)
    @boundscheck checkbounds(x, 3)
    @boundscheck checkbounds(y, 3)

    @inbounds for i in 1:3
        y[i] += a * x[i]
    end

    return y
end

@propagate_inbounds function axpby3!(a, x::AbstractVector, b, y::AbstractVector)
    @boundscheck checkbounds(x, 3)
    @boundscheck checkbounds(y, 3)

    @inbounds for i in 1:3
        y[i] = a * x[i] + b * y[i]
    end

    return y
end

@propagate_inbounds function lmul3!(a::Number, x::AbstractVector)
    @boundscheck checkbounds(x, 3)

    @inbounds for i in 1:3
        x[i] *= a
    end

    return x
end

@propagate_inbounds function lmul3!(a::Number, M::AbstractMatrix)
    @boundscheck checkbounds(M, 3, 3)

    @inbounds for j in 1:3, i in 1:3
        M[i,j] *= a
    end

    return M
end

@propagate_inbounds function ldiv3!(a::Number, x::AbstractVector)
    @boundscheck checkbounds(x, 3)

    @inbounds for i in 1:3
        x[i] /= a
    end

    return x
end

# ============================================================================
# mul - matrix multiplication
# ============================================================================

# ----------------------------------------------------------------------------
# General gemm: C = α*A*B + β*C (generated)
# ----------------------------------------------------------------------------

for n in (2, 3, 4)
    fname = Symbol(:mul, n, :!)
    @eval begin
        $fname(C::AbstractArray, A::AbstractArray, B::AbstractArray) = $fname(C, A, B, true, false)

        @propagate_inbounds function $fname(C::AbstractMatrix, A::AbstractMatrix, B::AbstractMatrix, α::Number, β::Number)
            @boundscheck checkbounds(C, $n, $n)
            @boundscheck checkbounds(A, $n, $n)
            @boundscheck checkbounds(B, $n, $n)

            T = promote_eltype(A, B)

            if iszero(β)
                @inbounds for i in 1:$n
                    for j in 1:$n
                        s = zero(T)

                        for k in 1:$n
                            s += A[i,k] * B[k,j]
                        end

                        C[i,j] = α * s
                    end
                end
            else
                @inbounds for i in 1:$n
                    for j in 1:$n
                        s = zero(T)

                        for k in 1:$n
                            s += A[i,k] * B[k,j]
                        end

                        C[i,j] = α * s + β * C[i,j]
                    end
                end
            end

            return C
        end
    end
end

# Matrix-vector: c = α*A*b + β*c (3×3 only)
@propagate_inbounds function mul3!(C::AbstractVector, A::AbstractMatrix, B::AbstractVector, α::Number, β::Number)
    @boundscheck checkbounds(C, 3)
    @boundscheck checkbounds(A, 3, 3)
    @boundscheck checkbounds(B, 3)

    T = promote_eltype(A, B)

    if iszero(β)
        @inbounds for i in 1:3
            s = zero(T)

            for k in 1:3
                s += A[i,k] * B[k]
            end

            C[i] = α * s
        end
    else
        @inbounds for i in 1:3
            s = zero(T)

            for k in 1:3
                s += A[i,k] * B[k]
            end

            C[i] = α * s + β * C[i]
        end
    end

    return C
end

# ----------------------------------------------------------------------------
# Matrix × Diagonal: C = α*A*D + β*C (generated)
# ----------------------------------------------------------------------------

for n in (2, 3, 4)
    fname = Symbol(:mul, n, :!)

    @eval @propagate_inbounds function $fname(C::AbstractMatrix, A::AbstractMatrix, D::Diagonal, α::Number, β::Number)
        d = D.diag

        @boundscheck checkbounds(C, $n, $n)
        @boundscheck checkbounds(A, $n, $n)
        @boundscheck checkbounds(d, $n)

        if iszero(β)
            @inbounds for j in 1:$n
                for i in 1:$n
                    C[i,j] = α * A[i,j] * d[j]
                end
            end
        else
            @inbounds for j in 1:$n
                for i in 1:$n
                    C[i,j] = α * A[i,j] * d[j] + β * C[i,j]
                end
            end
        end

        return C
    end
end

# ----------------------------------------------------------------------------
# Upper × Lower: C = U * L (inlined - sparsity pattern)
# ----------------------------------------------------------------------------

@propagate_inbounds function mul2!(C::AbstractMatrix, U::UpperTriangular, L::LowerTriangular)
    UP = parent(U)
    LP = parent(L)
    @boundscheck checkbounds(C, 2, 2)
    @boundscheck checkbounds(UP, 2, 2)
    @boundscheck checkbounds(LP, 2, 2)
    @inbounds U₁₁, U₁₂ = UP[1,1], UP[1,2]
    @inbounds U₂₂ = UP[2,2]
    @inbounds L₁₁, L₂₁, L₂₂ = LP[1,1], LP[2,1], LP[2,2]
    @inbounds C[1,1] = U₁₁ * L₁₁ + U₁₂ * L₂₁
    @inbounds C[1,2] = U₁₂ * L₂₂
    @inbounds C[2,1] = U₂₂ * L₂₁
    @inbounds C[2,2] = U₂₂ * L₂₂
    return C
end

@propagate_inbounds function mul3!(C::AbstractMatrix, U::UpperTriangular, L::LowerTriangular)
    UP = parent(U)
    LP = parent(L)
    @boundscheck checkbounds(C, 3, 3)
    @boundscheck checkbounds(UP, 3, 3)
    @boundscheck checkbounds(LP, 3, 3)
    @inbounds U₁₁, U₁₂, U₁₃ = UP[1,1], UP[1,2], UP[1,3]
    @inbounds U₂₂, U₂₃ = UP[2,2], UP[2,3]
    @inbounds U₃₃ = UP[3,3]
    @inbounds L₁₁, L₂₁, L₃₁ = LP[1,1], LP[2,1], LP[3,1]
    @inbounds L₂₂, L₃₂ = LP[2,2], LP[3,2]
    @inbounds L₃₃ = LP[3,3]
    @inbounds C[1,1] = U₁₁ * L₁₁ + U₁₂ * L₂₁ + U₁₃ * L₃₁
    @inbounds C[1,2] = U₁₂ * L₂₂ + U₁₃ * L₃₂
    @inbounds C[1,3] = U₁₃ * L₃₃
    @inbounds C[2,1] = U₂₂ * L₂₁ + U₂₃ * L₃₁
    @inbounds C[2,2] = U₂₂ * L₂₂ + U₂₃ * L₃₂
    @inbounds C[2,3] = U₂₃ * L₃₃
    @inbounds C[3,1] = U₃₃ * L₃₁
    @inbounds C[3,2] = U₃₃ * L₃₂
    @inbounds C[3,3] = U₃₃ * L₃₃
    return C
end

@propagate_inbounds function mul4!(C::AbstractMatrix, U::UpperTriangular, L::LowerTriangular)
    UP = parent(U)
    LP = parent(L)
    @boundscheck checkbounds(C, 4, 4)
    @boundscheck checkbounds(UP, 4, 4)
    @boundscheck checkbounds(LP, 4, 4)
    @inbounds U₁₁, U₁₂, U₁₃, U₁₄ = UP[1,1], UP[1,2], UP[1,3], UP[1,4]
    @inbounds U₂₂, U₂₃, U₂₄ = UP[2,2], UP[2,3], UP[2,4]
    @inbounds U₃₃, U₃₄ = UP[3,3], UP[3,4]
    @inbounds U₄₄ = UP[4,4]
    @inbounds L₁₁, L₂₁, L₃₁, L₄₁ = LP[1,1], LP[2,1], LP[3,1], LP[4,1]
    @inbounds L₂₂, L₃₂, L₄₂ = LP[2,2], LP[3,2], LP[4,2]
    @inbounds L₃₃, L₄₃ = LP[3,3], LP[4,3]
    @inbounds L₄₄ = LP[4,4]
    @inbounds C[1,1] = U₁₁ * L₁₁ + U₁₂ * L₂₁ + U₁₃ * L₃₁ + U₁₄ * L₄₁
    @inbounds C[1,2] = U₁₂ * L₂₂ + U₁₃ * L₃₂ + U₁₄ * L₄₂
    @inbounds C[1,3] = U₁₃ * L₃₃ + U₁₄ * L₄₃
    @inbounds C[1,4] = U₁₄ * L₄₄
    @inbounds C[2,1] = U₂₂ * L₂₁ + U₂₃ * L₃₁ + U₂₄ * L₄₁
    @inbounds C[2,2] = U₂₂ * L₂₂ + U₂₃ * L₃₂ + U₂₄ * L₄₂
    @inbounds C[2,3] = U₂₃ * L₃₃ + U₂₄ * L₄₃
    @inbounds C[2,4] = U₂₄ * L₄₄
    @inbounds C[3,1] = U₃₃ * L₃₁ + U₃₄ * L₄₁
    @inbounds C[3,2] = U₃₃ * L₃₂ + U₃₄ * L₄₂
    @inbounds C[3,3] = U₃₃ * L₃₃ + U₃₄ * L₄₃
    @inbounds C[3,4] = U₃₄ * L₄₄
    @inbounds C[4,1] = U₄₄ * L₄₁
    @inbounds C[4,2] = U₄₄ * L₄₂
    @inbounds C[4,3] = U₄₄ * L₄₃
    @inbounds C[4,4] = U₄₄ * L₄₄
    return C
end

# ----------------------------------------------------------------------------
# Rank-1 update: M = α*x*y' + β*M (3×3)
# ----------------------------------------------------------------------------

@propagate_inbounds function ger3!(M::AbstractMatrix, x::AbstractVector, y::AbstractVector, α::Number, β::Number)
    @boundscheck checkbounds(M, 3, 3)
    @boundscheck checkbounds(x, 3)
    @boundscheck checkbounds(y, 3)

    if iszero(β)
        @inbounds for j in 1:3, i in 1:3
            M[i,j] = α * x[i] * y[j]
        end
    else
        @inbounds for j in 1:3, i in 1:3
            M[i,j] = α * x[i] * y[j] + β * M[i,j]
        end
    end

    return M
end

# ============================================================================
# chol - Cholesky factorization (inlined, in-place, returns success)
# ============================================================================

@propagate_inbounds function chol2!(M::Symmetric{T}) where T
    A = parent(M)
    @boundscheck checkbounds(A, 2, 2)
    @inbounds A₁₁ = A[1,1]
    A₁₁ > zero(T) || return false
    @inbounds L₁₁ = sqrt(A₁₁)
    @inbounds L₂₁ = A[2,1] / L₁₁
    @inbounds A₂₂ = A[2,2] - L₂₁ * L₂₁
    A₂₂ > zero(T) || return false
    @inbounds A[1,1] = L₁₁
    @inbounds A[2,1] = L₂₁
    @inbounds A[2,2] = sqrt(A₂₂)
    return true
end

@propagate_inbounds function chol3!(M::Symmetric{T}) where T
    A = parent(M)
    @boundscheck checkbounds(A, 3, 3)
    @inbounds A₁₁ = A[1,1]
    A₁₁ > zero(T) || return false
    @inbounds L₁₁ = sqrt(A₁₁)
    @inbounds L₂₁ = A[2,1] / L₁₁
    @inbounds L₃₁ = A[3,1] / L₁₁
    @inbounds A₂₂ = A[2,2] - L₂₁ * L₂₁
    A₂₂ > zero(T) || return false
    @inbounds L₂₂ = sqrt(A₂₂)
    @inbounds L₃₂ = (A[3,2] - L₃₁ * L₂₁) / L₂₂
    @inbounds A₃₃ = A[3,3] - L₃₁ * L₃₁ - L₃₂ * L₃₂
    A₃₃ > zero(T) || return false
    @inbounds A[1,1] = L₁₁
    @inbounds A[2,1] = L₂₁
    @inbounds A[3,1] = L₃₁
    @inbounds A[2,2] = L₂₂
    @inbounds A[3,2] = L₃₂
    @inbounds A[3,3] = sqrt(A₃₃)
    return true
end

@propagate_inbounds function chol4!(M::Symmetric{T}) where T
    A = parent(M)
    @boundscheck checkbounds(A, 4, 4)
    @inbounds A₁₁ = A[1,1]
    A₁₁ > zero(T) || return false
    @inbounds L₁₁ = sqrt(A₁₁)
    @inbounds L₂₁ = A[2,1] / L₁₁
    @inbounds L₃₁ = A[3,1] / L₁₁
    @inbounds L₄₁ = A[4,1] / L₁₁
    @inbounds A₂₂ = A[2,2] - L₂₁ * L₂₁
    A₂₂ > zero(T) || return false
    @inbounds L₂₂ = sqrt(A₂₂)
    @inbounds L₃₂ = (A[3,2] - L₃₁ * L₂₁) / L₂₂
    @inbounds L₄₂ = (A[4,2] - L₄₁ * L₂₁) / L₂₂
    @inbounds A₃₃ = A[3,3] - L₃₁ * L₃₁ - L₃₂ * L₃₂
    A₃₃ > zero(T) || return false
    @inbounds L₃₃ = sqrt(A₃₃)
    @inbounds L₄₃ = (A[4,3] - L₄₁ * L₃₁ - L₄₂ * L₃₂) / L₃₃
    @inbounds A₄₄ = A[4,4] - L₄₁ * L₄₁ - L₄₂ * L₄₂ - L₄₃ * L₄₃
    A₄₄ > zero(T) || return false
    @inbounds A[1,1] = L₁₁
    @inbounds A[2,1] = L₂₁
    @inbounds A[3,1] = L₃₁
    @inbounds A[4,1] = L₄₁
    @inbounds A[2,2] = L₂₂
    @inbounds A[3,2] = L₃₂
    @inbounds A[4,2] = L₄₂
    @inbounds A[3,3] = L₃₃
    @inbounds A[4,3] = L₄₃
    @inbounds A[4,4] = sqrt(A₄₄)
    return true
end

# ============================================================================
# ldiv - left triangular solves (inlined)
# ============================================================================

# ----------------------------------------------------------------------------
# Vector solves: b = L⁻¹ b or b = U⁻¹ b
# ----------------------------------------------------------------------------

@propagate_inbounds function ldiv1!(L::LowerTriangular, b::AbstractVector)
    LP = parent(L)
    @boundscheck checkbounds(LP, 1, 1)
    @boundscheck checkbounds(b, 1)
    @inbounds b[1] = b[1] / LP[1,1]
    return b
end

@propagate_inbounds function ldiv1!(U::UpperTriangular, b::AbstractVector)
    UP = parent(U)
    @boundscheck checkbounds(UP, 1, 1)
    @boundscheck checkbounds(b, 1)
    @inbounds b[1] = b[1] / UP[1,1]
    return b
end

@propagate_inbounds function ldiv2!(L::LowerTriangular, b::AbstractVector)
    LP = parent(L)
    @boundscheck checkbounds(LP, 2, 2)
    @boundscheck checkbounds(b, 2)
    @inbounds b[1] = b[1] / LP[1,1]
    @inbounds b[2] = (b[2] - LP[2,1] * b[1]) / LP[2,2]
    return b
end

@propagate_inbounds function ldiv2!(U::UpperTriangular, b::AbstractVector)
    UP = parent(U)
    @boundscheck checkbounds(UP, 2, 2)
    @boundscheck checkbounds(b, 2)
    @inbounds b[2] = b[2] / UP[2,2]
    @inbounds b[1] = (b[1] - UP[1,2] * b[2]) / UP[1,1]
    return b
end

@propagate_inbounds function ldiv3!(L::LowerTriangular, b::AbstractVector)
    LP = parent(L)
    @boundscheck checkbounds(LP, 3, 3)
    @boundscheck checkbounds(b, 3)
    @inbounds b[1] = b[1] / LP[1,1]
    @inbounds b[2] = (b[2] - LP[2,1] * b[1]) / LP[2,2]
    @inbounds b[3] = (b[3] - LP[3,1] * b[1] - LP[3,2] * b[2]) / LP[3,3]
    return b
end

@propagate_inbounds function ldiv3!(U::UpperTriangular, b::AbstractVector)
    UP = parent(U)
    @boundscheck checkbounds(UP, 3, 3)
    @boundscheck checkbounds(b, 3)
    @inbounds b[3] = b[3] / UP[3,3]
    @inbounds b[2] = (b[2] - UP[2,3] * b[3]) / UP[2,2]
    @inbounds b[1] = (b[1] - UP[1,2] * b[2] - UP[1,3] * b[3]) / UP[1,1]
    return b
end

@propagate_inbounds function ldiv4!(L::LowerTriangular, b::AbstractVector)
    LP = parent(L)
    @boundscheck checkbounds(LP, 4, 4)
    @boundscheck checkbounds(b, 4)
    @inbounds b[1] = b[1] / LP[1,1]
    @inbounds b[2] = (b[2] - LP[2,1] * b[1]) / LP[2,2]
    @inbounds b[3] = (b[3] - LP[3,1] * b[1] - LP[3,2] * b[2]) / LP[3,3]
    @inbounds b[4] = (b[4] - LP[4,1] * b[1] - LP[4,2] * b[2] - LP[4,3] * b[3]) / LP[4,4]
    return b
end

@propagate_inbounds function ldiv4!(U::UpperTriangular, b::AbstractVector)
    UP = parent(U)
    @boundscheck checkbounds(UP, 4, 4)
    @boundscheck checkbounds(b, 4)
    @inbounds b[4] = b[4] / UP[4,4]
    @inbounds b[3] = (b[3] - UP[3,4] * b[4]) / UP[3,3]
    @inbounds b[2] = (b[2] - UP[2,3] * b[3] - UP[2,4] * b[4]) / UP[2,2]
    @inbounds b[1] = (b[1] - UP[1,2] * b[2] - UP[1,3] * b[3] - UP[1,4] * b[4]) / UP[1,1]
    return b
end

# ----------------------------------------------------------------------------
# Matrix solves with Symmetric input: M = L⁻¹ M (reads lower tri, writes full)
# ----------------------------------------------------------------------------

@propagate_inbounds function ldiv1!(L::LowerTriangular, M::Symmetric)
    LP = parent(L)
    MP = parent(M)
    @boundscheck checkbounds(LP, 1, 1)
    @boundscheck checkbounds(MP, 1, 1)
    @inbounds MP[1,1] = MP[1,1] / LP[1,1]
    return MP
end

@propagate_inbounds function ldiv2!(L::LowerTriangular, M::Symmetric)
    LP = parent(L)
    MP = parent(M)
    @boundscheck checkbounds(LP, 2, 2)
    @boundscheck checkbounds(MP, 2, 2)
    @inbounds L₁₁, L₂₁, L₂₂ = LP[1,1], LP[2,1], LP[2,2]
    @inbounds M₂₁ = MP[2,1]
    @inbounds MP[1,1] = MP[1,1] / L₁₁
    @inbounds MP[2,1] = (M₂₁ - L₂₁ * MP[1,1]) / L₂₂
    @inbounds MP[1,2] = M₂₁ / L₁₁
    @inbounds MP[2,2] = (MP[2,2] - L₂₁ * MP[1,2]) / L₂₂
    return MP
end

@propagate_inbounds function ldiv3!(L::LowerTriangular, M::Symmetric)
    LP = parent(L)
    MP = parent(M)
    @boundscheck checkbounds(LP, 3, 3)
    @boundscheck checkbounds(MP, 3, 3)
    @inbounds L₁₁, L₂₁, L₃₁ = LP[1,1], LP[2,1], LP[3,1]
    @inbounds L₂₂, L₃₂, L₃₃ = LP[2,2], LP[3,2], LP[3,3]
    @inbounds M₂₁, M₃₁, M₃₂ = MP[2,1], MP[3,1], MP[3,2]
    @inbounds MP[1,1] = MP[1,1] / L₁₁
    @inbounds MP[2,1] = (M₂₁ - L₂₁ * MP[1,1]) / L₂₂
    @inbounds MP[3,1] = (M₃₁ - L₃₁ * MP[1,1] - L₃₂ * MP[2,1]) / L₃₃
    @inbounds MP[1,2] = M₂₁ / L₁₁
    @inbounds MP[2,2] = (MP[2,2] - L₂₁ * MP[1,2]) / L₂₂
    @inbounds MP[3,2] = (M₃₂ - L₃₁ * MP[1,2] - L₃₂ * MP[2,2]) / L₃₃
    @inbounds MP[1,3] = M₃₁ / L₁₁
    @inbounds MP[2,3] = (M₃₂ - L₂₁ * MP[1,3]) / L₂₂
    @inbounds MP[3,3] = (MP[3,3] - L₃₁ * MP[1,3] - L₃₂ * MP[2,3]) / L₃₃
    return MP
end

@propagate_inbounds function ldiv4!(L::LowerTriangular, M::Symmetric)
    LP = parent(L)
    MP = parent(M)
    @boundscheck checkbounds(LP, 4, 4)
    @boundscheck checkbounds(MP, 4, 4)
    @inbounds L₁₁, L₂₁, L₃₁, L₄₁ = LP[1,1], LP[2,1], LP[3,1], LP[4,1]
    @inbounds L₂₂, L₃₂, L₄₂ = LP[2,2], LP[3,2], LP[4,2]
    @inbounds L₃₃, L₄₃ = LP[3,3], LP[4,3]
    @inbounds L₄₄ = LP[4,4]
    @inbounds M₂₁, M₃₁, M₄₁ = MP[2,1], MP[3,1], MP[4,1]
    @inbounds M₃₂, M₄₂ = MP[3,2], MP[4,2]
    @inbounds M₄₃ = MP[4,3]
    @inbounds MP[1,1] = MP[1,1] / L₁₁
    @inbounds MP[2,1] = (M₂₁ - L₂₁ * MP[1,1]) / L₂₂
    @inbounds MP[3,1] = (M₃₁ - L₃₁ * MP[1,1] - L₃₂ * MP[2,1]) / L₃₃
    @inbounds MP[4,1] = (M₄₁ - L₄₁ * MP[1,1] - L₄₂ * MP[2,1] - L₄₃ * MP[3,1]) / L₄₄
    @inbounds MP[1,2] = M₂₁ / L₁₁
    @inbounds MP[2,2] = (MP[2,2] - L₂₁ * MP[1,2]) / L₂₂
    @inbounds MP[3,2] = (M₃₂ - L₃₁ * MP[1,2] - L₃₂ * MP[2,2]) / L₃₃
    @inbounds MP[4,2] = (M₄₂ - L₄₁ * MP[1,2] - L₄₂ * MP[2,2] - L₄₃ * MP[3,2]) / L₄₄
    @inbounds MP[1,3] = M₃₁ / L₁₁
    @inbounds MP[2,3] = (M₃₂ - L₂₁ * MP[1,3]) / L₂₂
    @inbounds MP[3,3] = (MP[3,3] - L₃₁ * MP[1,3] - L₃₂ * MP[2,3]) / L₃₃
    @inbounds MP[4,3] = (M₄₃ - L₄₁ * MP[1,3] - L₄₂ * MP[2,3] - L₄₃ * MP[3,3]) / L₄₄
    @inbounds MP[1,4] = M₄₁ / L₁₁
    @inbounds MP[2,4] = (M₄₂ - L₂₁ * MP[1,4]) / L₂₂
    @inbounds MP[3,4] = (M₄₃ - L₃₁ * MP[1,4] - L₃₂ * MP[2,4]) / L₃₃
    @inbounds MP[4,4] = (MP[4,4] - L₄₁ * MP[1,4] - L₄₂ * MP[2,4] - L₄₃ * MP[3,4]) / L₄₄
    return MP
end

# ----------------------------------------------------------------------------
# Matrix solves with Symmetric input: M = U⁻¹ M (reads lower tri, writes full)
# ----------------------------------------------------------------------------

@propagate_inbounds function ldiv1!(U::UpperTriangular, M::Symmetric)
    UP = parent(U)
    MP = parent(M)
    @boundscheck checkbounds(UP, 1, 1)
    @boundscheck checkbounds(MP, 1, 1)
    @inbounds MP[1,1] = MP[1,1] / UP[1,1]
    return MP
end

@propagate_inbounds function ldiv2!(U::UpperTriangular, M::Symmetric)
    UP = parent(U)
    MP = parent(M)
    @boundscheck checkbounds(UP, 2, 2)
    @boundscheck checkbounds(MP, 2, 2)
    @inbounds U₁₁, U₁₂ = UP[1,1], UP[1,2]
    @inbounds U₂₂ = UP[2,2]
    @inbounds M₂₁ = MP[2,1]
    @inbounds MP[2,1] = MP[2,1] / U₂₂
    @inbounds MP[1,1] = (MP[1,1] - U₁₂ * MP[2,1]) / U₁₁
    @inbounds MP[2,2] = MP[2,2] / U₂₂
    @inbounds MP[1,2] = (M₂₁ - U₁₂ * MP[2,2]) / U₁₁
    return MP
end

@propagate_inbounds function ldiv3!(U::UpperTriangular, M::Symmetric)
    UP = parent(U)
    MP = parent(M)
    @boundscheck checkbounds(UP, 3, 3)
    @boundscheck checkbounds(MP, 3, 3)
    @inbounds U₁₁, U₁₂, U₁₃ = UP[1,1], UP[1,2], UP[1,3]
    @inbounds U₂₂, U₂₃ = UP[2,2], UP[2,3]
    @inbounds U₃₃ = UP[3,3]
    @inbounds M₂₁, M₃₁, M₃₂ = MP[2,1], MP[3,1], MP[3,2]
    @inbounds MP[3,1] = MP[3,1] / U₃₃
    @inbounds MP[2,1] = (M₂₁ - U₂₃ * MP[3,1]) / U₂₂
    @inbounds MP[1,1] = (MP[1,1] - U₁₂ * MP[2,1] - U₁₃ * MP[3,1]) / U₁₁
    @inbounds MP[3,2] = MP[3,2] / U₃₃
    @inbounds MP[2,2] = (MP[2,2] - U₂₃ * MP[3,2]) / U₂₂
    @inbounds MP[1,2] = (M₂₁ - U₁₂ * MP[2,2] - U₁₃ * MP[3,2]) / U₁₁
    @inbounds MP[3,3] = MP[3,3] / U₃₃
    @inbounds MP[2,3] = (M₃₂ - U₂₃ * MP[3,3]) / U₂₂
    @inbounds MP[1,3] = (M₃₁ - U₁₂ * MP[2,3] - U₁₃ * MP[3,3]) / U₁₁
    return MP
end

@propagate_inbounds function ldiv4!(U::UpperTriangular, M::Symmetric)
    UP = parent(U)
    MP = parent(M)
    @boundscheck checkbounds(UP, 4, 4)
    @boundscheck checkbounds(MP, 4, 4)
    @inbounds U₁₁, U₁₂, U₁₃, U₁₄ = UP[1,1], UP[1,2], UP[1,3], UP[1,4]
    @inbounds U₂₂, U₂₃, U₂₄ = UP[2,2], UP[2,3], UP[2,4]
    @inbounds U₃₃, U₃₄ = UP[3,3], UP[3,4]
    @inbounds U₄₄ = UP[4,4]
    @inbounds M₂₁, M₃₁, M₄₁ = MP[2,1], MP[3,1], MP[4,1]
    @inbounds M₃₂, M₄₂ = MP[3,2], MP[4,2]
    @inbounds M₄₃ = MP[4,3]
    @inbounds MP[4,1] = MP[4,1] / U₄₄
    @inbounds MP[3,1] = (M₃₁ - U₃₄ * MP[4,1]) / U₃₃
    @inbounds MP[2,1] = (M₂₁ - U₂₃ * MP[3,1] - U₂₄ * MP[4,1]) / U₂₂
    @inbounds MP[1,1] = (MP[1,1] - U₁₂ * MP[2,1] - U₁₃ * MP[3,1] - U₁₄ * MP[4,1]) / U₁₁
    @inbounds MP[4,2] = MP[4,2] / U₄₄
    @inbounds MP[3,2] = (M₃₂ - U₃₄ * MP[4,2]) / U₃₃
    @inbounds MP[2,2] = (MP[2,2] - U₂₃ * MP[3,2] - U₂₄ * MP[4,2]) / U₂₂
    @inbounds MP[1,2] = (M₂₁ - U₁₂ * MP[2,2] - U₁₃ * MP[3,2] - U₁₄ * MP[4,2]) / U₁₁
    @inbounds MP[4,3] = MP[4,3] / U₄₄
    @inbounds MP[3,3] = (MP[3,3] - U₃₄ * MP[4,3]) / U₃₃
    @inbounds MP[2,3] = (M₃₂ - U₂₃ * MP[3,3] - U₂₄ * MP[4,3]) / U₂₂
    @inbounds MP[1,3] = (M₃₁ - U₁₂ * MP[2,3] - U₁₃ * MP[3,3] - U₁₄ * MP[4,3]) / U₁₁
    @inbounds MP[4,4] = MP[4,4] / U₄₄
    @inbounds MP[3,4] = (M₄₃ - U₃₄ * MP[4,4]) / U₃₃
    @inbounds MP[2,4] = (M₄₂ - U₂₃ * MP[3,4] - U₂₄ * MP[4,4]) / U₂₂
    @inbounds MP[1,4] = (M₄₁ - U₁₂ * MP[2,4] - U₁₃ * MP[3,4] - U₁₄ * MP[4,4]) / U₁₁
    return MP
end

# ============================================================================
# rdiv - right triangular solves (inlined)
# ============================================================================

# ----------------------------------------------------------------------------
# Matrix right-solves: M = M * U⁻¹ (row-wise forward substitution)
# ----------------------------------------------------------------------------

@propagate_inbounds function rdiv1!(M::AbstractMatrix, U::UpperTriangular)
    UP = parent(U)
    @boundscheck checkbounds(UP, 1, 1)
    @boundscheck checkbounds(M, 1, 1)
    @inbounds M[1,1] = M[1,1] / UP[1,1]
    return M
end

@propagate_inbounds function rdiv2!(M::AbstractMatrix, U::UpperTriangular)
    UP = parent(U)
    @boundscheck checkbounds(UP, 2, 2)
    @boundscheck checkbounds(M, 2, 2)
    @inbounds U₁₁, U₁₂ = UP[1,1], UP[1,2]
    @inbounds U₂₂ = UP[2,2]
    @inbounds M[1,1] = M[1,1] / U₁₁
    @inbounds M[1,2] = (M[1,2] - U₁₂ * M[1,1]) / U₂₂
    @inbounds M[2,1] = M[2,1] / U₁₁
    @inbounds M[2,2] = (M[2,2] - U₁₂ * M[2,1]) / U₂₂
    return M
end

@propagate_inbounds function rdiv3!(M::AbstractMatrix, U::UpperTriangular)
    UP = parent(U)
    @boundscheck checkbounds(UP, 3, 3)
    @boundscheck checkbounds(M, 3, 3)
    @inbounds U₁₁, U₁₂, U₁₃ = UP[1,1], UP[1,2], UP[1,3]
    @inbounds U₂₂, U₂₃ = UP[2,2], UP[2,3]
    @inbounds U₃₃ = UP[3,3]
    @inbounds M[1,1] = M[1,1] / U₁₁
    @inbounds M[1,2] = (M[1,2] - U₁₂ * M[1,1]) / U₂₂
    @inbounds M[1,3] = (M[1,3] - U₁₃ * M[1,1] - U₂₃ * M[1,2]) / U₃₃
    @inbounds M[2,1] = M[2,1] / U₁₁
    @inbounds M[2,2] = (M[2,2] - U₁₂ * M[2,1]) / U₂₂
    @inbounds M[2,3] = (M[2,3] - U₁₃ * M[2,1] - U₂₃ * M[2,2]) / U₃₃
    @inbounds M[3,1] = M[3,1] / U₁₁
    @inbounds M[3,2] = (M[3,2] - U₁₂ * M[3,1]) / U₂₂
    @inbounds M[3,3] = (M[3,3] - U₁₃ * M[3,1] - U₂₃ * M[3,2]) / U₃₃
    return M
end

@propagate_inbounds function rdiv4!(M::AbstractMatrix, U::UpperTriangular)
    UP = parent(U)
    @boundscheck checkbounds(UP, 4, 4)
    @boundscheck checkbounds(M, 4, 4)
    @inbounds U₁₁, U₁₂, U₁₃, U₁₄ = UP[1,1], UP[1,2], UP[1,3], UP[1,4]
    @inbounds U₂₂, U₂₃, U₂₄ = UP[2,2], UP[2,3], UP[2,4]
    @inbounds U₃₃, U₃₄ = UP[3,3], UP[3,4]
    @inbounds U₄₄ = UP[4,4]
    @inbounds M[1,1] = M[1,1] / U₁₁
    @inbounds M[1,2] = (M[1,2] - U₁₂ * M[1,1]) / U₂₂
    @inbounds M[1,3] = (M[1,3] - U₁₃ * M[1,1] - U₂₃ * M[1,2]) / U₃₃
    @inbounds M[1,4] = (M[1,4] - U₁₄ * M[1,1] - U₂₄ * M[1,2] - U₃₄ * M[1,3]) / U₄₄
    @inbounds M[2,1] = M[2,1] / U₁₁
    @inbounds M[2,2] = (M[2,2] - U₁₂ * M[2,1]) / U₂₂
    @inbounds M[2,3] = (M[2,3] - U₁₃ * M[2,1] - U₂₃ * M[2,2]) / U₃₃
    @inbounds M[2,4] = (M[2,4] - U₁₄ * M[2,1] - U₂₄ * M[2,2] - U₃₄ * M[2,3]) / U₄₄
    @inbounds M[3,1] = M[3,1] / U₁₁
    @inbounds M[3,2] = (M[3,2] - U₁₂ * M[3,1]) / U₂₂
    @inbounds M[3,3] = (M[3,3] - U₁₃ * M[3,1] - U₂₃ * M[3,2]) / U₃₃
    @inbounds M[3,4] = (M[3,4] - U₁₄ * M[3,1] - U₂₄ * M[3,2] - U₃₄ * M[3,3]) / U₄₄
    @inbounds M[4,1] = M[4,1] / U₁₁
    @inbounds M[4,2] = (M[4,2] - U₁₂ * M[4,1]) / U₂₂
    @inbounds M[4,3] = (M[4,3] - U₁₃ * M[4,1] - U₂₃ * M[4,2]) / U₃₃
    @inbounds M[4,4] = (M[4,4] - U₁₄ * M[4,1] - U₂₄ * M[4,2] - U₃₄ * M[4,3]) / U₄₄
    return M
end

# ----------------------------------------------------------------------------
# Matrix right-solves: M = M * L⁻¹ (column-wise backward substitution)
# ----------------------------------------------------------------------------

@propagate_inbounds function rdiv1!(M::AbstractMatrix, L::LowerTriangular)
    LP = parent(L)
    @boundscheck checkbounds(LP, 1, 1)
    @boundscheck checkbounds(M, 1, 1)
    @inbounds M[1,1] = M[1,1] / LP[1,1]
    return M
end

@propagate_inbounds function rdiv2!(M::AbstractMatrix, L::LowerTriangular)
    LP = parent(L)
    @boundscheck checkbounds(LP, 2, 2)
    @boundscheck checkbounds(M, 2, 2)
    @inbounds L₁₁, L₂₁, L₂₂ = LP[1,1], LP[2,1], LP[2,2]
    @inbounds M[1,2] = M[1,2] / L₂₂
    @inbounds M[2,2] = M[2,2] / L₂₂
    @inbounds M[1,1] = (M[1,1] - M[1,2] * L₂₁) / L₁₁
    @inbounds M[2,1] = (M[2,1] - M[2,2] * L₂₁) / L₁₁
    return M
end

@propagate_inbounds function rdiv3!(M::AbstractMatrix, L::LowerTriangular)
    LP = parent(L)
    @boundscheck checkbounds(LP, 3, 3)
    @boundscheck checkbounds(M, 3, 3)
    @inbounds L₁₁, L₂₁, L₃₁ = LP[1,1], LP[2,1], LP[3,1]
    @inbounds L₂₂, L₃₂ = LP[2,2], LP[3,2]
    @inbounds L₃₃ = LP[3,3]
    @inbounds M[1,3] = M[1,3] / L₃₃
    @inbounds M[2,3] = M[2,3] / L₃₃
    @inbounds M[3,3] = M[3,3] / L₃₃
    @inbounds M[1,2] = (M[1,2] - M[1,3] * L₃₂) / L₂₂
    @inbounds M[2,2] = (M[2,2] - M[2,3] * L₃₂) / L₂₂
    @inbounds M[3,2] = (M[3,2] - M[3,3] * L₃₂) / L₂₂
    @inbounds M[1,1] = (M[1,1] - M[1,2] * L₂₁ - M[1,3] * L₃₁) / L₁₁
    @inbounds M[2,1] = (M[2,1] - M[2,2] * L₂₁ - M[2,3] * L₃₁) / L₁₁
    @inbounds M[3,1] = (M[3,1] - M[3,2] * L₂₁ - M[3,3] * L₃₁) / L₁₁
    return M
end

@propagate_inbounds function rdiv4!(M::AbstractMatrix, L::LowerTriangular)
    LP = parent(L)
    @boundscheck checkbounds(LP, 4, 4)
    @boundscheck checkbounds(M, 4, 4)
    @inbounds L₁₁, L₂₁, L₃₁, L₄₁ = LP[1,1], LP[2,1], LP[3,1], LP[4,1]
    @inbounds L₂₂, L₃₂, L₄₂ = LP[2,2], LP[3,2], LP[4,2]
    @inbounds L₃₃, L₄₃ = LP[3,3], LP[4,3]
    @inbounds L₄₄ = LP[4,4]
    @inbounds M[1,4] = M[1,4] / L₄₄
    @inbounds M[2,4] = M[2,4] / L₄₄
    @inbounds M[3,4] = M[3,4] / L₄₄
    @inbounds M[4,4] = M[4,4] / L₄₄
    @inbounds M[1,3] = (M[1,3] - M[1,4] * L₄₃) / L₃₃
    @inbounds M[2,3] = (M[2,3] - M[2,4] * L₄₃) / L₃₃
    @inbounds M[3,3] = (M[3,3] - M[3,4] * L₄₃) / L₃₃
    @inbounds M[4,3] = (M[4,3] - M[4,4] * L₄₃) / L₃₃
    @inbounds M[1,2] = (M[1,2] - M[1,3] * L₃₂ - M[1,4] * L₄₂) / L₂₂
    @inbounds M[2,2] = (M[2,2] - M[2,3] * L₃₂ - M[2,4] * L₄₂) / L₂₂
    @inbounds M[3,2] = (M[3,2] - M[3,3] * L₃₂ - M[3,4] * L₄₂) / L₂₂
    @inbounds M[4,2] = (M[4,2] - M[4,3] * L₃₂ - M[4,4] * L₄₂) / L₂₂
    @inbounds M[1,1] = (M[1,1] - M[1,2] * L₂₁ - M[1,3] * L₃₁ - M[1,4] * L₄₁) / L₁₁
    @inbounds M[2,1] = (M[2,1] - M[2,2] * L₂₁ - M[2,3] * L₃₁ - M[2,4] * L₄₁) / L₁₁
    @inbounds M[3,1] = (M[3,1] - M[3,2] * L₂₁ - M[3,3] * L₃₁ - M[3,4] * L₄₁) / L₁₁
    @inbounds M[4,1] = (M[4,1] - M[4,2] * L₂₁ - M[4,3] * L₃₁ - M[4,4] * L₄₁) / L₁₁
    return M
end

# ============================================================================
# Newton solvers
# ============================================================================

# Safeguarded Newton on a scalar root of f in [lo, hi],
# where f(lo) < 0 < f(hi) (increasing across bracket).
function rtsafe(jet, lo::T, hi::T, r0::T, fr0::T, fpr0::T,
                flo::T, fhi::T;
                tol::T = T(1e-12), maxit::Int = 60) where {T}
    @assert isfinite(flo) && isfinite(fhi)
    @assert flo <= 0 && fhi >= 0

    if iszero(flo)
        r = lo
    elseif iszero(fhi)
        r = hi
    else
        r = clamp(r0, lo, hi)
        wlo = flo
        whi = fhi
        side = 0

        for i in 1:maxit
            rp = r

            if i == 1 && r == r0
                fr, fpr = fr0, fpr0
            else
                fr, fpr = jet(r)
            end

            if fr >= 0
                hi, whi = r, fr

                if side == +1
                    wlo /= 2
                end

                side = +1
            elseif fr < 0
                lo, wlo = r, fr

                if side == -1
                    whi /= 2
                end

                side = -1
            else
                break
            end

            r = r - fr / fpr

            if !(lo < r < hi)
                r = lo + wlo * (hi - lo) / (wlo - whi)

                if !(lo < r < hi)
                    r = (lo + hi) / 2
                end
            end

            if abs(r - rp) < tol * (one(T) + abs(r))
                break
            end
        end
    end

    return r
end

# Newton–Fourier line search for the last τ in [0, hi] with a
# *certified* margin g(τ) ≥ noise(τ); g concave (g′ non-increasing).
# jet(τ) returns (g, gp, noise) with noise the absolute uncertainty
# of g at τ (0 if unknown).
function nflast(jet, hi::T; tol::T = eps(T), maxit::Int = 53) where {T}
    lo = zero(T)
    glo, gplo, nlo = jet(lo)
    #
    # τ-resolution below which g's sign is noise
    #
    τres(noise, gp) = (isfinite(gp) && !iszero(gp)) ? noise / abs(gp) : zero(T)
    #
    # the starting point must be inside with certifiable margin;
    # a start inside the noise band is treated as on the boundary.
    #
    if glo >= nlo && glo >= 0
        ghi, gphi, nhi = jet(hi)

        if ghi >= nhi
            lo = hi
        else
            wlo = glo
            whi = τres(nhi, gphi)

            for _ in 1:maxit
                lon = lo + wlo * (hi - lo) / (wlo - ghi)

                if !(lo < lon < hi)
                    lon = lo + (hi - lo) / 2
                end

                g, gp, n = jet(lon)

                if g >= n
                    lo, glo, gplo, wlo = lon, g, gp, g
                else
                    hi, ghi, gphi, whi = lon, g, gp, τres(n, gp)
                    wlo /= 2
                end

                if hi - lo < max(whi, tol * (one(T) + hi))
                    break
                end

                hin = hi - ghi / gphi

                if !(lo < hin < hi)
                    hin = hi
                end

                tlo = lo - glo / gplo

                if lo < tlo < hin
                    hin = tlo
                end

                if !(lo < hin < hi)
                    hin = lo + (hi - lo) / 2
                end

                g, gp, n = jet(hin)

                if g >= n
                    lo, glo, gplo, wlo = hin, g, gp, g
                else
                    hi, ghi, gphi, whi = hin, g, gp, τres(n, gp)
                end

                if hi - lo < max(whi, tol * (one(T) + hi))
                    break
                end
            end
        end
    end

    return lo
end
