#
# Small linear algebra helpers for IPM cone operations
#
# Organization:
#   1. 3×3 vector/matrix helpers (cross, dot, norm, copy, axpy, lmul, etc.)
#   2. mul - matrix multiplication (dense loops, triangular inlined)
#   3. chol - Cholesky factorization (inlined)
#   4. ldiv - left triangular solves (inlined)
#   5. rdiv - right triangular solves (inlined)
#

# ============================================================================
# 3×3 vector/matrix helpers
# ============================================================================

function cross3!(z::AbstractVector, x::AbstractVector, y::AbstractVector)
    z[1] = muladd(x[2], y[3], -x[3] * y[2])
    z[2] = muladd(x[3], y[1], -x[1] * y[3])
    z[3] = muladd(x[1], y[2], -x[2] * y[1])
    return z
end

function dot3(x::AbstractVector, y::AbstractVector)
    return muladd(x[1], y[1], muladd(x[2], y[2], x[3] * y[3]))
end

function norm3(x::AbstractVector)
    return sqrt(muladd(x[1], x[1], muladd(x[2], x[2], x[3]^2)))
end

function normsq(x::AbstractVector{T}) where {T}
    out = zero(T)

    @inbounds for i in eachindex(x)
        out = muladd(x[i], x[i], out)
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
        y[i] = muladd(a, x[i], y[i])
    end

    return y
end

@propagate_inbounds function axpby3!(a, x::AbstractVector, b, y::AbstractVector)
    @boundscheck checkbounds(x, 3)
    @boundscheck checkbounds(y, 3)

    @inbounds for i in 1:3
        y[i] = muladd(a, x[i], b * y[i])
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
                            s = muladd(A[i,k], B[k,j], s)
                        end

                        C[i,j] = α * s
                    end
                end
            else
                @inbounds for i in 1:$n
                    for j in 1:$n
                        s = zero(T)

                        for k in 1:$n
                            s = muladd(A[i,k], B[k,j], s)
                        end

                        C[i,j] = muladd(β, C[i,j], α * s)
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
                s = muladd(A[i,k], B[k], s)
            end

            C[i] = α * s
        end
    else
        @inbounds for i in 1:3
            s = zero(T)

            for k in 1:3
                s = muladd(A[i,k], B[k], s)
            end

            C[i] = muladd(β, C[i], α * s)
        end
    end

    return C
end

# ----------------------------------------------------------------------------
# Matrix × Diagonal: C = α*A*D + β*C (generated)
# ----------------------------------------------------------------------------

for n in (2, 3, 4)
    fname = Symbol(:mul, n, :!)

    # αdⱼ = α·d[j]; C[i,j] = A[i,j]·αdⱼ (+ β·C[i,j])
    scale = [:($(Symbol(:αd, j)) = α * d[$j]) for j in 1:n]
    zbody = [:(C[$i,$j] = A[$i,$j] * $(Symbol(:αd, j))) for j in 1:n for i in 1:n]
    ebody = [:(C[$i,$j] = muladd(β, C[$i,$j], A[$i,$j] * $(Symbol(:αd, j)))) for j in 1:n for i in 1:n]

    @eval @propagate_inbounds function $fname(C::AbstractMatrix, A::AbstractMatrix, D::Diagonal, α::Number, β::Number)
        d = D.diag

        @boundscheck checkbounds(C, $n, $n)
        @boundscheck checkbounds(A, $n, $n)
        @boundscheck checkbounds(d, $n)

        @inbounds begin
            $(scale...)

            if iszero(β)
                $(zbody...)
            else
                $(ebody...)
            end
        end

        return C
    end
end

# ----------------------------------------------------------------------------
# Upper × Lower: C = U * L (inlined - sparsity pattern)
# ----------------------------------------------------------------------------

# C[i,j] = Σ_{p ≥ max(i,j)} U[i,p]·L[p,j]
for n in (2, 3, 4)
    body = Expr[]

    for i in 1:n, j in 1:n
        acc = :(UP[$i,$n] * LP[$n,$j])

        for p in n-1:-1:max(i,j)
            acc = :(muladd(UP[$i,$p], LP[$p,$j], $acc))
        end

        push!(body, :(C[$i,$j] = $acc))
    end

    @eval @propagate_inbounds function $(Symbol(:mul, n, :!))(C::AbstractMatrix, U::UpperTriangular, L::LowerTriangular)
        UP = parent(U)
        LP = parent(L)

        @boundscheck checkbounds(C, $n, $n)
        @boundscheck checkbounds(UP, $n, $n)
        @boundscheck checkbounds(LP, $n, $n)

        @inbounds begin
            $(body...)
        end

        return C
    end
end

# ----------------------------------------------------------------------------
# Rank-1 update: M = α*x*y' + β*M (3×3)
# ----------------------------------------------------------------------------

@propagate_inbounds function ger3!(M::AbstractMatrix, x::AbstractVector, y::AbstractVector, α::Number, β::Number)
    @boundscheck checkbounds(M, 3, 3)
    @boundscheck checkbounds(x, 3)
    @boundscheck checkbounds(y, 3)

    @inbounds x₁, x₂, x₃ = α * x[1], α * x[2], α * x[3]
    @inbounds y₁, y₂, y₃ = y[1], y[2], y[3]

    if iszero(β)
        @inbounds M[1,1] = x₁ * y₁
        @inbounds M[2,1] = x₂ * y₁
        @inbounds M[3,1] = x₃ * y₁
        @inbounds M[1,2] = x₁ * y₂
        @inbounds M[2,2] = x₂ * y₂
        @inbounds M[3,2] = x₃ * y₂
        @inbounds M[1,3] = x₁ * y₃
        @inbounds M[2,3] = x₂ * y₃
        @inbounds M[3,3] = x₃ * y₃
    else
        @inbounds M[1,1] = muladd(β, M[1,1], x₁ * y₁)
        @inbounds M[2,1] = muladd(β, M[2,1], x₂ * y₁)
        @inbounds M[3,1] = muladd(β, M[3,1], x₃ * y₁)
        @inbounds M[1,2] = muladd(β, M[1,2], x₁ * y₂)
        @inbounds M[2,2] = muladd(β, M[2,2], x₂ * y₂)
        @inbounds M[3,2] = muladd(β, M[3,2], x₃ * y₂)
        @inbounds M[1,3] = muladd(β, M[1,3], x₁ * y₃)
        @inbounds M[2,3] = muladd(β, M[2,3], x₂ * y₃)
        @inbounds M[3,3] = muladd(β, M[3,3], x₃ * y₃)
    end

    return M
end

# ============================================================================
# chol - Cholesky factorization (inlined, in-place, returns success)
# ============================================================================

# In-place LLᵀ (writes lower tri of A):
#   L[j,j] = √(A[j,j] - Σ_{p<j} L[j,p]²)
#   L[i,j] = (A[i,j] - Σ_{p<j} L[i,p]·L[j,p]) / L[j,j]
for n in (2, 3, 4)
    body = Expr[]

    for j in 1:n
        djj = :(A[$j,$j])

        for p in 1:j-1
            djj = :(muladd(-A[$j,$p], A[$j,$p], $djj))
        end

        push!(body, :(djj = $djj))
        push!(body, :(djj > zero(T) || return false))
        push!(body, :(ljj = sqrt(djj)))
        push!(body, :(A[$j,$j] = ljj))

        for i in j+1:n
            num = :(A[$i,$j])

            for p in 1:j-1
                num = :(muladd(-A[$i,$p], A[$j,$p], $num))
            end

            push!(body, :(A[$i,$j] = ($num) / ljj))
        end
    end

    @eval @propagate_inbounds function $(Symbol(:chol, n, :!))(M::Symmetric{T}) where T
        A = parent(M)

        @boundscheck checkbounds(A, $n, $n)

        @inbounds begin
            $(body...)
        end

        return true
    end
end

# ============================================================================
# ldiv - left triangular solves (inlined)
# ============================================================================

# M = R⁻¹ M, reading the lower triangle of symmetric M (result written full).
# R lower ⇒ forward substitution; R upper ⇒ backward.
#   X[r,c] = (M[r,c] - Σ_p R[r,p]·X[p,c]) / R[r,r]
for n in (1, 2, 3, 4), (Rtype, lower) in ((:LowerTriangular, true), (:UpperTriangular, false))
    body = Expr[]

    # save the strict-lower originals of M (overwritten as we go)
    for j in 1:n, i in j+1:n
        push!(body, :($(Symbol(:msv, i, :_, j)) = MP[$i,$j]))
    end

    if lower
        rows = 1:n
    else
        rows = n:-1:1
    end

    for c in 1:n, r in rows
        # M[r,c] by symmetry: diagonal, saved strict-lower, or its mirror
        if r == c
            num = :(MP[$r,$r])
        elseif r > c
            num = Symbol(:msv, r, :_, c)
        else
            num = Symbol(:msv, c, :_, r)
        end

        if lower
            prange = 1:r-1
        else
            prange = r+1:n
        end

        for p in prange
            num = :(muladd(-RP[$r,$p], MP[$p,$c], $num))
        end

        push!(body, :(MP[$r,$c] = ($num) / RP[$r,$r]))
    end

    @eval @propagate_inbounds function $(Symbol(:ldiv, n, :!))(R::$Rtype, M::Symmetric)
        RP = parent(R)
        MP = parent(M)

        @boundscheck checkbounds(RP, $n, $n)
        @boundscheck checkbounds(MP, $n, $n)

        @inbounds begin
            $(body...)
        end

        return MP
    end
end

# ============================================================================
# rdiv - right triangular solves (inlined)
# ============================================================================

# M = M R⁻¹ (M a plain matrix, solved in place).
# R upper ⇒ forward over columns; R lower ⇒ backward.
#   X[i,c] = (M[i,c] - Σ_p X[i,p]·R[p,c]) / R[c,c]
for n in (1, 2, 3, 4), (Rtype, upper) in ((:UpperTriangular, true), (:LowerTriangular, false))
    body = Expr[]

    if upper
        cols = 1:n
    else
        cols = n:-1:1
    end

    for c in cols, i in 1:n
        num = :(M[$i,$c])

        if upper
            prange = 1:c-1
        else
            prange = c+1:n
        end

        for p in prange
            num = :(muladd(-M[$i,$p], RP[$p,$c], $num))
        end

        push!(body, :(M[$i,$c] = ($num) / RP[$c,$c]))
    end

    @eval @propagate_inbounds function $(Symbol(:rdiv, n, :!))(M::AbstractMatrix, R::$Rtype)
        RP = parent(R)

        @boundscheck checkbounds(RP, $n, $n)
        @boundscheck checkbounds(M, $n, $n)

        @inbounds begin
            $(body...)
        end

        return M
    end
end
