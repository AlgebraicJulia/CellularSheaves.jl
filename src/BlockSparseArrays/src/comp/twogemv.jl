############################################################################################
# twogemv_one!
############################################################################################

@inline function twogemv_one!(s::AbstractVector, cs::AbstractVector, tA::Val{TA}, A::AbstractVector, b::AbstractVector, bstrt::I, i::I, j::I) where {TA, I}
    @inbounds Ae = A[bstrt]

    if TA === :N
        @inbounds twoaxpy!(s, cs, i, Ae, b[j])
    else
        @inbounds twoaxpy!(s, cs, j, Ae, b[i])
    end

    return
end

############################################################################################
# twogemv_blk!
############################################################################################

# lda is the column stride of the (sub)block in A — equal to nrow for a whole block, or the
# parent block's row count when this is a strided sub-panel (used by the symv diagonal).

@inline function twogemv_blk!(s::AbstractVector, cs::AbstractVector, ::Val{:N}, A::AbstractVector, b::AbstractVector, bstrt::I, istrt::I, jstrt::I, nrow::I, ncol::I, lda::I) where {I}
    istop = istrt + nrow - one(I)
    jstop = jstrt + ncol - one(I)

    @inbounds for j in jstrt:jstop
        bj = b[j]
        b0 = bstrt + (j - jstrt) * lda - istrt

        @simd for i in istrt:istop
            twoaxpy!(s, cs, i, A[b0 + i], bj)
        end
    end

    return
end

@inline function twogemv_blk!(s::AbstractVector, cs::AbstractVector, ::Val{TA}, A::AbstractVector, b::AbstractVector, bstrt::I, istrt::I, jstrt::I, nrow::I, ncol::I, lda::I) where {TA, I}
    T = promote_eltype(A, b)
    V = Vec{4, T}

    istop = istrt + nrow - one(I)
    jstop = jstrt + ncol - one(I)
    j = jstrt

    @inbounds while j + three(I) <= jstop
        Δ = zero(V)
        δ = zero(V)
        b0 = bstrt + (j - jstrt) * lda - istrt

        for i in istrt:istop
            bi = V(b[i])
            a  = V((A[b0 + i], A[b0 + i + lda], A[b0 + i + two(I) * lda], A[b0 + i + three(I) * lda]))
            Δ, δ = twoacc(Δ, δ, a, bi)
        end

        twofold!(s, cs, j,            Δ[1], δ[1])
        twofold!(s, cs, j + one(I),   Δ[2], δ[2])
        twofold!(s, cs, j + two(I),   Δ[3], δ[3])
        twofold!(s, cs, j + three(I), Δ[4], δ[4])

        j += four(I)
    end

    @inbounds while j <= jstop
        Δ = zero(T)
        δ = zero(T)
        b0 = bstrt + (j - jstrt) * lda - istrt

        for i in istrt:istop
            Δ, δ = twoacc(Δ, δ, A[b0 + i], b[i])
        end

        twofold!(s, cs, j, Δ, δ)
        j += one(I)
    end

    return
end

############################################################################################
# twogemv_impl!
############################################################################################

function twogemv_impl!(tA::Val{TA}, A::BlockSparseMatrix{<:Any, I}, b::AbstractVector, s::AbstractVector, cs::AbstractVector) where {TA, I}
    fill!(s,  false)
    fill!(cs, false)

    jstrt = one(I)
    estrt = one(I)
    bstrt = one(I)

    @inbounds for v in vtxs(A)
        jstop = A.xcol[v + one(I)] - one(I)
        estop = A.xsrc[v + one(I)] - one(I)

        ncol = jstop - jstrt + one(I)

        for e in estrt:estop
            u = A.tgt[e]

            istrt = A.xrow[u]
            nrow = A.xrow[u + one(I)] - istrt

            if isone(nrow) && isone(ncol)
                twogemv_one!(s, cs, tA, A.val, b, bstrt, istrt, jstrt)
            else
                twogemv_blk!(s, cs, tA, A.val, b, bstrt, istrt, jstrt, nrow, ncol, nrow)
            end

            bstrt += ncol * nrow
        end

        jstrt = jstop + one(I)
        estrt = estop + one(I)
    end

    return s, cs
end
