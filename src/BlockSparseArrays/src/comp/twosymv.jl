const SYMV_MAXN = 8

############################################################################################
# twosymv_diag_one! / twosymv_offd_one!
############################################################################################

@inline function twosymv_diag_one!(s::AbstractVector, cs::AbstractVector, A::AbstractVector, b::AbstractVector, bstrt::I, i::I) where {I}
    @inbounds twoaxpy!(s, cs, i, A[bstrt], b[i])
    return
end

@inline function twosymv_offd_one!(s::AbstractVector, cs::AbstractVector, A::AbstractVector, b::AbstractVector, bstrt::I, i::I, j::I) where {I}
    @inbounds Ae = A[bstrt]
    @inbounds twoaxpy!(s, cs, j, Ae, b[i])
    @inbounds twoaxpy!(s, cs, i, Ae, b[j])
    return
end

############################################################################################
# twosymv_diag_small! / twosymv_offd_small!
############################################################################################

@inline function twosymv_diag_small!(s::AbstractVector, cs::AbstractVector, A::AbstractVector, b::AbstractVector, bstrt::I, istrt::I, jstrt::I, nrow::I, ncol::I, ::Val{UPLO}) where {UPLO, I}
    T = promote_eltype(A, b)
    istop = istrt + nrow - one(I)
    jstop = jstrt + ncol - one(I)

    @inbounds for j in jstrt:jstop
        jloc = j - jstrt
        bj = b[j]
        Δ = zero(T)
        δ = zero(T)

        if UPLO === :L
            p = bstrt + jloc * (one(I) + nrow)
            irng = j + one(I):istop
        else
            p = bstrt + jloc * nrow
            irng = istrt:j - one(I)
        end

        if UPLO === :L
            Δ, δ = twoacc(Δ, δ, A[p], bj)
            p += one(I)
        end

        for i in irng
            Aij = A[p]
            p += one(I)
            Δ, δ = twoacc(Δ, δ, Aij, b[i])
            twoaxpy!(s, cs, i, Aij, bj)
        end

        if UPLO === :U
            Δ, δ = twoacc(Δ, δ, A[p], bj)
            p += one(I)
        end

        twofold!(s, cs, j, Δ, δ)
    end

    return
end

@inline function twosymv_offd_small!(s::AbstractVector, cs::AbstractVector, A::AbstractVector, b::AbstractVector, bstrt::I, istrt::I, jstrt::I, nrow::I, ncol::I) where {I}
    T = promote_eltype(A, b)
    istop = istrt + nrow - one(I)
    jstop = jstrt + ncol - one(I)
    p = bstrt

    @inbounds for j in jstrt:jstop
        bj = b[j]
        Δ = zero(T)
        δ = zero(T)

        for i in istrt:istop
            Aij = A[p]
            p += one(I)
            Δ, δ = twoacc(Δ, δ, Aij, b[i])
            twoaxpy!(s, cs, i, Aij, bj)
        end

        twofold!(s, cs, j, Δ, δ)
    end

    return
end

############################################################################################
# twosymv_offd_large! / twosymv_diag_large!
############################################################################################

@inline function twosymv_offd_large!(s::AbstractVector, cs::AbstractVector, A::AbstractVector, b::AbstractVector, bstrt::I, istrt::I, jstrt::I, nrow::I, ncol::I) where {I}
    twogemv_blk!(s, cs, Val(:N), A, b, bstrt, istrt, jstrt, nrow, ncol, nrow)
    twogemv_blk!(s, cs, Val(:T), A, b, bstrt, istrt, jstrt, nrow, ncol, nrow)
    return
end

@inline function twosymv_diag_large!(s::AbstractVector, cs::AbstractVector, A::AbstractVector, b::AbstractVector, bstrt::I, istrt::I, jstrt::I, nrow::I, ncol::I, ::Val{UPLO}) where {UPLO, I}
    T = promote_eltype(A, b)
    istop = istrt + nrow - one(I)
    jstop = jstrt + ncol - one(I)
    j = jstrt

    @inbounds while j + three(I) <= jstop
        cstop = j + three(I)

        for l in j:cstop
            b0 = bstrt + (l - jstrt) * nrow - istrt
            bl = b[l]
            Δ = zero(T)
            δ = zero(T)

            if UPLO === :L
                Δ, δ = twoacc(Δ, δ, A[b0 + l], bl)
                irng = l + one(I):cstop
            else
                irng = j:l - one(I)
            end

            for i in irng
                Ail = A[b0 + i]
                Δ, δ = twoacc(Δ, δ, Ail, b[i])
                twoaxpy!(s, cs, i, Ail, bl)
            end

            if UPLO === :U
                Δ, δ = twoacc(Δ, δ, A[b0 + l], bl)
            end

            twofold!(s, cs, l, Δ, δ)
        end

        if UPLO === :L
            pstrt = cstop + one(I)
            pstop = istop
        else
            pstrt = istrt
            pstop = j - one(I)
        end

        if pstrt <= pstop
            pnrow = pstop - pstrt + one(I)
            bp = bstrt + (j - jstrt) * nrow + (pstrt - istrt)
            twogemv_blk!(s, cs, Val(:N), A, b, bp, pstrt, j, pnrow, four(I), nrow)
            twogemv_blk!(s, cs, Val(:T), A, b, bp, pstrt, j, pnrow, four(I), nrow)
        end

        j += four(I)
    end

    @inbounds while j <= jstop
        b0 = bstrt + (j - jstrt) * nrow - istrt
        bj = b[j]
        Δ = zero(T)
        δ = zero(T)

        if UPLO === :L
            Δ, δ = twoacc(Δ, δ, A[b0 + j], bj)
            irng = j + one(I):istop
        else
            irng = istrt:j - one(I)
        end

        for i in irng
            Aij = A[b0 + i]
            Δ, δ = twoacc(Δ, δ, Aij, b[i])
            twoaxpy!(s, cs, i, Aij, bj)
        end

        if UPLO === :U
            Δ, δ = twoacc(Δ, δ, A[b0 + j], bj)
        end

        twofold!(s, cs, j, Δ, δ)
        j += one(I)
    end

    return
end

############################################################################################
# twosymv_impl!
############################################################################################

function twosymv_impl!(tA::Val, uplo::Symbol, A::BlockSparseMatrix{<:Any, I}, b::AbstractVector, s::AbstractVector, cs::AbstractVector) where {I}
    if uplo === :L
        return twosymv_impl!(tA, Val(:L), A, b, s, cs)
    else
        return twosymv_impl!(tA, Val(:U), A, b, s, cs)
    end
end

function twosymv_impl!(tA::Val, uplo::Val{UPLO}, A::BlockSparseMatrix{<:Any, I}, b::AbstractVector, s::AbstractVector, cs::AbstractVector) where {UPLO, I}
    fill!(s,  false)
    fill!(cs, false)

    jstrt = one(I)
    estrt = one(I)

    @inbounds for v in vtxs(A)
        jstop = A.xcol[v + one(I)] - one(I)
        estop = A.xsrc[v + one(I)] - one(I)

        ncol = jstop - jstrt + one(I)

        for e in estrt:estop
            u = A.tgt[e]

            if intriangle(u, v, UPLO)
                bstrt = A.xblk[e]
                istrt = A.xrow[u]
                nrow = A.xrow[u + one(I)] - istrt

                if u == v
                    if isone(nrow) && isone(ncol)
                        twosymv_diag_one!(s, cs, A.val, b, bstrt, istrt)
                    elseif nrow <= SYMV_MAXN
                        twosymv_diag_small!(s, cs, A.val, b, bstrt, istrt, jstrt, nrow, ncol, uplo)
                    else
                        twosymv_diag_large!(s, cs, A.val, b, bstrt, istrt, jstrt, nrow, ncol, uplo)
                    end
                else
                    if isone(nrow) && isone(ncol)
                        twosymv_offd_one!(s, cs, A.val, b, bstrt, istrt, jstrt)
                    elseif nrow <= SYMV_MAXN && ncol <= SYMV_MAXN
                        twosymv_offd_small!(s, cs, A.val, b, bstrt, istrt, jstrt, nrow, ncol)
                    else
                        twosymv_offd_large!(s, cs, A.val, b, bstrt, istrt, jstrt, nrow, ncol)
                    end
                end
            end
        end

        jstrt = jstop + one(I)
        estrt = estop + one(I)
    end

    return s, cs
end
