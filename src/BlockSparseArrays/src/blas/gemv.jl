@inline function gemv_blk_unit!(C::AbstractVector, tA::Val{TA}, A::BlockSparseMatrix{<:Any, I}, B::AbstractVector, α::Number, bstrt::I, i::I, j::I) where {TA, I}
    @inbounds Ae = A.val[bstrt]

    if TA === :N
        @inbounds C[i] = axpy(α *      Ae,  B[j], C[i])
    elseif TA === :T
        @inbounds C[j] = axpy(α *      Ae,  B[i], C[j])
    else
        @inbounds C[j] = axpy(α * conj(Ae), B[i], C[j])
    end

    return
end

@inline function gemv_blk_small!(C::AbstractVector, tA::Val{TA}, A::BlockSparseMatrix{<:Any, I}, B::AbstractVector, α::Number, bstrt::I, istrt::I, istop::I, jstrt::I, jstop::I) where {TA, I}
    T = promote_eltype(A, B)
    b = bstrt

    if TA === :N
        @inbounds for j in jstrt:jstop
            αBj = α * B[j]

            for i in istrt:istop
                C[i] = axpy(A.val[b], αBj, C[i])
                b += one(I)
            end
        end
    elseif TA === :T
        @inbounds for j in jstrt:jstop
            Δj = zero(T)

            for i in istrt:istop
                Δj = axpy(A.val[b], B[i], Δj)
                b += one(I)
            end

            C[j] = axpy(α, Δj, C[j])
        end
    else
        @inbounds for j in jstrt:jstop
            Δj = zero(T)

            for i in istrt:istop
                Δj = axpy(conj(A.val[b]), B[i], Δj)
                b += one(I)
            end

            C[j] = axpy(α, Δj, C[j])
        end
    end

    return
end

@inline function gemv_blk_large!(C::AbstractVector, tA::Val{TA}, A::BlockSparseMatrix{<:Any, I}, B::AbstractVector, α::Number, bstrt::I, nrow::I, ncol::I, istrt::I, istop::I, jstrt::I, jstop::I) where {TA, I}
    Ae = @inbounds block_impl(A, nrow, ncol, bstrt)

    if TA === :N
        U = @inbounds view(C, istrt:istop)
        V = @inbounds view(B, jstrt:jstop)
        mul!(U,           Ae,  V, α, true)
    elseif TA === :T
        U = @inbounds view(B, istrt:istop)
        V = @inbounds view(C, jstrt:jstop)
        mul!(V, transpose(Ae), U, α, true)
    else
        U = @inbounds view(B, istrt:istop)
        V = @inbounds view(C, jstrt:jstop)
        mul!(V,   adjoint(Ae), U, α, true)
    end

    return
end

function gemv_impl!(C::AbstractVector, tA::Val{TA}, A::BlockSparseMatrix{<:Any, I}, B::AbstractVector, α::Number, β::Number) where {TA, I}
    if iszero(β)
        fill!(C, β)
    elseif !isone(β)
        rmul!(C, β)
    end

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
            istop = A.xrow[u + one(I)] - one(I)

            nrow = istop - istrt + one(I)

            if isone(ncol) && isone(nrow)
                gemv_blk_unit!(C, tA, A, B, α, bstrt, istrt, jstrt)
            elseif ncol < 24 && nrow < 24
                gemv_blk_small!(C, tA, A, B, α, bstrt, istrt, istop, jstrt, jstop)
            else
                gemv_blk_large!(C, tA, A, B, α, bstrt, nrow, ncol, istrt, istop, jstrt, jstop)
            end

            bstrt += ncol * nrow
        end

        jstrt = jstop + one(I)
        estrt = estop + one(I)
    end

    return C
end
