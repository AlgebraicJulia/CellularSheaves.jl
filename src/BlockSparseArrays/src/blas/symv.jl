@inline function symv_blk_unit_diag!(C::AbstractVector, tA::Val{TA}, A::BlockSparseMatrix{<:Any, I}, B::AbstractVector, α::Number, bstrt::I, i::I) where {TA, I}
    @inbounds Ae = A.val[bstrt]
    @inbounds C[i] = axpy(α * Ae, B[i], C[i])
    return
end

@inline function symv_blk_unit_offd!(C::AbstractVector, tA::Val{TA}, A::BlockSparseMatrix{<:Any, I}, B::AbstractVector, α::Number, bstrt::I, i::I, j::I) where {TA, I}
    @inbounds Aij = A.val[bstrt]

    if TA === :N
        Aji =      Aij
    else
        Aji = conj(Aij)
    end

    @inbounds C[j] = axpy(α * Aji, B[i], C[j])
    @inbounds C[i] = axpy(α * Aij, B[j], C[i])
    return
end

@inline function symv_blk_small_diag!(C::AbstractVector, tA::Val{TA}, uplo::Symbol, A::BlockSparseMatrix{<:Any, I}, B::AbstractVector, α::Number, bstrt::I, istrt::I, istop::I, jstrt::I, jstop::I) where {TA, I}
    T = promote_eltype(A, B)
    nrow = istop - istrt + one(I)

    @inbounds for j in jstrt:jstop
        jloc = j - jstrt

        Bj = B[j]; Δj = zero(T)

        if uplo === :L
            b = bstrt + jloc * (one(I) + nrow)
            ir = j + one(I):istop
        else
            b = bstrt + jloc * nrow
            ir = istrt:j - one(I)
        end

        if uplo === :L
            Aij = A.val[b]; b += one(I)

            if TA === :N
                Aji =      Aij
            else
                Aji = conj(Aij)
            end

            Δj = axpy(Aji, Bj, Δj)
        end


        for i in ir
            Aij = A.val[b]; b += one(I)

            if TA === :N
                Aji =      Aij
            else
                Aji = conj(Aij)
            end

             Δj  = axpy(    Aji, B[i], Δj  )
            C[i] = axpy(α * Aij, Bj,   C[i])
        end

        if uplo === :U
            Aij = A.val[b]; b += one(I)

            if TA === :N
                Aji =      Aij
            else
                Aji = conj(Aij)
            end

            Δj = axpy(Aji, Bj, Δj)
        end

        C[j] = axpy(α, Δj, C[j])
    end

    return
end

@inline function symv_blk_small_offd!(C::AbstractVector, tA::Val{TA}, A::BlockSparseMatrix{<:Any, I}, B::AbstractVector, α::Number, bstrt::I, istrt::I, istop::I, jstrt::I, jstop::I) where {TA, I}
    T = promote_eltype(A, B)
    b = bstrt

    @inbounds for j in jstrt:jstop
        Bj = B[j]; Δj = zero(T)

        for i in istrt:istop
            Aij = A.val[b]; b += one(I)

            if TA === :N
                Aji =      Aij
            else
                Aji = conj(Aij)
            end

            Δj   = axpy(    Aji, B[i], Δj)
            C[i] = axpy(α * Aij, Bj,   C[i])
        end

        C[j] = axpy(α, Δj, C[j])
    end

    return
end

@inline function symv_blk_large_diag!(C::AbstractVector, tA::Val{TA}, uplo::Symbol, A::BlockSparseMatrix{<:Any, I}, B::AbstractVector, α::Number, bstrt::I, nrow::I, ncol::I, istrt::I, istop::I, jstrt::I, jstop::I) where {TA, I}
    Ae = @inbounds block_impl(A, nrow, ncol, bstrt)
    U = @inbounds view(C, istrt:istop)
    V = @inbounds view(B, jstrt:jstop)
    mul!(U, wrapsym(Ae, tA, uplo), V, α, true)
    return
end

@inline function symv_blk_large_offd!(C::AbstractVector, tA::Val{TA}, A::BlockSparseMatrix{<:Any, I}, B::AbstractVector, α::Number, bstrt::I, nrow::I, ncol::I, istrt::I, istop::I, jstrt::I, jstop::I) where {TA, I}
    Ae = @inbounds block_impl(A, nrow, ncol, bstrt)

    BU = @inbounds view(B, istrt:istop)
    CU = @inbounds view(C, istrt:istop)
    BV = @inbounds view(B, jstrt:jstop)
    CV = @inbounds view(C, jstrt:jstop)

    if TA === :N
        mul!(CV, transpose(Ae), BU, α, true)
    else
        mul!(CV,   adjoint(Ae), BU, α, true)
    end

    mul!(CU, Ae, BV, α, true)
    return
end

function symv_impl!(C::AbstractVector, tA::Val{TA}, uplo::Symbol, A::BlockSparseMatrix{<:Any, I}, B::AbstractVector, α::Number, β::Number) where {TA, I}
    if iszero(β)
        fill!(C, β)
    elseif !isone(β)
        rmul!(C, β)
    end

    jstrt = one(I)
    estrt = one(I)

    @inbounds for v in vtxs(A)
        jstop = A.xcol[v + one(I)] - one(I)
        estop = A.xsrc[v + one(I)] - one(I)

        ncol = jstop - jstrt + one(I)

        for e in estrt:estop
            u = A.tgt[e]

            if intriangle(u, v, uplo)
                bstrt = A.xblk[e]
                istrt = A.xrow[u]
                istop = A.xrow[u + one(I)] - one(I)

                nrow = istop - istrt + one(I)

                if u == v
                    if isone(ncol) && isone(nrow)
                        symv_blk_unit_diag!(C, tA, A, B, α, bstrt, istrt)
                    elseif ncol < 24 && nrow < 24
                        symv_blk_small_diag!(C, tA, uplo, A, B, α, bstrt, istrt, istop, jstrt, jstop)
                    else
                        symv_blk_large_diag!(C, tA, uplo, A, B, α, bstrt, nrow, ncol, istrt, istop, jstrt, jstop)
                    end
                else
                    if isone(ncol) && isone(nrow)
                        symv_blk_unit_offd!(C, tA, A, B, α, bstrt, istrt, jstrt)
                    elseif ncol < 24 && nrow < 24
                        symv_blk_small_offd!(C, tA, A, B, α, bstrt, istrt, istop, jstrt, jstop)
                    else
                        symv_blk_large_offd!(C, tA, A, B, α, bstrt, nrow, ncol, istrt, istop, jstrt, jstop)
                    end
                end
            end
        end

        jstrt = jstop + one(I)
        estrt = estop + one(I)
    end

    return C
end
