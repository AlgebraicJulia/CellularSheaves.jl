function symm_impl!(C::AbstractMatrix, tA::Val{TA}, uplo::Symbol, A::BlockSparseMatrix, tB::Val{TB}, B::AbstractMatrix, α::Number, β::Number) where {TA, TB}
    if iszero(β)
        fill!(C, β)
    elseif !isone(β)
        rmul!(C, β)
    end

    @inbounds for v in vtxs(A)
        cols = colrange(A, v)
        CV = view(C, cols, :)

        if TB === :N
            BV =           view(B, cols, :)
        elseif TB === :T
            BV = transpose(view(B, :, cols))
        else
            BV =   adjoint(view(B, :, cols))
        end

        for e in srcrange(A, v)
            u = A.tgt[e]

            if intriangle(u, v, uplo)
                AE = block(A, u, v, e)

                if u == v
                    mul!(CV, wrapsym(AE, tA, uplo), BV, α, 1)
                else
                    rows = rowrange(A, u)
                    CU = view(C, rows, :)

                    if TB === :N
                        BU =           view(B, rows, :)
                    elseif TB === :T
                        BU = transpose(view(B, :, rows))
                    else
                        BU =   adjoint(view(B, :, rows))
                    end

                    if TA === :N
                        mul!(CV, transpose(AE), BU, α, 1)
                    else
                        mul!(CV,   adjoint(AE), BU, α, 1)
                    end

                    mul!(CU, AE, BV, α, 1)
                end
            end
        end
    end

    return C
end
