function gemm_impl!(C::AbstractMatrix, tA::Val{TA}, A::BlockSparseMatrix, tB::Val{TB}, B::AbstractMatrix, α::Number, β::Number) where {TA, TB}
    if iszero(β)
        fill!(C, β)
    elseif !isone(β)
        rmul!(C, β)
    end

    @inbounds for v in vtxs(A)
        cols = colrange(A, v)

        if TA !== :N
            V =           view(C, cols, :)
        elseif TB === :N
            V =           view(B, cols, :)
        elseif TB === :T
            V = transpose(view(B, :, cols))
        else
            V =   adjoint(view(B, :, cols))
        end

        for e in srcrange(A, v)
            u = A.tgt[e]
            rows = rowrange(A, u)

            if TA === :N
                U =           view(C, rows, :)
            elseif TB === :N
                U =           view(B, rows, :)
            elseif TB === :T
                U = transpose(view(B, :, rows))
            else
                U =   adjoint(view(B, :, rows))
            end

            if TA === :N
                mul!(U,           block(A, u, v, e),  V, α, true)
            elseif TA === :T
                mul!(V, transpose(block(A, u, v, e)), U, α, true)
            else
                mul!(V,   adjoint(block(A, u, v, e)), U, α, true)
            end
        end
    end

    return C
end

function gemm_impl!(C::BlockSparseMatrix{T, I}, ::Val{:N}, A::BlockSparseMatrix, ::Val{:N}, B::BlockSparseMatrix, α::Number, β::Number) where {T, I}
    if iszero(β)
        fill!(C, β)
    elseif !isone(β)
        rmul!(C, β)
    end

    if !iszero(α)
        marker = FVector{I}(undef, nouts(C))

        @inbounds for Cv in vtxs(C)
            for Ce in srcrange(C, Cv)
                Cu = C.tgt[Ce]
                marker[Cu] = Ce
            end

            for Be in srcrange(B, Cv)
                Bu = B.tgt[Be]
                Bx = block(B, Bu, Cv, Be)

                for Ae in srcrange(A, Bu)
                    Au = A.tgt[Ae]
                    Ax = block(A, Au, Bu, Ae)

                    Ce = marker[Au]
                    Cx = block(C, Au, Cv, Ce)

                    mul!(Cx, Ax, Bx, α, true)
                end
            end
        end
    end

    return C
end

function gemm_impl!(C::BlockSparseMatrix, tA::Val, A::BlockSparseMatrix, tB::Val, B::BlockSparseMatrix, α::Number, β::Number)
    An = convert(BlockSparseMatrix, wrapadj(A, tA))
    Bn = convert(BlockSparseMatrix, wrapadj(B, tB))
    return gemm_impl!(C, Val(:N), An, Val(:N), Bn, α, β)
end
