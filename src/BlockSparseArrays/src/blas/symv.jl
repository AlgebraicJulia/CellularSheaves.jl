const SYMV_MAXN = 8

############################################################################################
# symv_offd_mn!
############################################################################################

@generated function symv_offd_mn!(C::AbstractVector, tA::Val{TA}, A::AbstractVector, B::AbstractVector, α::Number, bstrt::I, istrt::I, jstrt::I, ::Val{M}, ::Val{N}) where {TA, M, N, I}
    fwd = accum_body(M, N, :istrt, :jstrt, k -> :((t - 1) * $M + $(k - 1)), :(Val(:N)))
    trn = accum_body(N, M, :jstrt, :istrt, k -> :($(k - 1) * $M + (t - 1)), :tA)

    return quote
        $(Expr(:meta, :inline))

        $fwd

        $trn

        return
    end
end

############################################################################################
# symv_diag_mn!
############################################################################################

@generated function symv_diag_mn!(C::AbstractVector, tA::Val{TA}, A::AbstractVector, B::AbstractVector, α::Number, bstrt::I, strt::I, ::Val{N}, ::Val{UPLO}) where {TA, N, UPLO, I}
    st = Expr[]

    for k in 1:N
        push!(st, :($(Symbol(:Δ, k)) = C[strt + $(k - 1)]))
    end

    for l in 1:N
        push!(st, :($(Symbol(:x, l)) = α * B[strt + $(l - 1)]))
    end

    for k in 1:N, l in 1:N
        intriangle(k, l, UPLO) || continue

        Akl = :(A[bstrt + $((l - 1) * N + (k - 1))])

        if k == l
            push!(st, :($(Symbol(:Δ, k)) = muladd($Akl, $(Symbol(:x, l)), $(Symbol(:Δ, k)))))
        else
            push!(st, :($(Symbol(:Δ, k)) = muladd($Akl, $(Symbol(:x, l)), $(Symbol(:Δ, k)))))
            push!(st, :($(Symbol(:Δ, l)) = muladd(wrapadj($Akl, tA), $(Symbol(:x, k)), $(Symbol(:Δ, l)))))
        end
    end

    for k in 1:N
        push!(st, :(C[strt + $(k - 1)] = $(Symbol(:Δ, k))))
    end

    return quote
        $(Expr(:meta, :inline))

        @inbounds begin
            $(st...)
        end

        return
    end
end

############################################################################################
# symv_jam4!
############################################################################################

@inline function symv_jam4!(C::AbstractVector, tA::Val{TA}, A::AbstractVector, B::AbstractVector, α::Number, b0::I, nrow::I, istrt::I, istop::I, j::I) where {TA, I}
    T = promote_eltype(A, B)

    Bj0 = B[j]
    Bj1 = B[j +   one(I)]
    Bj2 = B[j +   two(I)]
    Bj3 = B[j + three(I)]

    Δ0 = zero(T)
    Δ1 = zero(T)
    Δ2 = zero(T)
    Δ3 = zero(T)

    @inbounds @simd for i in istrt:istop
        Bi = B[i]

        a0 = A[b0 +                   i]
        a1 = A[b0 +            nrow + i]
        a2 = A[b0 +   two(I) * nrow + i]
        a3 = A[b0 + three(I) * nrow + i]

        Δ0 = muladd(wrapadj(a0, tA), Bi, Δ0)
        Δ1 = muladd(wrapadj(a1, tA), Bi, Δ1)
        Δ2 = muladd(wrapadj(a2, tA), Bi, Δ2)
        Δ3 = muladd(wrapadj(a3, tA), Bi, Δ3)

        Ci = C[i]
        Ci = muladd(α * a0, Bj0, Ci)
        Ci = muladd(α * a1, Bj1, Ci)
        Ci = muladd(α * a2, Bj2, Ci)
        Ci = muladd(α * a3, Bj3, Ci)
        C[i] = Ci
    end

    @inbounds begin
        C[j]            = muladd(α, Δ0, C[j           ])
        C[j +   one(I)] = muladd(α, Δ1, C[j +   one(I)])
        C[j +   two(I)] = muladd(α, Δ2, C[j +   two(I)])
        C[j + three(I)] = muladd(α, Δ3, C[j + three(I)])
    end

    return
end

############################################################################################
# symv_offd_blk!
############################################################################################

@generated function symv_offd_blk!(C::AbstractVector, tA::Val{TA}, A::AbstractVector, B::AbstractVector, α::Number, bstrt::I, istrt::I, jstrt::I, nrow::I, ncol::Val{NCOL}) where {TA, NCOL, I}
    body = :nothing

    for m in SYMV_MAXN:-1:1
        body = Expr(:if, :(nrow == $m), :(symv_offd_mn!(C, tA, A, B, α, bstrt, istrt, jstrt, Val($m), Val($NCOL))), body)
    end

    return quote
        $(Expr(:meta, :inline))
        @assert 1 <= nrow <= SYMV_MAXN

        $body

        return
    end
end

@inline function symv_offd_blk!(C::AbstractVector, tA::Val{TA}, A::AbstractVector, B::AbstractVector, α::Number, bstrt::I, istrt::I, jstrt::I, nrow::I, ncol::I) where {TA, I}
    T = promote_eltype(A, B)

    istop = istrt + nrow - one(I)
    jstop = jstrt + ncol - one(I)
    b = bstrt
    j = jstrt

    @inbounds while j + three(I) <= jstop
        b0 = b - istrt

        symv_jam4!(C, tA, A, B, α, b0, nrow, istrt, istop, j)

        b += four(I) * nrow
        j += four(I)
    end

    @inbounds while j <= jstop
        Bj = B[j]
        Δj = zero(T)

        @simd for i in istrt:istop
            Aij  = A[b]
            Δj   = muladd(wrapadj(Aij, tA), B[i], Δj)
            C[i] = muladd(α * Aij, Bj, C[i])
            b += one(I)
        end

        C[j] = muladd(α, Δj, C[j])
        j += one(I)
    end

    return
end

############################################################################################
# symv_diag_blk!
############################################################################################

@inline function symv_diag_blk!(C::AbstractVector, tA::Val{TA}, A::AbstractVector, B::AbstractVector, α::Number, bstrt::I, strt::I, ncol::Val{NCOL}, uplo::Val{UPLO}) where {TA, NCOL, UPLO, I}
    symv_diag_mn!(C, tA, A, B, α, bstrt, strt, ncol, uplo)
    return
end

@inline function symv_diag_blk!(C::AbstractVector, tA::Val{TA}, A::AbstractVector, B::AbstractVector, α::Number, bstrt::I, strt::I, nrow::I, uplo::Val{UPLO}) where {TA, UPLO, I}
    T = promote_eltype(A, B)

    istrt = strt
    istop = strt + nrow - one(I)
    jstrt = strt
    jstop = strt + nrow - one(I)
    j = jstrt

    @inbounds while j + three(I) <= jstop
        jloc = j - jstrt
        b0 = bstrt + jloc * nrow - istrt
        cstop = j + three(I)

        for l in j:cstop
            bcol = b0 + (l - j) * nrow
            Bl = B[l]

            C[l] = muladd(α * A[bcol + l], Bl, C[l])

            if UPLO === :L
                krng = l + one(I):cstop
            else
                krng = j:l - one(I)
            end

            for k in krng
                Akl = A[bcol + k]

                C[k] = muladd(α * Akl,              Bl,   C[k])
                C[l] = muladd(α * wrapadj(Akl, tA), B[k], C[l])
            end
        end

        if UPLO === :L
            ibstrt = j + four(I)
            ibstop = istop
        else
            ibstrt = istrt
            ibstop = j - one(I)
        end

        symv_jam4!(C, tA, A, B, α, b0, nrow, ibstrt, ibstop, j)

        j += four(I)
    end

    @inbounds while j <= jstop
        jloc = j - jstrt
        Bj = B[j]; Δj = zero(T)

        if UPLO === :L
            b = bstrt + jloc * (one(I) + nrow)
            irng = j + one(I):istop
        else
            b = bstrt + jloc * nrow
            irng = istrt:j - one(I)
        end

        if UPLO === :L
            Δj = muladd(wrapadj(A[b], tA), Bj, Δj)
            b += one(I)
        end

        @simd for i in irng
            Aij  = A[b]
            Δj   = muladd(wrapadj(Aij, tA), B[i], Δj)
            C[i] = muladd(α * Aij, Bj, C[i])
            b += one(I)
        end

        if UPLO === :U
            Δj = muladd(wrapadj(A[b], tA), Bj, Δj)
            b += one(I)
        end

        C[j] = muladd(α, Δj, C[j])
        j += one(I)
    end

    return
end

############################################################################################
# symv_vtx_impl!
############################################################################################

@inline function symv_vtx_impl!(C::AbstractVector, tA::Val{TA}, A::BlockSparseMatrix{<:Any, I}, B::AbstractVector, α::Number, v::I, estrt::I, estop::I, jstrt::I, ncol::Val{NCOL}, uplo::Val{UPLO}) where {TA, NCOL, UPLO, I}
    @inbounds for e in estrt:estop
        u = A.tgt[e]

        if intriangle(u, v, uplo)
            bstrt = A.xblk[e]
            istrt = A.xrow[u]
            nrow = A.xrow[u + one(I)] - istrt

            if u == v
                symv_diag_blk!(C, tA, A.val, B, α, bstrt, istrt, ncol, uplo)
            elseif nrow <= SYMV_MAXN
                symv_offd_blk!(C, tA, A.val, B, α, bstrt, istrt, jstrt, nrow, ncol)
            else
                symv_offd_blk!(C, tA, A.val, B, α, bstrt, istrt, jstrt, nrow, convert(I, NCOL))
            end
        end
    end

    return
end

@inline function symv_vtx_impl!(C::AbstractVector, tA::Val{TA}, A::BlockSparseMatrix{<:Any, I}, B::AbstractVector, α::Number, v::I, estrt::I, estop::I, jstrt::I, ncol::I, uplo::Val{UPLO}) where {TA, UPLO, I}
    @inbounds for e in estrt:estop
        u = A.tgt[e]

        if intriangle(u, v, uplo)
            bstrt = A.xblk[e]
            istrt = A.xrow[u]
            nrow = A.xrow[u + one(I)] - istrt

            if u == v
                symv_diag_blk!(C, tA, A.val, B, α, bstrt, istrt, nrow, uplo)
            else
                symv_offd_blk!(C, tA, A.val, B, α, bstrt, istrt, jstrt, nrow, ncol)
            end
        end
    end

    return
end

############################################################################################
# symv_vtx!
############################################################################################

@generated function symv_vtx!(C::AbstractVector, tA::Val{TA}, A::BlockSparseMatrix{<:Any, I}, B::AbstractVector, α::Number, v::I, estrt::I, estop::I, jstrt::I, ncol::I, uplo::Val{UPLO}) where {TA, UPLO, I}
    body = :(symv_vtx_impl!(C, tA, A, B, α, v, estrt, estop, jstrt, ncol, uplo))

    for n in SYMV_MAXN:-1:1
        body = Expr(:if, :(ncol == $n), :(symv_vtx_impl!(C, tA, A, B, α, v, estrt, estop, jstrt, Val($n), uplo)), body)
    end

    return Expr(:block, Expr(:meta, :inline), body)
end

############################################################################################
# symv_impl!
############################################################################################

function symv_impl!(C::AbstractVector, tA::Val{TA}, uplo::Symbol, A::BlockSparseMatrix{<:Any, I}, B::AbstractVector, α::Number, β::Number) where {TA, I}
    if uplo === :L
        return symv_impl!(C, tA, Val(:L), A, B, α, β)
    else
        return symv_impl!(C, tA, Val(:U), A, B, α, β)
    end
end

function symv_impl!(C::AbstractVector, tA::Val{TA}, uplo::Val{UPLO}, A::BlockSparseMatrix{<:Any, I}, B::AbstractVector, α::Number, β::Number) where {TA, UPLO, I}
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

        symv_vtx!(C, tA, A, B, α, v, estrt, estop, jstrt, jstop - jstrt + one(I), uplo)

        jstrt = jstop + one(I)
        estrt = estop + one(I)
    end

    return C
end
