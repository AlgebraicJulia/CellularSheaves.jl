const GEMV_MAXN = 8
const GEMV_FWD_TILE = 16

############################################################################################
# gemv_mn!
############################################################################################

@generated function gemv_mn!(C::AbstractVector, tA::Val{TA}, A::AbstractVector, B::AbstractVector, α::Number, bstrt::I, istrt::I, jstrt::I, ::Val{M}, ::Val{N}) where {TA, M, N, I}
    if TA === :N
        body = accum_body(M, N, :istrt, :jstrt, k -> :((t - 1) * $M + $(k - 1)), :tA)
    else
        body = accum_body(N, M, :jstrt, :istrt, k -> :($(k - 1) * $M + (t - 1)), :tA)
    end

    return quote
        $(Expr(:meta, :inline))

        $body

        return
    end
end

############################################################################################
# gemv_jam4!
############################################################################################

@inline function gemv_jam4!(C::AbstractVector, tA::Val{TA}, A::AbstractVector, B::AbstractVector, α::Number, b0::I, nrow::I, istrt::I, istop::I, j::I) where {TA, I}
    T = promote_eltype(A, B)

    Δ0 = zero(T)
    Δ1 = zero(T)
    Δ2 = zero(T)
    Δ3 = zero(T)

    @inbounds @simd for i in istrt:istop
        Bi = B[i]

        a0 = wrapadj(A[b0 +                   i], tA)
        a1 = wrapadj(A[b0 +            nrow + i], tA)
        a2 = wrapadj(A[b0 +   two(I) * nrow + i], tA)
        a3 = wrapadj(A[b0 + three(I) * nrow + i], tA)

        Δ0 = muladd(a0, Bi, Δ0)
        Δ1 = muladd(a1, Bi, Δ1)
        Δ2 = muladd(a2, Bi, Δ2)
        Δ3 = muladd(a3, Bi, Δ3)
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
# gemv_bwd!
############################################################################################

@generated function gemv_bwd!(C::AbstractVector, tA::Val{TA}, A::AbstractVector, B::AbstractVector, α::Number, bstrt::I, istrt::I, jstrt::I, nrow::Val{NROW}, ncol::I) where {TA, NROW, I}
    s    = [Symbol(:s, i) for i in 1:NROW]
    sset = [:($(s[i]) = α * B[istrt + $(i - 1)]) for i in 1:NROW]

    function chain()
        acc = :(C[jstrt + j - one($I)])

        for i in 1:NROW
            acc = :(muladd(wrapadj(A[b0 + $(i - 1)], tA), $(s[i]), $acc))
        end

        return acc
    end

    return quote
        $(Expr(:meta, :inline))
        $(sset...)

        @inbounds @simd for j in oneto(ncol)
            b0 = bstrt + (j - one($I)) * $NROW
            C[jstrt + j - one($I)] = $(chain())
        end

        return
    end
end

@inline function gemv_bwd!(C::AbstractVector, tA::Val{TA}, A::AbstractVector, B::AbstractVector, α::Number, bstrt::I, istrt::I, jstrt::I, nrow::I, ncol::I) where {TA, I}
    T = promote_eltype(A, B)

    istop = istrt + nrow - one(I)
    jstop = jstrt + ncol - one(I)
    b = bstrt
    j = jstrt

    @inbounds while j + three(I) <= jstop
        b0 = b - istrt

        gemv_jam4!(C, tA, A, B, α, b0, nrow, istrt, istop, j)

        b += four(I) * nrow
        j += four(I)
    end

    @inbounds while j <= jstop
        Δj = zero(T)

        @simd for i in istrt:istop
            Δj = muladd(wrapadj(A[b], tA), B[i], Δj)
            b += one(I)
        end

        C[j] = muladd(α, Δj, C[j])
        j += one(I)
    end

    return
end

############################################################################################
# gemv_fwd!
############################################################################################

@generated function gemv_fwd!(C::AbstractVector, A::AbstractVector, B::AbstractVector, α::Number, bstrt::I, istrt::I, jstrt::I, nrow::I, ncol::I, ::Val{TILE}) where {TILE, I}
    K = cld(8, TILE)
    Δ = [Symbol(:Δ, i, k) for i in 1:TILE, k in 1:K]

    function colupd(j, k)
        return quote
            Bj = B[jstrt + $j]
            $([:($(Δ[i, k]) = muladd(A[bstrt + $j * nrow + $(i - 1)], Bj, $(Δ[i, k]))) for i in 1:TILE]...)
        end
    end

    init  = [:($(Δ[i, k]) = zero(promote_eltype(A, B))) for i in 1:TILE for k in 1:K]
    store = [:(C[istrt + $(i - 1)] = muladd(α, $(Expr(:call, :+, Δ[i, :]...)), C[istrt + $(i - 1)])) for i in 1:TILE]
    group = [colupd(:(j + $(k - 1)), k) for k in 1:K]
    tail  = colupd(:j, 1)

    return quote
        $(Expr(:meta, :inline))
        @inbounds begin
            $(init...)

            j = zero($I)

            while j + $K <= ncol
                $(group...)
                j += $K
            end

            while j < ncol
                $tail
                j += one($I)
            end

            $(store...)
        end

        return
    end
end

@generated function gemv_fwd_rem!(C::AbstractVector, A::AbstractVector, B::AbstractVector, α::Number, bstrt::I, istrt::I, jstrt::I, nrow::I, ncol::I, i::I) where {I}
    body = :nothing

    for m in GEMV_FWD_TILE - 1:-1:1
        body = Expr(:if, :(rem == $m), :(gemv_fwd!(C, A, B, α, bstrt + i, istrt + i, jstrt, nrow, ncol, Val($m))), body)
    end

    return quote
        $(Expr(:meta, :inline))
        rem = nrow - i
        $body
        return
    end
end

############################################################################################
# gemv_blk!
############################################################################################

@generated function gemv_blk!(C::AbstractVector, tA::Val{TA}, A::AbstractVector, B::AbstractVector, α::Number, bstrt::I, istrt::I, jstrt::I, nrow::I, ncol::Val{NCOL}) where {TA, NCOL, I}
    body = :nothing

    for m in GEMV_MAXN:-1:1
        body = Expr(:if, :(nrow == $m), :(gemv_mn!(C, tA, A, B, α, bstrt, istrt, jstrt, Val($m), Val($NCOL))), body)
    end

    return quote
        $(Expr(:meta, :inline))
        @assert 1 <= nrow <= GEMV_MAXN

        $body

        return
    end
end

@inline function gemv_blk!(C::AbstractVector, ::Val{:N}, A::AbstractVector, B::AbstractVector, α::Number, bstrt::I, istrt::I, jstrt::I, nrow::I, ncol::I) where {I}
    i = zero(I)

    while i + GEMV_FWD_TILE <= nrow
        gemv_fwd!(C, A, B, α, bstrt + i, istrt + i, jstrt, nrow, ncol, Val(GEMV_FWD_TILE))
        i += GEMV_FWD_TILE
    end

    i < nrow && gemv_fwd_rem!(C, A, B, α, bstrt, istrt, jstrt, nrow, ncol, i)

    return
end

@generated function gemv_blk!(C::AbstractVector, tA::Val{TA}, A::AbstractVector, B::AbstractVector, α::Number, bstrt::I, istrt::I, jstrt::I, nrow::I, ncol::I) where {TA, I}
    body = :(gemv_bwd!(C, tA, A, B, α, bstrt, istrt, jstrt, nrow, ncol))

    for m in GEMV_MAXN:-1:1
        body = Expr(:if, :(nrow == $m), :(gemv_bwd!(C, tA, A, B, α, bstrt, istrt, jstrt, Val($m), ncol)), body)
    end

    return quote
        $(Expr(:meta, :inline))
        $body
        return
    end
end

############################################################################################
# gemv_vtx_impl!
############################################################################################

@generated function gemv_vtx_impl!(C::AbstractVector, tA::Val{TA}, A::BlockSparseMatrix{<:Any, I}, B::AbstractVector, α::Number, estrt::I, estop::I, bstrt::I, jstrt::I, ncol::Val{NCOL}) where {TA, NCOL, I}
    Δ     = [Symbol(:Δ, k) for k in 1:NCOL]
    init  = [:($(Δ[k]) = zero(promote_eltype(A, B))) for k in 1:NCOL]
    accum = [:($(Δ[k]) = muladd(wrapadj(A.val[b0 + $(k - 1) * nrow + i], tA), Bi, $(Δ[k]))) for k in 1:NCOL]
    store = [:(C[jstrt + $(k - 1)] = muladd(α, $(Δ[k]), C[jstrt + $(k - 1)])) for k in 1:NCOL]

    return quote
        $(Expr(:meta, :inline))

        $(init...)

        b = bstrt

        @inbounds for e in estrt:estop
            u = A.tgt[e]
            istrt = A.xrow[u]
            istop = A.xrow[u + one(I)] - one(I)
            nrow = istop - istrt + one(I)
            b0 = b - istrt

            @simd for i in istrt:istop
                Bi = B[i]

                $(accum...)
            end

            b += convert(I, NCOL) * nrow
        end

        @inbounds begin
            $(store...)
        end

        return b
    end
end

@inline function gemv_vtx_impl!(C::AbstractVector, tA::Val{:N}, A::BlockSparseMatrix{<:Any, I}, B::AbstractVector, α::Number, estrt::I, estop::I, bstrt::I, jstrt::I, ncol::Val{NCOL}) where {NCOL, I}
    @inbounds for e in estrt:estop
        u = A.tgt[e]
        istrt = A.xrow[u]
        nrow = A.xrow[u + one(I)] - istrt

        if nrow <= GEMV_MAXN
            gemv_blk!(C, tA, A.val, B, α, bstrt, istrt, jstrt, nrow, ncol)
        else
            gemv_blk!(C, tA, A.val, B, α, bstrt, istrt, jstrt, nrow, convert(I, NCOL))
        end

        bstrt += convert(I, NCOL) * nrow
    end

    return bstrt
end

@inline function gemv_vtx_impl!(C::AbstractVector, tA::Val{TA}, A::BlockSparseMatrix{<:Any, I}, B::AbstractVector, α::Number, estrt::I, estop::I, bstrt::I, jstrt::I, ncol::I) where {TA, I}
    @inbounds for e in estrt:estop
        u = A.tgt[e]
        istrt = A.xrow[u]
        nrow = A.xrow[u + one(I)] - istrt

        gemv_blk!(C, tA, A.val, B, α, bstrt, istrt, jstrt, nrow, ncol)

        bstrt += ncol * nrow
    end

    return bstrt
end

############################################################################################
# gemv_vtx!
############################################################################################

@generated function gemv_vtx!(C::AbstractVector, tA::Val{TA}, A::BlockSparseMatrix{<:Any, I}, B::AbstractVector, α::Number, estrt::I, estop::I, bstrt::I, jstrt::I, ncol::I) where {TA, I}
    body = :(gemv_vtx_impl!(C, tA, A, B, α, estrt, estop, bstrt, jstrt, ncol))

    for n in GEMV_MAXN:-1:1
        body = Expr(:if, :(ncol == $n), :(gemv_vtx_impl!(C, tA, A, B, α, estrt, estop, bstrt, jstrt, Val($n))), body)
    end

    return Expr(:block, Expr(:meta, :inline), body)
end

############################################################################################
# gemv_impl!
############################################################################################

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

        bstrt = gemv_vtx!(C, tA, A, B, α, estrt, estop, bstrt, jstrt, jstop - jstrt + one(I))

        jstrt = jstop + one(I)
        estrt = estop + one(I)
    end

    return C
end
