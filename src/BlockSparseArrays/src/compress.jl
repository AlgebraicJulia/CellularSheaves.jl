function twins(A::AdjOrTrans)
    return twins(copy(A))
end

function twins(A::SparseMatrixCSC)
    return twins(size(A)..., A.colptr, A.rowval)
end

function twins(nout::Integer, nvtx::Integer, xsrc::AbstractVector{I}, tgt::AbstractVector{I}) where {I}
    new  = FScalar{I}(undef)
    var  = FVector{I}(undef, nout + 1)
    svar = FVector{I}(undef, nout)
    flag = FVector{I}(undef, nout)
    card = FVector{I}(undef, nout)
    head = FVector{I}(undef, nout)
    prev = FVector{I}(undef, nout)
    next = FVector{I}(undef, nout)
    nout = convert(I, nout)
    nvtx = convert(I, nvtx)
    nsv = twins_impl!(new, var, svar, flag, card, head, prev, next, nout, nvtx, xsrc, tgt)
    return nsv, var, flag
end

function twins_impl!(
        new::AbstractScalar{I},
        var::AbstractVector{I},
        svar::AbstractVector{I},
        flag::AbstractVector{I},
        card::AbstractVector{I},
        head::AbstractVector{I},
        prev::AbstractVector{I},
        next::AbstractVector{I},
        nout::I,
        nvtx::I,
        xsrc::AbstractVector{I},
        tgt::AbstractVector{I},
    ) where {I}
    @assert nout < length(var)
    @assert nout <= length(svar)
    @assert nout <= length(flag)
    @assert nout <= length(card)
    @assert nout <= length(head)
    @assert nout <= length(prev)
    @assert nout <= length(next)

    nsv = zero(I)

    @inbounds for i in oneto(nout)
        flag[i] = zero(I)
        card[i] = zero(I)
        head[i] = zero(I)
    end

    new[] = zero(I); free = SinglyLinkedList(new, var)
    @inbounds prepend!(free, oneto(nout))

    function set(i::I)
        @inbounds h = view(head, i)
        return DoublyLinkedList(h, prev, next)
    end

    if ispositive(nout)
        s = popfirst!(free); nsv += one(I)
        @inbounds prepend!(set(s), oneto(nout)); card[s] = nout

        @inbounds for i in oneto(nout)
            svar[i] = s
        end
    end

    @inbounds for j in oneto(nvtx)
        # rows present in column j
        for e in xsrc[j]:(xsrc[j + one(I)] - one(I))
            i = tgt[e]
            # `s` is the old supervariable of row `i`
            s = svar[i]

            # first occurrence of `s` for column `j`
            if flag[s] < j
                flag[s] = j

                if card[s] > one(I)
                    delete!(set(s), i); card[s] -= one(I)
                    var[s] = i
                    ns = svar[i] = popfirst!(free); nsv += one(I)
                    pushfirst!(set(ns), i); card[ns] += one(I)
                end
            # second or later occurrence of `s` for column `j`
            else
                delete!(set(s), i); card[s] -= one(I)
                k = var[s]
                ns = svar[i] = svar[k]
                pushfirst!(set(ns), i); card[ns] += one(I)

                if isempty(set(s))
                    pushfirst!(free, s); nsv -= one(I)
                end
            end
        end
    end

    t = one(I); var[t] = p = one(I)

    @inbounds for s in oneto(nout)
        isempty(set(s)) && continue

        for i in set(s)
            svar[i] = t
            flag[p] = i; p += one(I)
        end

        t += one(I); var[t] = p
    end

    return nsv
end

function compress(A::AbstractMatrix)
    return compress(sparse(A))
end

function compress(A::SparseMatrixCSC{T, I}) where {T, I}
    nvtx, xcol, colprm = twins(transpose(A))
    nout, xrow, rowprm = twins(A)

    B = permute(A, rowprm, colprm)

    return rowprm, colprm, compress(B, nvtx, nout, xcol, xrow)
end

function compress(A::SparseMatrixCSC{T, I}, nvtx::I, nout::I, xcol::FVector{I}, xrow::FVector{I}) where {T, I}
    nrow = convert(I, size(A, 1))
    ncol = convert(I, size(A, 2))
    nbnz = convert(I, nnz(A))

     xsrc = FVector{I}(undef, nvtx + one(I))
    xblk = FVector{I}(undef, nbnz + one(I))
     row = FVector{I}(undef, nrow)
     tgt = FVector{I}(undef, nbnz)
     val = FVector{T}(undef, nbnz)

    @inbounds for u in oneto(nout)
        istrt = xrow[u]
        istop = xrow[u + one(I)] - one(I)

        for i in istrt:istop
            row[i] = u
        end
    end

    e = b = zero(I)

    @inbounds for v in oneto(nvtx)
        xsrc[v] = e + one(I)

        jstrt = xcol[v]
        jstop = xcol[v + one(I)] - one(I)

        t = zero(I)

        for f in nzrange(A, jstrt)
            i = rowvals(A)[f]
            u = row[i]

            if t < u
                t = u; e += one(I)

                 tgt[e] = u
                xblk[e] = b + one(I)
            end

            b += jstop - jstrt + one(I)
        end
    end

    narc = e
    xsrc[nvtx + one(I)] = narc + one(I)
    xblk[narc + one(I)] = nbnz + one(I)

    @inbounds for v in oneto(nvtx)
        estrt = xsrc[v]
        estop = xsrc[v + one(I)] - one(I)

        jstrt = xcol[v]
        jstop = xcol[v + one(I)] - one(I)

        for j in jstrt:jstop
            f = A.colptr[j]

            for e in estrt:estop
                u = tgt[e]

                istrt = xrow[u]
                istop = xrow[u + one(I)] - one(I)

                b = xblk[e] + (j - jstrt) * (istop - istrt + one(I))

                for i in istrt:istop
                    val[b] = nonzeros(A)[f]
                    b += one(I)
                    f += one(I)
                end
            end
        end
    end

    return BlockSparseMatrix{T, I}(nout, nvtx, narc, ncol, nrow, nbnz, xsrc, xcol, xrow, xblk, tgt, val)
end
