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

function twins(A::AdjOrTrans, At::AbstractMatrix, tau::Number)
    if isone(tau)
        return twins(A)
    else
        return twins(copy(A), At, tau)
    end
end

function twins(A::SparseMatrixCSC, At::AdjOrTrans, tau::Number)
    if isone(tau)
        return twins(A)
    else
        return twins(A, copy(At), tau)
    end
end

function twins(A::SparseMatrixCSC{T, I}, At::SparseMatrixCSC{T, I}, tau::Number) where {T, I}
    if isone(tau)
        return twins(A)
    else
        nout = convert(I, size(A, 1))
        nvtx = convert(I, size(A, 2))
        xrow = At.colptr; cadj = At.rowval
        xcol = A.colptr;  radj = A.rowval

        return twins(nout, nvtx, xrow, cadj, xcol, radj, tau)
    end
end

function twins(
        nout::I,
        nvtx::I,
        xrow::AbstractVector{I},
        cadj::AbstractVector{I},
        xcol::AbstractVector{I},
        radj::AbstractVector{I},
        tau::Number,
    ) where {I}
    head   = FScalar{I}(undef)
    prev   = FVector{I}(undef, nout)
    next   = FVector{I}(undef, nout)
    adjmap = FVector{I}(undef, nout)
    cosine = FVector{I}(undef, nout)
    degree = FVector{I}(undef, nout)
    xptr   = FVector{I}(undef, nout + one(I))
    perm   = FVector{I}(undef, nout)
    nb = twins_impl!(head, prev, next, adjmap, cosine, degree, xptr, perm,
                     nout, nvtx, xrow, cadj, xcol, radj, convert(Float64, tau))
    return nb, xptr, perm
end

function twins_impl!(
        head::AbstractScalar{I},
        prev::AbstractVector{I},
        next::AbstractVector{I},
        adjmap::AbstractVector{I},
        cosine::AbstractVector{I},
        degree::AbstractVector{I},
        xptr::AbstractVector{I},
        perm::AbstractVector{I},
        nout::I,
        nvtx::I,
        xrow::AbstractVector{I},
        cadj::AbstractVector{I},
        xcol::AbstractVector{I},
        radj::AbstractVector{I},
        tau::T,
    ) where {I, T}
    @assert zero(T) < tau <= one(T)

    list = DoublyLinkedList(head, prev, next)

    @inbounds for i in oneto(nout)
        adjmap[i] = zero(I)
        cosine[i] = zero(I)
        degree[i] = xrow[i + one(I)] - xrow[i]
    end

    head[] = zero(I)
    nb = zero(I)

    @inbounds for i in oneto(nout)
        iszero(adjmap[i]) || continue

        # start a new set
        adjmap[i] = nb += one(I)

        # for each entry j in row i...
        for e in xrow[i]:xrow[i + one(I)] - one(I)
            j = cadj[e]

            for f in (xcol[j + one(I)] - one(I)):-one(I):xcol[j]
                k = radj[f]

                # both rows i and k have an entry in column j;
                # the column is sorted, so stop when k reaches i
                k > i || break

                # k has not yet been added to a part
                if iszero(adjmap[k])
                    # increase partial dot product
                    iszero(cosine[k]) && pushfirst!(list, k)
                    cosine[k] += one(I)
                end
            end
        end

        for k in list
            cos = convert(T, cosine[k])
            nzi = convert(T, degree[i])
            nzk = convert(T, degree[k])

            # test similarity of row patterns
            if cos * cos >= tau * tau * nzi * nzk
                adjmap[k] = nb
            end

            delete!(list, k)
            cosine[k] = zero(I)
        end
    end

    @inbounds for t in oneto(nb + one(I))
        xptr[t] = zero(I)
    end

    @inbounds for i in oneto(nout)
        xptr[adjmap[i] + one(I)] += one(I)
    end

    xptr[one(I)] = one(I)

    @inbounds for t in oneto(nb)
        xptr[t + one(I)] += xptr[t]
    end

    @inbounds for t in oneto(nb)
        cosine[t] = xptr[t]
    end

    @inbounds for i in oneto(nout)
        t = adjmap[i]
        perm[cosine[t]] = i
        cosine[t] += one(I)
    end

    return nb
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

function compress2(A::SparseMatrixCSC{T, I}, nvtx::I, nout::I, xcol::FVector{I}, xrow::FVector{I}) where {T, I}
    nrow = convert(I, size(A, 1))
    ncol = convert(I, size(A, 2))
    marc = convert(I, nnz(A))

    xsrc = FVector{I}(undef, nvtx + one(I))
    xblk = FVector{I}(undef, marc + one(I))
     row = FVector{I}(undef, nrow)
     col = FVector{I}(undef, ncol)
     tgt = FVector{I}(undef, marc)

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

        vncol = jstop - jstrt + one(I)

        t = zero(I); i = typemax(I)

        for j in jstrt:jstop
            f = col[j] = A.colptr[j]

            if f < A.colptr[j + one(I)]
                i = min(i, rowvals(A)[f])
            end
        end

        while i < typemax(I)
            u = row[i]

            istrt = xrow[u]
            istop = xrow[u + one(I)] - one(I)

            unrow = istop - istrt + one(I)

            if t < u
                t = u; e += one(I)

                 tgt[e] = u
                xblk[e] = b + one(I)

                b += unrow * vncol
            end

            inext = typemax(I)

            for j in jstrt:jstop
                f = col[j]
                fstop = A.colptr[j + one(I)] - one(I)

                if f ≤ fstop
                    if rowvals(A)[f] == i
                        f = col[j] = f + one(I)
                    end

                    if f ≤ fstop
                        inext = min(inext, rowvals(A)[f])
                    end
                end
            end

            i = inext
        end
    end

    narc = e
    nbnz = b

    xsrc[nvtx + one(I)] = narc + one(I)
    xblk[narc + one(I)] = nbnz + one(I)

    val = FVector{T}(undef, nbnz)

    @inbounds for v in oneto(nvtx)
        estrt = xsrc[v]
        estop = xsrc[v + one(I)] - one(I)

        jstrt = xcol[v]
        jstop = xcol[v + one(I)] - one(I)

        for j in jstrt:jstop
            col[j] = A.colptr[j]
        end

        for e in estrt:estop
            u =  tgt[e]
            b = xblk[e]

            istrt = xrow[u]
            istop = xrow[u + one(I)] - one(I)

            for j in jstrt:jstop
                f =  col[j]
                fstop = A.colptr[j + one(I)] - one(I)

                for i in istrt:istop
                    if f ≤ fstop && rowvals(A)[f] == i
                        Aij = nonzeros(A)[f]; f += one(I)
                    else
                        Aij = zero(T)
                    end

                    val[b] = Aij; b += one(I)
                end

                col[j] = f
            end
        end
    end

    return BlockSparseMatrix{T, I}(nout, nvtx, narc, ncol, nrow, nbnz, xsrc, xcol, xrow, xblk, tgt, val)
end
