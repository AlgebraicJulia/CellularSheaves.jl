"""
    IPMProblem{T, I, C}

A primal-dual pair of conic quadratic programs.

```math
\\begin{aligned}
\\text{(P)} \\quad & \\min_p \\quad \\tfrac{1}{2} p^\\top Q p - f^\\top p \\\\
& \\text{s.t.} \\quad B p = g, \\quad p \\in K
\\end{aligned}
\\qquad
\\begin{aligned}
\\text{(D)} \\quad & \\max_{p, y} \\quad -\\tfrac{1}{2} p^\\top Q p + g^\\top y \\\\
& \\text{s.t.} \\quad Q p - f - B^\\top y = d, \\quad d \\in K^*
\\end{aligned}
```

Solve using [`solve`](@ref).
"""

struct IPMProblem{T, I, V <: AbstractCone}
    Q::BlockSparseMatrix{T, I}
    B::BlockSparseMatrix{T, I}
    f::FVector{T}
    g::FVector{T}
    K::FVector{V}
    R::FPermutation{I}
    C::FPermutation{I}

    function IPMProblem{T, I, V}(Q::BlockSparseMatrix, B::BlockSparseMatrix, f::FVector, g::FVector, K::FVector, R::FPermutation, C::FPermutation) where {T, I, V <: AbstractCone}
        @assert nrows(B) == length(g)
        @assert ncols(B) == ncols(Q) == length(f)
        @assert nvtxs(B) == nvtxs(Q) == length(K)
        @assert length(R.perm) == nrows(B)
        @assert length(C.perm) == ncols(B)

        for v in vtxs(B)
            @assert ncols(B, v) == ncols(Q, v)
        end

        return new{T, I, V}(Q, B, f, g, K, R, C)
    end
end

"""
    IPMProblem(Q, B, f, g, K)

Construct an [`IPMProblem`](@ref).
"""
function IPMProblem(Q::BlockSparseMatrix{T, I}, B::BlockSparseMatrix{T, I}, f::AbstractVector{T}, g::AbstractVector{T}, K::AbstractVector{V}) where {T, I, V <: AbstractCone}
    R = FPermutation{I}(rows(B))
    C = FPermutation{I}(cols(B))
    return IPMProblem{T, I, V}(Q, B, f, g, K, R, C)
end

function IPMProblem(Q::BlockSparseMatrix{T, I}, B::BlockSparseMatrix{T, I}, f::AbstractVector{T}, g::AbstractVector{T}, K::AbstractVector{V}, R::FPermutation{I}, C::FPermutation{I}) where {T, I, V <: AbstractCone}
    return IPMProblem{T, I, V}(Q, B, f, g, K, R, C)
end

function IPMProblem{T, I, V}(Q::BlockSparseMatrix, B::BlockSparseMatrix, f::AbstractVector, g::AbstractVector, K::AbstractVector, R::FPermutation, C::FPermutation) where {T, I, V <: AbstractCone}
    if !(f isa FVector{T})
        f = FVector{T}(f)
    end

    if !(g isa FVector{T})
        g = FVector{T}(g)
    end

    if !(K isa FVector{V})
        K = FVector{V}(K)
    end

    return IPMProblem{T, I, V}(Q, B, f, g, K, R, C)
end

function symbkkt(prob::IPMProblem, alg::EliminationAlgorithm)
    return symbkkt(prob.B, prob.Q, prob.f, prob.g, prob.K, prob.C, prob.R, alg)
end

"""
    IPMProblem(Q, B, f, g, K, s)

Construct an [`IPMProblem`](@ref).
"""
function IPMProblem(
        Q::SparseMatrixCSC{T},
        B::SparseMatrixCSC{T},
        f::AbstractVector,
        g::AbstractVector,
        K::AbstractVector,
        s::AbstractVector;
        tau::Real = 1.0,
    ) where {T}
    m, n = size(B); k = length(K)

    @assert size(Q, 1) == n
    @assert size(Q, 2) == n
    @assert length(f)  == n
    @assert length(g)  == m
    @assert length(s)  == k
    @assert    sum(s)  == n
    @assert  0 <  tau  <= 1

    Q = dropoffdz(Q)
    B = dropzeros(B)

    nvtx, xcol, colperm, K = colpartition(B, Q, K, s, tau)
    C = colcompress(B, nvtx, xcol, colperm)
    nout, xrow, rowperm = twins(C, transpose(C), tau)

    Bp = permute(B, rowperm, colperm)
    Qp = permute(Q, colperm, colperm)

    Bb = compress2(Bp, nvtx, nout, xcol, xrow)
    Qb = compress2(Qp, nvtx, nvtx, xcol, xcol)

    cp = FVector{T}(undef, n)
    gp = FVector{T}(undef, m)

    for j in 1:n
        cp[j] = f[colperm[j]]
    end

    for i in 1:m
        gp[i] = g[rowperm[i]]
    end

    return IPMProblem(Qb, Bb, cp, gp, K, Permutation(rowperm), Permutation(colperm))
end

function dropoffdz(Q::SparseMatrixCSC{T, I}) where {T, I}
    n = size(Q, 1)
    k = n

    for j in 1:n
        for p in nzrange(Q, j)
            i = Q.rowval[p]
            x = Q.nzval[p]

            if i != j && !iszero(x)
                k += 1
            end
        end
    end

    colptr = Vector{I}(undef, n + 1)
    rowval = Vector{I}(undef, k)
     nzval = Vector{T}(undef, k)

    q = 0

    for j in 1:n
        colptr[j] = q + 1

        iprv = 0

        for p in nzrange(Q, j)
            i = Q.rowval[p]
            x = Q.nzval[p]

            if !iszero(x) && iprv < j < i
                q += 1
                rowval[q] = iprv = j
                nzval[q] = zero(T)
            end

            if !iszero(x) || i == j
                q += 1
                rowval[q] = iprv = i
                nzval[q] = x
            end
        end

        if iprv < j
            q += 1
            rowval[q] = j
            nzval[q] = zero(T)
        end
    end

    colptr[n + 1] = k + 1

    return SparseMatrixCSC{T, I}(n, n, colptr, rowval, nzval)
end

function colcompress(A::SparseMatrixCSC{T, I}, nvtx::I, xcol::AbstractVector{I}, colperm::AbstractVector{I}) where {T, I}
    m = size(A, 1)
    n = nnz(A)

    colptr = Vector{I}(undef, nvtx + one(I))
    rowval = Vector{I}(undef, n)
    curptr = Vector{I}(undef, size(A, 2))

    e = zero(I)

    @inbounds for v in oneto(nvtx)
        colptr[v] = e + one(I)

        jstrt = xcol[v]
        jstop = xcol[v + one(I)] - one(I)

        i = typemax(I)

        for j in jstrt:jstop
            k = colperm[j]
            f = curptr[k] = A.colptr[k]

            if f < A.colptr[k + one(I)]
                i = min(i, A.rowval[f])
            end
        end

        while i < typemax(I)
            e += one(I)
            rowval[e] = i

            inext = typemax(I)

            for j in jstrt:jstop
                k = colperm[j]
                f = curptr[k]
                fstop = A.colptr[k + one(I)] - one(I)

                if f ≤ fstop
                    if A.rowval[f] == i
                        f = curptr[k] = f + one(I)
                    end

                    if f ≤ fstop
                        inext = min(inext, A.rowval[f])
                    end
                end
            end

            i = inext
        end
    end

    colptr[nvtx + one(I)] = e + one(I)
    resize!(rowval, e)

    return SparseMatrixCSC{T, I}(m, nvtx, colptr, rowval, ones(T, e))
end

function colpartition(B::SparseMatrixCSC{T, I}, Q::SparseMatrixCSC{T, I}, K::AbstractVector, s::AbstractVector, tau::Real) where {T, I}
    ncol = convert(I, size(B, 2))

    perm = FVector{I}(undef, ncol)
    work = FVector{I}(undef, ncol)

    ncolP = ncolC = nvtxR = zero(I)

    for (k, n) in zip(K, s)
        m = convert(I, n)

        if k isa PositiveCone
            ncolP += m
        elseif k isa CofreeCone
            ncolC += m
        else
            nvtxR += one(I)
        end
    end

    p = one(I)
    c = one(I)
    r = one(I)
    j = one(I)

    for (k, n) in zip(K, s)
        m = convert(I, n)

        if k isa PositiveCone
            for _ in oneto(m)
                work[                p] = j
                p += one(I)
                j += one(I)
            end
        elseif k isa CofreeCone
            for _ in oneto(m)
                work[ncolP         + c] = j
                c += one(I)
                j += one(I)
            end
        else
            for _ in oneto(m)
                work[ncolP + ncolC + r] = j
                r += one(I)
                j += one(I)
            end
        end
    end

    indP = view(work,         one(I):ncolP)
    indC = view(work, ncolP + one(I):ncolP + ncolC)

    XP = [Q[:, indP]; B[:, indP]]
    XC = [Q[:, indC]; B[:, indC]]

    nvtxP, xcolP, permP = twins(transpose(XP), XP, tau)
    nvtxC, xcolC, permC = twins(transpose(XC), XC, tau)

    nvtx = nvtxP + nvtxC + nvtxR

    xcol = FVector{I}(undef, nvtx + one(I))
    cone = FVector{AbstractCone}(undef, nvtx)

    for j in oneto(ncolP)
        perm[j] = work[permP[j]]
    end

    for v in oneto(nvtxP)
        xcol[v] = xcolP[v]
        cone[v] = PositiveCone()
    end

    for j in oneto(ncolC)
        perm[ncolP + j] = work[ncolP + permC[j]]
    end

    for v in oneto(nvtxC)
        xcol[nvtxP + v] = ncolP + xcolC[v]
        cone[nvtxP + v] = CofreeCone()
    end

    for j in ncolP + ncolC + one(I):ncol
        perm[j] = work[j]
    end

    r = one(I)
    j = one(I)

    for (k, n) in zip(K, s)
        m = convert(I, n)

        if !(k isa PositiveCone || k isa CofreeCone)
            xcol[nvtxP + nvtxC + r] = ncolP + ncolC + j
            cone[nvtxP + nvtxC + r] = k
            r += one(I)
            j += m
        end
    end

    xcol[nvtx + one(I)] = ncol + one(I)

    return nvtx, xcol, perm, cone
end

function showproblem(io::IO, prob::IPMProblem; indent::Integer=0)
    pad = " "^indent
    println(io, pad, "(P)  min   ½ pᵀ Q p - fᵀ p")
    println(io, pad, "     s.t.  B p = g,  p ∈ K")
    println(io)
    println(io, pad, "(D)  max  -½ pᵀ Q p + gᵀ y")
    println(io, pad, "     s.t.  Q p - f - Bᵀ y ∈ K*")
    println(io)
    @printf(io, "%sQ: %6d × %-6d  f: %d\n", pad, size(prob.Q, 1), size(prob.Q, 2), length(prob.f))
    @printf(io, "%sB: %6d × %-6d  g: %d\n", pad, size(prob.B, 1), size(prob.B, 2), length(prob.g))
    println(io)

    nnoc = npos = nsoc = nsdp = nexp = npow = 0

    for k in prob.K
        if k isa CofreeCone
            nnoc += 1
        elseif k isa PositiveCone
            npos += 1
        elseif k isa SecondOrderCone
            nsoc += 1
        elseif k isa SemidefiniteCone
            nsdp += 1
        elseif k isa ExponentialCone
            nexp += 1
        elseif k isa PowerCone
            npow += 1
        end
    end

    println(io, pad, "K:")
    nnoc > 0 && @printf(io, "%s    CofreeCone:       %d\n", pad, nnoc)
    npos > 0 && @printf(io, "%s    PositiveCone:     %d\n", pad, npos)
    nsoc > 0 && @printf(io, "%s    SecondOrderCone:  %d\n", pad, nsoc)
    nsdp > 0 && @printf(io, "%s    SemidefiniteCone: %d\n", pad, nsdp)
    nexp > 0 && @printf(io, "%s    ExponentialCone:  %d\n", pad, nexp)
    npow > 0 && @printf(io, "%s    PowerCone:        %d\n", pad, npow)
    return
end

function Base.show(io::IO, ::MIME"text/plain", prob::T) where {T <: IPMProblem}
    println(io, T, ":")
    return showproblem(io, prob; indent=2)
end

function CommonSolve.solve(prob::IPMProblem{T}, settings::AbstractSettings{T}) where {T}
    return solve!(init(prob, settings))
end
