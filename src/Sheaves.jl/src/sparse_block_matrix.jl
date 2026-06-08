############################################################################################
# Construction
############################################################################################

"""
    SparseBlockMatrix{T, I} <: AbstractMatrix{T}

A sparse block matrix.
"""
struct SparseBlockMatrix{T, I} <: AbstractMatrix{T}
    #
    # We represent sparse block matrices as diagrams D: I → Set:
    #
    #   ┌───────────────────────────┐
    #   │       R           C       │
    #   │   row ↓ tgt   src ↓ col   │
    #   │       U  ←  E  →  V       │
    #   │             ↑ blk         │
    #   │             B             │
    #   │             ↓ val         │
    #   │             T             │
    #   └───────────────────────────┘
    #
    # where
    #
    #   - U is a set of out vertices
    #   - V is a set of vertices
    #   - E is a set of arcs
    #   - R is a set of rows
    #   - C is a set of columns
    #   - B is a set of nonzero indices
    #
    # and T is the element type. Several of these functions are stored
    # directly as arrays, and several are stored indirectly, via their preimages,
    # which are prefixed with the letter "x".
    #
    # # Attributes
    #
    # The number of out vertices:
    #
    #    | U |
    #
    nout::I
    #
    # The number of vertices:
    #
    #    | V |
    #
    nvtx::I
    #
    # The number of arcs:
    #
    #    | E |
    #
    narc::I
    #
    # The number of columns:
    #
    #    | C |
    #
    ncol::I
    #
    # The number of rows:
    #
    #    | R |
    #
    nrow::I
    #
    # The number of nonzero entries:
    #
    #    | B |
    #
    nblk::I
    #
    # Each vertex v ∈ V is incident to the arcs
    #
    #    xsrc[v] ... xsrc[v + 1] - 1 ∈ E
    #
    xsrc::FVector{I}
    #
    # Each vertex v ∈ V is incident to the columns
    #
    #    xcol[v] ... xcol[v + 1] - 1 ∈ C
    #
    xcol::FVector{I}
    #
    # Each out vertex u ∈ U is incident to the rows
    #
    #    xrow[u] ... xrow[u + 1] - 1 ∈ R
    #
    xrow::FVector{I}
    #
    # Each arc e ∈ E is incident to the nonzero indices
    #
    #    xblk[e] ... xblk[e + 1] - 1 ∈ B
    #
    xblk::FVector{I}
    #
    # Each arc e ∈ E targets an out vertex
    #
    #    tgt[e] ∈ U
    #
    tgt::FVector{I}
    #
    # The nonzero entries:
    #
    #    val[b] ∈ T
    #
    val::FVector{T}

    function SparseBlockMatrix{T, I}(
            nout::I,
            nvtx::I,
            narc::I,
            ncol::I,
            nrow::I,
            nblk::I,
            xsrc::FVector{I},
            xcol::FVector{I},
            xrow::FVector{I},
            xblk::FVector{I},
            tgt::FVector{I},
            val::FVector{T},
        ) where {T, I}
        @assert nvtx < length(xsrc)
        @assert nvtx < length(xcol)
        @assert nout < length(xrow)
        @assert narc < length(xblk)
        @assert narc ≤ length(tgt)
        @assert nblk ≤ length(val)

        return new{T, I}(nout, nvtx, narc, ncol, nrow, nblk, xsrc, xcol, xrow, xblk, tgt, val)
    end
end

const AdjOrTransSparseBlockMatrix{T, I} = AdjOrTrans{T, SparseBlockMatrix{T, I}}
const MaybeAdjOrTransSparseBlockMatrix{T, I} = Union{SparseBlockMatrix{T, I}, AdjOrTransSparseBlockMatrix{T, I}}
const HermOrSymSparseBlockMatrix{T, I} = HermOrSym{T, SparseBlockMatrix{T, I}}

function SparseBlockMatrix{T, I}(nout::Integer, nvtx::Integer, narc::Integer, ncol::Integer, nrow::Integer, nblk::Integer) where {T, I}
    xsrc = FVector{I}(undef, nvtx + one(I))
    xcol = FVector{I}(undef, nvtx + one(I))
    xrow = FVector{I}(undef, nout + one(I))
    xblk = FVector{I}(undef, narc + one(I))
    tgt = FVector{I}(undef, narc)
    val = FVector{T}(undef, nblk)

    return SparseBlockMatrix{T, I}(nout, nvtx, narc, ncol, nrow, nblk, xsrc, xcol, xrow, xblk, tgt, val)
end

function SparseBlockMatrix(A::SparseMatrixCSC{Mat, I}) where {T, I, Mat <: AbstractMatrix{T}}
    return SparseBlockMatrix{T, I}(A)
end

function SparseBlockMatrix{T, I}(A::SparseMatrixCSC{Mat, I}) where {T, I, Mat <: AbstractMatrix{T}}
    nout = convert(I, size(A, 1))
    nvtx = convert(I, size(A, 2))
    narc = convert(I, nnz(A))
    nblk = convert(I, sum(length, nonzeros(A)))

    xsrc = FVector{I}(undef, nvtx + one(I))
    xcol = FVector{I}(undef, nvtx + one(I)); fill!(xcol, zero(I))
    xrow = FVector{I}(undef, nout + one(I)); fill!(xrow, zero(I))
    xblk = FVector{I}(undef, narc + one(I))
     tgt = FVector{I}(undef, narc)
     val = FVector{T}(undef, nblk)

    b = zero(I)

    for v in oneto(nvtx)
        xsrc[v] = A.colptr[v]

        for e in nzrange(A, v)
            u = tgt[e] = rowvals(A)[e]
            Ae         = nonzeros(A)[e]

            xcol[v + one(I)] = size(Ae, 2)
            xrow[u + one(I)] = size(Ae, 1)
            xblk[e]          = b + one(I)

            for jloc in axes(Ae, 2)
                for iloc in axes(Ae, 1)
                    b += one(I); val[b] = Ae[iloc, jloc]
                end
            end
        end
    end

    ncol = nrow = zero(I)

    for v in oneto(nvtx)
        xcol[v] = ncol + one(I)
        ncol += xcol[v + one(I)]
    end

    for u in oneto(nout)
        xrow[u] = nrow + one(I)
        nrow += xrow[u + one(I)]
    end

    xsrc[nvtx + one(I)] = narc + one(I)
    xcol[nvtx + one(I)] = ncol + one(I)
    xrow[nout + one(I)] = nrow + one(I)
    xblk[narc + one(I)] = nblk + one(I)

    return SparseBlockMatrix{T, I}(nout, nvtx, narc, ncol, nrow, nblk, xsrc, xcol, xrow, xblk, tgt, val)
end

function sparseblock(args...)
    return SparseBlockMatrix(sparse(args...))
end

function SparseArrays.sparse(A::SparseBlockMatrix{T, I}) where {T, I}
    nrow = nrows(A)
    ncol = ncols(A)
    nblk = nblks(A)

    colptr = Vector{I}(undef, ncol + one(I))
    rowval = Vector{I}(undef, nblk)
    nzval = Vector{T}(undef, nblk)

    p = zero(I)

    for v in vtxs(A)
        cols = colrange(A, v)

        for j in cols
            colptr[j] = p + one(I)
            jloc = j - first(cols) + one(I)

            for e in srcrange(A, v)
                u = A.tgt[e]
                Ae = block(A, u, v, e)
                rows = rowrange(A, u)

                for i in rows
                    iloc = i - first(rows) + one(I)
                    p += one(I)
                    rowval[p] = i
                    nzval[p] = Ae[iloc, jloc]
                end
            end
        end
    end

    colptr[ncol + one(I)] = p + one(I)
    return SparseMatrixCSC(nrow, ncol, colptr, rowval, nzval)
end

function Base.Matrix{T}(A::SparseBlockMatrix) where {T}
    nrow = nrows(A)
    ncol = ncols(A)

    M = Matrix{T}(undef, nrow, ncol)
    fill!(M, zero(T))

    for v in vtxs(A)
        cols = colrange(A, v)

        for e in srcrange(A, v)
            u = A.tgt[e]
            Ae = block(A, u, v, e)
            rows = rowrange(A, u)

            for jloc in axes(Ae, 2)
                for iloc in axes(Ae, 1)
                    M[rows[iloc], cols[jloc]] = Ae[iloc, jloc]
                end
            end
        end
    end

    return M
end

function Base.Matrix(A::SparseBlockMatrix{T}) where {T}
    return Matrix{T}(A)
end

############################################################################################
# AbstractArray Interface
############################################################################################

function Base.size(A::SparseBlockMatrix)
    return (convert(Int, nrows(A)), convert(Int, ncols(A)))
end

function Base.getindex(A::SparseBlockMatrix{T}, i::Integer, j::Integer) where {T}
    @boundscheck checkbounds(A, i, j)

    b = findblk(A, i, j)
 
    if !iszero(b)
        x = A.val[b]
    else
        x = zero(T)
    end

    return x
end

function Base.setindex!(A::SparseBlockMatrix, x, i::Integer, j::Integer)
    @boundscheck checkbounds(A, i, j)

    b = findblk(A, i, j)

    if !iszero(b)
        A.val[b] = x
    elseif !iszero(x)
        e = ArgumentError("cannot set out-of-pattern entry ($i, $j) to nonzero value ($x)")
        throw(e)
    end

    return A
end

############################################################################################
# Accessors
############################################################################################

# Get the number of vertices
#
#   nvtxs(A) = | V |
#
function nvtxs(A::SparseBlockMatrix)
    return A.nvtx
end

# Get the number of arcs
#
#   narcs(A) = | E |
#
function narcs(A::SparseBlockMatrix)
    return A.narc
end

# Get the number of out vertices
#
#   nouts(A) = | U |
#
function nouts(A::SparseBlockMatrix)
    return A.nout
end

# Get the number of columns
#
#   ncols(A) = | C |
#
function ncols(A::SparseBlockMatrix)
    return A.ncol
end

# Get the number of rows
#
#   nrows(A) = | R |
#
function nrows(A::SparseBlockMatrix)
    return A.nrow
end

# Get the number of nonzero indices
#
#   nblks(A) = | B |
#
function nblks(A::SparseBlockMatrix)
    return A.nblk
end

# Get the set of vertices
#
#   vtxs(A) = V
#
function vtxs(A::SparseBlockMatrix)
    return oneto(nvtxs(A))
end

# Get the set of arcs
#
#   arcs(A) = E
#
function arcs(A::SparseBlockMatrix)
    return oneto(narcs(A))
end

# Get the set of out vertices
#
#   outs(A) = U
#
function outs(A::SparseBlockMatrix)
    return oneto(nouts(A))
end

# Get the set of columns
#
#   cols(A) = C
#
function cols(A::SparseBlockMatrix)
    return oneto(ncols(A))
end

# Get the set of rows
#
#   rows(A) = R
#
function rows(A::SparseBlockMatrix)
    return oneto(nrows(A))
end

# Get the set of nonzero indices
#
#   blks(A) = B
#
function blks(A::SparseBlockMatrix)
    return oneto(nblks(A))
end

# Get the number
#
#   narcs(A, v) = | srcrange(A, v) |
#
# of arcs incident to a vertex v ∈ V.
@propagate_inbounds function narcs(A::SparseBlockMatrix, v::J) where {J}
    @boundscheck checkbounds(vtxs(A), v)
    estrt = @inbounds A.xsrc[v]
    estop = @inbounds A.xsrc[v + one(J)]
    return estop - estrt
end

# Get the number
#
#   ncols(A, v) = | colrange(A, v) |
#
# of columns incident to a vertex v ∈ V.
@propagate_inbounds function ncols(A::SparseBlockMatrix, v::J) where {J}
    @boundscheck checkbounds(vtxs(A), v)
    jstrt = @inbounds A.xcol[v]
    jstop = @inbounds A.xcol[v + one(J)]
    return jstop - jstrt
end

# Get the number
#
#   nrows(A, u) = | rowrange(A, u) |
#
# of rows incident to an out vertex u ∈ U.
@propagate_inbounds function nrows(A::SparseBlockMatrix, u::J) where {J}
    @boundscheck checkbounds(outs(A), u)
    istrt = @inbounds A.xrow[u]
    istop = @inbounds A.xrow[u + one(J)]
    return istop - istrt
end

# Get the number
#
#   nblks(A, e) = | blkrange(A, e) |
#
# of nonzero indices incident to an arc e ∈ E.
@propagate_inbounds function nblks(A::SparseBlockMatrix, e::J) where {J}
    @boundscheck checkbounds(arcs(A), e)
    bstrt = @inbounds A.xblk[e]
    bstop = @inbounds A.xblk[e + one(J)]
    return bstop - bstrt
end

# Get the arcs
#
#   srcrange(A, v) ⊆ E
#
# incident to a vertex v ∈ V.
@propagate_inbounds function srcrange(A::SparseBlockMatrix{<:Any, I}, v::J) where {I, J}
    @boundscheck checkbounds(vtxs(A), v)
    estrt = @inbounds A.xsrc[v]
    estop = @inbounds A.xsrc[v + one(J)]
    return estrt:estop - one(I)
end

# Get the columns
#
#   colrange(A, v) ⊆ C
#
# incident to a vertex v ∈ V.
@propagate_inbounds function colrange(A::SparseBlockMatrix{<:Any, I}, v::J) where {I, J}
    @boundscheck checkbounds(vtxs(A), v)
    jstrt = @inbounds A.xcol[v]
    jstop = @inbounds A.xcol[v + one(J)]
    return jstrt:jstop - one(I)
end

# Get the rows
#
#   rowrange(A, u) ⊆ R
#
# incident to an out vertex u ∈ U.
@propagate_inbounds function rowrange(A::SparseBlockMatrix{<:Any, I}, u::J) where {I, J}
    @boundscheck checkbounds(outs(A), u)
    istrt = @inbounds A.xrow[u]
    istop = @inbounds A.xrow[u + one(J)]
    return istrt:istop - one(I)
end

# Get the nonzero indices
#
#   blkrange(A, e) ⊆ B
#
# incident to an arc e ∈ E.
@propagate_inbounds function blkrange(A::SparseBlockMatrix{<:Any, I}, e::J) where {I, J}
    @boundscheck checkbounds(arcs(A), e)
    bstrt = @inbounds A.xblk[e]
    bstop = @inbounds A.xblk[e + one(J)]
    return bstrt:bstop - one(I)
end

# Get the block of A corresponding to the arc e = (v, u):
#
#          v
#   Ae = [ A ] u
#
@propagate_inbounds function block(A::SparseBlockMatrix, u::Integer, v::Integer, e::Integer)
    @boundscheck checkbounds(outs(A), u)
    @boundscheck checkbounds(vtxs(A), v)
    @boundscheck checkbounds(arcs(A), e)
    return @inbounds reshape(view(A.val, blkrange(A, e)), nrows(A, u), ncols(A, v))
end

# Find the nonzero index
#
#   b ∈ B
#
# such that blk[b] = e, where e ∈ E is the arc
#
#   e = (row[i], col[j]).
#
# Returns zero if no such block exists.
function findblk(A::SparseBlockMatrix{<:Any, I}, i::Integer, j::Integer) where {I}
    b = zero(I)

    if i in rows(A) && j in cols(A)
        u = binarysearchlast(A.xrow, i, one(I), nouts(A))
        v = binarysearchlast(A.xcol, j, one(I), nvtxs(A))

        strt = A.xsrc[v]
        stop = A.xsrc[v + one(I)] - one(I)

        e = binarysearchlast(A.tgt, u, strt, stop)

        if strt <= e && A.tgt[e] == u
            iloc = convert(I, i) - A.xrow[u]
            jloc = convert(I, j) - A.xcol[v]

            b = A.xblk[e] + iloc + jloc * nrows(A, u)
        end
    end

    return b
end

############################################################################################
# Validation
############################################################################################

# Ensure that a sparse block matrix is well-formed.
function validate(A::SparseBlockMatrix{T, I}) where {T, I}
    nout = nouts(A)
    nvtx = nvtxs(A)
    narc = narcs(A)
    ncol = ncols(A)
    nrow = nrows(A)
    nblk = nblks(A)

    @assert zero(I) < nout
    @assert zero(I) < nvtx
    @assert zero(I) < narc
    @assert zero(I) < ncol
    @assert zero(I) < nrow
    @assert zero(I) < nblk

    @assert nvtx < length(A.xsrc)
    @assert nvtx < length(A.xcol)
    @assert nout < length(A.xrow)
    @assert narc < length(A.xblk)
    @assert narc ≤ length(A.tgt)
    @assert nblk ≤ length(A.val)

    @assert A.xsrc[one(I)] == one(I)
    @assert A.xcol[one(I)] == one(I)
    @assert A.xrow[one(I)] == one(I)
    @assert A.xblk[one(I)] == one(I)

    @assert A.xsrc[nvtx + one(I)] == narc + one(I)
    @assert A.xcol[nvtx + one(I)] == ncol + one(I)
    @assert A.xrow[nout + one(I)] == nrow + one(I)
    @assert A.xblk[narc + one(I)] == nblk + one(I)

    for u in outs(A)
        @assert A.xrow[u] ≤ A.xrow[u + one(I)]
    end

    for v in vtxs(A)
        @assert A.xsrc[v] ≤ A.xsrc[v + one(I)]
        @assert A.xcol[v] ≤ A.xcol[v + one(I)]
    end

    for e in arcs(A)
        @assert A.xblk[e] ≤ A.xblk[e + one(I)]
    end

    for v in vtxs(A)
        p = zero(I)

        for e in srcrange(A, v)
            u = A.tgt[e]
            @assert p < u ≤ nout
            p = u
        end
    end

    return
end

############################################################################################
# Operators
############################################################################################

# Form a submatrix by specifying subsets
#
#   outmask: U → 2
#   vtxmask: V → 2
#
# of out vertices and vertices to preserve.
function submatrix(A::SparseBlockMatrix{T, I}, outmask::AbstractVector{Bool}, vtxmask::AbstractVector{Bool}) where {T, I}
    nBout = zero(I)
    nBvtx = zero(I)
    nBarc = zero(I)
    nBcol = zero(I)
    nBrow = zero(I)
    nBblk = zero(I)

    for u in outs(A)
        if outmask[u]
            nBout += one(I)
            nBrow += nrows(A, u)
        end
    end

    for v in vtxs(A)
        if vtxmask[v]
            nBvtx += one(I)
            nBcol += ncols(A, v)

            for e in srcrange(A, v)
                u = A.tgt[e]

                if outmask[u]
                    nBarc += one(I)
                    nBblk += nblks(A, e)
                end
            end
        end
    end

    B = SparseBlockMatrix{T, I}(nBout, nBvtx, nBarc, nBcol, nBrow, nBblk)

    return submatrix!(B, A, outmask, vtxmask)
end

function submatrix!(B::SparseBlockMatrix{T, I}, A::SparseBlockMatrix{T, I}, outmask::AbstractVector{Bool}, vtxmask::AbstractVector{Bool}) where {T, I}
    outmap = FVector{I}(undef, nouts(A))

    Bu = zero(I)
    Br = zero(I)

    for Au in outs(A)
        if outmask[Au]
            outmap[Au] = Bu += one(I)
            B.xrow[Bu] = Br +  one(I)
            Br += nrows(A, Au)
        end
    end

    Bv = zero(I)
    Bc = zero(I)
    Be = zero(I)
    Bb = zero(I)

    for Av in vtxs(A)
        if vtxmask[Av]
            Bv += one(I)
            B.xsrc[Bv] = Be + one(I)
            B.xcol[Bv] = Bc + one(I)
            Bc += ncols(A, Av)

            for Ae in srcrange(A, Av)
                Au = A.tgt[Ae]

                if outmask[Au]
                    Be += one(I)
                    B.tgt[Be] = outmap[Au]
                    B.xblk[Be] = Bb + one(I)

                    for Ab in blkrange(A, Ae)
                        Bb += one(I)
                        B.val[Bb] = A.val[Ab]
                    end
                end
            end
        end
    end

    B.xrow[nouts(B) + one(I)] = nrows(B) + one(I)
    B.xsrc[nvtxs(B) + one(I)] = narcs(B) + one(I)
    B.xcol[nvtxs(B) + one(I)] = ncols(B) + one(I)
    B.xblk[narcs(B) + one(I)] = nblks(B) + one(I)

    return B
end

############################################################################################
# Matrix Multiplication
############################################################################################

function LinearAlgebra.mul!(C::AbstractVector, A::MaybeAdjOrTransSparseBlockMatrix{T, I}, B::AbstractVector, α::Number, β::Number) where {T, I}
    @assert length(C) == size(A, 1)
    @assert length(B) == size(A, 2)
    uA, tA = unwrapadj(A)
    return gemv_impl!(C, tA, uA, B, α, β)
end

function LinearAlgebra.mul!(C::AbstractMatrix, A::MaybeAdjOrTransSparseBlockMatrix{T, I}, B::AbstractMatrix, α::Number, β::Number) where {T, I}
    @assert size(C, 1) == size(A, 1)
    @assert size(C, 2) == size(B, 2)
    @assert size(A, 2) == size(B, 1)
    uA, tA = unwrapadj(A)
    uB, tB = unwrapadj(B)
    return gemm_impl!(C, tA, uA, tB, uB, α, β)
end

function LinearAlgebra.mul!(C::AbstractVector, A::HermOrSymSparseBlockMatrix{T, I}, B::AbstractVector, α::Number, β::Number) where {T, I}
    @assert length(C) == size(A, 1)
    @assert length(B) == size(A, 2)
    uA, tA = unwrapsym(A)
    return symv_impl!(C, tA, Symbol(A.uplo), uA, B, α, β)
end

function LinearAlgebra.mul!(C::AbstractMatrix, A::HermOrSymSparseBlockMatrix{T, I}, B::AbstractMatrix, α::Number, β::Number) where {T, I}
    @assert size(C, 1) == size(A, 1)
    @assert size(C, 2) == size(B, 2)
    @assert size(A, 2) == size(B, 1)
    uA, tA = unwrapsym(A)
    uB, tB = unwrapadj(B)
    return symm_impl!(C, tA, Symbol(A.uplo), uA, tB, uB, α, β)
end

function gemv_impl!(C::AbstractVector, tA::Val{TA}, A::SparseBlockMatrix, B::AbstractVector, α::Number, β::Number) where {TA}
    if iszero(β)
        fill!(C, β)
    elseif !isone(β)
        rmul!(C, β)
    end

    @inbounds for v in vtxs(A)
        cols = colrange(A, v)

        if TA !== :N
            V = view(C, cols)
        else
            V = view(B, cols)
        end

        for e in srcrange(A, v)
            u = A.tgt[e]
            rows = rowrange(A, u)

            if TA === :N
                U = view(C, rows)
            else
                U = view(B, rows)
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

function gemm_impl!(C::AbstractMatrix, tA::Val{TA}, A::SparseBlockMatrix, tB::Val{TB}, B::AbstractMatrix, α::Number, β::Number) where {TA, TB}
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

function symv_impl!(C::AbstractVector, tA::Val{TA}, uplo::Symbol, A::SparseBlockMatrix, B::AbstractVector, α::Number, β::Number) where {TA}
    if iszero(β)
        fill!(C, β)
    elseif !isone(β)
        rmul!(C, β)
    end

    @inbounds for v in vtxs(A)
        BV = view(B, colrange(A, v))
        CV = view(C, colrange(A, v))

        for e in srcrange(A, v)
            u = A.tgt[e]
            BU = view(B, rowrange(A, u))
            CU = view(C, rowrange(A, u))
            AE = block(A, u, v, e)

            if u == v
                if TA === :N
                    SE = Symmetric(AE, uplo)
                else
                    SE = Hermitian(AE, uplo)
                end

                mul!(CU, SE, BV, α, true)
            else
                if TA === :N
                    mul!(CV, transpose(AE), BU, α, true)
                else
                    mul!(CV,   adjoint(AE), BU, α, true)
                end

                mul!(CU, AE, BV, α, true)
            end
        end
    end

    return C
end

function symm_impl!(C::AbstractMatrix, tA::Val{TA}, uplo::Symbol, A::SparseBlockMatrix, tB::Val{TB}, B::AbstractMatrix, α::Number, β::Number) where {TA, TB}
    if iszero(β)
        fill!(C, β)
    elseif !isone(β)
        rmul!(C, β)
    end

    @inbounds for v in vtxs(A)
        cols = colrange(A, v)

        if TB === :N
            BV =           view(B, cols, :)
        elseif TB === :T
            BV = transpose(view(B, :, cols))
        else
            BV =   adjoint(view(B, :, cols))
        end

        CV = view(C, cols, :)

        for e in srcrange(A, v)
            u = A.tgt[e]
            rows = rowrange(A, u)

            if TB === :N
                BU =           view(B, rows, :)
            elseif TB === :T
                BU = transpose(view(B, :, rows))
            else
                BU =   adjoint(view(B, :, rows))
            end

            CU = view(C, rows, :)
            AE = block(A, u, v, e)

            if u == v
                if TA === :N
                    SE = Symmetric(AE, uplo)
                else
                    SE = Hermitian(AE, uplo)
                end

                mul!(CU, SE, BV, α, true)
            else
                if TA === :N
                    mul!(CV, transpose(AE), BU, α, true)
                else
                    mul!(CV,   adjoint(AE), BU, α, true)
                end

                mul!(CU, AE, BV, α, true)
            end
        end
    end

    return C
end

############################################################################################
# Printing
############################################################################################

function Base.show(io::IO, A::Mat) where {Mat <: SparseBlockMatrix}
    nrow = nrows(A)
    ncol = ncols(A)
    nblk = nblks(A)
    print(io, "$nrow×$ncol $Mat with $nblk stored entries")
    return
end

function Base.show(io::IO, ::MIME"text/plain", A::Mat) where {Mat <: SparseBlockMatrix}
    nrow = nrows(A)
    ncol = ncols(A)
    nblk = nblks(A)
    println(io, "$nrow×$ncol $Mat with $nblk stored entries:")

    if ncol < 16
        show_matrix_dense(io, A)
    else
        show_matrix_sparse(io, A)
    end

    return
end

function show_matrix_sparse(io::IO, A::SparseBlockMatrix)
    nrow = nrows(A)
    ncol = ncols(A)
    grid, rowscale, colscale = braille_grid(io, nrow, ncol)
    grow, gcol = size(grid)

    for v in vtxs(A)
        for e in srcrange(A, v)
            u = A.tgt[e]

            for i in rowrange(A, u)
                for j in colrange(A, v)
                    setbraille!(grid, i, j, rowscale, colscale)
                end
            end
        end
    end

    for j in 1:gcol - 1
        for i in 1:grow
            print(io, grid[i, j] |> Char)
        end
    end

    for i in 1:grow - 1
        print(io, grid[i, gcol] |> Char)
    end

    return
end

function show_matrix_dense(io::IO, A::SparseBlockMatrix)
    print_matrix(io, A)
    return
end
