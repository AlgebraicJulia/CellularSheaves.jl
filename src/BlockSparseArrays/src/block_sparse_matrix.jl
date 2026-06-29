############################################################################################
# Construction
############################################################################################

"""
    BlockSparseMatrix{T, I} <: AbstractMatrix{T}

A sparse block matrix.
"""
struct BlockSparseMatrix{T, I} <: AbstractMatrix{T}
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
    nbnz::I
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

    function BlockSparseMatrix{T, I}(
            nout::Integer,
            nvtx::Integer,
            narc::Integer,
            ncol::Integer,
            nrow::Integer,
            nbnz::Integer,
            xsrc::AbstractVector,
            xcol::AbstractVector,
            xrow::AbstractVector,
            xblk::AbstractVector,
            tgt::AbstractVector,
            val::AbstractVector,
        ) where {T, I}
        @assert nvtx < length(xsrc)
        @assert nvtx < length(xcol)
        @assert nout < length(xrow)
        @assert narc < length(xblk)
        @assert narc ≤ length(tgt)
        @assert nbnz ≤ length(val)

        return new{T, I}(nout, nvtx, narc, ncol, nrow, nbnz, xsrc, xcol, xrow, xblk, tgt, val)
    end
end

const AdjOrTransBlockSparseMatrix{T, I} = AdjOrTrans{T, BlockSparseMatrix{T, I}}
const  HermOrSymBlockSparseMatrix{T, I} = HermOrSym{T, BlockSparseMatrix{T, I}}

const MaybeAdjOrTransBlockSparseMatrix{T, I} = Union{
              BlockSparseMatrix{T, I},
    AdjOrTransBlockSparseMatrix{T, I},
}

function BlockSparseMatrix{T, I}(nout::Integer, nvtx::Integer, narc::Integer, ncol::Integer, nrow::Integer, nbnz::Integer) where {T, I}
    xsrc = FVector{I}(undef, nvtx + one(I))
    xcol = FVector{I}(undef, nvtx + one(I))
    xrow = FVector{I}(undef, nout + one(I))
    xblk = FVector{I}(undef, narc + one(I))
    tgt = FVector{I}(undef, narc)
    val = FVector{T}(undef, nbnz)

    return BlockSparseMatrix{T, I}(nout, nvtx, narc, ncol, nrow, nbnz, xsrc, xcol, xrow, xblk, tgt, val)
end

function BlockSparseMatrix(A::SparseMatrixCSC{Mat, I}, R::AbstractVector{I}, C::AbstractVector{I}) where {T, I, Mat <: AbstractMatrix{T}}
    return BlockSparseMatrix{T, I}(A, R, C)
end

function BlockSparseMatrix{T, I}(A::SparseMatrixCSC{Mat, I}, R::AbstractVector{I}, C::AbstractVector{I}) where {T, I, Mat <: AbstractMatrix{T}}
    nout = convert(I, size(A, 1))
    nvtx = convert(I, size(A, 2))
    narc = convert(I, nnz(A))

    ncol = zero(I)
    nrow = zero(I)
    nbnz = zero(I)

    for v in oneto(nvtx)
        ncol += C[v]
    end

    for u in oneto(nout)
        nrow += R[u]
    end

    for Ae in nonzeros(A)
        nbnz += convert(I, length(Ae))
    end

    xsrc = FVector{I}(undef, nvtx + one(I))
    xcol = FVector{I}(undef, nvtx + one(I))
    xrow = FVector{I}(undef, nout + one(I))
    xblk = FVector{I}(undef, narc + one(I))
     tgt = FVector{I}(undef, narc)
     val = FVector{T}(undef, nbnz)

    i = zero(I)
    j = zero(I)
    b = zero(I)

    for v in oneto(nvtx)
        xcol[v] = j + one(I); j += C[v]
        xsrc[v] = A.colptr[v]

        for e in nzrange(A, v)
             tgt[e] = rowvals(A)[e]
            xblk[e] = b + one(I)

            for x in nonzeros(A)[e]
                b += one(I); val[b] = x
            end
        end
    end

    for u in oneto(nout)
        xrow[u] = i + one(I); i += R[u]
    end

    xsrc[nvtx + one(I)] = narc + one(I)
    xcol[nvtx + one(I)] = ncol + one(I)
    xrow[nout + one(I)] = nrow + one(I)
    xblk[narc + one(I)] = nbnz + one(I)

    return BlockSparseMatrix{T, I}(nout, nvtx, narc, ncol, nrow, nbnz, xsrc, xcol, xrow, xblk, tgt, val)
end

function blocksparse(I::AbstractVector{K}, J::AbstractVector{K}, V::AbstractVector) where {K}
    A = sparse(I, J, V)

    R = FVector{K}(undef, size(A, 1))
    C = FVector{K}(undef, size(A, 2))

    fill!(R, zero(K))

    for v in axes(A, 2)
        Cv = 0

        for e in nzrange(A, v)
            u =   rowvals(A)[e]
            Ae = nonzeros(A)[e]

            R[u] = size(Ae, 1)
            Cv   = size(Ae, 2)
        end

        C[v] = Cv
    end

    return BlockSparseMatrix(A, R, C)
end

function blocksparse(I::AbstractVector, J::AbstractVector, V::AbstractVector, R::AbstractVector, C::AbstractVector)
    m = length(R)
    n = length(C)
    return BlockSparseMatrix(sparse(I, J, V, m, n), R, C)
end

function Base.similar(A::BlockSparseMatrix{T, I}, ::Type{U}=T, ::Type{J}=I) where {T, U, I, J}
    return BlockSparseMatrix{U, J}(A.nout, A.nvtx, A.narc, A.ncol, A.nrow, A.nbnz)
end

function Base.copyto!(B::BlockSparseMatrix, A::BlockSparseMatrix)
    copyto!(B.xsrc, 1, A.xsrc, 1, A.nvtx + 1)
    copyto!(B.xcol, 1, A.xcol, 1, A.nvtx + 1)
    copyto!(B.xrow, 1, A.xrow, 1, A.nout + 1)
    copyto!(B.xblk, 1, A.xblk, 1, A.narc + 1)
    copyto!(B.tgt,  1, A.tgt,  1, A.narc)
    copyto!(B.val,  1, A.val,  1, A.nbnz)

    return B
end

function Base.copy(A::BlockSparseMatrix)
    return copyto!(similar(A), A)
end

function SparseArrays.sparse(A::BlockSparseMatrix{T, I}) where {T, I}
    nrow = nrows(A)
    ncol = ncols(A)
    nbnz = nbnzs(A)

    colptr = Vector{I}(undef, ncol + one(I))
    rowval = Vector{I}(undef, nbnz)
    nzval = Vector{T}(undef, nbnz)

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

function Base.Matrix{T}(A::BlockSparseMatrix) where {T}
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

function Base.Matrix(A::BlockSparseMatrix{T}) where {T}
    return Matrix{T}(A)
end

############################################################################################
# AbstractArray Interface
############################################################################################

function Base.size(A::BlockSparseMatrix)
    return (convert(Int, nrows(A)), convert(Int, ncols(A)))
end

function Base.getindex(A::BlockSparseMatrix{T}, i::Integer, j::Integer) where {T}
    @boundscheck checkbounds(A, i, j)

    b = findbnz(A, i, j)

    if !iszero(b)
        x = A.val[b]
    else
        x = zero(T)
    end

    return x
end

function Base.setindex!(A::BlockSparseMatrix, x, i::Integer, j::Integer)
    @boundscheck checkbounds(A, i, j)

    b = findbnz(A, i, j)

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
function nvtxs(A::BlockSparseMatrix)
    return A.nvtx
end

# Get the number of arcs
#
#   narcs(A) = | E |
#
function narcs(A::BlockSparseMatrix)
    return A.narc
end

# Get the number of out vertices
#
#   nouts(A) = | U |
#
function nouts(A::BlockSparseMatrix)
    return A.nout
end

# Get the number of columns
#
#   ncols(A) = | C |
#
function ncols(A::BlockSparseMatrix)
    return A.ncol
end

# Get the number of rows
#
#   nrows(A) = | R |
#
function nrows(A::BlockSparseMatrix)
    return A.nrow
end

# Get the number of nonzero indices
#
#   nbnzs(A) = | B |
#
function nbnzs(A::BlockSparseMatrix)
    return A.nbnz
end

# Get the set of vertices
#
#   vtxs(A) = V
#
function vtxs(A::BlockSparseMatrix)
    return oneto(nvtxs(A))
end

# Get the set of arcs
#
#   arcs(A) = E
#
function arcs(A::BlockSparseMatrix)
    return oneto(narcs(A))
end

# Get the set of out vertices
#
#   outs(A) = U
#
function outs(A::BlockSparseMatrix)
    return oneto(nouts(A))
end

# Get the set of columns
#
#   cols(A) = C
#
function cols(A::BlockSparseMatrix)
    return oneto(ncols(A))
end

# Get the set of rows
#
#   rows(A) = R
#
function rows(A::BlockSparseMatrix)
    return oneto(nrows(A))
end

# Get the set of nonzero indices
#
#   bnzs(A) = B
#
function bnzs(A::BlockSparseMatrix)
    return oneto(nbnzs(A))
end

# Get the number
#
#   narcs(A, v) = | srcrange(A, v) |
#
# of arcs incident to a vertex v ∈ V.
@propagate_inbounds function narcs(A::BlockSparseMatrix, v::J) where {J}
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
@propagate_inbounds function ncols(A::BlockSparseMatrix, v::J) where {J}
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
@propagate_inbounds function nrows(A::BlockSparseMatrix, u::J) where {J}
    @boundscheck checkbounds(outs(A), u)
    istrt = @inbounds A.xrow[u]
    istop = @inbounds A.xrow[u + one(J)]
    return istop - istrt
end

# Get the number
#
#   nbnzs(A, e) = | blkrange(A, e) |
#
# of nonzero indices incident to an arc e ∈ E.
@propagate_inbounds function nbnzs(A::BlockSparseMatrix, e::J) where {J}
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
@propagate_inbounds function srcrange(A::BlockSparseMatrix{<:Any, I}, v::J) where {I, J}
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
@propagate_inbounds function colrange(A::BlockSparseMatrix{<:Any, I}, v::J) where {I, J}
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
@propagate_inbounds function rowrange(A::BlockSparseMatrix{<:Any, I}, u::J) where {I, J}
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
@propagate_inbounds function blkrange(A::BlockSparseMatrix{<:Any, I}, e::J) where {I, J}
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
@propagate_inbounds function block(A::BlockSparseMatrix, u::Integer, v::Integer, e::Integer)
    @boundscheck checkbounds(outs(A), u)
    @boundscheck checkbounds(vtxs(A), v)
    @boundscheck checkbounds(arcs(A), e)
    return @inbounds block_impl(A, nrows(A, u), ncols(A, v), A.xblk[e])
end

@propagate_inbounds function block_impl(A::BlockSparseMatrix, nrow::I, ncol::I, bstrt::I) where {I}
    bstop = bstrt + nrow * ncol - one(I)
    @boundscheck checkbounds(bnzs(A), bstrt)
    @boundscheck checkbounds(bnzs(A), bstop)
    return @inbounds reshape(view(A.val, bstrt:bstop), nrow, ncol)
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
function findbnz(A::BlockSparseMatrix{<:Any, I}, i::Integer, j::Integer) where {I}
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
function validate(A::BlockSparseMatrix{T, I}) where {T, I}
    nout = nouts(A)
    nvtx = nvtxs(A)
    narc = narcs(A)
    ncol = ncols(A)
    nrow = nrows(A)
    nbnz = nbnzs(A)

    @assert zero(I) < nout
    @assert zero(I) < nvtx
    @assert zero(I) < narc
    @assert zero(I) < ncol
    @assert zero(I) < nrow
    @assert zero(I) < nbnz

    @assert nvtx < length(A.xsrc)
    @assert nvtx < length(A.xcol)
    @assert nout < length(A.xrow)
    @assert narc < length(A.xblk)
    @assert narc ≤ length(A.tgt)
    @assert nbnz ≤ length(A.val)

    @assert A.xsrc[one(I)] == one(I)
    @assert A.xcol[one(I)] == one(I)
    @assert A.xrow[one(I)] == one(I)
    @assert A.xblk[one(I)] == one(I)

    @assert A.xsrc[nvtx + one(I)] == narc + one(I)
    @assert A.xcol[nvtx + one(I)] == ncol + one(I)
    @assert A.xrow[nout + one(I)] == nrow + one(I)
    @assert A.xblk[narc + one(I)] == nbnz + one(I)

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
# Matrix Multiplication
############################################################################################

function LinearAlgebra.mul!(C::AbstractVector, A::MaybeAdjOrTransBlockSparseMatrix{T, I}, B::AbstractVector, α::Number, β::Number) where {T, I}
    @assert length(C) == size(A, 1)
    @assert length(B) == size(A, 2)
    uA, tA = unwrapadj(A)
    return gemv_impl!(C, tA, uA, B, α, β)
end

function LinearAlgebra.mul!(C::AbstractMatrix, A::MaybeAdjOrTransBlockSparseMatrix{T, I}, B::AbstractMatrix, α::Number, β::Number) where {T, I}
    @assert size(C, 1) == size(A, 1)
    @assert size(C, 2) == size(B, 2)
    @assert size(A, 2) == size(B, 1)
    uA, tA = unwrapadj(A)
    uB, tB = unwrapadj(B)
    return gemm_impl!(C, tA, uA, tB, uB, α, β)
end

function LinearAlgebra.mul!(C::AbstractVector, A::HermOrSymBlockSparseMatrix{T, I}, B::AbstractVector, α::Number, β::Number) where {T, I}
    @assert length(C) == size(A, 1)
    @assert length(B) == size(A, 2)
    uA, tA = unwrapsym(A)

    if A.uplo == 'L'
        uplo = :L
    else
        uplo = :U
    end

    return symv_impl!(C, tA, uplo, uA, B, α, β)
end

function LinearAlgebra.mul!(C::AbstractMatrix, A::HermOrSymBlockSparseMatrix{T, I}, B::AbstractMatrix, α::Number, β::Number) where {T, I}
    @assert size(C, 1) == size(A, 1)
    @assert size(C, 2) == size(B, 2)
    @assert size(A, 2) == size(B, 1)
    uA, tA = unwrapsym(A)
    uB, tB = unwrapadj(B)
    return symm_impl!(C, tA, Symbol(A.uplo), uA, tB, uB, α, β)
end

function LinearAlgebra.mul!(C::BlockSparseMatrix, A::MaybeAdjOrTransBlockSparseMatrix, B::MaybeAdjOrTransBlockSparseMatrix, α::Number, β::Number)
    @assert size(C, 1) == size(A, 1)
    @assert size(C, 2) == size(B, 2)
    @assert size(A, 2) == size(B, 1)
    uA, tA = unwrapadj(A)
    uB, tB = unwrapadj(B)
    return gemm_impl!(C, tA, uA, tB, uB, α, β)
end

function Base.:*(A::MaybeAdjOrTransBlockSparseMatrix, B::MaybeAdjOrTransBlockSparseMatrix)
    An = convert(BlockSparseMatrix, A)
    Bn = convert(BlockSparseMatrix, B)
    T = promote_eltype(A, B)
    C = matprod_dest(An, Bn, T)
    return gemm_impl!(C, Val(:N), An, Val(:N), Bn, one(T), zero(T))
end

############################################################################################
# Operators
############################################################################################

function Base.copy(A::AdjOrTransBlockSparseMatrix)
    return copyto!(similar(A), A)
end

function Base.convert(::Type{BlockSparseMatrix}, A::AdjOrTransBlockSparseMatrix)
    return copy(A)
end

function Base.similar(A::AdjOrTransBlockSparseMatrix{T, I}, ::Type{U}=T, ::Type{J}=I) where {T, U, I, J}
    P = parent(A)
    return BlockSparseMatrix{U, J}(P.nvtx, P.nout, P.narc, P.nrow, P.ncol, P.nbnz)
end

function Base.copyto!(A::BlockSparseMatrix, B::AdjOrTransBlockSparseMatrix)
    uB, tB = unwrapadj(B)
    return halfselectvtxs_impl!(A, uB, tB, vtxs(uB))
end

function LinearAlgebra.adjoint!(A::BlockSparseMatrix, B::BlockSparseMatrix)
    return copyto!(A, adjoint(B))
end

function LinearAlgebra.transpose!(A::BlockSparseMatrix, B::BlockSparseMatrix)
    return copyto!(A, transpose(B))
end

function halfselectvtxs(A::BlockSparseMatrix{T, I}, sel::AbstractVector, tA::Val = Val(:T)) where {T, I}
    nBout = length(sel)
    nBvtx = nouts(A)
    nBcol = nrows(A)

    nBrow = zero(I)
    nBarc = zero(I)
    nBbnz = zero(I)

    @inbounds for Av in sel
        nBrow += ncols(A, Av)

        for Ae in srcrange(A, Av)
            nBarc += one(I)
            nBbnz += nbnzs(A, Ae)
        end
    end

    B = BlockSparseMatrix{T, I}(nBout, nBvtx, nBarc, nBcol, nBrow, nBbnz)
    return halfselectvtxs_impl!(B, A, tA, sel)
end

function halfselectvtxs_impl!(A::BlockSparseMatrix{T, I}, B::BlockSparseMatrix, tB::Val, sel::AbstractVector) where {T, I}
    @assert nouts(A) == length(sel)
    @assert nvtxs(A) == nouts(B)

    @inbounds for Av in vtxs(A)
        A.xsrc[Av + one(I)] = zero(I)
        A.xcol[Av + one(I)] = zero(I)
    end

    @inbounds for Au in outs(A)
        Bv = sel[Au]

        for Be in srcrange(B, Bv)
            Bu = B.tgt[Be]

            if Bu < nouts(B)
                A.xsrc[Bu + two(I)] += one(I)
                A.xcol[Bu + two(I)] += nbnzs(B, Be)
            end
        end
    end

    A.xsrc[one(I)] = Ae = one(I)
    A.xcol[one(I)] = Ab = one(I)

    @inbounds for Av in vtxs(A)
        Ae = A.xsrc[Av + one(I)] += Ae
        Ab = A.xcol[Av + one(I)] += Ab
    end

    @inbounds for Au in outs(A)
        Bv = sel[Au]

        for Be in srcrange(B, Bv)
            Bu = B.tgt[Be]

            Ae = A.xsrc[Bu + one(I)]
            Ab = A.xcol[Bu + one(I)]

            A.tgt[Ae] = Au
            A.xblk[Ae] = Ab

            for Bij in wrapadj(block(B, Bu, Bv, Be), tB)
                A.val[Ab] = Bij
                Ab += one(I)
            end

            A.xsrc[Bu + one(I)] += one(I)
            A.xcol[Bu + one(I)] += nbnzs(B, Be)
        end
    end

    Ar = one(I)

    @inbounds for Au in outs(A)
        Bv = sel[Au]
        A.xrow[Au] = Ar
        Ar += ncols(B, Bv)
    end

    A.xrow[nouts(A) + one(I)] = nrows(A) + one(I)
    A.xblk[narcs(A) + one(I)] = nbnzs(A) + one(I)
    copyto!(A.xcol, 1, B.xrow, 1, nvtxs(A) + one(I))

    return A
end

function Base.fill!(A::BlockSparseMatrix, α::Number)
    fill!(view(A.val, bnzs(A)), α)
    return A
end

function LinearAlgebra.lmul!(α::Number, A::BlockSparseMatrix)
    lmul!(α, view(A.val, bnzs(A)))
    return A
end

function LinearAlgebra.rmul!(A::BlockSparseMatrix, α::Number)
    return lmul!(α, A)
end

function LinearAlgebra.ldiv!(α::Number, A::BlockSparseMatrix)
    ldiv!(α, view(A.val, bnzs(A)))
    return A
end

function LinearAlgebra.rdiv!(A::BlockSparseMatrix, α::Number)
    return ldiv!(α, A)
end

function LinearAlgebra.axpy!(α::Number, J::UniformScaling, A::BlockSparseMatrix{T, I}) where {T, I}
    λ = α * J.λ

    @inbounds for v in vtxs(A)
        cols = colrange(A, v)

        for e in srcrange(A, v)
            u = A.tgt[e]
            rows = rowrange(A, u)
            both = rows ∩ cols

            if !isempty(both)
                Ae = block(A, u, v, e)

                for i in both
                    iloc = i - first(rows) + one(I)
                    jloc = i - first(cols) + one(I)
                    Ae[iloc, jloc] += λ
                end
            end
        end
    end

    return A
end

function LinearAlgebra.axpby!(α::Number, J::UniformScaling, β::Number, A::BlockSparseMatrix)
    if iszero(β)
        fill!(A, β)
    elseif !isone(β)
        rmul!(A, β)
    end

    return axpy!(α, J, A)
end

function LinearAlgebra.norm(A::MaybeAdjOrTransBlockSparseMatrix, p::Real=2)
    uA, tA = unwrapadj(A)
    return norm_impl(uA, p)
end

function LinearAlgebra.norm(A::HermOrSymBlockSparseMatrix, p::Real=2)
    uA, tA = unwrapsym(A)
    return symnorm_impl(tA, Symbol(A.uplo), uA, p)
end

function norm_impl(A::BlockSparseMatrix, p::Real)
    v = view(A.val, bnzs(A))
    return norm(v, p)
end

function symnorm_impl(tA::Val{TA}, uplo::Symbol, A::BlockSparseMatrix{T}, p::Real) where {T, TA}
    out = zero(real(T))

    @inbounds for v in vtxs(A)
        for e in srcrange(A, v)
            u = A.tgt[e]
            AE = block(A, u, v, e)

            if u == v
                if TA === :N
                    SE = Symmetric(AE, uplo)
                else
                    SE = Hermitian(AE, uplo)
                end

                nxt = norm(SE, p)
                out = norm((out, nxt), p)
            else
                nxt = norm(AE, p)
                out = norm((out, nxt, nxt), p)
            end
        end
    end

    return out
end

function selectvtxs(A::BlockSparseMatrix{T, I}, sel::AbstractVector{<:Integer}) where {T, I}
    nBout = nouts(A)
    nBrow = nrows(A)

    nBvtx = zero(I)
    nBcol = zero(I)
    nBarc = zero(I)
    nBbnz = zero(I)

    @inbounds for Av in sel
        nBvtx += one(I)
        nBcol += ncols(A, Av)

        for Ae in srcrange(A, Av)
            nBarc += one(I)
            nBbnz += nbnzs(A, Ae)
        end
    end

    B = BlockSparseMatrix{T, I}(nBout, nBvtx, nBarc, nBcol, nBrow, nBbnz)
    return selectvtxs!(B, A, sel)
end

function selectvtxs!(B::BlockSparseMatrix{T, I}, A::BlockSparseMatrix, sel::AbstractVector{<:Integer}) where {T, I}
    @assert nouts(B) == nouts(A)
    @assert nvtxs(B) == length(sel)

    copyto!(B.xrow, 1, A.xrow, 1, nouts(A) + 1)

    Bc = zero(I)
    Be = zero(I)
    Bb = zero(I)

    @inbounds for Bv in vtxs(B)
        Av = sel[Bv]
        B.xcol[Bv] = Bc + one(I)
        B.xsrc[Bv] = Be + one(I)
        Bc += ncols(A, Av)

        for Ae in srcrange(A, Av)
            Au = A.tgt[Ae]

            Be += one(I)
            B.tgt[Be] = Au
            B.xblk[Be] = Bb + one(I)

            for Ab in blkrange(A, Ae)
                Bb += one(I)
                B.val[Bb] = A.val[Ab]
            end
        end
    end

    B.xcol[nvtxs(B) + one(I)] = Bc + one(I)
    B.xsrc[nvtxs(B) + one(I)] = Be + one(I)
    B.xblk[narcs(B) + one(I)] = Bb + one(I)

    return B
end

function LinearAlgebra.matprod_dest(A::BlockSparseMatrix{<:Any, I}, B::BlockSparseMatrix, ::Type{T}) where {T, I}
    @assert nvtxs(A) == nouts(B)

    nCout = nouts(A)
    nCvtx = nvtxs(B)
    nCcol = ncols(B)
    nCrow = nrows(A)

    marker = FVector{I}(undef, nCout)
    fill!(marker, zero(I))

    nCarc = zero(I)
    nCblk = zero(I)

    @inbounds for Bv in vtxs(B)
        for Be in srcrange(B, Bv)
            Bu = B.tgt[Be]

            for Ae in srcrange(A, Bu)
                Au = A.tgt[Ae]

                if marker[Au] < Bv
                    marker[Au] = Bv
                    nCarc += one(I)
                    nCblk += nrows(A, Au) * ncols(B, Bv)
                end
            end
        end
    end

    C = BlockSparseMatrix{T, I}(nCout, nCvtx, nCarc, nCcol, nCrow, nCblk)

    matprod_dest!(C, A, B)

    return C
end

function matprod_dest!(C::BlockSparseMatrix{<:Any, I}, A::BlockSparseMatrix, B::BlockSparseMatrix) where {I}
    copyto!(C.xrow, 1, A.xrow, 1, nouts(A) + 1)
    copyto!(C.xcol, 1, B.xcol, 1, nvtxs(B) + 1)

    marker = FVector{I}(undef, nouts(C))
    fill!(marker, zero(I))

    Ce = zero(I)

    @inbounds for Bv in vtxs(B)
        C.xsrc[Bv] = Cestrt = Ce + one(I)

        for Be in srcrange(B, Bv)
            Bu = B.tgt[Be]

            for Ae in srcrange(A, Bu)
                Au = A.tgt[Ae]

                if marker[Au] < Bv
                    marker[Au] = Bv
                    Ce += one(I)
                    C.tgt[Ce] = Au
                end
            end
        end

        sort!(view(C.tgt, Cestrt:Ce))
    end

    C.xsrc[nvtxs(B) + one(I)] = Ce + one(I)

    Cb = zero(I)

    @inbounds for Cv in vtxs(C)
        for Ce in srcrange(C, Cv)
            Cu = C.tgt[Ce]
            C.xblk[Ce] = Cb + one(I)
            Cb += nrows(C, Cu) * ncols(C, Cv)
        end
    end

    C.xblk[narcs(C) + one(I)] = Cb + one(I)

    return C
end

############################################################################################
# Printing
############################################################################################

function Base.show(io::IO, A::Mat) where {Mat <: BlockSparseMatrix}
    nrow = nrows(A)
    ncol = ncols(A)
    nbnz = nbnzs(A)
    print(io, "$nrow×$ncol $Mat with $nbnz stored entries")
    return
end

function Base.show(io::IO, ::MIME"text/plain", A::Mat) where {Mat <: BlockSparseMatrix}
    nrow = nrows(A)
    ncol = ncols(A)
    nbnz = nbnzs(A)
    println(io, "$nrow×$ncol $Mat with $nbnz stored entries:")

    if ncol < 16
        show_matrix_dense(io, A)
    else
        show_matrix_sparse(io, A)
    end

    return
end

function show_matrix_sparse(io::IO, A::BlockSparseMatrix)
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

function show_matrix_dense(io::IO, A::BlockSparseMatrix)
    print_matrix(io, A)
    return
end
