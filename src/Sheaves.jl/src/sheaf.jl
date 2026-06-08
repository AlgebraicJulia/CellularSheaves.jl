############################################################################################
# Construction
############################################################################################

"""
    Sheaf{T, I}

A cellular sheaf. Construct one using [`sheaf`](@ref).
"""
struct Sheaf{T, I} <: AbstractGraph{I}
    #
    # We represent cellular sheaves as diagrams D: I → Set:
    #
    #   ┌───────────────────────┐
    #   │       R       C       │
    #   │   row ↓  src  ↓ col   │
    #   │       E   ⇉   V       │
    #   │   blk ↑  tgt          │
    #   │       B               │
    #   │   val ↓               │
    #   │       T               │
    #   └───────────────────────┘
    #
    # where
    #
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
    # Each arc e ∈ E is incident to the rows
    #
    #    xrow[e] ... xrow[e + 1] - 1 ∈ R
    #
    xrow::FVector{I}
    #
    # Each arc e ∈ E is incident to the nonzero indices
    #
    #    xblk[e] ... xblk[e + 1] - 1 ∈ B
    #
    xblk::FVector{I}
    #
    # Each arc e ∈ E targets a vertex
    #
    #    tgt[e] ∈ V
    #
    tgt::FVector{I}
    #
    # The nonzero entries:
    #
    #    val[b] ∈ T
    #
    val::FVector{T}

    function Sheaf{T, I}(
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
        @assert narc < length(xrow)
        @assert narc < length(xblk)
        @assert narc ≤ length(tgt)
        @assert nblk ≤ length(val)

        return new{T, I}(nvtx, narc, ncol, nrow, nblk, xsrc, xcol, xrow, xblk, tgt, val)
    end
end

function Sheaf(A::SparseMatrixCSC{Mat, I}) where {T, I, Mat <: AbstractMatrix{T}}
    return Sheaf{T, I}(A)
end

function Sheaf{T, I}(A::SparseMatrixCSC{Mat, I}) where {T, I, Mat <: AbstractMatrix{T}}
    nvtx = convert(I, size(A, 2))
    narc = convert(I, nnz(A))
    nblk = convert(I, sum(length, nonzeros(A)))

    xsrc = FVector{I}(undef, nvtx + one(I))
    xcol = FVector{I}(undef, nvtx + one(I)); fill!(xcol, zero(I))
    xrow = FVector{I}(undef, narc + one(I))
    xblk = FVector{I}(undef, narc + one(I))
     tgt = FVector{I}(undef, narc)
     val = FVector{T}(undef, nblk)

    r = b = zero(I)

    for v in oneto(nvtx)
        xsrc[v] = A.colptr[v]

        for e in nzrange(A, v)
            u = tgt[e] = rowvals(A)[e]
            Ae         = nonzeros(A)[e]

            xcol[v + one(I)] = size(Ae, 2)
            xrow[e]          = r + one(I)
            xblk[e]          = b + one(I)

            r += convert(I, size(Ae, 1))

            for jloc in axes(Ae, 2)
                for iloc in axes(Ae, 1)
                    b += one(I); val[b] = Ae[iloc, jloc]
                end
            end
        end
    end

    ncol = zero(I)

    for v in oneto(nvtx)
        xcol[v] = ncol + one(I)
        ncol += xcol[v + one(I)]
    end

    nrow = r

    xsrc[nvtx + one(I)] = narc + one(I)
    xcol[nvtx + one(I)] = ncol + one(I)
    xrow[narc + one(I)] = nrow + one(I)
    xblk[narc + one(I)] = nblk + one(I)

    return Sheaf{T, I}(nvtx, narc, ncol, nrow, nblk, xsrc, xcol, xrow, xblk, tgt, val)
end

"""
    sheaf(I, J, V[, n])

Construct a cellular sheaf by providing sequences
`I`, `J`, and `V` of

- `I`: target vertices (integers)
- `J`: source vertices (integers)
- `V`: restriction maps (matrices)

Each triple `(u, v, A) ∈ zip(I, J, V)` represents an
oriented edge ``(u, v)`` and a restriction map ``A``
from ``u`` to ``v``. Each edge ``\\{u, v\\}`` should
appear twice: once as ``(u, v)`` and once as ``(v, u)``.
Duplicates are ignored.

The optional fourth argument `n` specifies the number
of vertices in the sheaf.
"""
function sheaf(I::AbstractVector, J::AbstractVector, V::AbstractVector)
    return Sheaf(sparse(I, J, V))
end

function sheaf(I::AbstractVector, J::AbstractVector, V::AbstractVector, n::Integer)
    return Sheaf(sparse(I, J, V, n, n))
end

############################################################################################
# AbstractGraph Interface
############################################################################################

function Graphs.nv(S::Sheaf)
    return nvtxs(S)
end

function Graphs.ne(S::Sheaf)
    return half(narcs(S))
end

function Graphs.vertices(S::Sheaf)
    return vtxs(S)
end

function Graphs.is_directed(::Type{<:Sheaf})
    return false
end

function Graphs.is_directed(S::Sheaf)
    return false
end

function Graphs.has_vertex(S::Sheaf, v::Integer)
    return v in vtxs(S)
end

function Graphs.outneighbors(S::Sheaf, v::Integer)
    return view(S.tgt, srcrange(S, v))
end

function Graphs.inneighbors(S::Sheaf, v::Integer)
    return Graphs.outneighbors(S, v)
end

function Graphs.has_edge(S::Sheaf{T, I}, u::Integer, v::Integer) where {T, I}
    strt = S.xsrc[u]
    stop = S.xsrc[u + one(I)] - one(I)
    e = binarysearchlast(S.tgt, v, strt, stop)
    return strt ≤ e && S.tgt[e] == v
end

function Graphs.outdegree(S::Sheaf, v::Integer)
    return narcs(S, v)
end

############################################################################################
# Accessors
############################################################################################

# Get the number of vertices
#
#   nvtxs(S) = | V |
#
function nvtxs(S::Sheaf)
    return S.nvtx
end

# Get the number of arcs
#
#   narcs(S) = | E |
#
function narcs(S::Sheaf)
    return S.narc
end

# Get the number of columns
#
#   ncols(S) = | C |
#
function ncols(S::Sheaf)
    return S.ncol
end

# Get the number of rows
#
#   nrows(S) = | R |
#
function nrows(S::Sheaf)
    return S.nrow
end

# Get the number of nonzero indices
#
#   nblks(S) = | B |
#
function nblks(S::Sheaf)
    return S.nblk
end

# Get the set of vertices
#
#   vtxs(S) = V
#
function vtxs(S::Sheaf)
    return oneto(nvtxs(S))
end

# Get the set of arcs
#
#   arcs(S) = E
#
function arcs(S::Sheaf)
    return oneto(narcs(S))
end

# Get the set of columns
#
#   cols(S) = C
#
function cols(S::Sheaf)
    return oneto(ncols(S))
end

# Get the set of rows
#
#   rows(S) = R
#
function rows(S::Sheaf)
    return oneto(nrows(S))
end

# Get the set of nonzero indices
#
#   blks(S) = B
#
function blks(S::Sheaf)
    return oneto(nblks(S))
end

# Get the number
#
#   narcs(S, v) = | srcrange(S, v) |
#
# of arcs incident to a vertex v ∈ V.
@propagate_inbounds function narcs(S::Sheaf, v::J) where {J}
    @boundscheck checkbounds(vtxs(S), v)
    estrt = @inbounds S.xsrc[v]
    estop = @inbounds S.xsrc[v + one(J)]
    return estop - estrt
end

# Get the number
#
#   ncols(S, v) = | colrange(S, v) |
#
# of columns incident to a vertex v ∈ V.
@propagate_inbounds function ncols(S::Sheaf, v::J) where {J}
    @boundscheck checkbounds(vtxs(S), v)
    jstrt = @inbounds S.xcol[v]
    jstop = @inbounds S.xcol[v + one(J)]
    return jstop - jstrt
end

# Get the number
#
#   nrows(S, e) = | rowrange(S, e) |
#
# of rows incident to an arc e ∈ E.
@propagate_inbounds function nrows(S::Sheaf, e::J) where {J}
    @boundscheck checkbounds(arcs(S), e)
    istrt = @inbounds S.xrow[e]
    istop = @inbounds S.xrow[e + one(J)]
    return istop - istrt
end

# Get the number
#
#   nblks(S, e) = | blkrange(S, e) |
#
# of nonzero indices incident to an arc e ∈ E.
@propagate_inbounds function nblks(S::Sheaf, e::J) where {J}
    @boundscheck checkbounds(arcs(S), e)
    bstrt = @inbounds S.xblk[e]
    bstop = @inbounds S.xblk[e + one(J)]
    return bstop - bstrt
end

# Get the arcs
#
#   srcrange(S, v) ⊆ E
#
# incident to a vertex v ∈ V.
@propagate_inbounds function srcrange(S::Sheaf{<:Any, I}, v::J) where {I, J}
    @boundscheck checkbounds(vtxs(S), v)
    estrt = @inbounds S.xsrc[v]
    estop = @inbounds S.xsrc[v + one(J)]
    return estrt:estop - one(I)
end

# Get the columns
#
#   colrange(S, v) ⊆ C
#
# incident to a vertex v ∈ V.
@propagate_inbounds function colrange(S::Sheaf{<:Any, I}, v::J) where {I, J}
    @boundscheck checkbounds(vtxs(S), v)
    jstrt = @inbounds S.xcol[v]
    jstop = @inbounds S.xcol[v + one(J)]
    return jstrt:jstop - one(I)
end

# Get the rows
#
#   rowrange(S, e) ⊆ R
#
# incident to an arc e ∈ E.
@propagate_inbounds function rowrange(S::Sheaf{<:Any, I}, e::J) where {I, J}
    @boundscheck checkbounds(arcs(S), e)
    istrt = @inbounds S.xrow[e]
    istop = @inbounds S.xrow[e + one(J)]
    return istrt:istop - one(I)
end

# Get the nonzero indices
#
#   blkrange(S, e) ⊆ B
#
# incident to an arc e ∈ E.
@propagate_inbounds function blkrange(S::Sheaf{<:Any, I}, e::J) where {I, J}
    @boundscheck checkbounds(arcs(S), e)
    bstrt = @inbounds S.xblk[e]
    bstop = @inbounds S.xblk[e + one(J)]
    return bstrt:bstop - one(I)
end

# Get the restriction map of S corresponding to the arc e = (v, u):
#
#          v
#   Se = [ S ] e
#
@propagate_inbounds function block(S::Sheaf, v::Integer, e::Integer)
    @boundscheck checkbounds(vtxs(S), v)
    @boundscheck checkbounds(arcs(S), e)
    return @inbounds reshape(view(S.val, blkrange(S, e)), nrows(S, e), ncols(S, v))
end

############################################################################################
# Validation
############################################################################################

# Ensure that a sheaf is well-formed.
function validate(S::Sheaf{T, I}) where {T, I}
    nvtx = nvtxs(S)
    narc = narcs(S)
    ncol = ncols(S)
    nrow = nrows(S)
    nblk = nblks(S)

    @assert zero(I) < nvtx
    @assert zero(I) < narc
    @assert zero(I) < ncol
    @assert zero(I) < nrow
    @assert zero(I) < nblk

    @assert nvtx < length(S.xsrc)
    @assert nvtx < length(S.xcol)
    @assert narc < length(S.xrow)
    @assert narc < length(S.xblk)
    @assert narc ≤ length(S.tgt)
    @assert nblk ≤ length(S.val)

    @assert S.xsrc[one(I)] == one(I)
    @assert S.xcol[one(I)] == one(I)
    @assert S.xrow[one(I)] == one(I)
    @assert S.xblk[one(I)] == one(I) 

    @assert S.xsrc[nvtx + one(I)] == narc + one(I)
    @assert S.xcol[nvtx + one(I)] == ncol + one(I)
    @assert S.xrow[narc + one(I)] == nrow + one(I)
    @assert S.xblk[narc + one(I)] == nblk + one(I) 

    for v in vtxs(S)
        @assert S.xsrc[v] ≤ S.xsrc[v + one(I)]
        @assert S.xcol[v] ≤ S.xcol[v + one(I)]
    end

    for e in arcs(S)
        @assert S.xrow[e] ≤ S.xrow[e + one(I)]
        @assert S.xblk[e] ≤ S.xblk[e + one(I)]
    end

    for v in vtxs(S)
        p = zero(I)

        for e in srcrange(S, v)
            u = S.tgt[e]
            @assert p < u ≤ nvtx
            p = u
        end
    end

    return
end

############################################################################################
# Operators
############################################################################################

"""
    coboundary(S::Sheaf)

Form the coboundary map of a sheaf. Returns a sparse
block matrix.
"""
function coboundary(S::Sheaf{T, I}) where {T, I}
    A = coboundary_alloc(S)
    return coboundary!(A, S)
end

function coboundary_alloc(S::Sheaf{T, I}) where {T, I}
    nAvtx = nvtxs(S)
    nAarc = narcs(S)
    nAcol = ncols(S)
    nAblk = nblks(S)

    nAout = half(nAarc)
    nArow = half(nrows(S))

    return SparseBlockMatrix{T, I}(nAout, nAvtx, nAarc, nAcol, nArow, nAblk)
end

function coboundary!(A::SparseBlockMatrix{T, I}, S::Sheaf{T, I}) where {T, I}
    nAout = nouts(A)
    nAvtx = nvtxs(A)
    nAarc = narcs(A)
    nAcol = ncols(A)
    nArow = nrows(A)
    nAblk = nblks(A)

    A.xsrc[one(I)] = one(I)

    for v in vtxs(S)
        A.xsrc[v + one(I)] = S.xsrc[v]
    end

    Au = Ar = Ab = zero(I)

    for v in vtxs(S)
        A.xsrc[v] = S.xsrc[v]
        A.xcol[v] = S.xcol[v]

        for e in srcrange(S, v)
            Su = S.tgt[e]

            A.xblk[e] = Ab + one(I)

            if Su > v
                A.tgt[A.xsrc[Su + one(I)]] = A.tgt[e] = Au += one(I)
                      A.xsrc[Su + one(I)]                  += one(I)

                A.xrow[Au] = Ar + one(I)
                Ar += nrows(S, e)

                α = -one(T)
            else
                α =  one(T)
            end

            for Sb in blkrange(S, e)
                Ab += one(I); A.val[Ab] = α * S.val[Sb]
            end
        end
    end

    A.xsrc[nAvtx + one(I)] = nAarc + one(I)
    A.xcol[nAvtx + one(I)] = nAcol + one(I)
    A.xrow[nAout + one(I)] = nArow + one(I)
    A.xblk[nAarc + one(I)] = nAblk + one(I)

    return A
end

"""
    trilaplacian(S::Sheaf[, uplo::Symbol])

Form one triangle of the Laplacian. Returns a sparse
block matrix.
"""
function trilaplacian(S::Sheaf{T, I}, uplo::Symbol=:U) where {T, I}
    A = trilaplacian_alloc(S)
    return trilaplacian!(A, S, uplo)
end

function trilaplacian_alloc(S::Sheaf{T, I}) where {T, I}
    nAvtx = nvtxs(S)
    nSarc = narcs(S)
    nAcol = ncols(S)
    nArow = ncols(S)

    nAout = nAvtx
    nAarc = nAvtx + half(nSarc)

    nAblk = zero(I)

    for v in vtxs(S)
        nv = ncols(S, v)

        for e in srcrange(S, v)
            u = S.tgt[e]

            if u > v
                nAblk += nv * ncols(S, u)
            end
        end

        nAblk += nv^2
    end

    return SparseBlockMatrix{T, I}(nAout, nAvtx, nAarc, nAcol, nArow, nAblk)
end

function trilaplacian!(A::SparseBlockMatrix{T, I}, S::Sheaf{T, I}, uplo::Symbol) where {T, I}
    nAvtx = nvtxs(A)
    nAarc = narcs(A)
    nAcol = ncols(A)
    nAblk = nblks(A)

    Ae = zero(I)
    Ab = zero(I)

    A.xsrc[one(I)] = one(I)

    for v in vtxs(S)
        A.xcol[v] = S.xcol[v]
        A.xrow[v] = S.xcol[v]

        nv = ncols(S, v)

        if uplo == :L
            Ae += one(I)
            A.xblk[Ae] = Ab + one(I)
            Ab += nv * nv
        end

        A.xsrc[v + one(I)] = Ae + one(I)

        for e in srcrange(S, v)
            u = S.tgt[e]

            if intriangle(u, v, uplo)
                nu = ncols(S, u); Ae += one(I)
                A.xblk[Ae] = Ab + one(I)
                Ab += nu * nv
            end
        end

        if uplo == :U
            Ae += one(I)
            A.xblk[Ae] = Ab + one(I)
            Ab += nv * nv
        end
    end

    for v in vtxs(S)
        for Se in srcrange(S, v)
            u = S.tgt[Se]

            if !intriangle(u, v, uplo)
                A.tgt[A.xsrc[u + one(I)]] = Se
                      A.xsrc[u + one(I)] += one(I)
            end
        end

        if uplo == :U
            A.xsrc[v + one(I)] += one(I)
        end
    end

    A.xsrc[nAvtx + one(I)] = nAarc + one(I)
    A.xcol[nAvtx + one(I)] = nAcol + one(I)
    A.xrow[nAvtx + one(I)] = nAcol + one(I)
    A.xblk[nAarc + one(I)] = nAblk + one(I)

    Ae = zero(I)

    for v in vtxs(S)
        if uplo == :L
            Ae += one(I); A.tgt[Ae] = v
        end

        if uplo == :U
            Ad = A.xsrc[v + one(I)] - one(I)
        else
            Ad = Ae
        end

        D = block(A, v, v, Ad); fill!(D, zero(T))

        for Se in srcrange(S, v)
            u = S.tgt[Se]

            V = block(S, v, Se)
            mul!(D, adjoint(V), V, one(T), one(T))

            if intriangle(u, v, uplo)
                Ae += one(I)
                U = block(S, u, A.tgt[Ae])
                A.tgt[Ae] = u

                W = block(A, u, v, Ae)
                mul!(W, adjoint(U), V, -one(T), zero(T))
            end
        end

        if uplo == :U
            Ae += one(I); A.tgt[Ae] = v
        end
    end

    return A
end

"""
    laplacian(S::Sheaf)

Form the Laplacian of a sheaf. Returns a sparse
block matrix.
"""
function laplacian(S::Sheaf{T, I}) where {T, I}
    A = laplacian_alloc(S)
    return laplacian!(A, S)
end

function laplacian_alloc(S::Sheaf{T, I}) where {T, I}
    nAvtx = nvtxs(S)
    nSarc = narcs(S)
    nAcol = ncols(S)
    nArow = ncols(S)

    nAout = nAvtx
    nAarc = nAvtx + nSarc

    nAblk = zero(I)

    for v in vtxs(S)
        nv = ncols(S, v)

        for e in srcrange(S, v)
            u = S.tgt[e]
            nAblk += nv * ncols(S, u)
        end

        nAblk += nv^2
    end

    return SparseBlockMatrix{T, I}(nAout, nAvtx, nAarc, nAcol, nArow, nAblk)
end

function laplacian!(A::SparseBlockMatrix{T, I}, S::Sheaf{T, I}) where {T, I}
    nAvtx = nvtxs(A)
    nAarc = narcs(A)
    nAcol = ncols(A)
    nAblk = nblks(A)

    Ae = zero(I)
    Ab = zero(I)

    A.xsrc[one(I)] = one(I)

    for v in vtxs(S)
        A.xcol[v] = S.xcol[v]
        A.xrow[v] = S.xcol[v]

        nv = ncols(S, v)

        A.xsrc[v + one(I)] = Ae + one(I)

        flag = true

        for e in srcrange(S, v)
            u = S.tgt[e]

            if flag && u > v
                flag = false

                Ae += one(I)
                A.xblk[Ae] = Ab + one(I)
                Ab += nv * nv
            end

            nu = ncols(S, u); Ae += one(I)
            A.xblk[Ae] = Ab + one(I)
            Ab += nu * nv
        end

        if flag
            Ae += one(I)
            A.xblk[Ae] = Ab + one(I)
            Ab += nv * nv
        end
    end

    for v in vtxs(S)
        for Se in srcrange(S, v)
            u = S.tgt[Se]

            if u > v
                A.tgt[A.xsrc[u + one(I)]] = Se
                      A.xsrc[u + one(I)] += one(I)
            end
        end
    end    

    A.xcol[nAvtx + one(I)] = nAcol + one(I)
    A.xrow[nAvtx + one(I)] = nAcol + one(I)
    A.xblk[nAarc + one(I)] = nAblk + one(I)

    Ae = zero(I)

    for v in vtxs(S)
        A.xsrc[v] = Ae + one(I)

        Ad = A.xsrc[v + one(I)]
        D = block(A, v, v, Ad); fill!(D, zero(T))

        flag = true

        for Se in srcrange(S, v)
            u = S.tgt[Se]

            if flag && u > v
                flag = false
                Ae += one(I); A.tgt[Ae] = v
            end

            Ae += one(I)

            V = block(S, v, Se)
            mul!(D, adjoint(V), V, one(T), one(T))

            if u < v
                Sr = A.tgt[Ae]
                Ar = Sr + u

                U = block(S, u, Sr)
                W = block(A, u, v, Ae)
                mul!(W, adjoint(U), V, -one(T), zero(T))

                X = block(A, v, u, Ar)
                adjoint!(X, W)
            end

            A.tgt[Ae] = u
        end

        if flag
            Ae += one(I); A.tgt[Ae] = v
        end
    end

    A.xsrc[nAvtx + one(I)] = nAarc + one(I)

    return A
end

############################################################################################
# Printing
############################################################################################

function Base.show(io::IO, S::Shf) where {Shf <: Sheaf}
    nvtx = nvtxs(S)
    narc = narcs(S)
    print(io, "{$nvtx, $narc} $Shf")
    return
end

function Base.show(io::IO, ::MIME"text/plain", S::Shf) where {Shf <: Sheaf}
    nvtx = nvtxs(S)
    narc = narcs(S)
    println(io, "{$nvtx, $narc} $Shf:")

    if nvtxs(S) < 16
        show_sheaf_dense(io, S)
    else
        show_sheaf_sparse(io, S)
    end

    return
end

function show_sheaf_dense(io::IO, S::Sheaf)
    nvtx = nvtxs(S)

    M = FMatrix{Int}(undef, nvtx, nvtx); fill!(M, 0)

    for v in vtxs(S)
        M[v, v] = ncols(S, v)

        for e in srcrange(S, v)
            u = S.tgt[e]
            M[u, v] = nrows(S, e)
        end
    end

    print_matrix(io, M)
    return
end

function show_sheaf_sparse(io::IO, S::Sheaf)
    nvtx = nvtxs(S)
    grid, rowscale, colscale = braille_grid(io, nvtx, nvtx)
    grow, gcol = size(grid)

    for v in vtxs(S)
        for e in srcrange(S, v)
            u = S.tgt[e]
            setbraille!(grid, u, v, rowscale, colscale)
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
