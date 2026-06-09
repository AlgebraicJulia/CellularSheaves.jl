module Sheaves

export SparseBlockMatrix, sparseblock, submatrix, nouts, nvtxs, narcs, ncols, nrows, nblks, outs, vtxs, arcs, cols, rows, blks, srcrange, colrange, rowrange, blkrange, block
export Sheaf, SheafEdgeIter, sheaf, coboundary, trilaplacian, laplacian

using Base: oneto, print_matrix, @propagate_inbounds, @boundscheck
using CliqueTrees: FMatrix, FVector, half
using Graphs: AbstractGraph, AbstractEdgeIter, SimpleEdge, ne, has_edge

import CliqueTrees.Multifrontal
import Graphs

using LinearAlgebra
using LinearAlgebra: AdjOrTrans, HermOrSym
using SparseArrays

include("utils.jl")
include("sparse_block_matrix.jl")
include("sheaf.jl")
include("sheaf_edge_iter.jl")

function Multifrontal.ChordalLDLt(A::SparseBlockMatrix; kw...)
    return Multifrontal.ChordalLDLt(sparse(A); kw...)
end

function Multifrontal.ChordalLDLt(A::HermOrSymSparseBlockMatrix; kw...)
    return chordalldlt(unwrapsym(A)..., Symbol(A.uplo); kw...)
end

function chordalldlt(A::SparseBlockMatrix, tA::Val{TA}, uplo::Symbol; kw...) where {TA}
    uB = sparse(A)

    if TA === :N
        B = Symmetric(uB, uplo)
    else
        B = Hermitian(uB, uplo)
    end

    return Multifrontal.ChordalLDLt(B; kw...)
end

end
