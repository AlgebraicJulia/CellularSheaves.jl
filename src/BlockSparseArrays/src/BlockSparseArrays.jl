module BlockSparseArrays

export BlockSparseMatrix, blocksparse, compress, extmul!, twomul!

using Base: oneto, print_matrix, @propagate_inbounds, @boundscheck, promote_eltype
using CliqueTrees.Multifrontal: ChordalFactorization
using FixedSizeArrays: FixedSizeArrayDefault

const FArray{T, N} = FixedSizeArrayDefault{T, N}
const FMatrix{T} = FArray{T, 2}
const FVector{T} = FArray{T, 1}
const FScalar{T} = FArray{T, 0}

using LinearAlgebra
using LinearAlgebra: AdjOrTrans, HermOrSym, matprod_dest
using SparseArrays
using SIMD: Vec

import CliqueTrees.Multifrontal

include("utils.jl")
include("abstract_linked_lists/abstract_linked_lists.jl")
include("compress.jl")
include("block_sparse_matrix.jl")
include("blas/blas.jl")
include("comp/comp.jl")

function Multifrontal.ChordalFactorization{DIAG, UPLO}(A::BlockSparseMatrix; kw...) where {DIAG, UPLO}
    return ChordalFactorization{DIAG, UPLO}(sparse(A); kw...)
end

end # module BlockSparseArrays
