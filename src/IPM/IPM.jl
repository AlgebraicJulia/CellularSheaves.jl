module IPM

using LinearAlgebra
using LinearAlgebra: chkstride1, BlasFloat, BlasInt, LowerTriangular, Adjoint
using Printf: @sprintf, @printf
using LinearAlgebra.BLAS: @blasfunc, libblastrampoline
using LinearAlgebra.LAPACK: chklapackerror
using Base: require_one_based_indexing, ReshapedArray
using FixedSizeArrays: FixedSizeArrayDefault
using SparseArrays

const FArray{T, N} = FixedSizeArrayDefault{T, N}
const FMatrix{T} = FArray{T, 2}
const FVector{T} = FArray{T, 1}
const FScalar{T} = FArray{T, 0}
const Scalar{T} = Array{T, 0}

const FScalarView{T} = SubArray{T, 0, FVector{T}, Tuple{Int64}, true}
const FVectorView{T} = SubArray{T, 1, FVector{T}, Tuple{UnitRange{Int64}}, true}
const FMatrixView{T} = ReshapedArray{T, 2, FVectorView{T}, Tuple{}}

using CliqueTrees: BipartiteGraph, linegraph
using CliqueTrees.Multifrontal: cholesky!, ChordalTriangular, FChordalTriangular, triangular, fronts, diagblock, offdblock,
                                 DivisionWorkspace, AbstractDivisionWorkspace, FactorizationWorkspace, symbolic, FPermutation
using Krylov: cg!, CgWorkspace, cr!, CrWorkspace
using LinearOperators: LinearOperator
using ..BlockSparseArrays: BlockSparseMatrix, block, colrange, rowrange, srcrange, nvtxs, vtxs, ncols, nrows, nouts, outs, nbnzs, narcs, blocksparse, selectvtxs, halfselectvtxs, rows, cols
using CommonSolve: init, solve!, solve
using Core.Compiler: tmerge

import CommonSolve

include("utils.jl")
include("cone/cone.jl")
include("kkt/kkt.jl")
include("scaling.jl")
include("ipm/ipm.jl")

export IPMProblem, IPMSettings, IPMSolver, IPMResult, IPMHistory, IPMHistoryRow, IPMStatus, OPTIMAL, NEAR_OPTIMAL, STALLED, NUMERICAL_FAILURE, ITERATION_LIMIT
export solve, solve!, step!
export AbstractCone, SemidefiniteCone, PositiveCone, SecondOrderCone, CofreeCone, ExponentialCone
export UzawaSettings

end
