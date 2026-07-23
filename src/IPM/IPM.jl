module IPM

using LinearAlgebra
using LinearAlgebra: chkstride1, BlasFloat, BlasInt, LowerTriangular, Adjoint
using Random: Xoshiro
using Printf: @sprintf, @printf
using LinearAlgebra.BLAS: @blasfunc, libblastrampoline
using LinearAlgebra.LAPACK: chklapackerror
using Base: require_one_based_indexing, ReshapedArray, @propagate_inbounds, promote_eltype, oneto
using FixedSizeArrays: FixedSizeArrayDefault
using SparseArrays
using TimerOutputs

const FArray{T, N} = FixedSizeArrayDefault{T, N}
const FMatrix{T} = FArray{T, 2}
const FVector{T} = FArray{T, 1}
const FScalar{T} = FArray{T, 0}
const Scalar{T} = Array{T, 0}

const FScalarView{T} = SubArray{T, 0, FVector{T}, Tuple{Int64}, true}
const FVectorView{T} = SubArray{T, 1, FVector{T}, Tuple{UnitRange{Int64}}, true}
const FMatrixView{T} = ReshapedArray{T, 2, FVectorView{T}, Tuple{}}

using CliqueTrees: BipartiteGraph, linegraph, EliminationAlgorithm, DEFAULT_ELIMINATION_ALGORITHM
using CliqueTrees.Multifrontal: cholesky!, ChordalTriangular, FChordalTriangular, triangular, fronts, diagblock, offdblock,
                                 DivisionWorkspace, FactorizationWorkspace, symbolic, FPermutation
using Krylov: craig!, CraigWorkspace
using LinearOperators: LinearOperator
using ..BlockSparseArrays: BlockSparseMatrix, block, colrange, rowrange, srcrange, nvtxs, vtxs, ncols, nrows, nouts, outs, nbnzs, narcs, blocksparse, selectvtxs, halfselectvtxs, rows, cols
using CommonSolve: init, solve!, solve
using Core.Compiler: tmerge
import CommonSolve

include("src/utils.jl")
include("src/cone/cone.jl")
include("src/kkt/kkt.jl")
include("src/scaling/scaling.jl")
include("src/ipm/ipm.jl")

export IPMProblem, IPMSettings, IPMSolver, IPMResult, IPMHistory, IPMHistoryRow, IPMStatus
export HSDSettings, HSDSolver, HSDResult, HSDHistory, HSDHistoryRow
export OPTIMAL, NEAR_OPTIMAL, STALLED, NUMERICAL_FAILURE, ITERATION_LIMIT, PRIMAL_INFEASIBLE, DUAL_INFEASIBLE, ILL_POSED
export solve, solve!, step!, init, reinit!
export AbstractCone, SemidefiniteCone, PositiveCone, SecondOrderCone, CofreeCone, ExponentialCone
export UzawaSettings
export print_timers

end
