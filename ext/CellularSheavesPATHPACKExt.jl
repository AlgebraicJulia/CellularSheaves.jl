# The coarsest-grid solver of the HJB multigrid (Multigrid(coarsest = :pathpack))
# with PathPack's chordal semiring LU. Policy evaluation A V = b with
# A = D - N an M-matrix is the fixpoint V = P V + c, P = D⁻¹N, c = D⁻¹b,
# whose least solution is the closure V = P* c over the (+, ×) semiring:
# PathPack factorizes P* = U* L* and applies it with lmul!. The symbolic
# analysis (elimination order, clique tree) is done once for the pattern of P;
# each policy copies its values in and refactorizes numerically.
module CellularSheavesPATHPACKExt

using CellularSheaves
using LinearAlgebra
using SparseArrays
using PATHPACK
using PATHPACK.CPU: PlusProd, ChordalSLU, TDWorkspace, sgetrs!
import CliqueTrees

const DIH = CellularSheaves.ControlSheaves.DoubleIntegratorHJB

struct PathPackCoarsest{F,T,W}
    factor::F
    trans::T                     # how lmul! applies the factor (from CliqueTrees' unwrap)
    workspace::W                 # the triangular solves' workspace, allocated once
end

function DIH._coarsest_solver(::Val{:pathpack}, P::SparseMatrixCSC)
    F = ChordalSLU(PlusProd(), P)                                   # symbolic analysis, once
    _, trans = CliqueTrees.Multifrontal.unwrap(F)
    return PathPackCoarsest(F, trans, TDWorkspace(F, 1))
end

function DIH._refactor!(s::PathPackCoarsest, P::SparseMatrixCSC)
    copyto!(s.factor, P)                                            # new values, same pattern
    lu!(s.factor)
    return s
end

function DIH._solve!(s::PathPackCoarsest, x::Vector{Float64}, c::Vector{Float64})
    copyto!(x, c)
    sgetrs!(s.factor, Val(:L), s.trans, x, s.workspace)             # x = P* c
    return x
end

end # module
