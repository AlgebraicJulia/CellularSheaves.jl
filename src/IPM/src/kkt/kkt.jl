abstract type KKTSolver{T} end

# Exit status of a KKT solve / one refinement step, mirroring CraigStatus.
#   KKT_CONTINUE   — in-progress sentinel (a refinement step corrected and the loop continues)
#   KKT_SOLVED     — the residual met the requested (ptol, ytol)
#   KKT_STAGNATED  — stopped short of tol at the achievable limit: residual at its roundoff floor,
#                    or the refinement stalled (not improving). The direction is as accurate as
#                    arithmetic allows, but the requested tol was NOT reached.
#   KKT_ITMAX      — hit the refinement pass cap still short of tol/floor (truncated, not at the limit)
@enum KKTStatus KKT_CONTINUE KKT_SOLVED KKT_STAGNATED KKT_ITMAX

#
# compute the KKT matrix-vector product
#
#   [ u ] = [ A  -Bᵀ ] [ x ]
#   [ v ]   [ B   0  ] [ y ]
#
function mulkkt!(
        u::AbstractVector{T},
        v::AbstractVector{T},
        A::AbstractMatrix{T},
        B::BlockSparseMatrix{T},
        x::AbstractVector{T},
        y::AbstractVector{T},
    ) where {T}
    mul!(u, Symmetric(A, :L), x)
    mul!(u, B', y, -one(T), one(T))
    mul!(v, B, x)
    return
end

include("craig.jl")
include("uzawa.jl")
include("bordered.jl")
