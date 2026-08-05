abstract type AbstractSolver{T} end

include("ipm.jl")
include("hsd.jl")
include("utils.jl")

# The ρ-shift floor: u·σ²max(B), the size of the perturbation the shift corrects, with σ²max bounded
# above by the Frobenius norm ‖B‖² (the augmentation's nB) and u = eps/2 the unit roundoff. It scales
# with the problem; per rgmin_note the old fixed 1e-9 overshot it by a problem-dependent 6-8 orders.
resolve_rgmin(B) = eps(eltype(B)) / 2 * norm(B)^2

function augmin(α::T, pfres::T, tol::T, state::Bool, pbase::Int, npass::Int, bdg::T) where {T}
    αmin = T(NaN)

    if state && npass <= 1 && isfinite(pfres) && pfres > zero(T)
        αmin = α * pfres / tol / exp10(bdg - max(pbase - 1, 0) / T(20))
    end

    return αmin
end

function augmax(α::T, pdres::T, tol::T, state::Bool, pbase::Int, gap::T) where {T}
    αmax = T(NaN)

    if state && pbase <= 3 && isfinite(pdres) && pdres > zero(T)
        αmax = α * tol / pdres * exp10(gap)
    end

    return αmax
end

# NaN-tolerant max: a role whose residual wasn't recorded (NaN) drops out of the aggregate; all-NaN → NaN.
_nanmax(a::T, b::T) where {T} = isnan(a) ? b : isnan(b) ? a : max(a, b)
_nanmax(a::T, b::T, c::T) where {T} = _nanmax(_nanmax(a, b), c)

# Tier 1 (policies_1_and_2.md §4, Policy 1) window boundaries, mutatis mutandis the aggregation across the
# step's 2–3 solves: ℓ_p = α·max_roles‖r0‖, ℓ_d = max_roles(s0)/α, giving the intersection of the per-solve
# windows [ℓ_p/ε, ε/ℓ_d] — the band an α must sit in to cost one application for EVERY solve. No fudge
# factors, no validity gate; a non-finite/nonpositive input drops that boundary to NaN (getaug then holds).
function augwindow1(α::T, fres::T, dres::T, tol::T) where {T}
    αmin = (isfinite(fres) && fres > zero(T)) ? α * fres / tol : T(NaN)   # ℓ_p / ε
    αmax = (isfinite(dres) && dres > zero(T)) ? α * tol / dres : T(NaN)   # ε / ℓ_d
    return αmin, αmax
end

function setaug!(s::AbstractSolver{T}, cap::T) where {T}
    s.settings.fix_alpha && return   # ORACLE ONLY: leave the caller-injected α in place
    if isempty(s.hist)
        s.α[] = s.settings.aaug + s.settings.raug * norm(Symmetric(s.H, :L)) / s.nB[]^2
    else
        s.α[] = getaug(s.hist, cap, s.settings.policy)
    end

    return
end

function showsolver(io::IO, s::AbstractSolver; indent::Integer=0)
    return showsettings(io, s.settings; indent)
end

function Base.show(io::IO, ::MIME"text/plain", s::T) where {T <: AbstractSolver}
    println(io, T, ":")
    return showsolver(io, s; indent=2)
end

function print_timers(s::AbstractSolver)
    display(s.timers)
end

function isstalled(s::AbstractSolver)
    return isstalled(s.hist)
end

function isnearoptimal(s::AbstractSolver)
    return isnearoptimal(s.hist; feas_tol=s.settings.feas_tol, gap_tol=s.settings.gap_tol, near_factor=s.settings.near_factor)
end

function isnumfail(s::AbstractSolver)
    return isnumfail(s.hist; threshold=s.settings.stall_tol)
end

function atfloor(s::AbstractSolver)
    return atfloor(s.hist; patience=s.settings.floor_patience)
end

function initkkt!(s::AbstractSolver{T}) where {T}
    # s.ρ[] is the running ρ-shift floor (state): it feeds in as the ladder's starting rung and rises
    # monotonically so a later solve never re-climbs from below a shift an earlier one already needed.
    # `ρ` is the shift ACTUALLY applied this solve (0 if the unshifted factorization succeeded) — that
    # is the diagnostic value the caller records in the history, distinct from the floor.
    flag, ρ = initkkt!(s.kkt, s.H; α=s.α[], rgmin=s.ρ[])
    s.ρ[] = max(s.ρ[], ρ)
    return flag, ρ
end

function scale!(
        cone::AbstractCone,
        v::Integer,
        H::BlockSparseMatrix,
        caches::Caches,
        p::AbstractVector,
        d::AbstractVector,
        B::BlockSparseMatrix,
        conewrk::ConeWorkspace,
    )
    r = colrange(B, v)
    Hv = block(H, v, v, v)
    cv = cache(caches, v, cone)
    return scale!(Hv, view(p, r), view(d, r), cv, conewrk)
end

function initcorrector!(
        cone::AbstractCone,
        v::Integer,
        f::AbstractVector,
        caches::Caches,
        p::AbstractVector,
        d::AbstractVector,
        Δp::AbstractVector,
        Δd::AbstractVector,
        σμ::Real,
        B::BlockSparseMatrix,
        conewrk::ConeWorkspace,
    )
    r = colrange(B, v)
    cv = cache(caches, v, cone)
    corr!(view(f, r), view(p, r), view(d, r), view(Δp, r), view(Δd, r), σμ, cv, conewrk)
    return
end

function maxsteps(
        cone::AbstractCone,
        v::Integer,
        p::AbstractVector{T},
        d::AbstractVector{T},
        Δp::AbstractVector{T},
        Δd::AbstractVector{T},
        caches::Caches{T},
        B::BlockSparseMatrix{T},
        conewrk::ConeWorkspace{T};
        step_frac::T = one(T),
    ) where {T}
    r = colrange(B, v)
    τp, τd = maxsteps(view(p, r), view(Δp, r), view(d, r), view(Δd, r), cache(caches, v, cone), conewrk)

    if !(cone isa CofreeCone)
        τp *= step_frac
        τd *= step_frac
    end

    return τp, τd
end

#
# compute the Hessian
#
#   f''(w)
#
# of the primal barrier function f
# at the Nestorov-Todd scaling point w
#
# for non-symmetric cones, no such point
# exists, so the Hessian is replaced
# by a Tuncel scaling matrix
#
function scale!(s::AbstractSolver)
    flag = true

    for v in vtxs(s.B)
        flag = flag && scale!(s.K[v], v, s.H, s.caches, s.p, s.d, s.B, s.conewrk)
    end

    return flag
end

function maxsteps(s::AbstractSolver{T}, v::Integer, Δp::AbstractVector{T}, Δd::AbstractVector{T}; step_frac::T=one(T)) where {T}
    return maxsteps(s.K[v], v, s.p, s.d, Δp, Δd, s.caches, s.B, s.conewrk; step_frac)
end

function maxsteps(s::AbstractSolver{T}, Δp::AbstractVector{T}, Δd::AbstractVector{T}; step_frac::T=one(T)) where {T}
    τp = one(T)
    τd = one(T)

    for v in vtxs(s.B)
        τpv, τdv = maxsteps(s, v, Δp, Δd; step_frac)
        τp = min(τp, τpv)
        τd = min(τd, τdv)
    end

    return τp, τd
end

function solve_impl!(s::AbstractSolver, io::IO)
    status = CONTINUE; i = 0

    if s.settings.verbose > 0
        showsettings(io, s.settings)
        println(io)
        showtop(io, s.hist)
    end

    while status == CONTINUE
        if i ≥ s.settings.itmax
            status = ITERATION_LIMIT
        else
            status = step!(s); i += 1

            if s.settings.verbose > 0
                showrow(io, i, s.hist[i])
            end
        end

        if status != CONTINUE
            break
        end
    end

    if s.settings.verbose > 0
        showbot(io, s.hist)
    end

    return status
end

function CommonSolve.solve!(s::AbstractSolver)
    status = @timeit s.timers "solve" solve_impl!(s, stdout)
    return result(s, status)
end
