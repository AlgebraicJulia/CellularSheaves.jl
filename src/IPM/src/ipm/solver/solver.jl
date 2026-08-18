abstract type AbstractSolver{T} end

const FORCING_FRAC = 0.1
const FORCING_CEIL = 0.3

include("ipm.jl")
include("hsd.jl")
include("utils.jl")

function initreg(nB::T) where {T}
    return eps(T) / 2 * nB^2
end

function isoptimal(s::AbstractSolver, pobj, dobj, pres, dres)
    return isoptimal(pobj, dobj, pres, dres; gap_tol=s.settings.gap_tol, feas_tol=s.settings.feas_tol)
end

function isoptimal(pobj::T, dobj::T, pres::T, dres::T; gap_tol::T, feas_tol::T) where {T}
    return pobj - dobj < gap_tol * max(one(T), min(abs(pobj), abs(dobj))) &&
        pres < feas_tol && dres < feas_tol
end

function setaug!(s::AbstractSolver{T}) where {T}
    if iszero(s.nB[])
        s.δ[] = one(T)
    elseif isempty(s.hist)
        s.δ[] = s.settings.aug_tol * s.nB[]^2 / norm(s.H)
    else
        s.δ[] = getaug(s.hist)
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
    return isstalled(s.hist, s.settings.stall_tol)
end

function isnearoptimal(s::AbstractSolver)
    return isnearoptimal(s.hist; feas_tol=s.settings.feas_tol, gap_tol=s.settings.gap_tol, near_factor=s.settings.near_factor)
end

function initkkt!(s::AbstractSolver{T}) where {T}
    flag, ρ = initkkt!(s.kkt, s.H; δ=s.δ[], rgmin=s.ρ[])
    s.ρ[] = max(s.ρ[], ρ)
    return flag, ρ
end

function scale!(
        cone::AbstractCone,
        v::Integer,
        H::BlockSparseMatrix,
        Q::BlockSparseMatrix,
        caches::Caches,
        p::AbstractVector,
        d::AbstractVector,
        B::BlockSparseMatrix,
        conewrk::ConeWorkspace,
    )
    rv = colrange(B, v)
    cv = cache(caches, v, cone)

    for e in srcrange(H, v)
        if H.tgt[e] == v
            Hv = block(H, v, v, e)
            Qv = block(Q, v, v, e)

            pv = view(p, rv)
            dv = view(d, rv)

            copyto!(Hv, Qv)
            return scale!(Hv, pv, dv, cv, conewrk)
        end
    end

    error()
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

function maxsteps(
        sched::ConeSchedule{T, I},
        K::AbstractVector,
        p::AbstractVector{T},
        d::AbstractVector{T},
        Δp::AbstractVector{T},
        Δd::AbstractVector{T},
        caches::Caches,
        B::BlockSparseMatrix{T},
        step::AbstractVector{T};
        step_frac::T,
    ) where {T, I}
    @inbounds for v in vtxs(B)
        if ncols(B, v) > SMALL_CONE_THRESHOLD
            a, b = maxsteps(K[v], v, p, d, Δp, Δd, caches, B, sched.large; step_frac)
            step[v] = min(a, b)
        end
    end

    @threads for s in oneto(sched.nsmll)
        ws = sched.small[s]

        sstrt = sched.xsmll[s]
        sstop = sched.xsmll[s + one(I)] - one(I)

        @inbounds for v in sstrt:sstop
            if ncols(B, v) <= SMALL_CONE_THRESHOLD
                a, b = maxsteps(K[v], v, p, d, Δp, Δd, caches, B, ws; step_frac)
                step[v] = min(a, b)
            end
        end
    end

    return min(one(T), minimum(step))
end

function maxsteps(s::AbstractSolver{T}, Δp::AbstractVector{T}, Δd::AbstractVector{T}; step_frac::T=one(T)) where {T}
    return maxsteps(s.sched, s.K, s.p, s.d, Δp, Δd, s.caches, s.B, s.wrk.step; step_frac)
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
    return scale!(s.sched, s.K, s.H, s.Q, s.caches, s.p, s.d, s.B, s.wrk.flag)
end

function scale!(
        sched::ConeSchedule{<:Any, I},
        K::AbstractVector,
        H::BlockSparseMatrix,
        Q::BlockSparseMatrix,
        caches::Caches,
        p::AbstractVector,
        d::AbstractVector,
        B::BlockSparseMatrix,
        flags::AbstractVector,
    ) where {I}
    @inbounds for v in vtxs(B)
        if ncols(B, v) > SMALL_CONE_THRESHOLD
            flags[v] = scale!(K[v], v, H, Q, caches, p, d, B, sched.large)
        end
    end

    @threads for s in oneto(sched.nsmll)
        ws = sched.small[s]

        sstrt = sched.xsmll[s]
        sstop = sched.xsmll[s + one(I)] - one(I)

        @inbounds for v in sstrt:sstop
            if ncols(B, v) <= SMALL_CONE_THRESHOLD
                flags[v] = scale!(K[v], v, H, Q, caches, p, d, B, ws)
            end
        end
    end

    return all(flags)
end

function initcorrector!(
        sched::ConeSchedule{<:Any, I},
        K::AbstractVector,
        f::AbstractVector,
        caches::Caches,
        p::AbstractVector,
        d::AbstractVector,
        Δp::AbstractVector,
        Δd::AbstractVector,
        σμ::Real,
        B::BlockSparseMatrix,
    ) where {I}
    @inbounds for v in vtxs(B)
        if ncols(B, v) > SMALL_CONE_THRESHOLD
            initcorrector!(K[v], v, f, caches, p, d, Δp, Δd, σμ, B, sched.large)
        end
    end

    @threads for s in oneto(sched.nsmll)
        ws = sched.small[s]

        sstrt = sched.xsmll[s]
        sstop = sched.xsmll[s + one(I)] - one(I)

        @inbounds for v in sstrt:sstop
            if ncols(B, v) <= SMALL_CONE_THRESHOLD
                initcorrector!(K[v], v, f, caches, p, d, Δp, Δd, σμ, B, ws)
            end
        end
    end

    return
end

function solve_impl!(s::AbstractSolver, io::IO)
    status = CONTINUE; i = 0

    if s.settings.verbose > 0
        showsettings(io, s.settings)
        println(io)
        showtop(io, s.hist)
    end

    while status == CONTINUE
        if i ≥ s.settings.max_iter
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
