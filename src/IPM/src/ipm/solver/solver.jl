abstract type AbstractSolver{T} end

# ── solve-state snapshot hook (snapshots spec) ──────────────────────────────────────────
# SNAP is nothing in production (the guard below is a single === test — no behavior change).
# When a capture driver arms it (SNAP[] = NamedTuple[]), each base solve records its pre-solve
# input state — the symmetric first block A (= H completed from its lower triangle, the operator
# the solver factors), the constraint block B, and the rhs pair f/g plus warm start y0.
const SNAP = Base.RefValue{Any}(nothing)
function _snap!(role::Symbol, H, B, f, g, y0)
    S = sparse(H); L = tril(S)
    A = L + permutedims(L) - spdiagm(0 => diag(L))     # full symmetric operator, storage-agnostic
    push!(SNAP[]::Vector, (; role, A = A, Bmat = sparse(B),
        f = Vector(f), g = Vector(g), y0 = y0 === nothing ? nothing : Vector(y0)))
    return
end

include("ipm.jl")
include("hsd.jl")
include("utils.jl")
include("oracle.jl")

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

# Median and mid-cluster fraction of the (pre-Q) barrier-Hessian diagonal over cone-constrained
# coordinates. Free-cone columns contribute zero barrier curvature and are excluded. Must be called
# AFTER scale!(s) but BEFORE Q is folded into s.H: the diagonal hⱼ = dⱼ/pⱼ is the complementarity
# (degeneracy) signature, which a nonzero Q would mask. Diagnostic-only (observation, no side effects).
function barrier_hdiag_stats(s::AbstractSolver{T}) where {T}
    hs = T[]
    for v in vtxs(s.B)
        s.K[v] isa CofreeCone && continue
        Hvv = block(s.H, v, v, v)
        @inbounds for k in 1:size(Hvv, 1)
            push!(hs, Hvv[k, k])
        end
    end
    isempty(hs) && return T(NaN), T(NaN)
    sort!(hs)
    n = length(hs)
    med = isodd(n) ? hs[(n + 1) ÷ 2] : (hs[n ÷ 2] + hs[n ÷ 2 + 1]) / 2
    mid = count(h -> h > zero(T) && abs(log10(h)) ≤ one(T), hs) / n
    return med, T(mid)
end

function initkkt!(s::AbstractSolver{T}) where {T}
    flag, s.ρ[] = initkkt!(s.kkt, s.H; α=s.α[])
    return flag
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
