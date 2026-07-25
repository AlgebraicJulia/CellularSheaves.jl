abstract type AbstractSolver{T} end

include("ipm.jl")
include("hsd.jl")
include("utils.jl")

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
    flag, s.ρ[] = initkkt!(s.kkt, s.H; α=s.α[])
    return flag
end

############################################################################################
# updateaug! — the augmentation controller (rebuild_spec.md §3)
############################################################################################
#
# TRANSCRIBED from the measured incumbent (docs/augmentation/). These constants are the incumbent's
# mechanics, verbatim in behavior; the rebuild does not resize, add to, or "fix" them. Improvements
# require new measurements first — which is why they are module consts, not settings: changing one is
# a code edit with a rationale, not a config tweak. There are no α_min/α_max rails: the emergency arms
# are the working bounds (descent raises base until the jump fires; climb stops at the secant/refine
# arms), and no measured run ever needed a clamp.
#
const AUG_INC          = 4.0   # climb factor per quiet iter; also the floor of the emergency jump
const AUG_CRAIG_MIN    = 2     # climb gate: don't climb when summed base is already at the absolute floor
const AUG_CRAIG_HI     = 8     # emergency tripwire; ~2.5× the in-valley summed base (2–4 per iteration)
const AUG_CRAIG_TARGET = 2     # jump law  α *= (craig1/craig_target)²   (from base ∝ α^{-1/2})
const AUG_JUMP_MAX     = 1e4   # single-jump cap; a library-capped base solve saturates it too (badlo recoveries)
const AUG_DESC         = 4.0   # refine-arm cut (the band's ÷4, proven across every race)
const AUG_REFINE_HI    = 2     # descend trigger: refinement passes ≥ AUG_REFINE_HI

#
# The controller reads ONE base observable, `nbase` = Σ base CRAIG over this iteration's solves (the
# archived `craig1_total` / history `ncraig1` column). This is the observable the incumbent's constants
# were tuned against — the archived tables' "in-valley base ≈ 3" and "13.8 at low raug" were per-iteration
# sums — so AUG_CRAIG_HI, AUG_CRAIG_MIN, and the jump law's (nbase/target)² all consume it, verbatim.
# (A prose transcription earlier compressed this to a per-solve max; archive-wins per §3, see README.md.)
# For HSD the sum includes the warm-started woodbury solve — it was in `d.solves` and the incumbent's
# constants were earned with it counted; the cold-only exclusion is an unshipped refinement, docs-only.
# `refine_passes` stays the worst (max) refinement-pass count over the solves — the shipped descend arm's
# observable, distinct from the base sum.
#
# Mutate s.α for the next iteration. Priority order:
#
#   1. Emergency-up: if nbase ≥ AUG_CRAIG_HI, the model jump from base ∝ α^{-1/2}:
#      α *= clamp((nbase/AUG_CRAIG_TARGET)², AUG_INC, AUG_JUMP_MAX). A base solve that runs to Krylov's
#      dimension cap (a hang / badly-low α) returns a huge nbase that trips this and saturates AUG_JUMP_MAX
#      — so the old explicit truncation escape is subsumed here, with no CRAIG budget of our own.
#   2. Descend (the refine arm): if worst refinement passes ≥ AUG_REFINE_HI, α /= AUG_DESC. (Count-
#      based, no stall-status clause — the status-triggered variant is the refuted stall guard.)
#   3. Marginal climb (single-strike secant): if all solves quiet (refine == 0) and nbase above the
#      floor, compare the last two iterations — if α was climbed and nbase did not drop, hold (at the
#      floor); otherwise α *= AUG_INC. The secant self-re-arms when the summed base later rises.
#      (KNOWN ISSUE — this one-row lookback ratchets α upward on a degenerate endgame; see
#      examples/e16/AUG_CLIMB_RATCHET.md for the diagnosis and a proposed latched fix, deferred
#      pending a proper cross-suite test harness.)
#
function updateaug!(s::AbstractSolver{T}) where {T}
    n = length(s.hist)

    if n > 0
        # observables derived from the just-pushed row; nbase/npass dispatch on the history type, so this
        # stays solver-agnostic (HSD's nbase sums the woodbury solve, IPM's does not).
        nb = nbase(s.hist, n)
        np = npass(s.hist, n)
        α  = s.α[]

        if nb ≥ AUG_CRAIG_HI
            jump = clamp((nb / AUG_CRAIG_TARGET)^2, AUG_INC, AUG_JUMP_MAX)
            s.α[] = α * T(jump)
        elseif np ≥ AUG_REFINE_HI
            s.α[] = α / T(AUG_DESC)
        elseif np == 0 && nb > AUG_CRAIG_MIN
            # single-strike secant from the last two rows: hold only if α climbed and nbase did not drop.
            climbed = n ≥ 2 && s.hist.α[n] > s.hist.α[n-1]
            held    = climbed && nb ≥ nbase(s.hist, n-1)

            if !held
                s.α[] = α * T(AUG_INC)
            end
        end
    end

    return
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
