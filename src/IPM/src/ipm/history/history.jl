include("ipm.jl")
include("hsd.jl")

const AbstractHistory{T} = Union{IPMHistory{T}, HSDHistory{T}}

function Base.size(hist::AbstractHistory)
    return size(hist.μ)
end

# Predictor-only α-controller (the one shipped policy). αmin/αmax are the predictor's min-cost window
# (kktwindow!, stored as pamin/pamax); the trigger is the OBSERVED predictor refinement (ppass):
#   • predictor refined last step (ppass > 0) → we sat above the refinement ceiling, so drop one decade
#     below min(α used, ceiling αmax) — anchor on whichever is lower so a mis-high α still descends;
#   • else → geo-mean the predictor window (even if inverted).
# Observed refinement is the trigger, not the sign of the computed window — robust to a mis-estimated αmax
# (which fails to invert the window when the true ceiling collapses in a single step). Only the predictor:
# corrector refinement is inherent to the endgame stagnation regime and fires at ANY α, so triggering on it
# would make α free-fall into a Woodbury blowup.
function getaug(hist::AbstractHistory{T}) where {T}
    α = hist.α[end]
    αmin = hist.pamin[end]
    αmax = hist.pamax[end]

    if hist.ppass[end] > 0
        anchor = isfinite(αmax) ? min(α, αmax) : α
        return anchor / T(10)
    end
    (isfinite(αmin) && isfinite(αmax)) && return sqrt(αmin) * sqrt(αmax)
    return α
end

function isstalled(hist::AbstractHistory{T}, μ::T, stats::KKTStatus...; window=6, threshold=0.5) where {T}
    #
    # endgame stagnation (fast): a KKT solve is floor-limited (STAGNATED — cannot reach tol) AND μ has
    # stopped improving this step (gain < 1/threshold vs the previous iterate). The requested tol is then
    # below the arithmetic floor, so grinding on wastes iterations. Fires at the endgame ONSET; it cannot
    # trip mid-convergence (the solves stay SOLVED while μ still drops fast). Mirrors the oracle's
    # all-stagnated stop, which the slow μ-window below only reproduces ~`window` iterations later.
    #
    if any(==(KKT_STAGNATED), stats) && !isempty(hist) && μ > threshold * hist.μ[end]
        return true
    end

    n = length(hist); flag = false

    if n >= 2 * window
        floor = √eps(T)

        pstrt = n - 2 * window + 1
        pstop = n - window
        cstrt = n - window + 1
        cstop = n

        μprev = minimum(view(hist.μ, pstrt:pstop))
        μcurr = minimum(view(hist.μ, cstrt:cstop))

        pprev = minimum(view(hist.pres, pstrt:pstop))
        pcurr = minimum(view(hist.pres, cstrt:cstop))

        dprev = minimum(view(hist.dres, pstrt:pstop))
        dcurr = minimum(view(hist.dres, cstrt:cstop))

        μflag = μcurr < threshold * μprev
        pflag = pprev > floor && pcurr < threshold * pprev
        dflag = dprev > floor && dcurr < threshold * dprev

        flag = !μflag && !pflag && !dflag
    end

    return flag
end

function showhistory(io::IO, hist::AbstractHistory; nrows::Int=10, indent::Integer=0)
    showtop(io, hist; indent)

    n = length(hist)

    stop = min(    nrows ÷ 2, n   )
    strt = max(n - nrows ÷ 2, stop) + 1

    for i in 1:stop
        showrow(io, i, hist[i]; indent)
    end

    if stop + 1 < strt
        showmid(io, hist; indent)
    end

    for i in strt:n
        showrow(io, i, hist[i]; indent)
    end

    showbot(io, hist; indent)
    return
end
