include("ipm.jl")
include("hsd.jl")

const AbstractHistory{T} = Union{IPMHistory{T}, HSDHistory{T}}

function Base.size(hist::AbstractHistory)
    return size(hist.μ)
end

function getaug(hist::AbstractHistory{T}, cap::T, policy::Int = 0) where {T}
    α = hist.α[end]
    αmin = hist.αmin[end]
    αmax = hist.αmax[end]

    if policy == 1
        #
        # Policy 1 + closed-window descent (policies_1_and_2.md §4–§5). αmin/αmax are the aggregated
        # window boundaries stored by step! under this policy:
        #   • window open (a0 ≤ ac) → geometric mean (nothing to optimise inside, maximin midpoint);
        #   • window closed (a0 > ac) → one decade below the ceiling a_c, so α DESCENDS with the
        #     collapsing ceiling rather than freezing. Holding froze e04-hsd two decades too high and
        #     stalled it; the one-decade drop is a stand-in for Policy 2's measured level spacing until
        #     the correction-residual history is wired.
        # If a boundary can't be priced (NaN) there is nothing to descend toward, so hold.
        #
        if isfinite(αmin) && isfinite(αmax)
            α = αmin ≤ αmax ? sqrt(αmin) * sqrt(αmax) : αmax / exp10(one(T))
        end
    elseif isnan(αmin)
        α /= 100
    elseif isnan(αmax)
        α = αmin * exp10(T(1.5))
    else
        α = min(sqrt(αmin) * sqrt(αmax), αmax / exp10(cap))
    end

    return α
end

function isstalled(hist::AbstractHistory{T}; window=6, threshold=0.5) where {T}
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
