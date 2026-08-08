include("ipm.jl")
include("hsd.jl")

const AbstractHistory{T} = Union{IPMHistory{T}, HSDHistory{T}}

function Base.size(hist::AbstractHistory)
    return size(hist.μ)
end

function atfloor(hist::AbstractHistory; patience::Int=3)
    n = length(hist)
    return n ≥ patience && all(hist.cstat[i] != REACHED_FORCE for i in n-patience+1:n)
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
