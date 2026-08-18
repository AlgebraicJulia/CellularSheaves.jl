include("ipm.jl")
include("hsd.jl")

const AbstractHistory{T} = Union{IPMHistory{T}, HSDHistory{T}}

function Base.size(hist::AbstractHistory)
    return size(hist.μ)
end

function getaug(hist::AbstractHistory{T}) where {T}
    δ = hist.δ[end]
    dmin = hist.dmin[end]
    dmax = hist.dmax[end]

    if hist.ppass[end] > 0
        if isfinite(dmin)
            δ = max(δ, dmin)
        end

        δ *= 10
    elseif isfinite(dmin) && isfinite(dmax)
        δ = sqrt(dmin) * sqrt(dmax)
    end

    return δ
end

function isstalled(hist::AbstractHistory{T}, tol::T; window::Int=5) where {T}
    flag = false; n = length(hist)

    if n > window
        flag = true

        for i in n - window + 1:n
            μok = hist.μ[i]    < (1 - tol) * hist.μ[i - 1]
            pok = hist.pres[i] < (1 - tol) * hist.pres[i - 1]
            dok = hist.dres[i] < (1 - tol) * hist.dres[i - 1]

            if μok || pok || dok
                flag = false
                break
            end
        end
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
