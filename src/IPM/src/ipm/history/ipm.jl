const IPMHistoryRow{T} = @NamedTuple{μ::T, pstep::T, dstep::T, pres::T, dres::T, npred::Int, ncorr::Int, refstat::RefStatus}

struct IPMHistory{T} <: AbstractVector{IPMHistoryRow{T}}
    μ::Vector{T}
    pstep::Vector{T}
    dstep::Vector{T}
    pres::Vector{T}
    dres::Vector{T}
    npred::Vector{Int}
    ncorr::Vector{Int}
    refstat::Vector{RefStatus}
end

function IPMHistory{T}() where {T}
    return IPMHistory{T}(T[], T[], T[], T[], T[], Int[], Int[], RefStatus[])
end

function Base.getindex(hist::IPMHistory, i::Int)
    μ       = hist.μ[i]
    pstep   = hist.pstep[i]
    dstep   = hist.dstep[i]
    pres    = hist.pres[i]
    dres    = hist.dres[i]
    npred   = hist.npred[i]
    ncorr   = hist.ncorr[i]
    refstat = hist.refstat[i]
    return (; μ, pstep, dstep, pres, dres, npred, ncorr, refstat)
end

function Base.push!(hist::IPMHistory, row::NamedTuple)
    push!(hist.μ,       row.μ)
    push!(hist.pstep,   row.pstep)
    push!(hist.dstep,   row.dstep)
    push!(hist.pres,    row.pres)
    push!(hist.dres,    row.dres)
    push!(hist.npred,   row.npred)
    push!(hist.ncorr,   row.ncorr)
    push!(hist.refstat, row.refstat)
    return hist
end

function Base.empty!(hist::IPMHistory)
    empty!(hist.μ)
    empty!(hist.pstep)
    empty!(hist.dstep)
    empty!(hist.pres)
    empty!(hist.dres)
    empty!(hist.npred)
    empty!(hist.ncorr)
    empty!(hist.refstat)
    return hist
end

function showtop(io::IO, ::IPMHistory; indent::Integer=0)
    pad = " "^indent
    println(io, pad, "┌──────┬──────────┬──────────┬────────┬────────┬───────┬───────┬──────────┐")
    println(io, pad, "│ iter │   pres   │   dres   │ pstep  │ dstep  │ npred │ ncorr │    μ     │")
    return
end

function showbot(io::IO, ::IPMHistory; indent::Integer=0)
    pad = " "^indent
    println(io, pad, "└──────┴──────────┴──────────┴────────┴────────┴───────┴───────┴──────────┘")
    return
end

function showmid(io::IO, ::IPMHistory; indent::Integer=0)
    pad = " "^indent
    println(io, pad, "├──────┼──────────┼──────────┼────────┼────────┼───────┼───────┼──────────┤")
    println(io, pad, "│    ⋮ │        ⋮ │        ⋮ │      ⋮ │      ⋮ │     ⋮ │     ⋮ │        ⋮ │")
    return
end

function showrow(io::IO, i::Integer, row::IPMHistoryRow; indent::Integer=0)
    pad = " "^indent
    println(io, pad, "├──────┼──────────┼──────────┼────────┼────────┼───────┼───────┼──────────┤")
    print(io, pad)
    @printf(io, "│ %4d │ %8.2e │ %8.2e │ %6.4f │ %6.4f │ %5d │ %5d │ %8.2e │\n",
            i, row.pres, row.dres, row.pstep, row.dstep, row.npred, row.ncorr, row.μ)
    return
end

function Base.show(io::IO, ::MIME"text/plain", hist::T) where {T <: IPMHistory}
    println(io, T, ":")
    return showhistory(io, hist; indent=2)
end

function isnearoptimal(hist::IPMHistory; feas_tol, gap_tol, near_factor)
    flag = false

    if !isempty(hist)
        μ  = hist.μ[end]
        rp = hist.pres[end]
        rd = hist.dres[end]

        flag = rp < near_factor * feas_tol && rd < near_factor * feas_tol && μ < near_factor * gap_tol
    end

    return flag
end

function isnumfail(hist::IPMHistory; window=3, threshold=1e-6)
    flag = false

    if length(hist.pstep) >= window
        τavg = sum(hist.pstep[end-window+1:end]) / window
        τavg = min(τavg, sum(hist.dstep[end-window+1:end]) / window)

        if τavg <= threshold
            if length(hist.pres) < window + 1
                flag = true
            else
                flag = hist.pres[end] > 0.9 * hist.pres[end - window] || hist.dres[end] > 0.9 * hist.dres[end - window]
            end
        end
    end

    return flag
end
