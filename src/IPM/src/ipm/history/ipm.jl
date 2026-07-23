const IPMHistoryRow{T} = @NamedTuple{μ::T, step::T, pres::T, dres::T, npred::Int, ncorr::Int, rpred::RefStatus, rcorr::RefStatus}

struct IPMHistory{T} <: AbstractVector{IPMHistoryRow{T}}
    μ::Vector{T}
    step::Vector{T}
    pres::Vector{T}
    dres::Vector{T}
    npred::Vector{Int}
    ncorr::Vector{Int}
    rpred::Vector{RefStatus}
    rcorr::Vector{RefStatus}
end

function IPMHistory{T}() where {T}
    return IPMHistory{T}(T[], T[], T[], T[], Int[], Int[], RefStatus[], RefStatus[])
end

function Base.getindex(hist::IPMHistory, i::Int)
    μ       = hist.μ[i]
    step    = hist.step[i]
    pres    = hist.pres[i]
    dres    = hist.dres[i]
    npred   = hist.npred[i]
    ncorr   = hist.ncorr[i]
    rpred   = hist.rpred[i]
    rcorr   = hist.rcorr[i]
    return (; μ, step, pres, dres, npred, ncorr, rpred, rcorr)
end

function Base.push!(hist::IPMHistory, row::NamedTuple)
    push!(hist.μ,       row.μ)
    push!(hist.step,    row.step)
    push!(hist.pres,    row.pres)
    push!(hist.dres,    row.dres)
    push!(hist.npred,   row.npred)
    push!(hist.ncorr,   row.ncorr)
    push!(hist.rpred,   row.rpred)
    push!(hist.rcorr,   row.rcorr)
    return hist
end

function Base.empty!(hist::IPMHistory)
    empty!(hist.μ)
    empty!(hist.step)
    empty!(hist.pres)
    empty!(hist.dres)
    empty!(hist.npred)
    empty!(hist.ncorr)
    empty!(hist.rpred)
    empty!(hist.rcorr)
    return hist
end

function showtop(io::IO, ::IPMHistory; indent::Integer=0)
    pad = " "^indent
    println(io, pad, "┌──────┬──────────┬──────────┬────────┬───────┬───────┬──────────┐")
    println(io, pad, "│ iter │   pres   │   dres   │  step  │ npred │ ncorr │    μ     │")
    return
end

function showbot(io::IO, ::IPMHistory; indent::Integer=0)
    pad = " "^indent
    println(io, pad, "└──────┴──────────┴──────────┴────────┴───────┴───────┴──────────┘")
    return
end

function showmid(io::IO, ::IPMHistory; indent::Integer=0)
    pad = " "^indent
    println(io, pad, "├──────┼──────────┼──────────┼────────┼───────┼───────┼──────────┤")
    println(io, pad, "│    ⋮ │        ⋮ │        ⋮ │      ⋮ │     ⋮ │     ⋮ │        ⋮ │")
    return
end

function showrow(io::IO, i::Integer, row::IPMHistoryRow; indent::Integer=0)
    pad = " "^indent
    println(io, pad, "├──────┼──────────┼──────────┼────────┼───────┼───────┼──────────┤")
    print(io, pad)
    @printf(io, "│ %4d │ %8.2e │ %8.2e │ %6.4f │ %5d │ %5d │ %8.2e │\n",
            i, row.pres, row.dres, row.step, row.npred, row.ncorr, row.μ)
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

    if length(hist.step) >= window
        τavg = sum(hist.step[end-window+1:end]) / window

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
