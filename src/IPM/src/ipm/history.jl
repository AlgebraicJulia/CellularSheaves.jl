const IPMHistoryRow{T} = @NamedTuple{μ::T, pstep::T, dstep::T, pres::T, dres::T, npred::Int, ncorr::Int}

struct IPMHistory{T} <: AbstractVector{IPMHistoryRow{T}}
    μ::Vector{T}
    pstep::Vector{T}
    dstep::Vector{T}
    pres::Vector{T}
    dres::Vector{T}
    npred::Vector{Int}
    ncorr::Vector{Int}
end

function IPMHistory{T}() where {T}
    return IPMHistory{T}(T[], T[], T[], T[], T[], Int[], Int[])
end

function Base.size(hist::IPMHistory)
    return size(hist.μ)
end

function Base.getindex(hist::IPMHistory, i::Int)
    μ     = hist.μ[i]
    pstep = hist.pstep[i]
    dstep = hist.dstep[i]
    pres  = hist.pres[i]
    dres  = hist.dres[i]
    npred = hist.npred[i]
    ncorr = hist.ncorr[i]
    return (; μ, pstep, dstep, pres, dres, npred, ncorr)
end

function Base.push!(hist::IPMHistory, row::NamedTuple)
    push!(hist.μ,     row.μ)
    push!(hist.pstep, row.pstep)
    push!(hist.dstep, row.dstep)
    push!(hist.pres,  row.pres)
    push!(hist.dres,  row.dres)
    push!(hist.npred, row.npred)
    push!(hist.ncorr, row.ncorr)
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
    return hist
end

function showtop(io::IO; indent::Integer=0)
    pad = " "^indent
    println(io, pad, "┌──────┬──────────┬──────────┬──────────┬────────┬────────┬───────┬───────┐")
    println(io, pad, "│ iter │    μ     │   pres   │   dres   │ pstep  │ dstep  │ npred │ ncorr │")
    return
end

function showbot(io::IO; indent::Integer=0)
    pad = " "^indent
    println(io, pad, "└──────┴──────────┴──────────┴──────────┴────────┴────────┴───────┴───────┘")
    return
end

function showmid(io::IO; indent::Integer=0)
    pad = " "^indent
    println(io, pad, "├──────┼──────────┼──────────┼──────────┼────────┼────────┼───────┼───────┤")
    println(io, pad, "│    ⋮ │        ⋮ │        ⋮ │        ⋮ │      ⋮ │      ⋮ │     ⋮ │     ⋮ │")
    return
end

function showhistory(io::IO, hist::IPMHistory; nrows::Int=10, indent::Integer=0)
    showtop(io; indent)

    n = length(hist)

    stop = min(    nrows ÷ 2, n   )
    strt = max(n - nrows ÷ 2, stop) + 1

    for i in 1:stop
        showrow(io, i, hist[i]; indent)
    end

    if stop + 1 < strt
        showmid(io; indent)
    end

    for i in strt:n
        showrow(io, i, hist[i]; indent)
    end

    showbot(io; indent)
    return
end

function showrow(io::IO, i::Integer, row::IPMHistoryRow; indent::Integer=0)
    pad = " "^indent
    println(io, pad, "├──────┼──────────┼──────────┼──────────┼────────┼────────┼───────┼───────┤")
    print(io, pad)
    @printf(io, "│ %4d │ %8.2e │ %8.2e │ %8.2e │ %6.4f │ %6.4f │ %5d │ %5d │\n",
            i, row.μ, row.pres, row.dres, row.pstep, row.dstep, row.npred, row.ncorr)
    return
end

function Base.show(io::IO, ::MIME"text/plain", hist::T) where {T <: IPMHistory}
    println(io, T, ":")
    return showhistory(io, hist; indent=2)
end
