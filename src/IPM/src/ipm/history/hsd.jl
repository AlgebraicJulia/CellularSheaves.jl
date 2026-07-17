const HSDHistoryRow{T} = @NamedTuple{μ::T, step::T, pres::T, dres::T, npred::Int, ncorr::Int, nwood::Int, τ::T, κ::T, refstat::RefStatus}

struct HSDHistory{T} <: AbstractVector{HSDHistoryRow{T}}
    μ::Vector{T}
    step::Vector{T}
    pres::Vector{T}
    dres::Vector{T}
    npred::Vector{Int}
    ncorr::Vector{Int}
    nwood::Vector{Int}
    τ::Vector{T}
    κ::Vector{T}
    refstat::Vector{RefStatus}
end

function HSDHistory{T}() where {T}
    return HSDHistory{T}(T[], T[], T[], T[], Int[], Int[], Int[], T[], T[], RefStatus[])
end

function Base.getindex(hist::HSDHistory, i::Int)
    μ       = hist.μ[i]
    step    = hist.step[i]
    pres    = hist.pres[i]
    dres    = hist.dres[i]
    npred   = hist.npred[i]
    ncorr   = hist.ncorr[i]
    nwood   = hist.nwood[i]
    τ       = hist.τ[i]
    κ       = hist.κ[i]
    refstat = hist.refstat[i]
    return (; μ, step, pres, dres, npred, ncorr, nwood, τ, κ, refstat)
end

function Base.push!(hist::HSDHistory, row::NamedTuple)
    push!(hist.μ,       row.μ)
    push!(hist.step,    row.step)
    push!(hist.pres,    row.pres)
    push!(hist.dres,    row.dres)
    push!(hist.npred,   row.npred)
    push!(hist.ncorr,   row.ncorr)
    push!(hist.nwood,   row.nwood)
    push!(hist.τ,       row.τ)
    push!(hist.κ,       row.κ)
    push!(hist.refstat, row.refstat)
    return hist
end

function Base.empty!(hist::HSDHistory)
    empty!(hist.μ)
    empty!(hist.step)
    empty!(hist.pres)
    empty!(hist.dres)
    empty!(hist.npred)
    empty!(hist.ncorr)
    empty!(hist.nwood)
    empty!(hist.τ)
    empty!(hist.κ)
    empty!(hist.refstat)
    return hist
end

function showtop(io::IO, ::HSDHistory; indent::Integer=0)
    pad = " "^indent
    println(io, pad, "┌──────┬──────────┬──────────┬────────┬───────┬───────┬──────────┬──────────┬──────────┐")
    println(io, pad, "│ iter │   pres   │   dres   │  step  │ npred │ ncorr │    μ     │    τ     │    κ     │")
    return
end

function showbot(io::IO, ::HSDHistory; indent::Integer=0)
    pad = " "^indent
    println(io, pad, "└──────┴──────────┴──────────┴────────┴───────┴───────┴──────────┴──────────┴──────────┘")
    return
end

function showmid(io::IO, ::HSDHistory; indent::Integer=0)
    pad = " "^indent
    println(io, pad, "├──────┼──────────┼──────────┼────────┼───────┼───────┼──────────┼──────────┼──────────┤")
    println(io, pad, "│    ⋮ │        ⋮ │        ⋮ │      ⋮ │     ⋮ │     ⋮ │        ⋮ │        ⋮ │        ⋮ │")
    return
end

function showrow(io::IO, i::Integer, row::HSDHistoryRow; indent::Integer=0)
    pad = " "^indent
    println(io, pad, "├──────┼──────────┼──────────┼────────┼───────┼───────┼──────────┼──────────┼──────────┤")
    print(io, pad)
    @printf(io, "│ %4d │ %8.2e │ %8.2e │ %6.4f │ %5d │ %5d │ %8.2e │ %8.2e │ %8.2e │\n",
            i, row.pres, row.dres, row.step, row.npred, row.ncorr, row.μ, row.τ, row.κ)
    return
end

function Base.show(io::IO, ::MIME"text/plain", hist::T) where {T <: HSDHistory}
    println(io, T, ":")
    return showhistory(io, hist; indent=2)
end

function isnearoptimal(hist::HSDHistory; feas_tol, gap_tol, near_factor)
    flag = false

    if !isempty(hist)
        τ = hist.τ[end]
        κ = hist.κ[end]
        μ = hist.μ[end]
        rp = hist.pres[end]
        rd = hist.dres[end]
        gap = μ / τ^2

        flag = τ > near_factor * κ && rp < near_factor * feas_tol && rd < near_factor * feas_tol && gap < near_factor * gap_tol
    end

    return flag
end

function isnumfail(hist::HSDHistory; window=3, threshold=1e-6)
    flag = false

    if length(hist.step) >= window
        αavg = sum(hist.step[end-window+1:end]) / window

        if αavg <= threshold
            if length(hist.pres) < window + 1
                flag = true
            else
                flag = hist.pres[end] > 0.9 * hist.pres[end - window] || hist.dres[end] > 0.9 * hist.dres[end - window]
            end
        end
    end

    return flag
end

function isillposed(hist::HSDHistory; tol=1e-10)
    flag = false

    if !isempty(hist)
        flag = max(hist.τ[end], hist.κ[end]) <= tol
    end

    return flag
end

