const IPMHistoryRow{T} = @NamedTuple{μ::T, step::T, pres::T, dres::T, α::T, ρ::T,
    piter::Int, ppass::Int, pstat::KKTStatus,
    citer::Int, cpass::Int, cstat::KKTStatus,
    pamin::T, pamax::T, camin::T, camax::T}

struct IPMHistory{T} <: AbstractVector{IPMHistoryRow{T}}
    μ::Vector{T}
    step::Vector{T}
    pres::Vector{T}
    dres::Vector{T}
    α::Vector{T}
    ρ::Vector{T}
    piter::Vector{Int}
    ppass::Vector{Int}
    pstat::Vector{KKTStatus}
    citer::Vector{Int}
    cpass::Vector{Int}
    cstat::Vector{KKTStatus}
    pamin::Vector{T}
    pamax::Vector{T}
    camin::Vector{T}
    camax::Vector{T}
end

function IPMHistory{T}() where {T}
    return IPMHistory{T}(T[], T[], T[], T[], T[], T[],
        Int[], Int[], KKTStatus[],
        Int[], Int[], KKTStatus[],
        T[], T[], T[], T[])
end

function Base.getindex(hist::IPMHistory, i::Int)
    μ       = hist.μ[i]
    step    = hist.step[i]
    pres    = hist.pres[i]
    dres    = hist.dres[i]
    α       = hist.α[i]
    ρ       = hist.ρ[i]
    piter   = hist.piter[i]; ppass = hist.ppass[i]; pstat = hist.pstat[i]
    citer   = hist.citer[i]; cpass = hist.cpass[i]; cstat = hist.cstat[i]
    pamin   = hist.pamin[i]; pamax = hist.pamax[i]; camin = hist.camin[i]; camax = hist.camax[i]
    return (; μ, step, pres, dres, α, ρ, piter, ppass, pstat, citer, cpass, cstat,
        pamin, pamax, camin, camax)
end

function Base.push!(hist::IPMHistory, row::NamedTuple)
    push!(hist.μ,       row.μ)
    push!(hist.step,    row.step)
    push!(hist.pres,    row.pres)
    push!(hist.dres,    row.dres)
    push!(hist.α,       row.α)
    push!(hist.ρ,       row.ρ)
    push!(hist.piter,   row.piter); push!(hist.ppass, row.ppass); push!(hist.pstat, row.pstat)
    push!(hist.citer,   row.citer); push!(hist.cpass, row.cpass); push!(hist.cstat, row.cstat)
    push!(hist.pamin, row.pamin); push!(hist.pamax, row.pamax)
    push!(hist.camin, row.camin); push!(hist.camax, row.camax)
    return hist
end

function Base.empty!(hist::IPMHistory)
    empty!(hist.μ)
    empty!(hist.step)
    empty!(hist.pres)
    empty!(hist.dres)
    empty!(hist.α)
    empty!(hist.ρ)
    empty!(hist.piter); empty!(hist.ppass); empty!(hist.pstat)
    empty!(hist.citer); empty!(hist.cpass); empty!(hist.cstat)
    empty!(hist.pamin); empty!(hist.pamax); empty!(hist.camin); empty!(hist.camax)
    return hist
end

function showtop(io::IO, ::IPMHistory; indent::Integer=0)
    pad = " "^indent
    println(io, pad, "┌──────┬──────────┬──────────┬────────┬───────┬──────────┬──────────┬──────────┐")
    println(io, pad, "│ iter │   pres   │   dres   │  step  │ solve │    α     │    ρ     │    μ     │")
    return
end

function showbot(io::IO, ::IPMHistory; indent::Integer=0)
    pad = " "^indent
    println(io, pad, "└──────┴──────────┴──────────┴────────┴───────┴──────────┴──────────┴──────────┘")
    return
end

function showmid(io::IO, ::IPMHistory; indent::Integer=0)
    pad = " "^indent
    println(io, pad, "├──────┼──────────┼──────────┼────────┼───────┼──────────┼──────────┼──────────┤")
    println(io, pad, "│    ⋮ │        ⋮ │        ⋮ │      ⋮ │     ⋮ │        ⋮ │        ⋮ │        ⋮ │")
    return
end

function showrow(io::IO, i::Integer, row::IPMHistoryRow; indent::Integer=0)
    pad = " "^indent
    println(io, pad, "├──────┼──────────┼──────────┼────────┼───────┼──────────┼──────────┼──────────┤")
    print(io, pad)
    # solve = total solve (CRAIG) iterations this step (base + refinement, both solves); α = penalty, ρ = regularization.
    solve = row.piter + row.citer
    @printf(io, "│ %4d │ %8.2e │ %8.2e │ %6.4f │ %5d │ %8.2e │ %8.2e │ %8.2e │\n",
            i, row.pres, row.dres, row.step, solve, row.α, row.ρ, row.μ)
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
