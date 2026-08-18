const IPMHistoryRow{T} = @NamedTuple{μ::T, step::T, pres::T, dres::T, pobj::T, dobj::T, δ::T, ρ::T,
    piter::Int, ppass::Int, pstat::KKTStatus,
    citer::Int, cpass::Int, cstat::KKTStatus,
    dmin::T, dmax::T}

struct IPMHistory{T} <: AbstractVector{IPMHistoryRow{T}}
    μ::Vector{T}
    step::Vector{T}
    pres::Vector{T}
    dres::Vector{T}
    pobj::Vector{T}
    dobj::Vector{T}
    δ::Vector{T}
    ρ::Vector{T}
    piter::Vector{Int}
    ppass::Vector{Int}
    pstat::Vector{KKTStatus}
    citer::Vector{Int}
    cpass::Vector{Int}
    cstat::Vector{KKTStatus}
    dmin::Vector{T}
    dmax::Vector{T}
end

function IPMHistory{T}() where {T}
    return IPMHistory{T}(T[], T[], T[], T[], T[], T[], T[], T[],
        Int[], Int[], KKTStatus[],
        Int[], Int[], KKTStatus[],
        T[], T[])
end

function Base.getindex(hist::IPMHistory, i::Int)
    μ       = hist.μ[i]
    step    = hist.step[i]
    pres    = hist.pres[i]
    dres    = hist.dres[i]
    pobj    = hist.pobj[i]
    dobj    = hist.dobj[i]
    δ       = hist.δ[i]
    ρ       = hist.ρ[i]
    piter   = hist.piter[i]; ppass = hist.ppass[i]; pstat = hist.pstat[i]
    citer   = hist.citer[i]; cpass = hist.cpass[i]; cstat = hist.cstat[i]
    dmin   = hist.dmin[i]; dmax = hist.dmax[i]
    return (; μ, step, pres, dres, pobj, dobj, δ, ρ, piter, ppass, pstat, citer, cpass, cstat,
        dmin, dmax)
end

function Base.push!(hist::IPMHistory, row::NamedTuple)
    push!(hist.μ,       row.μ)
    push!(hist.step,    row.step)
    push!(hist.pres,    row.pres)
    push!(hist.dres,    row.dres)
    push!(hist.pobj,    row.pobj)
    push!(hist.dobj,    row.dobj)
    push!(hist.δ,       row.δ)
    push!(hist.ρ,       row.ρ)
    push!(hist.piter,   row.piter); push!(hist.ppass, row.ppass); push!(hist.pstat, row.pstat)
    push!(hist.citer,   row.citer); push!(hist.cpass, row.cpass); push!(hist.cstat, row.cstat)
    push!(hist.dmin, row.dmin); push!(hist.dmax, row.dmax)
    return hist
end

function Base.empty!(hist::IPMHistory)
    empty!(hist.μ)
    empty!(hist.step)
    empty!(hist.pres)
    empty!(hist.dres)
    empty!(hist.pobj)
    empty!(hist.dobj)
    empty!(hist.δ)
    empty!(hist.ρ)
    empty!(hist.piter); empty!(hist.ppass); empty!(hist.pstat)
    empty!(hist.citer); empty!(hist.cpass); empty!(hist.cstat)
    empty!(hist.dmin); empty!(hist.dmax)
    return hist
end

function showtop(io::IO, ::IPMHistory; indent::Integer=0)
    pad = " "^indent
    println(io, pad, "┌──────┬──────────┬──────────┬───────────┬───────────┬────────┬───────┬──────────┬──────────┬──────────┐")
    println(io, pad, "│ iter │   pres   │   dres   │   pobj    │   dobj    │  step  │ solve │    δ     │    ρ     │    μ     │")
    return
end

function showbot(io::IO, ::IPMHistory; indent::Integer=0)
    pad = " "^indent
    println(io, pad, "└──────┴──────────┴──────────┴───────────┴───────────┴────────┴───────┴──────────┴──────────┴──────────┘")
    return
end

function showmid(io::IO, ::IPMHistory; indent::Integer=0)
    pad = " "^indent
    println(io, pad, "├──────┼──────────┼──────────┼───────────┼───────────┼────────┼───────┼──────────┼──────────┼──────────┤")
    println(io, pad, "│    ⋮ │        ⋮ │        ⋮ │         ⋮ │         ⋮ │      ⋮ │     ⋮ │        ⋮ │        ⋮ │        ⋮ │")
    return
end

function showrow(io::IO, i::Integer, row::IPMHistoryRow; indent::Integer=0)
    pad = " "^indent
    println(io, pad, "├──────┼──────────┼──────────┼───────────┼───────────┼────────┼───────┼──────────┼──────────┼──────────┤")
    print(io, pad)
    # solve = total solve (CG) iterations this step (base + refinement, both solves); δ = penalty, ρ = regularization.
    solve = row.piter + row.citer
    @printf(io, "│ %4d │ %8.2e │ %8.2e │ %9.2e │ %9.2e │ %6.4f │ %5d │ %8.2e │ %8.2e │ %8.2e │\n",
            i, row.pres, row.dres, row.pobj, row.dobj, row.step, solve, row.δ, row.ρ, row.μ)
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
