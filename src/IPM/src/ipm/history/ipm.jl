const IPMHistoryRow{T} = @NamedTuple{μ::T, step::T, pres::T, dres::T, α::T, ρ::T,
    pbase::Int, prefn::Int, ppass::Int, pstat::RefStatus,
    cbase::Int, crefn::Int, cpass::Int, cstat::RefStatus,
    pres0::T, pres1::T, cres0::T, cres1::T, r0_p::T, r0_c::T, craig_p::String, craig_c::String}

struct IPMHistory{T} <: AbstractVector{IPMHistoryRow{T}}
    μ::Vector{T}
    step::Vector{T}
    pres::Vector{T}
    dres::Vector{T}
    α::Vector{T}
    ρ::Vector{T}
    pbase::Vector{Int}
    prefn::Vector{Int}
    ppass::Vector{Int}
    pstat::Vector{RefStatus}
    cbase::Vector{Int}
    crefn::Vector{Int}
    cpass::Vector{Int}
    cstat::Vector{RefStatus}
    pres0::Vector{T}
    pres1::Vector{T}
    cres0::Vector{T}
    cres1::Vector{T}
    r0_p::Vector{T}
    r0_c::Vector{T}
    craig_p::Vector{String}
    craig_c::Vector{String}
end

function IPMHistory{T}() where {T}
    return IPMHistory{T}(T[], T[], T[], T[], T[], T[],
        Int[], Int[], Int[], RefStatus[],
        Int[], Int[], Int[], RefStatus[],
        T[], T[], T[], T[],
        T[], T[],
        String[], String[])
end

function Base.getindex(hist::IPMHistory, i::Int)
    μ       = hist.μ[i]
    step    = hist.step[i]
    pres    = hist.pres[i]
    dres    = hist.dres[i]
    α       = hist.α[i]
    ρ       = hist.ρ[i]
    pbase   = hist.pbase[i]
    prefn   = hist.prefn[i]
    ppass   = hist.ppass[i]
    pstat   = hist.pstat[i]
    cbase   = hist.cbase[i]
    crefn   = hist.crefn[i]
    cpass   = hist.cpass[i]
    cstat   = hist.cstat[i]
    pres0   = hist.pres0[i]
    pres1   = hist.pres1[i]
    cres0   = hist.cres0[i]
    cres1   = hist.cres1[i]
    r0_p    = hist.r0_p[i]
    r0_c    = hist.r0_c[i]
    craig_p = hist.craig_p[i]
    craig_c = hist.craig_c[i]
    return (; μ, step, pres, dres, α, ρ, pbase, prefn, ppass, pstat, cbase, crefn, cpass, cstat, pres0, pres1, cres0, cres1, r0_p, r0_c, craig_p, craig_c)
end

function Base.push!(hist::IPMHistory, row::NamedTuple)
    push!(hist.μ,       row.μ)
    push!(hist.step,    row.step)
    push!(hist.pres,    row.pres)
    push!(hist.dres,    row.dres)
    push!(hist.α,       row.α)
    push!(hist.ρ,       row.ρ)
    push!(hist.pbase,   row.pbase)
    push!(hist.prefn,   row.prefn)
    push!(hist.ppass,   row.ppass)
    push!(hist.pstat,   row.pstat)
    push!(hist.cbase,   row.cbase)
    push!(hist.crefn,   row.crefn)
    push!(hist.cpass,   row.cpass)
    push!(hist.cstat,   row.cstat)
    push!(hist.pres0,   row.pres0)
    push!(hist.pres1,   row.pres1)
    push!(hist.cres0,   row.cres0)
    push!(hist.cres1,   row.cres1)
    push!(hist.r0_p,    row.r0_p)
    push!(hist.r0_c,    row.r0_c)
    push!(hist.craig_p, row.craig_p)
    push!(hist.craig_c, row.craig_c)
    return hist
end

function Base.empty!(hist::IPMHistory)
    empty!(hist.μ)
    empty!(hist.step)
    empty!(hist.pres)
    empty!(hist.dres)
    empty!(hist.α)
    empty!(hist.ρ)
    empty!(hist.pbase)
    empty!(hist.prefn)
    empty!(hist.ppass)
    empty!(hist.pstat)
    empty!(hist.cbase)
    empty!(hist.crefn)
    empty!(hist.cpass)
    empty!(hist.cstat)
    empty!(hist.pres0)
    empty!(hist.pres1)
    empty!(hist.cres0)
    empty!(hist.cres1)
    empty!(hist.r0_p)
    empty!(hist.r0_c)
    empty!(hist.craig_p)
    empty!(hist.craig_c)
    return hist
end

function nbase(hist::IPMHistory, i::Integer)
    return hist.pbase[i] + hist.cbase[i]
end

function npass(hist::IPMHistory, i::Integer)
    return max(hist.ppass[i], hist.cpass[i])
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
    solve = row.pbase + row.prefn + row.cbase + row.crefn
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
