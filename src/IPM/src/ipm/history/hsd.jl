const HSDHistoryRow{T} = @NamedTuple{μ::T, step::T, pres::T, dres::T, gap::T, α::T, ρ::T, τ::T, κ::T,
    pbase::Int, prefn::Int, ppass::Int, pstat::KKTStatus,
    cbase::Int, crefn::Int, cpass::Int, cstat::KKTStatus,
    wbase::Int, wrefn::Int, wpass::Int, wstat::KKTStatus,
    pfres::T, ppres::T, pdres::T, cfres::T, cpres::T, cdres::T, wfres::T, wpres::T, wdres::T, αmin::T, αmax::T}

struct HSDHistory{T} <: AbstractVector{HSDHistoryRow{T}}
    μ::Vector{T}
    step::Vector{T}
    pres::Vector{T}
    dres::Vector{T}
    gap::Vector{T}
    α::Vector{T}
    ρ::Vector{T}
    τ::Vector{T}
    κ::Vector{T}
    pbase::Vector{Int}
    prefn::Vector{Int}
    ppass::Vector{Int}
    pstat::Vector{KKTStatus}
    cbase::Vector{Int}
    crefn::Vector{Int}
    cpass::Vector{Int}
    cstat::Vector{KKTStatus}
    wbase::Vector{Int}
    wrefn::Vector{Int}
    wpass::Vector{Int}
    wstat::Vector{KKTStatus}
    pfres::Vector{T}
    ppres::Vector{T}
    pdres::Vector{T}
    cfres::Vector{T}
    cpres::Vector{T}
    cdres::Vector{T}
    wfres::Vector{T}
    wpres::Vector{T}
    wdres::Vector{T}
    αmin::Vector{T}
    αmax::Vector{T}
end

function HSDHistory{T}() where {T}
    return HSDHistory{T}(T[], T[], T[], T[], T[], T[], T[], T[], T[],
        Int[], Int[], Int[], KKTStatus[],
        Int[], Int[], Int[], KKTStatus[],
        Int[], Int[], Int[], KKTStatus[],
        T[], T[], T[], T[], T[], T[], T[], T[], T[], T[], T[])
end

function Base.getindex(hist::HSDHistory, i::Int)
    μ       = hist.μ[i]
    step    = hist.step[i]
    pres    = hist.pres[i]
    dres    = hist.dres[i]
    gap     = hist.gap[i]
    α       = hist.α[i]
    ρ       = hist.ρ[i]
    τ       = hist.τ[i]
    κ       = hist.κ[i]
    pbase   = hist.pbase[i]; prefn = hist.prefn[i]; ppass = hist.ppass[i]; pstat = hist.pstat[i]
    cbase   = hist.cbase[i]; crefn = hist.crefn[i]; cpass = hist.cpass[i]; cstat = hist.cstat[i]
    wbase   = hist.wbase[i]; wrefn = hist.wrefn[i]; wpass = hist.wpass[i]; wstat = hist.wstat[i]
    pfres   = hist.pfres[i]; ppres = hist.ppres[i]; pdres = hist.pdres[i]
    cfres   = hist.cfres[i]; cpres = hist.cpres[i]; cdres = hist.cdres[i]
    wfres   = hist.wfres[i]; wpres = hist.wpres[i]; wdres = hist.wdres[i]
    αmin    = hist.αmin[i];  αmax  = hist.αmax[i]
    return (; μ, step, pres, dres, gap, α, ρ, τ, κ,
        pbase, prefn, ppass, pstat, cbase, crefn, cpass, cstat, wbase, wrefn, wpass, wstat,
        pfres, ppres, pdres, cfres, cpres, cdres, wfres, wpres, wdres, αmin, αmax)
end

function Base.push!(hist::HSDHistory, row::NamedTuple)
    push!(hist.μ,       row.μ)
    push!(hist.step,    row.step)
    push!(hist.pres,    row.pres)
    push!(hist.dres,    row.dres)
    push!(hist.gap,     row.gap)
    push!(hist.α,       row.α)
    push!(hist.ρ,       row.ρ)
    push!(hist.τ,       row.τ)
    push!(hist.κ,       row.κ)
    push!(hist.pbase,   row.pbase); push!(hist.prefn, row.prefn); push!(hist.ppass, row.ppass); push!(hist.pstat, row.pstat)
    push!(hist.cbase,   row.cbase); push!(hist.crefn, row.crefn); push!(hist.cpass, row.cpass); push!(hist.cstat, row.cstat)
    push!(hist.wbase,   row.wbase); push!(hist.wrefn, row.wrefn); push!(hist.wpass, row.wpass); push!(hist.wstat, row.wstat)
    push!(hist.pfres, row.pfres); push!(hist.ppres, row.ppres); push!(hist.pdres, row.pdres)
    push!(hist.cfres, row.cfres); push!(hist.cpres, row.cpres); push!(hist.cdres, row.cdres)
    push!(hist.wfres, row.wfres); push!(hist.wpres, row.wpres); push!(hist.wdres, row.wdres)
    push!(hist.αmin, row.αmin); push!(hist.αmax, row.αmax)
    return hist
end

function Base.empty!(hist::HSDHistory)
    empty!(hist.μ)
    empty!(hist.step)
    empty!(hist.pres)
    empty!(hist.dres)
    empty!(hist.gap)
    empty!(hist.α)
    empty!(hist.ρ)
    empty!(hist.τ)
    empty!(hist.κ)
    empty!(hist.pbase); empty!(hist.prefn); empty!(hist.ppass); empty!(hist.pstat)
    empty!(hist.cbase); empty!(hist.crefn); empty!(hist.cpass); empty!(hist.cstat)
    empty!(hist.wbase); empty!(hist.wrefn); empty!(hist.wpass); empty!(hist.wstat)
    return hist
end

function nbase(hist::HSDHistory, i::Integer)
    return hist.pbase[i] + hist.cbase[i] + hist.wbase[i]
end

function npass(hist::HSDHistory, i::Integer)
    return max(hist.ppass[i], hist.cpass[i], hist.wpass[i])
end

function showtop(io::IO, ::HSDHistory; indent::Integer=0)
    pad = " "^indent
    println(io, pad, "┌──────┬──────────┬──────────┬────────┬───────┬──────────┬──────────┬──────────┬──────────┬──────────┐")
    println(io, pad, "│ iter │   pres   │   dres   │  step  │ solve │    α     │    ρ     │    μ     │    τ     │    κ     │")
    return
end

function showbot(io::IO, ::HSDHistory; indent::Integer=0)
    pad = " "^indent
    println(io, pad, "└──────┴──────────┴──────────┴────────┴───────┴──────────┴──────────┴──────────┴──────────┴──────────┘")
    return
end

function showmid(io::IO, ::HSDHistory; indent::Integer=0)
    pad = " "^indent
    println(io, pad, "├──────┼──────────┼──────────┼────────┼───────┼──────────┼──────────┼──────────┼──────────┼──────────┤")
    println(io, pad, "│    ⋮ │        ⋮ │        ⋮ │      ⋮ │     ⋮ │        ⋮ │        ⋮ │        ⋮ │        ⋮ │        ⋮ │")
    return
end

function showrow(io::IO, i::Integer, row::HSDHistoryRow; indent::Integer=0)
    pad = " "^indent
    println(io, pad, "├──────┼──────────┼──────────┼────────┼───────┼──────────┼──────────┼──────────┼──────────┼──────────┤")
    print(io, pad)
    # solve = total solve (CRAIG) iterations this step (base + refinement, all three solves); α = penalty, ρ = regularization.
    solve = row.pbase + row.prefn + row.cbase + row.crefn + row.wbase + row.wrefn
    @printf(io, "│ %4d │ %8.2e │ %8.2e │ %6.4f │ %5d │ %8.2e │ %8.2e │ %8.2e │ %8.2e │ %8.2e │\n",
            i, row.pres, row.dres, row.step, solve, row.α, row.ρ, row.μ, row.τ, row.κ)
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
