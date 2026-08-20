const HSDHistoryRow{T} = @NamedTuple{μ::T, step::T, pres::T, dres::T, pobj::T, dobj::T, δ::T, ρ::T, τ::T, κ::T,
    piter::Int, ppass::Int, pstat::KKTStatus,
    citer::Int, cpass::Int, cstat::KKTStatus,
    witer::Int, wpass::Int,
    dmin::T, dmax::T}

struct HSDHistory{T} <: AbstractVector{HSDHistoryRow{T}}
    μ::Vector{T}
    step::Vector{T}
    pres::Vector{T}
    dres::Vector{T}
    pobj::Vector{T}
    dobj::Vector{T}
    δ::Vector{T}
    ρ::Vector{T}
    τ::Vector{T}
    κ::Vector{T}
    piter::Vector{Int}
    ppass::Vector{Int}
    pstat::Vector{KKTStatus}
    citer::Vector{Int}
    cpass::Vector{Int}
    cstat::Vector{KKTStatus}
    witer::Vector{Int}
    wpass::Vector{Int}
    dmin::Vector{T}
    dmax::Vector{T}
end

function HSDHistory{T}() where {T}
    return HSDHistory{T}(T[], T[], T[], T[], T[], T[], T[], T[], T[], T[],
        Int[], Int[], KKTStatus[],
        Int[], Int[], KKTStatus[],
        Int[], Int[],
        T[], T[])
end

function Base.getindex(hist::HSDHistory, i::Int)
    μ       = hist.μ[i]
    step    = hist.step[i]
    pres    = hist.pres[i]
    dres    = hist.dres[i]
    pobj    = hist.pobj[i]
    dobj    = hist.dobj[i]
    δ       = hist.δ[i]
    ρ       = hist.ρ[i]
    τ       = hist.τ[i]
    κ       = hist.κ[i]
    piter   = hist.piter[i]; ppass = hist.ppass[i]; pstat = hist.pstat[i]
    citer   = hist.citer[i]; cpass = hist.cpass[i]; cstat = hist.cstat[i]
    witer   = hist.witer[i]; wpass = hist.wpass[i]
    dmin   = hist.dmin[i]; dmax = hist.dmax[i]
    return (; μ, step, pres, dres, pobj, dobj, δ, ρ, τ, κ,
        piter, ppass, pstat, citer, cpass, cstat, witer, wpass,
        dmin, dmax)
end

function Base.push!(hist::HSDHistory, row::NamedTuple)
    push!(hist.μ,       row.μ)
    push!(hist.step,    row.step)
    push!(hist.pres,    row.pres)
    push!(hist.dres,    row.dres)
    push!(hist.pobj,    row.pobj)
    push!(hist.dobj,    row.dobj)
    push!(hist.δ,       row.δ)
    push!(hist.ρ,       row.ρ)
    push!(hist.τ,       row.τ)
    push!(hist.κ,       row.κ)
    push!(hist.piter,   row.piter); push!(hist.ppass, row.ppass); push!(hist.pstat, row.pstat)
    push!(hist.citer,   row.citer); push!(hist.cpass, row.cpass); push!(hist.cstat, row.cstat)
    push!(hist.witer,   row.witer); push!(hist.wpass, row.wpass)
    push!(hist.dmin, row.dmin); push!(hist.dmax, row.dmax)
    return hist
end

function Base.empty!(hist::HSDHistory)
    empty!(hist.μ)
    empty!(hist.step)
    empty!(hist.pres)
    empty!(hist.dres)
    empty!(hist.pobj)
    empty!(hist.dobj)
    empty!(hist.δ)
    empty!(hist.ρ)
    empty!(hist.τ)
    empty!(hist.κ)
    empty!(hist.piter); empty!(hist.ppass); empty!(hist.pstat)
    empty!(hist.citer); empty!(hist.cpass); empty!(hist.cstat)
    empty!(hist.witer); empty!(hist.wpass)
    empty!(hist.dmin); empty!(hist.dmax)
    return hist
end

function showtop(io::IO, ::HSDHistory; indent::Integer=0)
    pad = " "^indent
    println(io, pad, "┌──────┬──────────┬──────────┬───────────┬───────────┬────────┬───────┬──────────┬──────────┬──────────┬──────────┬──────────┐")
    println(io, pad, "│ iter │   pres   │   dres   │   pobj    │   dobj    │  step  │ solve │    δ     │    ρ     │    μ     │    τ     │    κ     │")
    return
end

function showbot(io::IO, ::HSDHistory; indent::Integer=0)
    pad = " "^indent
    println(io, pad, "└──────┴──────────┴──────────┴───────────┴───────────┴────────┴───────┴──────────┴──────────┴──────────┴──────────┴──────────┘")
    return
end

function showmid(io::IO, ::HSDHistory; indent::Integer=0)
    pad = " "^indent
    println(io, pad, "├──────┼──────────┼──────────┼───────────┼───────────┼────────┼───────┼──────────┼──────────┼──────────┼──────────┼──────────┤")
    println(io, pad, "│    ⋮ │        ⋮ │        ⋮ │         ⋮ │         ⋮ │      ⋮ │     ⋮ │        ⋮ │        ⋮ │        ⋮ │        ⋮ │        ⋮ │")
    return
end

function showrow(io::IO, i::Integer, row::HSDHistoryRow; indent::Integer=0)
    pad = " "^indent
    println(io, pad, "├──────┼──────────┼──────────┼───────────┼───────────┼────────┼───────┼──────────┼──────────┼──────────┼──────────┼──────────┤")
    print(io, pad)
    # solve = total triangular solve-pairs this step (predictor + corrector directions + woodbury column,
    # each counting CG iterations and refinement passes); δ = penalty, ρ = regularization.
    solve = row.piter + row.ppass + row.citer + row.cpass + row.witer + row.wpass
    @printf(io, "│ %4d │ %8.2e │ %8.2e │ %9.2e │ %9.2e │ %6.4f │ %5d │ %8.2e │ %8.2e │ %8.2e │ %8.2e │ %8.2e │\n",
            i, row.pres, row.dres, row.pobj, row.dobj, row.step, solve, row.δ, row.ρ, row.μ, row.τ, row.κ)
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

function isillposed(hist::HSDHistory; tol=1e-10)
    flag = false

    if !isempty(hist)
        flag = max(hist.τ[end], hist.κ[end]) <= tol
    end

    return flag
end
