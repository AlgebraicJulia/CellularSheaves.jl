const HSDHistoryRow{T} = @NamedTuple{μ::T, step::T, pres::T, dres::T, gap::T, α::T, ρ::T, τ::T, κ::T,
    pbase::Int, prefn::Int, ppass::Int, pstat::RefStatus,
    cbase::Int, crefn::Int, cpass::Int, cstat::RefStatus,
    wbase::Int, wrefn::Int, wpass::Int, wstat::RefStatus,
    pres0::T, pres1::T, cres0::T, cres1::T, wres0::T, wres1::T, r0_p::T, r0_c::T, r0_w::T, craig_p::String, craig_c::String, craig_w::String,
    r1_p::T, r1_c::T, r1_w::T, pres0_d::T, pres0_p::T, cres0_d::T, cres0_p::T, wres0_d::T, wres0_p::T,
    pres_exit::T, cres_exit::T, wres_exit::T, bar_hdiag_med::T, bar_hdiag_frac_mid::T, s2min_p::T, s2max_p::T,
    ritz_beta_p::T, omm_p::T, ritz_θ_p::Vector{T}, ritz_w_p::Vector{T}}

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
    pstat::Vector{RefStatus}
    cbase::Vector{Int}
    crefn::Vector{Int}
    cpass::Vector{Int}
    cstat::Vector{RefStatus}
    wbase::Vector{Int}
    wrefn::Vector{Int}
    wpass::Vector{Int}
    wstat::Vector{RefStatus}
    pres0::Vector{T}
    pres1::Vector{T}
    cres0::Vector{T}
    cres1::Vector{T}
    wres0::Vector{T}
    wres1::Vector{T}
    r0_p::Vector{T}
    r0_c::Vector{T}
    r0_w::Vector{T}
    craig_p::Vector{String}
    craig_c::Vector{String}
    craig_w::Vector{String}
    r1_p::Vector{T}
    r1_c::Vector{T}
    r1_w::Vector{T}
    pres0_d::Vector{T}
    pres0_p::Vector{T}
    cres0_d::Vector{T}
    cres0_p::Vector{T}
    wres0_d::Vector{T}
    wres0_p::Vector{T}
    pres_exit::Vector{T}
    cres_exit::Vector{T}
    wres_exit::Vector{T}
    bar_hdiag_med::Vector{T}
    bar_hdiag_frac_mid::Vector{T}
    s2min_p::Vector{T}
    s2max_p::Vector{T}
    ritz_beta_p::Vector{T}
    omm_p::Vector{T}
    ritz_θ_p::Vector{Vector{T}}
    ritz_w_p::Vector{Vector{T}}
end

function HSDHistory{T}() where {T}
    return HSDHistory{T}(T[], T[], T[], T[], T[], T[], T[], T[], T[],
        Int[], Int[], Int[], RefStatus[],
        Int[], Int[], Int[], RefStatus[],
        Int[], Int[], Int[], RefStatus[],
        T[], T[], T[], T[], T[], T[],
        T[], T[], T[],
        String[], String[], String[],
        T[], T[], T[], T[], T[], T[], T[], T[], T[],
        T[], T[], T[], T[], T[], T[], T[],
        T[], T[], Vector{T}[], Vector{T}[])
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
    pres0 = hist.pres0[i]; pres1 = hist.pres1[i]; cres0 = hist.cres0[i]; cres1 = hist.cres1[i]; wres0 = hist.wres0[i]; wres1 = hist.wres1[i]
    r0_p = hist.r0_p[i]; r0_c = hist.r0_c[i]; r0_w = hist.r0_w[i]
    craig_p = hist.craig_p[i]; craig_c = hist.craig_c[i]; craig_w = hist.craig_w[i]
    r1_p = hist.r1_p[i]; r1_c = hist.r1_c[i]; r1_w = hist.r1_w[i]
    pres0_d = hist.pres0_d[i]; pres0_p = hist.pres0_p[i]; cres0_d = hist.cres0_d[i]; cres0_p = hist.cres0_p[i]; wres0_d = hist.wres0_d[i]; wres0_p = hist.wres0_p[i]
    pres_exit = hist.pres_exit[i]; cres_exit = hist.cres_exit[i]; wres_exit = hist.wres_exit[i]
    bar_hdiag_med = hist.bar_hdiag_med[i]; bar_hdiag_frac_mid = hist.bar_hdiag_frac_mid[i]
    s2min_p = hist.s2min_p[i]; s2max_p = hist.s2max_p[i]
    ritz_beta_p = hist.ritz_beta_p[i]; omm_p = hist.omm_p[i]; ritz_θ_p = hist.ritz_θ_p[i]; ritz_w_p = hist.ritz_w_p[i]
    return (; μ, step, pres, dres, gap, α, ρ, τ, κ,
        pbase, prefn, ppass, pstat, cbase, crefn, cpass, cstat, wbase, wrefn, wpass, wstat, pres0, pres1, cres0, cres1, wres0, wres1, r0_p, r0_c, r0_w, craig_p, craig_c, craig_w,
        r1_p, r1_c, r1_w, pres0_d, pres0_p, cres0_d, cres0_p, wres0_d, wres0_p,
        pres_exit, cres_exit, wres_exit, bar_hdiag_med, bar_hdiag_frac_mid, s2min_p, s2max_p, ritz_beta_p, omm_p, ritz_θ_p, ritz_w_p)
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
    push!(hist.pres0, row.pres0); push!(hist.pres1, row.pres1); push!(hist.cres0, row.cres0); push!(hist.cres1, row.cres1); push!(hist.wres0, row.wres0); push!(hist.wres1, row.wres1)
    push!(hist.r0_p, row.r0_p); push!(hist.r0_c, row.r0_c); push!(hist.r0_w, row.r0_w)
    push!(hist.craig_p, row.craig_p); push!(hist.craig_c, row.craig_c); push!(hist.craig_w, row.craig_w)
    push!(hist.r1_p, row.r1_p); push!(hist.r1_c, row.r1_c); push!(hist.r1_w, row.r1_w)
    push!(hist.pres0_d, row.pres0_d); push!(hist.pres0_p, row.pres0_p); push!(hist.cres0_d, row.cres0_d); push!(hist.cres0_p, row.cres0_p); push!(hist.wres0_d, row.wres0_d); push!(hist.wres0_p, row.wres0_p)
    push!(hist.pres_exit, row.pres_exit); push!(hist.cres_exit, row.cres_exit); push!(hist.wres_exit, row.wres_exit)
    push!(hist.bar_hdiag_med, row.bar_hdiag_med); push!(hist.bar_hdiag_frac_mid, row.bar_hdiag_frac_mid)
    push!(hist.s2min_p, row.s2min_p); push!(hist.s2max_p, row.s2max_p)
    push!(hist.ritz_beta_p, row.ritz_beta_p); push!(hist.omm_p, row.omm_p)
    push!(hist.ritz_θ_p, row.ritz_θ_p); push!(hist.ritz_w_p, row.ritz_w_p)
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
    empty!(hist.pres0); empty!(hist.pres1); empty!(hist.cres0); empty!(hist.cres1); empty!(hist.wres0); empty!(hist.wres1)
    empty!(hist.r0_p); empty!(hist.r0_c); empty!(hist.r0_w)
    empty!(hist.craig_p); empty!(hist.craig_c); empty!(hist.craig_w)
    empty!(hist.r1_p); empty!(hist.r1_c); empty!(hist.r1_w)
    empty!(hist.pres0_d); empty!(hist.pres0_p); empty!(hist.cres0_d); empty!(hist.cres0_p); empty!(hist.wres0_d); empty!(hist.wres0_p)
    empty!(hist.pres_exit); empty!(hist.cres_exit); empty!(hist.wres_exit)
    empty!(hist.bar_hdiag_med); empty!(hist.bar_hdiag_frac_mid)
    empty!(hist.s2min_p); empty!(hist.s2max_p)
    empty!(hist.ritz_beta_p); empty!(hist.omm_p); empty!(hist.ritz_θ_p); empty!(hist.ritz_w_p)
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
