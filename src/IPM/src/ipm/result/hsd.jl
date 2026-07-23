struct HSDResult{T} <: AbstractResult{T}
    p::Vector{T}
    d::Vector{T}
    y::Vector{T}
    status::IPMStatus
    niter::Int
    npred::Int
    ncorr::Int
    τ::T
    κ::T
    history::HSDHistory{T}
    timers::TimerOutput
    # C5 (addendum): the answer at the returned point, in USER frame (no frame suffixes — the
    # container disambiguates: history columns are embedding-frame, result fields user-frame).
    # Every field is evaluated at the returned point; nothing is a rescaled history scalar.
    mu::T     # dot(p,d)/(ν·τ²) — the affine solver's μ at the returned pair
    pres::T   # primal residual, original data and units, freshly recomputed (norms don't survive
              # equilibration reweighting, so the stored embedding scalar is not this quantity)
    dres::T   # dual residual, original data and units
    pobj::T   # primal objective ½pᵀQp − cᵀp in problem.jl's documented convention
    dobj::T   # dual objective, same convention. (The user's gap is pobj − dobj.)
    # Certificate carve-out: for infeasibility statuses mu/pres/dres are NaN and the objective
    # slots carry the certificate quantities — pobj ← gᵀy (primal infeasible), dobj ← cᵀp (dual).
end

# C5 (addendum): non-stored relgap convenience accessor.
relgap(r::HSDResult) = abs(r.pobj - r.dobj) / (1 + max(abs(r.pobj), abs(r.dobj)))

function showresult(io::IO, result::HSDResult; indent::Integer=0)
    pad = " "^indent
    println(io, pad, "status: ", result.status)
    println(io, pad, "niter: ",  result.niter)
    println(io, pad, "npred: ",  result.npred)
    println(io, pad, "ncorr: ",  result.ncorr)
    println(io, pad, "τ: ",      result.τ)
    println(io, pad, "κ: ",      result.κ)
    # C5 (addendum): the answer at the returned point, in user frame.
    println(io, pad, "mu: ",     result.mu)
    println(io, pad, "pres: ",   result.pres)
    println(io, pad, "dres: ",   result.dres)
    println(io, pad, "pobj: ",   result.pobj)
    println(io, pad, "dobj: ",   result.dobj)

    return
end
