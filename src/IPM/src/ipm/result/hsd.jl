struct HSDResult{T} <: AbstractResult{T}
    p::Vector{T}
    d::Vector{T}
    y::Vector{T}
    status::IPMStatus
    niter::Int
    nsolve::Int
    τ::T
    κ::T
    history::HSDHistory{T}
    timers::TimerOutput
    mu::T
    pres::T
    dres::T
    pobj::T
    dobj::T
end

function showresult(io::IO, result::HSDResult; indent::Integer=0)
    pad = " "^indent
    println(io, pad, "status: ", result.status)
    println(io, pad, "niter: ",  result.niter)
    println(io, pad, "nsolve: ", result.nsolve)
    println(io, pad, "τ: ",      result.τ)
    println(io, pad, "κ: ",      result.κ)
    println(io, pad, "mu: ",     result.mu)
    println(io, pad, "pres: ",   result.pres)
    println(io, pad, "dres: ",   result.dres)
    println(io, pad, "pobj: ",   result.pobj)
    println(io, pad, "dobj: ",   result.dobj)

    return
end
