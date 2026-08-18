struct IPMResult{T} <: AbstractResult{T}
    p::Vector{T}
    d::Vector{T}
    y::Vector{T}
    status::IPMStatus
    niter::Int
    nsolve::Int
    history::IPMHistory{T}
    timers::TimerOutput
    mu::T
    pres::T
    dres::T
    pobj::T
    dobj::T
end

function showresult(io::IO, result::IPMResult; indent::Integer=0)
    pad = " "^indent
    println(io, pad, "status: ", result.status)
    println(io, pad, "niter: ",  result.niter)
    println(io, pad, "nsolve: ", result.nsolve)
    println(io, pad, "mu: ",     result.mu)
    println(io, pad, "pres: ",   result.pres)
    println(io, pad, "dres: ",   result.dres)
    println(io, pad, "pobj: ",   result.pobj)
    println(io, pad, "dobj: ",   result.dobj)

    return
end
