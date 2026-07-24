struct IPMResult{T} <: AbstractResult{T}
    p::Vector{T}
    d::Vector{T}
    y::Vector{T}
    status::IPMStatus
    niter::Int
    nsolve::Int          # total solve (CRAIG) iterations (base + refinement, all solves, all iterations)
    history::IPMHistory{T}
    timers::TimerOutput
end

function showresult(io::IO, result::IPMResult; indent::Integer=0)
    pad = " "^indent
    println(io, pad, "status: ", result.status)
    println(io, pad, "niter: ",  result.niter)
    println(io, pad, "nsolve: ", result.nsolve)

    if !isempty(result.history)
        println(io, pad, "α: ",    result.history.α[end])
        println(io, pad, "μ: ",    result.history.μ[end])
        println(io, pad, "pres: ", result.history.pres[end])
        println(io, pad, "dres: ", result.history.dres[end])
    end

    return
end
