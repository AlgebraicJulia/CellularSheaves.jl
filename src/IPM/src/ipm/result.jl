struct IPMResult{T}
    p::Vector{T}
    d::Vector{T}
    y::Vector{T}
    status::IPMStatus
    niter::Int
    npred::Int
    ncorr::Int
    history::IPMHistory{T}
    timers::TimerOutput
end

function showresult(io::IO, result::IPMResult; indent::Integer=0)
    pad = " "^indent
    println(io, pad, "status: ", result.status)
    println(io, pad, "niter: ",  result.niter)
    println(io, pad, "npred: ",  result.npred)
    println(io, pad, "ncorr: ",  result.ncorr)

    if !isempty(result.history)
        println(io, pad, "μ: ",    result.history.μ[end])
        println(io, pad, "pres: ", result.history.pres[end])
        println(io, pad, "dres: ", result.history.dres[end])
    end

    return
end

function Base.show(io::IO, ::MIME"text/plain", result::T) where {T <: IPMResult}
    println(io, T, ":")
    return showresult(io, result; indent=2)
end

function print_timers(r::IPMResult)
    display(r.timers)
end

struct HSDResult{T}
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
end

function showresult(io::IO, result::HSDResult; indent::Integer=0)
    pad = " "^indent
    println(io, pad, "status: ", result.status)
    println(io, pad, "niter: ",  result.niter)
    println(io, pad, "npred: ",  result.npred)
    println(io, pad, "ncorr: ",  result.ncorr)
    println(io, pad, "τ: ",      result.τ)
    println(io, pad, "κ: ",      result.κ)

    if !isempty(result.history)
        println(io, pad, "μ: ",    result.history.μ[end])
        println(io, pad, "pres: ", result.history.pres[end])
        println(io, pad, "dres: ", result.history.dres[end])
    end

    return
end

function Base.show(io::IO, ::MIME"text/plain", result::T) where {T <: HSDResult}
    println(io, T, ":")
    return showresult(io, result; indent=2)
end

function print_timers(r::HSDResult)
    display(r.timers)
end
