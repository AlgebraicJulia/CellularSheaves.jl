abstract type AbstractResult{T} end

include("ipm.jl")
include("hsd.jl")

function Base.show(io::IO, ::MIME"text/plain", result::AbstractResult)
    println(io, typeof(result), ":")
    return showresult(io, result; indent=2)
end

print_timers(r::AbstractResult) = display(r.timers)
