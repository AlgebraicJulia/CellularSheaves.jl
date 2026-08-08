abstract type AbstractSettings{T} end

include("ipm.jl")
include("hsd.jl")

function Base.show(io::IO, ::MIME"text/plain", set::T) where {T <: AbstractSettings}
    println(io, T, ":")
    return showsettings(io, set; indent=2)
end
