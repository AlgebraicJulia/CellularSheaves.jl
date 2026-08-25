struct IPMScaling{T} <: AbstractScaling{T}
    pscl::FVector{T}   # length n: per-column scaling (constant on each cone block)
    yscl::FVector{T}   # length m: per-B-row scaling
end

"""
    IPMScaling{T}(n, m)

Construct an identity (trivial) scaling with all entries equal to 1.
"""
function IPMScaling{T}(n::Int, m::Int) where {T}
    pscl  = FVector{T}(undef, n)
    yscl  = FVector{T}(undef, m)

    fill!(pscl,  one(T))
    fill!(yscl,  one(T))

    return IPMScaling{T}(pscl, yscl)
end
