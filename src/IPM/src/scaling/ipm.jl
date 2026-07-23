struct IPMScaling{T} <: AbstractScaling{T}
    cscl::FVector{T}   # length n: per-column scaling (constant on each cone block)
    rscl::FVector{T}   # length m: per-row scaling
end

"""
    IPMScaling{T}(n, m)

Construct an identity (trivial) scaling with all entries equal to 1.
"""
function IPMScaling{T}(n::Int, m::Int) where {T}
    cscl = FVector{T}(undef, n)
    rscl = FVector{T}(undef, m)

    fill!(cscl, one(T))
    fill!(rscl, one(T))

    return IPMScaling{T}(cscl, rscl)
end
