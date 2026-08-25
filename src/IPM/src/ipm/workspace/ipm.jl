struct IPMWorkspace{T} <: AbstractWorkspace{T}
    #
    # residuals: Δf dual (n), Δg B-primal (m)
    #
    Δf::FVector{T}
    Δg::FVector{T}
    #
    # Newton right-hand-side
    #
    f::FVector{T}
    #
    # predictor directions
    #
    Δpa::FVector{T}
    Δya::FVector{T}
    Δda::FVector{T}
    #
    # corrector directions
    #
    Δp::FVector{T}
    Δy::FVector{T}
    Δd::FVector{T}
    #
    # step-length and status
    #
    step::FVector{T}
    flag::FVector{Bool}
end

function IPMWorkspace{T}(m::Integer, n::Integer, nv::Integer) where {T}
    return IPMWorkspace{T}(
        FVector{T}(undef, n),  # Δf
        FVector{T}(undef, m),  # Δg
        FVector{T}(undef, n),  # f
        FVector{T}(undef, n),  # Δpa
        FVector{T}(undef, m),  # Δya
        FVector{T}(undef, n),  # Δda
        FVector{T}(undef, n),  # Δp
        FVector{T}(undef, m),  # Δy
        FVector{T}(undef, n),  # Δd
        FVector{T}(undef, nv),     # step
        FVector{Bool}(undef, nv),  # flag
    )
end
