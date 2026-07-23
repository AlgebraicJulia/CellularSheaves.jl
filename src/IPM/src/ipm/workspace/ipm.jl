struct IPMWorkspace{T} <: AbstractWorkspace{T}
    #
    # residuals
    #
    rp::FVector{T}
    rd::FVector{T}
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
    # iterative refinement workspace
    #
    sp::FVector{T}
    sy::FVector{T}
    dp::FVector{T}
    dy::FVector{T}
end

function IPMWorkspace{T}(m::Integer, n::Integer) where {T}
    return IPMWorkspace{T}(
        FVector{T}(undef, m),  # rp
        FVector{T}(undef, n),  # rd
        FVector{T}(undef, n),  # f
        FVector{T}(undef, n),  # Δpa
        FVector{T}(undef, m),  # Δya
        FVector{T}(undef, n),  # Δda
        FVector{T}(undef, n),  # Δp
        FVector{T}(undef, m),  # Δy
        FVector{T}(undef, n),  # Δd
        FVector{T}(undef, n),  # sp
        FVector{T}(undef, m),  # sy
        FVector{T}(undef, n),  # dp
        FVector{T}(undef, m),  # dy
    )
end
