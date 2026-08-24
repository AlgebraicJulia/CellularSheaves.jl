struct HSDWorkspace{T} <: AbstractWorkspace{T}
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
    # HSD embedding workspace
    #
    Δp2::FVector{T}  # border direction primal
    Δy2::FVector{T}  # border direction dual
    aτ::FVector{T}   # border row: f - 2Qp/τ
    Qp::FVector{T}   # cached Q*p
    #
    # step-length and status
    #
    step::FVector{T}
    flag::FVector{Bool}
end

function HSDWorkspace{T}(m::Integer, n::Integer, nv::Integer) where {T}
    return HSDWorkspace{T}(
        FVector{T}(undef, n),  # Δf
        FVector{T}(undef, m),  # Δg
        FVector{T}(undef, n),  # f
        FVector{T}(undef, n),  # Δpa
        FVector{T}(undef, m),  # Δya
        FVector{T}(undef, n),  # Δda
        FVector{T}(undef, n),  # Δp
        FVector{T}(undef, m),  # Δy
        FVector{T}(undef, n),  # Δd
        FVector{T}(undef, n),  # Δp2
        FVector{T}(undef, m),  # Δy2
        FVector{T}(undef, n),  # aτ
        FVector{T}(undef, n),  # Qp
        FVector{T}(undef, nv),     # step
        FVector{Bool}(undef, nv),  # flag
    )
end
