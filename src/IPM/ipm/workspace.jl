struct IPMWorkspace{T}
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

struct IPMSolver{T, I, W, C}
    Q::BlockSparseMatrix{T, I}
    H::BlockSparseMatrix{T, I}
    B::BlockSparseMatrix{T, I}
    c::FVector{T}
    g::FVector{T}
    p::FVector{T}
    d::FVector{T}
    y::FVector{T}
    K::FVector{C}
    scaling::Scaling{T}
    P::FPermutation{I}
    wrk::IPMWorkspace{T}
    caches::Caches{T, I}
    conewrk::ConeWorkspace{T}
    kkt::W
    hist::IPMHistory{T}
    ν::Int
    settings::IPMSettings{T}
end

function showsolver(io::IO, s::IPMSolver; indent::Integer=0)
    return showsettings(io, s.settings; indent)
end

function Base.show(io::IO, ::MIME"text/plain", s::T) where {T <: IPMSolver}
    println(io, T, ":")
    return showsolver(io, s; indent=2)
end
