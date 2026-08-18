const SMALL_CONE_THRESHOLD = 3

"""
    AbstractCone

A convex cone.
"""
abstract type AbstractCone end

abstract type AbstractCache{C <: AbstractCone} end

struct Caches{T, I}
    #
    # The ith cache corresponds to the columns
    #
    #   xcol[i] ... xcol[i + 1] - 1
    #
    xcol::FVector{I}
    #
    # The ith cache corresponds to the slots
    #
    #   xblk[i] ... xblk[i + 1] - 1
    #
    xblk::FVector{I}
    #
    # The value
    #
    #   val(b)
    #
    # at slot b.
    #
    val::FVector{T}
end

struct ConeWorkspace{T}
    data::FVector{T}
    work::Vector{T}
    iwork::Vector{BlasInt}
    rng::Xoshiro
end

struct ConeSchedule{T, I}
    nsmll::I
    xsmll::FVector{I}
    small::FVector{ConeWorkspace{T}}
    large::ConeWorkspace{T}
end

function Caches(cones::AbstractVector, B::BlockSparseMatrix{T, I}) where {T, I}
    xcol = FVector{I}(undef, nvtxs(B) + one(I))
    xblk = FVector{I}(undef, nvtxs(B) + one(I))

    c = zero(I)
    b = zero(I)

    for v in vtxs(B)
        ncol = ncols(B, v)
        xcol[v] = c + one(I); c += ncol
        xblk[v] = b + one(I); b += cachesize(cones[v], ncol)
    end

    val = FVector{T}(undef, b)

    xcol[nvtxs(B) + one(I)] = c + one(I)
    xblk[nvtxs(B) + one(I)] = b + one(I)

    return Caches(xcol, xblk, val)
end

function ConeWorkspace{T}(m::Integer, n::Integer) where {T}
    data = FVector{T}(undef, m)
    work = Vector{T}(undef, 1)
    iwork = Vector{BlasInt}(undef, 8triroot(n))
    rng = Xoshiro(0x2654435761)
    return ConeWorkspace{T}(data, work, iwork, rng)
end

function ConeSchedule{T}(cones::AbstractVector, B::BlockSparseMatrix{T, I}, tdmax::I) where {T, I}
    rsmll = zero(I)
    ssmll = zero(I)

    rlarg = zero(I)
    slarg = zero(I)

    nsmll = zero(I)
    tsmll = zero(I)

    for v in vtxs(B)
        ncol = ncols(B, v)

        rlarg = max(rlarg, workspacesize(cones[v], ncol))

        if cones[v] isa SemidefiniteCone
            slarg = max(slarg, ncol)
        end

        if ncol <= SMALL_CONE_THRESHOLD
            rsmll = max(rsmll, workspacesize(cones[v], ncol))

            if cones[v] isa SemidefiniteCone
                ssmll = max(ssmll, ncol)
            end

            nsmll += one(I)
            tsmll += ncol
        end
    end

    nsmll = min(nsmll, tdmax)

    xsmll = FVector{I}(undef, nsmll + one(I))
    small = FVector{ConeWorkspace{T}}(undef, nsmll)
    large = ConeWorkspace{T}(rlarg, slarg)

    v = zero(I)
    Δ = zero(I)

    xsmll[one(I)] = one(I)

    for s in oneto(nsmll)
        small[s] = ConeWorkspace{T}(rsmll, ssmll)

        while v < nvtxs(B) && nsmll * Δ < tsmll
            v += one(I); ncol = ncols(B, v)

            if ncol ≤ SMALL_CONE_THRESHOLD
                Δ += ncol
            end
        end

        xsmll[s + one(I)] = v + one(I)
        Δ = zero(I)
    end

    return ConeSchedule{T, I}(nsmll, xsmll, small, large)
end

"""
    degree(cone::AbstractCone, n::Integer)

Get the rank of a cone with embedding
dimension `n`.
"""
degree(cone::AbstractCone)

"""
    identity!(x::AbstractVector, cone::AbstractCone)

Set x to the fixed point -f'(e) = e of the barrier.
"""
identity!(x::AbstractVector, cone::AbstractCone)

"""
    scale!(H, p, d, cache)

Set H to the Tuncel scaling matrix. If p
and d are elements of a symmetric cone, this
is the Hessian f''(w) of the barrier at the
Nesterov-Todd scaling point w.
"""
scale!(H::AbstractMatrix, p::AbstractVector, d::AbstractVector, cache::AbstractCache)

"""
    corr!(r, p, d, Δp, Δd, σμ, cache)

Set r to the Mehrotra corrector term
r = -d - σμ f'(p) - η, where η is the third-order
correction η = -½ f'''(p)[Δp, f''(p)⁻¹ Δd]. If
p and d are elements of a symmetric cone, this
formula simplifies to r = -d + (σμ e - Δp ∘ Δd) / p,
where Δp ∘ Δd is the Jordan product of Δp and Δd.
"""
corr!(r::AbstractVector, p::AbstractVector, d::AbstractVector, Δp::AbstractVector, Δd::AbstractVector, σμ::Number, cache::AbstractCache)

"""
    maxsteps(p, Δp, d, Δd, cache)

Compute the largest numbers 0 < τp, τd ≤ 1 such that
p + τp Δp and d + τd Δd lie in the interior of their
respective cones
"""
maxsteps(p::AbstractVector, Δp::AbstractVector, d::AbstractVector, Δd::AbstractVector, cache::AbstractCache)

"""
    dualshadow!(sd, p, cache)

Set sd to the dual "shadow" iterate d* = -f'(p).
"""
dualshadow!(sd::AbstractVector, p::AbstractVector, cache::AbstractCache)

"""
    primalshadow!(sp, d, cache)

Set sp to the primal "shadow" iterate p*, solving -f'(p*) = d.
"""
primalshadow!(sp::AbstractVector, d::AbstractVector, cache::AbstractCache)

"""
    cachesize(cone, n)

Return the number of cache slots needed for a cone
with embedding dimension n.
"""
function cachesize(cone::C, n::Integer) where {C <: AbstractCone}
    return cachesize(C, n)
end

"""
    workspacesize(cone, n)

Return the number of workspace floats needed for a cone
with embedding dimension n.
"""
function workspacesize(cone::C, n::Integer) where {C <: AbstractCone}
    return workspacesize(C, n)
end

function workspacesize(::Type{<:AbstractCone}, ::Integer)
    return 0
end

"""
    initcache!(cache)

Initialise a cache.
"""
function initcache!(c::AbstractCache)
    return c
end

"""
    cache(caches, i, cone)

Get the ith cache.
"""
cache(caches::Caches, i::Integer, cone::AbstractCone)

function cachedata(c::Caches, i::Integer)
    return view(c.val, c.xblk[i]:c.xblk[i + 1] - 1)
end

include("eig.jl")
include("sdp.jl")
include("tdc/tdc.jl")
include("pos.jl")
include("soc.jl")
include("noc.jl")
include("utils.jl")
