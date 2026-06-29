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
    cachesize(cone, n)

Return the number of cache slots needed for a cone
with embedding dimension n.
"""
cachesize(cone::AbstractCone, n::Integer)

"""
    workspacesize(cone, n)

Return the number of workspace floats needed for a cone
with embedding dimension n.
"""
workspacesize(::AbstractCone, ::Integer) = 0

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

struct ConeWorkspace{T}
    data::FVector{T}
    work::Vector{T}
    iwork::Vector{BlasInt}
end

function ConeWorkspace{T}(cones::AbstractVector, B::BlockSparseMatrix) where {T}
    max_size = 0
    max_iwork = 0

    for v in vtxs(B)
        cone = cones[v]
        n = ncols(B, v)
        max_size = max(max_size, workspacesize(cone, n))
        if cone isa SemidefiniteCone
            max_iwork = max(max_iwork, 8 * triroot(n))
        end
    end

    return ConeWorkspace{T}(
        FVector{T}(undef, max_size),
        Vector{T}(undef, 1),
        Vector{BlasInt}(undef, max_iwork),
    )
end

include("sdp.jl")
include("tdc/tdc.jl")
include("pos.jl")
include("soc.jl")
include("noc.jl")
include("utils.jl")
