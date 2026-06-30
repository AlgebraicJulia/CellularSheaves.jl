"""
    CofreeCone <: AbstractCone

The cone of all n-dimensional Euclidean
vectors.
"""
struct CofreeCone <: AbstractCone end

struct CofreeConeCache <: AbstractCache{CofreeCone}
    cone::CofreeCone
end

function CofreeConeCache()
    return CofreeConeCache(CofreeCone())
end

#
# AbstractCone Interface
#

function degree(::CofreeCone, n::Integer)
    return 0
end

function cachesize(::CofreeCone, n::Integer)
    return 0
end

function cache(::Caches, ::Integer, c::CofreeCone)
    return CofreeConeCache(c)
end

function identity!(x::AbstractVector{T}, ::CofreeCone) where {T}
    fill!(x, zero(T))
    return x
end

function scale!(H::AbstractMatrix{T}, ::AbstractVector{T}, ::AbstractVector{T}, ::CofreeConeCache, ::ConeWorkspace) where {T}
    fill!(H, zero(T))
    return H
end

function corr!(
        r::AbstractVector{T},
        ::AbstractVector{T},
        ::AbstractVector{T},
        ::AbstractVector{T},
        ::AbstractVector{T},
        ::Real,
        ::CofreeConeCache,
        ::ConeWorkspace,
    ) where {T}
    fill!(r, zero(T))
    return r
end

function maxsteps(::AbstractVector{T}, ::AbstractVector{T}, ::AbstractVector{T}, ::AbstractVector{T}, ::CofreeConeCache, ::ConeWorkspace) where {T}
    return one(T), one(T)
end
