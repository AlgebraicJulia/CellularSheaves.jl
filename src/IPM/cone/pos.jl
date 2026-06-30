"""
    PositiveCone <: AbstractCone

An n-dimensional positive orthant.
"""
struct PositiveCone <: AbstractCone end

struct PositiveConeCache <: AbstractCache{PositiveCone}
    cone::PositiveCone
end

function PositiveConeCache()
    return PositiveConeCache(PositiveCone())
end

# construct the ones vector
#
#   e = (1, …, 1)
#
function posid!(x::AbstractVector)
    fill!(x, true)
    return x
end

# construct the diagonal scaling matrix
#
#   H = diag(w)⁻²
#
# where
#
#   wᵢ = √pᵢ / √dᵢ
#
# is the Nesterov-Todd scaling point.
#
function posscale!(H::AbstractMatrix, p::AbstractVector, d::AbstractVector)
    fill!(H, false)

    for i in eachindex(p)
        H[i, i] = d[i] / p[i]
    end

    return H
end

# Compute the corrector term
#
#   rᵢ = (σμ - Δpᵢ Δdᵢ) / pᵢ - dᵢ.
#
function poscorr!(
        r::AbstractVector,
        p::AbstractVector,
        d::AbstractVector,
        Δp::AbstractVector,
        Δd::AbstractVector,
        σμ::Real
    )
    for i in eachindex(r)
        r[i] = (σμ - Δp[i] * Δd[i]) / p[i] - d[i]
    end

    return r
end

# Find the largest number 0 < τ ≤ 1 such that
#
#   x + τ Δx ≥ 0.
#
# This is precisely τ = min {1, κ}, where
#
#   κ = min { -xᵢ / Δxᵢ : Δxᵢ < 0 }.
#
function posmaxstep(x::AbstractVector{T}, Δx::AbstractVector{T}) where {T}
    τ = one(T)

    for i in eachindex(x)
        Δxi = Δx[i]

        if Δxi < 0
            τ = min(τ, -x[i] / Δxi)
        end
    end

    return τ
end

#
# AbstractCone Interface
#

function degree(::PositiveCone, n::Integer)
    return n
end

function cachesize(::PositiveCone, n::Integer)
    return 0
end

function cache(::Caches, ::Integer, c::PositiveCone)
    return PositiveConeCache(c)
end

function identity!(x::AbstractVector, ::PositiveCone)
    return posid!(x)
end

function scale!(H::AbstractMatrix, p::AbstractVector, d::AbstractVector, ::PositiveConeCache, ::ConeWorkspace)
    return posscale!(H, p, d)
end

function corr!(
        r::AbstractVector,
        p::AbstractVector,
        d::AbstractVector,
        Δp::AbstractVector,
        Δd::AbstractVector,
        σμ::Real,
        ::PositiveConeCache,
        ::ConeWorkspace,
    )
    return poscorr!(r, p, d, Δp, Δd, σμ)
end

function maxsteps(p::AbstractVector, Δp::AbstractVector, d::AbstractVector, Δd::AbstractVector, ::PositiveConeCache, ::ConeWorkspace)
    return posmaxstep(p, Δp), posmaxstep(d, Δd)
end
