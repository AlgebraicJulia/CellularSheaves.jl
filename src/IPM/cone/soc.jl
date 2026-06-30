"""
    SecondOrderCone <: AbstractCone

The n + 1-dimensional second-order cone,
consisting of all pairs (x, y) such that
y ≥ ‖x‖.
"""
struct SecondOrderCone <: AbstractCone end

struct SecondOrderConeCache{T} <: AbstractCache{SecondOrderCone}
    cone::SecondOrderCone
    #
    # the square root
    #
    #   β = √det(ω)
    #
    # of the deminant of the
    # Nesterov-Todd scaling point.
    #
    β::FScalarView{T}
    #
    # the normalized Nesterov-Todd scaling
    # point:
    #
    #   w = ω / β
    #
    w::FVectorView{T}
end

# compute the J-inner product
#
#   xᵀ J y = x₁y₁ - x₂y₂ - … - xₙyₙ
#
# using compensated arithmetic.
function socdot(x::AbstractVector{T}, y::AbstractVector{T}) where {T}
    @assert length(x) == length(y)
    n = length(x)

    s, c = twoprod(x[1], y[1])

    @inbounds for i in 2:n
        p, e = twoprod(x[i], y[i])
        s, e2 = twosum(s, -p)
        c += -e + e2
    end

    return s + c
end

# compute the Jordan determinant
#
#   xᵀ J x = x₁x₁ - x₂x₂ - … - xₙxₙ
#
# using compensated arithmetic.
function socdet(x::AbstractVector)
    return socdot(x, x)
end

# compute the Jordan product
#
#   x ∘ y = (xᵀy, x₁y₂ + y₁x₂, ..., x₁yₙ + y₁xₙ)
#
function socmul!(out::AbstractVector{T}, x::AbstractVector{T}, y::AbstractVector{T}) where {T}
    n = length(x)

    x1 = x[1]
    y1 = y[1]

    out[1] = dot(x, y)

    for i in 2:n
        out[i] = y1 * x[i] + x1 * y[i]
    end

    return out
end

# solve for z in
#
#   y ∘ z = b,
#
# where y ∘ z is the Jordan product
# of y and z.
function socdiv!(z::AbstractVector{T}, b::AbstractVector{T}) where {T}
    n = length(z)
    δ = socdet(z)

    z1 = z[1]
    b1 = b[1] * z1

    for i in 2:n
        b1 -= z[i] * b[i]
    end

    b[1] = b1 /= δ

    for i in 2:n
        b[i] = (b[i] - b1 * z[i]) / z1
    end

    return b
end

# If flag = false, evaluate the product
#
#   y = √H x.
#
# If flag = true, solve for y:
#
#   √H y = x
#
# where H is the Hessian of the primal barrier
# function at the Nesterov-Todd scaling point.
function socroot!(x::AbstractVector{T}, w::AbstractVector{T}, β::T, flag::Bool) where {T}
    n = length(x)

    if !flag
        σ = -one(T)
        α =  inv(β)
    else
        σ =  one(T)
        α =      β
    end

    w1p1 = w[1] + 1

    wpdx = w1p1 * x[1]

    for i in 2:n
        wpdx += σ * w[i] * x[i]
    end

    wpdx /= w1p1

    x[1] = α * (wpdx * w1p1 - x[1])

    for i in 2:n
        x[i] = α * (σ * wpdx * w[i] + x[i])
    end

    return x
end

# Compute the Nesterov-Todd scaling point
# ω in factored form:
#
#   β = √det(ω)
#   w = ω / β
#
# Then assemble the Hessian of the primal
# barrier function evaluated at ω:
#
#   H = (2 (Jw)(Jw)ᵀ - J) / β²
#
function socscale!(
        H::AbstractMatrix{T},
        w::AbstractVector{T},
        p::AbstractVector{T},
        d::AbstractVector{T}
    ) where {T}
    n = length(p)
    #
    # compute the scaling point ω
    #
    pdet = socdet(p)
    ddet = socdet(d)

    pdot = cdot(p, d)

    spdet = sqrt(pdet)
    sddet = sqrt(ddet)

    κ = sqrt(two(T) * (spdet * sddet + pdot))
    β = sqrt(spdet / sddet)

    w[1] = (p[1] / β + d[1] * β) / κ

    for i in 2:n
        w[i] = (p[i] / β - d[i] * β) / κ
    end
    #
    # assemble the Hessian
    #
    #   H = (2 (Jw)(Jw)ᵀ - J) / β²
    #
    η = inv(β^2)

    w1 = w[1]; H[1, 1] = 2η * w1^2 - η

    for i in 2:n
        H[i, 1] = -2η * w[i] * w1
    end

    for j in 2:n
        wj = w[j]; H[j, j] = 2η * wj^2 + η

        for i in j + 1:n
            H[i, j] = 2η * w[i] * wj
        end
    end

    return β
end

function soccorr!(
        r::AbstractVector{T},
        w::AbstractVector{T},
        β::T,
        p::AbstractVector{T},
        Δp::AbstractVector{T},
        Δd::AbstractVector{T},
        σμ::Real
    ) where {T}
    n = length(p)
    pdet = socdet(p)

    socroot!(Δp, w, β, false)
    socroot!(Δd, w, β, true)

    socmul!(r, Δp, Δd)

    copyto!(Δp, p)
    socroot!(Δp, w, β, false)

    socdiv!(Δp, r)

    axpy!(one(T), Δp, r)
    socroot!(r, w, β, false)

    r[1] = 2σμ * p[1] / pdet - r[1]

    for i in 2:n
        r[i] = -2σμ * p[i] / pdet - r[i]
    end

    return r
end

# construct the identity element
#
#   e = (√2, 0, …, 0)
#
function socid!(x::AbstractVector{T}) where {T}
    fill!(x, zero(T))
    x[1] = sqrt(T(2))
    return x
end

function socmaxstep(x::AbstractVector{T}, Δx::AbstractVector{T}) where {T}
    n = length(x)

    # compute scalars
    #
    #   a = Δx J Δx
    #   b = 2x J Δx
    #   c =  x J  x
    #
    # such that that
    #
    #   det(x + Δx τ) = (x + Δx τ)ᵀ J (x + Δx τ)
    #                 = aτ² + bτ + c
    #
    a = socdet(Δx)
    c = socdet( x)

    b = x[1] * Δx[1]

    for i in 2:n
        b -= x[i] * Δx[i]
    end

    b *= 2

    # find the largest number τ ≤ 1 such that
    #
    #   1. a τ² +   b τ + c  ≥ 0
    #   2.        Δx₁ τ + x₁ ≥ 0
    #
    τ = one(T)
    #
    # ensure that
    #
    #   a τ² + b τ + c  ≥ 0
    #
    if abs(a) > eps(T)
        #
        # d is the discriminant
        #
        #   d = b² - 4ac
        #
        d = b^2 - 4*a*c

        if d > -eps(T)
            #
            # - a > 0: roots have the same sign
            #          and τ1 is the smaller one
            #
            # - a < 0: roots have opposite signs
            #          and τ1 is the positive one
            #
            s = sqrt(max(d, zero(T)))
            q = -(b + copysign(s, b)) / 2

            if b ≥ 0
                τ1 = q / a
            else
                τ1 = c / q
            end

            if τ1 > 0
                τ = min(τ, τ1)
            end
        end
    elseif b < -eps(T)
        τ = min(τ, -c / b)
    end
    #
    # ensure that
    #
    #   Δx₁ τ + x₁ ≥ 0
    #
    if Δx[1] < 0
        τ = min(τ, -x[1] / Δx[1])
    end

    return τ
end

#
# AbstractCone Interface
#

function degree(::SecondOrderCone, n::Integer)
    return 2
end

function cachesize(::SecondOrderCone, n::Integer)
    return 1 + n
end

function cache(c::Caches, i::Integer, cone::SecondOrderCone)
    data = cachedata(c, i)
    β = view(data, 1)
    w = view(data, 2:length(data))
    SecondOrderConeCache(cone, β, w)
end

function identity!(x::AbstractVector, ::SecondOrderCone)
    return socid!(x)
end

function scale!(H::AbstractMatrix, p::AbstractVector, d::AbstractVector, cache::SecondOrderConeCache, ::ConeWorkspace)
    cache.β[] = socscale!(H, cache.w, p, d)
    return H
end

function corr!(
        r::AbstractVector,
        p::AbstractVector,
        ::AbstractVector,
        Δp::AbstractVector,
        Δd::AbstractVector,
        σμ::Real,
        cache::SecondOrderConeCache,
        ::ConeWorkspace,
    )
    soccorr!(r, cache.w, cache.β[], p, Δp, Δd, σμ)
end

function maxsteps(p::AbstractVector, Δp::AbstractVector, d::AbstractVector, Δd::AbstractVector, ::SecondOrderConeCache, ::ConeWorkspace)
    return socmaxstep(p, Δp), socmaxstep(d, Δd)
end
