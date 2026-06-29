"""
    ExponentialCone <: AbstractTDCone

The exponential cone, consisting of all triples (x, y, z)
such that x > 0, y > 0, and y log(x/y) ≥ z.
"""
struct ExponentialCone <: AbstractTDCone end

const ExponentialConeCache{T} = AbstractTDConeCache{ExponentialCone, T}

function cache(c::Caches{T}, i::Integer, cone::ExponentialCone) where T
    data = cachedata(c, i)
    L    = reshape(view(data, 1:9), 3, 3)
    sd   = view(data, 10:12)
    seed = view(data, 13)
    AbstractTDConeCache(cone, L, sd, seed)
end

# compute the fixed point
#
#   f'(e) = -e
#
function expid!(x::AbstractVector)
    x[1] =  1.2909282315382298
    x[2] =  0.8051015526498357
    x[3] = -0.8278379086082098
    return x
end

function identity!(x::AbstractVector, ::ExponentialCone)
    return expid!(x)
end

function initcache!(cache::ExponentialConeCache)
    cache.seed[] = 0.8051015526498357
    return cache
end

# the barrier argument
#
#   ψ(x) = x₂ log(x₁/x₂) - x₃
#
function exppsi(x::AbstractVector)
    return x[2] * log(x[1] / x[2]) - x[3]
end

# the gradient of ψ
#
#   ψ'(x) = (x₂/x₁, log(x₁/x₂) - 1, -1)
#
function exppsigrad!(g::AbstractVector{T}, x::AbstractVector{T}) where {T}
    g[1] = x[2] / x[1]
    g[2] = log(x[1] / x[2]) - one(T)
    g[3] = -one(T)
    return g
end

# The gradient
#
#   f'(x) = -ψ'(x)/ψ(x) - (1/x₁, 1/x₂, 0)
#
# of the barrier function.
function expbarrgrad!(g::AbstractVector{T}, x::AbstractVector{T}) where {T}
    ψ = exppsi(x)
    exppsigrad!(g, x)

    g[1] = -g[1] / ψ - inv(x[1])
    g[2] = -g[2] / ψ - inv(x[2])
    g[3] = -g[3] / ψ

    return g
end

# Compute the analytic factor R
#
#   f''(p) = R Rᵀ
#
# of the barrier function f.
function expfact!(R::AbstractMatrix{T}, x::AbstractVector{T}) where {T}
    x1, x2 = x[1], x[2]

    ψ = exppsi(x)

    σ = sqrt(one(T) + 2x2 / ψ)

    R[1,1] = (one(T) - σ) / 2x1
    R[1,2] = (one(T) + σ) / 2x1
    R[2,1] = (one(T) + σ) / 2x2
    R[2,2] = (one(T) - σ) / 2x2

    R[1,3] = x2 / x1 / ψ
    R[2,3] = (log(x1 / x2) - one(T)) / ψ
    R[3,3] = -inv(ψ)

    R[3,1] = zero(T)
    R[3,2] = zero(T)

    return R
end

# Compute the third-order directional derivative
#
#   f'''(x)[u]
#
# as a 3x3 matrix.
function expbarrthird!(D::AbstractMatrix{T}, x::AbstractVector{T}, u::AbstractVector{T}) where {T}
    ψ = exppsi(x)

    x1, x2     = x[1], x[2]
    u1, u2, u3 = u[1], u[2], u[3]

    ψg1 = x2 / x1
    ψg2 = log(x1 / x2) - one(T)
    ψgu = ψg1 * u1 + ψg2 * u2 - u3

    ψH11 = -x2 / x1^2
    ψH21 =  inv(x1)
    ψH22 = -inv(x2)

    ψHu1 = ψH11 * u1 + ψH21 * u2
    ψHu2 = ψH21 * u1 + ψH22 * u2

    ψ2 = ψ^2
    ψ3 = ψ^3

    α = -2ψgu / ψ3

    D[1,1] =  α * ψg1^2 + (2ψg1 * ψHu1 + ψH11 * ψgu) / ψ2 - (2x2 * u1 / x1^3 - u2 / x1^2) / ψ - 2u1 / x1^3
    D[2,1] =  α * ψg1 * ψg2 + (ψg1 * ψHu2 + ψg2 * ψHu1 + ψH21 * ψgu) / ψ2 + u1 / x1^2 / ψ
    D[3,1] = -α * ψg1 - ψHu1 / ψ2
    D[2,2] =  α * ψg2^2 + (2ψg2 * ψHu2 + ψH22 * ψgu) / ψ2 - u2 / x2^2 / ψ - 2u2 / x2^3
    D[3,2] = -α * ψg2 - ψHu2 / ψ2
    D[3,3] =  α

    D[1,2] = D[2,1]
    D[1,3] = D[3,1]
    D[2,3] = D[3,2]

    return D
end

# Determine if x is in the exponential cone.
function expincone(x::AbstractVector)
    x[1] > 0 && x[2] > 0 && exppsi(x) > 0
end

# Determine if z is in the dual exponential cone.
function expindual(z::AbstractVector)
    z[1] > 0 && z[3] < 0 && ℯ * z[1] >= -z[3] * exp(z[2] / z[3])
end

# compute the "shadow" primal, solving
#
#   f'(p*) = -d
#
# using a 1-D scalar root-find on the function
#
#   h(p₂*) = d₃ (log(p₁* / p₂*) - 1) - 1 / p₂* + d₂
#
# with
#
#   p₁* = (1 - p₂* d₃) / d₁
#   p₃* = p₂* log(p₁* / p₂*) + 1 / d₃
#
function expdualgrad!(sp::AbstractVector{T}, seed::T, d::AbstractVector{T}) where {T}
    d1, d2, d3 = d[1], d[2], d[3]

    function h(p2)
        d3 * (log((one(T) - p2 * d3) / (d1 * p2)) - one(T)) - inv(p2) + d2
    end

    function hp(p2)
        -d3^2 / (one(T) - p2 * d3) - d3 / p2 + inv(p2^2)
    end
    #
    # bracket the function h, finding points lo and hi
    # such that
    #
    #   h(lo) < 0 < h(hi)
    #
    if h(seed) < 0
        lo = seed
        hi = seed * 2

        while h(hi) < 0
            hi *= 2
        end
    else
        hi = seed
        lo = seed / 2

        while h(lo) > 0
            lo /= 2
        end
    end
    #
    # find p₂* such that
    #
    #   h(p₂*) = 0
    #
    p2 = rtsafe(h, hp, lo, hi, seed)
    #
    # recover p* as
    #
    #   p₁* = (1 - p₂* d₃) / d₁
    #   p₃* = p₂* log(p₁* / p₂*) + 1 / d₃
    #
    sp[1] = p1 = (one(T) - p2 * d3) / d1
    sp[2] = p2
    sp[3] = p2 * log(p1 / p2) + inv(d3)

    return p2
end

#
# AbstractTDCone Interface
#

function tdfact!(L::AbstractMatrix, p::AbstractVector, ::ExponentialConeCache)
    return expfact!(L, p)
end

function tdbarrgrad!(g::AbstractVector, p::AbstractVector, ::ExponentialConeCache)
    return expbarrgrad!(g, p)
end

function tdbarrthird!(D::AbstractMatrix, p::AbstractVector, u::AbstractVector, ::ExponentialConeCache)
    return expbarrthird!(D, p, u)
end

function tddualgrad!(sp::AbstractVector, seed, d::AbstractVector, ::ExponentialConeCache)
    return expdualgrad!(sp, seed, d)
end

function tdincone(p::AbstractVector, ::ExponentialConeCache)
    return expincone(p)
end

function tdindual(d::AbstractVector, ::ExponentialConeCache)
    return expindual(d)
end

# apply the Hessian
#
#   f''(p) v = R Rᵀ v
#
function tdhessmul!(u::AbstractVector{T}, L::AbstractMatrix{T}, v::AbstractVector{T}, ::ExponentialConeCache) where {T}
    t1 = L[1,1] * v[1] + L[2,1] * v[2] + L[3,1] * v[3]
    t2 = L[1,2] * v[1] + L[2,2] * v[2] + L[3,2] * v[3]
    t3 = L[1,3] * v[1] + L[2,3] * v[2] + L[3,3] * v[3]

    u[1] = L[1,1] * t1 + L[1,2] * t2 + L[1,3] * t3
    u[2] = L[2,1] * t1 + L[2,2] * t2 + L[2,3] * t3
    u[3] = L[3,1] * t1 + L[3,2] * t2 + L[3,3] * t3
    return u
end

# materialize the scaled Hessian
#
#   H = μ f''(p) = μ R Rᵀ
#
function tdgram!(H::AbstractMatrix, L::AbstractMatrix, μ, ::ExponentialConeCache)
    mul3!(H, L, L', μ, 0)
    return H
end

# solve for v in
#
#   f''(p) v = b
#
# using the factorization f''(p) = R Rᵀ.
#
function tdhessldiv!(L::AbstractMatrix{T}, v::AbstractVector{T}, ::ExponentialConeCache) where {T}
    ldiv3!(L,  v)
    ldiv3!(L', v)
    return v
end
