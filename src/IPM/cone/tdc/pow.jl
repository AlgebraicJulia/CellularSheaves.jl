"""
    PowerCone{T} <: AbstractTDCone

The three-dimensional power cone with parameter α ∈ (0, 1),
consisting of all triples (x₁, x₂, x₃) such that
x₁ ≥ 0, x₂ ≥ 0, and x₁^α x₂^(1-α) ≥ |x₃|.
"""
struct PowerCone{T} <: AbstractTDCone
    α::T

    function PowerCone{T}(α) where {T}
        @assert 0 < α < 1
        return new{T}(α)
    end
end

function PowerCone(α::T) where {T}
    return PowerCone{T}(α)
end

const PowerConeCache{T} = AbstractTDConeCache{PowerCone{T}, T}

function cache(c::Caches{T}, i::Integer, cone::PowerCone{T}) where {T}
    data = cachedata(c, i)
    L    = reshape(view(data, 1:9), 3, 3)
    sd   = view(data, 10:12)
    seed = view(data, 13)
    AbstractTDConeCache(cone, L, sd, seed)
end

# construct the identity element
#
#   e = (√(1 + α), √(2 - α), 0)
#
function identity!(x::AbstractVector, cone::PowerCone)
    α = cone.α
    x[1] = sqrt(1 + α)
    x[2] = sqrt(2 - α)
    x[3] = false
    return x
end

function initcache!(cache::PowerConeCache{T}) where {T}
    cache.seed[] = true
    return cache
end

# evaluate the barrier argument
#
#   φ(p) = p₁ᵃ p₂ᵇ - p₃²
#
# where
#
#   - a = 2α
#   - b = 2 - 2α
#
function powphi(x::AbstractVector, α)
    a = 2α
    b = 2 - 2α
    return x[1]^a * x[2]^b - x[3]^2
end

# evaluate the gradient
#
#   f'(p)
#
# of the barrier function f.
function powbarrgrad!(g::AbstractVector, x::AbstractVector, α)
    a = 2α
    b = 2 - 2α
    p = x[1]^a * x[2]^b
    φ = p - x[3] * x[3]
    ρ = p / φ

    g[1] = -(a * ρ + 1 - α) / x[1]
    g[2] = -(b * ρ     + α) / x[2]
    g[3] =  2x[3] / φ

    return g
end

# solve for x in
#
#   f''(p) x = b
#
# using a pre-computed factorization
#
#   P f''(p) Pᵀ = L Lᵀ
#
function powldiv3!(L::AbstractMatrix{T}, b::AbstractVector{T}) where {T}
    b1, b2, b3 = b[1], b[2], b[3]

    y1 =  b3                              / L[1,1]
    y2 = (b1 - L[2,1] * y1)               / L[2,2]
    y3 = (b2 - L[3,1] * y1 - L[3,2] * y2) / L[3,3]

    z3 = y3                               / L[3,3]
    z2 = (y2 - L[3,2] * z3)               / L[2,2]
    z1 = (y1 - L[2,1] * z2 - L[3,1] * z3) / L[1,1]

    b[1] = z2
    b[2] = z3
    b[3] = z1
    return b
end

# factorize the Hessian
#
#   P f''(p) Pᵀ = L Lᵀ
#
# of the barrier function f, using the pivot order
# (3, 1, 2). The factor is built directly from p;
# the Hessian is never assembled.
function powfact!(L::AbstractMatrix, x::AbstractVector, α)
    x1, x2, x3 = x[1], x[2], x[3]

    a = 2α
    b = 2 - 2α
    p = x1^a * x2^b
    φ = p - x3 * x3

    ρ = p / φ
    w = p / (φ * φ)
    s = p + x3 * x3

    l1 = a / x1
    l2 = b / x2

    d1 = (2ρ * a + b) / 2x1^2
    d2 = (2ρ * b + a) / 2x2^2
    c  = p * x3^2 / (φ * s)

    r1 =  d1      - c *  l1^2
    r2 = (d1 * d2 - c * (l2^2 * d1 + l1^2 * d2)) / r1

    L[1,1] = sqrt(2s) / φ            # √H₃₃
    L[2,1] = -2l1 * x3 * w / L[1,1]  # H₁₃ / L₁₁
    L[3,1] = -2l2 * x3 * w / L[1,1]  # H₂₃ / L₁₁
    L[2,2] = sqrt(r1)
    L[3,3] = sqrt(r2)
    L[3,2] = -c * l1 * l2 / L[2,2]

    L[1,2] = false
    L[1,3] = false
    L[2,3] = false

    return L
end

# apply the Hessian
#
#   f''(p) v = Pᵀ L Lᵀ P v
#
# using a pre-computed factorization
#
#   P f''(p) Pᵀ = L Lᵀ.
#
function powmul3!(out::AbstractVector{T}, L::AbstractMatrix{T}, v::AbstractVector{T}) where {T}
    v1, v2, v3 = v[1], v[2], v[3]

    r1 = L[1,1] * v3 + L[2,1] * v1 + L[3,1] * v2
    r2 =               L[2,2] * v1 + L[3,2] * v2
    r3 =                             L[3,3] * v2

    t1 = L[1,1] * r1
    t2 = L[2,1] * r1 + L[2,2] * r2
    t3 = L[3,1] * r1 + L[3,2] * r2 + L[3,3] * r3

    out[1] = t2
    out[2] = t3
    out[3] = t1

    return out
end

# materialize the scaled Hessian
#
#   H = μ f''(p) = μ Pᵀ L Lᵀ P
#
# from a pre-computed factorization
#
#   P f''(p) Pᵀ = L Lᵀ.
#
function powgram!(H::AbstractMatrix{T}, L::AbstractMatrix{T}, μ::T) where {T}
    L11 = L[1,1]
    L21 = L[2,1]; L22 = L[2,2]
    L31 = L[3,1]; L32 = L[3,2]; L33 = L[3,3]

    H[3,3] = μ * L11^2
    H[1,1] = μ * (L21^2 + L22^2)
    H[2,2] = μ * (L31^2 + L32^2 + L33^2)

    H[1,3] = H[3,1] = μ * L21 * L11
    H[2,3] = H[3,2] = μ * L31 * L11
    H[1,2] = H[2,1] = μ * (L31 * L21 + L32 * L22)

    return H
end

# compute the third-order directional derivative
#
#   f'''(p)[u]
#
# as a 3x3 matrix
function powbarrthird!(D::AbstractMatrix, x::AbstractVector, u::AbstractVector, α)
    x1, x2, x3 = x[1], x[2], x[3]
    u1, u2, u3 = u[1], u[2], u[3]

    a = 2α
    b = 2 - 2α
    p = x1^a * x2^b
    φ = p - x3 * x3

    l1 = a / x1
    l2 = b / x2
    m1 = u1 / x1
    m2 = u2 / x2

    φp1 = p * l1
    φp2 = p * l2
    φp3 = -2x3

    φdot = φp1 * u1 + φp2 * u2 + φp3 * u3

    φpp11 = p * l1 * (l1 - inv(x1))
    φpp22 = p * l2 * (l2 - inv(x2))
    φpp12 = p * l1 * l2

    φdotp1 = φpp11 * u1 + φpp12 * u2
    φdotp2 = φpp12 * u1 + φpp22 * u2
    φdotp3 = -2u3

    φ2 = φ * φ
    φ3 = φ * φ * φ

    c2 =  φdot / φ2
    c3 = 2φdot / φ3

    K  = p / φ * a * b * (a - 1) * (m2 - m1)

    D[1,1] = 2φdotp1 * φp1                 / φ2 - c3 * φp1 * φp1 - K / (x1 * x1) +  c2 * φpp11 - b * m1 / (x1 * x1)
    D[2,2] = 2φdotp2 * φp2                 / φ2 - c3 * φp2 * φp2 - K / (x2 * x2) +  c2 * φpp22 - a * m2 / (x2 * x2)
    D[2,1] = (φdotp2 * φp1 + φp2 * φdotp1) / φ2 - c3 * φp2 * φp1 + K / (x1 * x2) +  c2 * φpp12
    D[3,3] = 2φdotp3 * φp3                 / φ2 - c3 * φp3 * φp3                 - 2c2
    D[3,1] = (φdotp3 * φp1 + φp3 * φdotp1) / φ2 - c3 * φp3 * φp1
    D[3,2] = (φdotp3 * φp2 + φp3 * φdotp2) / φ2 - c3 * φp3 * φp2

    D[1,2] = D[2,1]
    D[1,3] = D[3,1]
    D[2,3] = D[3,2]

    return D
end

# Determine if x is in the power cone.
function powincone(x::AbstractVector{T}, α::T) where {T}
    x[1] > 0 && x[2] > 0 && powphi(x, α) > 0
end

# Determine if s is in the dual power cone.
function powindual(s::AbstractVector{T}, α::T) where {T}
    s[1] > 0 && s[2] > 0 && (s[1] / α)^α * (s[2] / (one(T) - α))^(one(T) - α) > abs(s[3])
end

# compute the primal "shadow" iterate p*, solving
#
#   -f'(p*) = d,
#
# using a 1-D scalar root-find on the function
#
#   g(ρ) = ρ(ρ - 1) - ¼ d₃² ((aρ + ½b) / d₁)ᵃ ((bρ + ½a) / d₂)ᵇ
#
# with
#
#   a = 2α
#   b = 2 - 2α
#
# and
#
#   p₁* = ((aρ + ½b) / d₁)
#   p₂* = ((bρ + ½a) / d₂)
#   p₃* = -d₃ / 2ρ (p₁*)ᵃ (p₂*)ᵇ
#
function powdualgrad!(sp::AbstractVector{T}, seed::T, d::AbstractVector{T}, α::T) where {T}
    d1, d2, d3 = d[1], d[2], d[3]

    a =     2α
    b = 2 - 2α

    function g(ρ)
        return ρ * (ρ - 1) - d3^2 / 4 * ((a * ρ + 1 - α) / d1)^a * ((b * ρ + α) / d2)^b
    end

    function gp(ρ)
        return (2ρ - 1) - d3^2 / 4 * ((a * ρ + 1 - α) / d1)^a * ((b * ρ + α) / d2)^b * (a^2 / (a * ρ + 1 - α) + b^2 / (b * ρ + α))
    end

    hi = 2max(seed, 1)

    while g(hi) ≤ 0
        hi *= 2
    end

    ρ = rtsafe(g, gp, one(T), hi, seed)

    sp[1] = p1 = (a * ρ + 1 - α) / d1
    sp[2] = p2 = (b * ρ +     α) / d2
    sp[3] = -d3 * (p1^a * p2^b / ρ) / 2

    return ρ
end

#
# AbstractTDCone Interface
#

function tdfact!(L::AbstractMatrix, p::AbstractVector, cache::PowerConeCache{T}) where {T}
    return powfact!(L, p, cache.cone.α)
end

function tdbarrgrad!(g::AbstractVector, p::AbstractVector, cache::PowerConeCache{T}) where {T}
    return powbarrgrad!(g, p, cache.cone.α)
end

function tdbarrthird!(D::AbstractMatrix, p::AbstractVector, u::AbstractVector, cache::PowerConeCache{T}) where {T}
    return powbarrthird!(D, p, u, cache.cone.α)
end

function tddualgrad!(sp::AbstractVector, seed, d::AbstractVector, cache::PowerConeCache{T}) where {T}
    return powdualgrad!(sp, seed, d, cache.cone.α)
end

function tdincone(p::AbstractVector, cache::PowerConeCache{T}) where {T}
    return powincone(p, cache.cone.α)
end

function tdindual(d::AbstractVector, cache::PowerConeCache{T}) where {T}
    return powindual(d, cache.cone.α)
end

# apply the Hessian
#
#   f''(p) v = Pᵀ L Lᵀ P v
#
function tdhessmul!(u::AbstractVector{T}, L::AbstractMatrix{T}, v::AbstractVector{T}, ::PowerConeCache) where {T}
    return powmul3!(u, L, v)
end

# materialize the scaled Hessian
#
#   H = μ f''(p) = μ Pᵀ L Lᵀ P
#
function tdgram!(H::AbstractMatrix{T}, L::AbstractMatrix{T}, μ::T, ::PowerConeCache) where {T}
    return powgram!(H, L, μ)
end

# solve for v in
#
#   f''(p) v = b
#
# using the factorization P f''(p) Pᵀ = L Lᵀ.
#
function tdhessldiv!(L::AbstractMatrix{T}, v::AbstractVector{T}, ::PowerConeCache) where {T}
    return powldiv3!(L, v)
end
