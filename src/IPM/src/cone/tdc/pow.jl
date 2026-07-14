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

function powboundprim(p::AbstractVector{T}, Δp::AbstractVector{T}) where {T}
    @inbounds  p1,  p2 =  p[1],  p[2]
    @inbounds Δp1, Δp2 = Δp[1], Δp[2]

    hi = one(T)

    if Δp1 < 0
        hi = min(hi, -p1 / Δp1)
    end

    if Δp2 < 0
        hi = min(hi, -p2 / Δp2)
    end

    return hi
end

function powbounddual(d::AbstractVector{T}, Δd::AbstractVector{T}) where {T}
    @inbounds  d1,  d2 =  d[1],  d[2]
    @inbounds Δd1, Δd2 = Δd[1], Δd[2]

    hi = one(T)

    if Δd1 < 0
        hi = min(hi, -d1 / Δd1)
    end

    if Δd2 < 0
        hi = min(hi, -d2 / Δd2)
    end

    return hi
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

# factorize the Hessian
#
#   P f''(p) Pᵀ = L Lᵀ
#
# of the barrier function f, using the pivot order
# (3, 1, 2). The factor is built directly from p;
# the Hessian is never assembled.
#
# Returns true on success, false if the Schur complements
# are non-positive (p too close to the boundary).
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

    flag = r1 > 0 && r2 > 0 && isfinite(r2)

    if flag
        L[1,1] = sqrt(2s) / φ            # √H₃₃
        L[2,1] = -2l1 * x3 * w / L[1,1]  # H₁₃ / L₁₁
        L[3,1] = -2l2 * x3 * w / L[1,1]  # H₂₃ / L₁₁
        L[2,2] = sqrt(r1)
        L[3,3] = sqrt(r2)
        L[3,2] = -c * l1 * l2 / L[2,2]

        L[1,2] = false
        L[1,3] = false
        L[2,3] = false
    end

    return flag
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
# PRECONDITION: powdualmargin(d, α) > sqrt(eps(T)).
# Callers must gate on the margin (see tdscale!); below the
# threshold the root ρ ≈ 1/(2m) degenerates and the shadow
# primal is numerically meaningless.
#
function powdualgrad!(sp::AbstractVector{T}, seed::T, d::AbstractVector{T}, α::T) where {T}
    @inbounds d1, d2, d3 = d[1], d[2], d[3]

    a =     2α
    b = 2 - 2α
    c = d3^2 / 4

    seed = max(seed, one(T))

    function jet(ρ)
        t1 = (a * ρ + 1 - α) / d1
        t2 = (b * ρ +     α) / d2
        w  = c * t1^a * t2^b
        g  = ρ * (ρ - 1) - w
        gp = (2ρ - 1) - w * (a^2 / (a * ρ + 1 - α) + b^2 / (b * ρ + α))
        return g, gp
    end

    g0, gp0 = jet(seed)

    # bracket: g(lo) < 0 < g(hi), with lo = 1
    lo = one(T)
    glo, _ = jet(lo)
    hi = 2seed
    ghi, gphi = jet(hi)

    for _ in 1:200
        ghi > 0 && break

        if isfinite(gphi) && gphi > 0
            hi = clamp(hi - ghi / gphi, 2hi, 16hi)
        else
            hi = 2hi
        end

        ghi, gphi = jet(hi)
    end

    @assert isfinite(hi)

    ρ = rtsafe(jet, lo, hi, seed, g0, gp0, glo, ghi)

    @inbounds sp[1] = (a * ρ + 1 - α) / d1
    @inbounds sp[2] = (b * ρ +     α) / d2
    @inbounds sp[3] = iszero(d3) ? zero(T) : -2(ρ - one(T)) / d3

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

function tdboundprim(p::AbstractVector{T}, Δp::AbstractVector{T}, ::PowerConeCache{T}) where {T}
    return powboundprim(p, Δp)
end

function tdbounddual(d::AbstractVector{T}, Δd::AbstractVector{T}, ::PowerConeCache{T}) where {T}
    return powbounddual(d, Δd)
end

# evaluate the jet (g, gp, noise), where
#
#   g     = ψ(p + τΔp)
#   gp    = ⟨ψ'(p + τΔp), Δp⟩
#   noise = absolute uncertainty of g at this point; 0 if unknown
#
# ψ is the degree-1 barrier argument
#
#   ψ(x) = x₁^α x₂^β - |x₃|,  β = 1 - α
#
function powjetprim(τ::T, p::AbstractVector{T}, Δp::AbstractVector{T}, α::T) where {T}
    g     = -T(Inf)
    gp    =  T(NaN)
    noise =  zero(T)

    @inbounds Δp1 = Δp[1]
    @inbounds Δp2 = Δp[2]
    @inbounds Δp3 = Δp[3]

    @inbounds x1 = muladd(τ, Δp1, p[1])
    @inbounds x2 = muladd(τ, Δp2, p[2])
    @inbounds x3 = muladd(τ, Δp3, p[3])

    if x1 > 0 && x2 > 0
        β  = one(T) - α
        t1 = x1^α * x2^β
        t2 = abs(x3)

        g     = t1 - t2
        gp    = t1 * (α * Δp1 / x1 + β * Δp2 / x2) - sign(x3) * Δp3
        noise = 64eps(T) * (t1 + t2)
    elseif (x1 == 0 && x2 >= 0) || (x1 >= 0 && x2 == 0)
        g  = -abs(x3)
        gp = -sign(x3) * Δp3
    end

    return g, gp, noise
end

# evaluate the jet (g, gp, noise), where
#
#   g     = φ(d + τΔd)
#   gp    = ⟨φ'(d + τΔd), Δd⟩
#   noise = absolute uncertainty of g at this point; 0 if unknown
#
# φ is the degree-1 dual barrier argument
#
#   φ(s) = (s₁/α)^α (s₂/β)^β - |s₃|,  β = 1 - α
#
function powjetdual(τ::T, d::AbstractVector{T}, Δd::AbstractVector{T}, α::T) where {T}
    g     = -T(Inf)
    gp    =  T(NaN)
    noise =  zero(T)

    @inbounds Δd1 = Δd[1]
    @inbounds Δd2 = Δd[2]
    @inbounds Δd3 = Δd[3]

    @inbounds s1 = muladd(τ, Δd1, d[1])
    @inbounds s2 = muladd(τ, Δd2, d[2])
    @inbounds s3 = muladd(τ, Δd3, d[3])

    if s1 > 0 && s2 > 0
        β  = one(T) - α
        c  = max(s1, s2)
        t1 = c * (s1 / (α * c))^α * (s2 / (β * c))^β
        t2 = abs(s3)

        g     = t1 - t2
        gp    = t1 * (α * Δd1 / s1 + β * Δd2 / s2) - sign(s3) * Δd3
        noise = 64eps(T) * (t1 + t2)
    elseif (s1 == 0 && s2 >= 0) || (s1 >= 0 && s2 == 0)
        g  = -abs(s3)
        gp = -sign(s3) * Δd3
    end

    return g, gp, noise
end

function tdjetprim(τ::T, p::AbstractVector{T}, Δp::AbstractVector{T}, cache::PowerConeCache{T}) where {T}
    return powjetprim(τ, p, Δp, cache.cone.α)
end

function tdjetdual(τ::T, d::AbstractVector{T}, Δd::AbstractVector{T}, cache::PowerConeCache{T}) where {T}
    return powjetdual(τ, d, Δd, cache.cone.α)
end
