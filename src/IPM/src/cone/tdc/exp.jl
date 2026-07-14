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
    @inbounds return x[2] * log(x[1] / x[2]) - x[3]
end

function expboundprim(p::AbstractVector{T}, Δp::AbstractVector{T}) where {T}
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

function expbounddual(d::AbstractVector{T}, Δd::AbstractVector{T}) where {T}
    @inbounds  d1,  d3 =  d[1],  d[3]
    @inbounds Δd1, Δd3 = Δd[1], Δd[3]

    hi = one(T)

    if Δd1 < 0
        hi = min(hi, -d1 / Δd1)
    end

    if Δd3 > 0
        hi = min(hi, -d3 / Δd3)
    end

    return hi
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

# Compute the pivoted Cholesky factor L
#
#   P f''(p) Pᵀ = L Lᵀ.
#
# Returns true on success, false if x is too close to the boundary.
function expfact!(L::AbstractMatrix{T}, x::AbstractVector{T}) where {T}
    x1, x2 = x[1], x[2]
    ψ = exppsi(x)

    r1 = (ψ + x2) / ψ
    r2 = (ψ + 2x2) / (ψ + x2)

    flag = ψ > 0 && r1 > 0 && r2 > 0 && isfinite(r2)

    if flag
        L[1,1] = inv(ψ)
        L[2,1] = -x2 / (x1 * ψ)
        L[3,1] = -(log(x1 / x2) - one(T)) / ψ
        L[2,2] = sqrt(r1) / x1
        L[3,2] = -inv(sqrt(ψ * (ψ + x2)))
        L[3,3] = sqrt(r2) / x2

        L[1,2] = zero(T)
        L[1,3] = zero(T)
        L[2,3] = zero(T)
    end

    return flag
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
# PRECONDITION: expdualmargin(d) > sqrt(eps(T)).
# Callers must gate on the margin (see tdscale!); below the
# threshold the shadow primal is numerically meaningless.
#
function expdualgrad!(sp::AbstractVector{T}, seed::T, d::AbstractVector{T}) where {T}
    @inbounds d1, d2, d3 = d[1], d[2], d[3]

    seed = max(seed, floatmin(T))

    function jet(p2)
        w  = one(T) - p2 * d3
        r  = (inv(p2) - d3) / d1
        h  = d3 * (log(r) - one(T)) - inv(p2) + d2
        hp = (2 - inv(w)) / p2^2
        return h, hp
    end

    h0, hp0 = jet(seed)
    #
    # bracket the function h, finding points lo and hi
    # such that
    #
    #   h(lo) < 0 < h(hi)
    #
    #
    if h0 < 0
        lo, hlo = seed, h0
        hi = seed * 2
        hhi, hphi = jet(hi)

        for _ in 1:200
            hhi >= 0 && break

            if isfinite(hphi) && hphi > 0
                hi = clamp(hi - hhi / hphi, 2hi, 16hi)
            else
                hi = 2hi
            end

            hhi, hphi = jet(hi)
        end
    else
        hi, hhi = seed, h0
        lo = seed / 2
        hlo, hplo = jet(lo)

        for _ in 1:200
            hlo < 0 && break

            if isfinite(hplo) && hplo > 0
                lo = clamp(lo - hlo / hplo, lo / 16, lo / 2)
            else
                lo = lo / 2
            end

            hlo, hplo = jet(lo)
        end
    end
    #
    # find p₂* such that
    #
    #   h(p₂*) = 0
    #
    p2 = rtsafe(jet, lo, hi, seed, h0, hp0, hlo, hhi)
    #
    # recover p* as
    #
    #   p₁* = (1 - p₂* d₃) / d₁
    #   p₃* = p₂* log(p₁* / p₂*) + 1 / d₃
    #
    s = (inv(p2) - d3) / d1

    @inbounds sp[1] = s * p2
    @inbounds sp[2] = p2
    @inbounds sp[3] = p2 * log(s) + inv(d3)

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

function tdboundprim(p::AbstractVector{T}, Δp::AbstractVector{T}, ::ExponentialConeCache) where {T}
    return expboundprim(p, Δp)
end

function tdbounddual(d::AbstractVector{T}, Δd::AbstractVector{T}, ::ExponentialConeCache) where {T}
    return expbounddual(d, Δd)
end

# evaluate the jet (g, gp, noise), where
#
#   g     = ψ(p + τΔp)
#   gp    = ⟨ψ'(p + τΔp), Δp⟩
#   noise = absolute uncertainty of g at this point; 0 if unknown
#
# ψ is the barrier argument
#
#   ψ(x) = x₂ log(x₁/x₂) - x₃
#
function expjetprim(τ::T, p::AbstractVector{T}, Δp::AbstractVector{T}) where {T}
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
        d = x1 - x2

        if abs(d) < x2 / 2
            s = log1p(d / x2)
        else
            s = log(x1 / x2)
        end

        t1 = x2 * s
        t2 = x3

        g     = t1 - t2
        gp    = Δp1 * x2 / x1 + Δp2 * (s - one(T)) - Δp3
        noise = 64eps(T) * (abs(t1) + abs(t2))
    elseif x1 >= 0 && x2 == 0
        g  = -x3
        gp = -Δp3
    end

    return g, gp, noise
end

# evaluate the jet (g, gp, noise), where
#
#   g     = φ(d + τ Δd)
#   gp    = ⟨φ'(d + τΔd), Δd⟩
#   noise = absolute uncertainty of g at this point; 0 if unknown
#
# φ is the dual barrier argument
#
#   φ(z) = e z₁ + z₃ exp(z₂/z₃)
#
function expjetdual(τ::T, d::AbstractVector{T}, Δd::AbstractVector{T}) where {T}
    g     = -T(Inf)
    gp    =  T(NaN)
    noise =  zero(T)

    @inbounds Δd1 = Δd[1]
    @inbounds Δd2 = Δd[2]
    @inbounds Δd3 = Δd[3]

    @inbounds z1 = muladd(τ, Δd1, d[1])
    @inbounds z2 = muladd(τ, Δd2, d[2])
    @inbounds z3 = muladd(τ, Δd3, d[3])

    if z1 > 0 && z3 < 0
        r = z2 / z3
        s = exp(r)

        t1 = ℯ * z1
        t2 = z3 * s

        g     = t1 + t2
        gp    = ℯ * Δd1 + Δd2 * s + Δd3 * s * (one(T) - r)
        noise = 64eps(T) * (abs(t1) + abs(t2) * (one(T) + abs(r)))
    elseif z1 >= 0 && z2 >= 0 && z3 == 0
        g  = ℯ * z1
        gp = ℯ * Δd1
    end

    return g, gp, noise
end

function tdjetprim(τ::T, p::AbstractVector{T}, Δp::AbstractVector{T}, ::ExponentialConeCache) where {T}
    return expjetprim(τ, p, Δp)
end

function tdjetdual(τ::T, d::AbstractVector{T}, Δd::AbstractVector{T}, ::ExponentialConeCache) where {T}
    return expjetdual(τ, d, Δd)
end
