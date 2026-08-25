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

# Compute the "determinant"
#
#   det(x) = x₁²ᵃ x₂²⁻²ᵃ - x₃²
#
# in double precision. This function
# satisfies
#
#   f(x) = -log(det(x)) - (1 - α) log(x₁) - α log(x₂).
#
# Returns a triple
#
#   (x₁ᵃx₂ᵇ, x₁ᵃx₂ᵇ / det(x), 1 / det(x)).
#
function powdet(x1::T, x2::T, x3::T, α::T) where {T <: AbstractFloat}
    a = 2α                                        # exact; uses the stored α
    bh, bl = twosum(oftype(a, 2), -a)             # b exact in two words
    l1h, l1l = twolog(x1)
    l2h, l2l = twolog(x2)
    p1h, e1 = twoprod(l1h, a);   e1 = muladd(l1l, a, e1)
    p2h, e2 = twoprod(l2h, bh);  e2 = muladd(l2l, bh, e2)
    e2 = muladd(l2h, bl, e2)                      # fold the lo word of b
    sh, se = twosum(p1h, p2h)
    Lh, Ll = twosum(sh, se + e1 + e2)             # L = a·log x₁ + b·log x₂
    Ph, Pl = twoexp(Lh, Ll)
    qh, ql = twoprod(x3, x3)
    fh, fe = twosum(Ph, -qh)
    dh, dl = twosum(fh, fe + (Pl - ql))           # det, two words
    ih, il = twodiv(one(T), zero(T), dh, dl)
    ρh, ρe = twoprod(Ph, ih)                      # ρ = P·invdet
    ρe = muladd(Ph, il, muladd(Pl, ih, ρe))
    ρh, ρl = twosum(ρh, ρe)
    return T(Ph + Pl), T(ρh + ρl), T(ih + il)
end

# evaluate the gradient
#
#   f'(p)
#
# of the barrier function f, from the compensated determinant chain
# (ρ, invdet).
function powbarrgrad!(g::AbstractVector, x::AbstractVector, α, ρ, invdet)
    a = 2α
    b = 2 - 2α

    g[1] = -(a * ρ + 1 - α) / x[1]
    g[2] = -(b * ρ     + α) / x[2]
    g[3] =  2x[3] * invdet

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
function powfact!(L::AbstractMatrix, x::AbstractVector, α, p, ρ, invdet)
    x1, x2, x3 = x[1], x[2], x[3]

    a = 2α
    b = 2 - 2α

    w = ρ * invdet
    s = p + x3 * x3

    l1 = a / x1
    l2 = b / x2

    d1 = (2ρ * a + b) / 2x1^2
    d2 = (2ρ * b + a) / 2x2^2
    c  = ρ * x3^2 / s

    r1 =  d1      - c *  l1^2
    r2 = (d1 * d2 - c * (l2^2 * d1 + l1^2 * d2)) / r1

    flag = r1 > 0 && r2 > 0 && isfinite(r2)

    if flag
        L[1,1] = sqrt(2s) * invdet       # √H₃₃
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
function powbarrthird!(D::AbstractMatrix, x::AbstractVector, u::AbstractVector, α, p, ρ, invdet)
    x1, x2, x3 = x[1], x[2], x[3]
    u1, u2, u3 = u[1], u[2], u[3]

    a = 2α
    b = 2 - 2α

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

    invdet2 = invdet * invdet
    invdet3 = invdet2 * invdet

    c2 =  φdot * invdet2
    c3 = 2φdot * invdet3

    K  = ρ * a * b * (a - 1) * (m2 - m1)

    D[1,1] = 2φdotp1 * φp1                 * invdet2 - c3 * φp1 * φp1 - K / (x1 * x1) +  c2 * φpp11 - b * m1 / (x1 * x1)
    D[2,2] = 2φdotp2 * φp2                 * invdet2 - c3 * φp2 * φp2 - K / (x2 * x2) +  c2 * φpp22 - a * m2 / (x2 * x2)
    D[2,1] = (φdotp2 * φp1 + φp2 * φdotp1) * invdet2 - c3 * φp2 * φp1 + K / (x1 * x2) +  c2 * φpp12
    D[3,3] = 2φdotp3 * φp3                 * invdet2 - c3 * φp3 * φp3                 - 2c2
    D[3,1] = (φdotp3 * φp1 + φp3 * φdotp1) * invdet2 - c3 * φp3 * φp1
    D[3,2] = (φdotp3 * φp2 + φp3 * φdotp2) * invdet2 - c3 * φp3 * φp2

    D[1,2] = D[2,1]
    D[1,3] = D[3,1]
    D[2,3] = D[3,2]

    return D
end

# Compute the largest step to the domain boundary,
#
#   hi = min(1, first affine clamp of each strictly-signed
#             domain coordinate),
#
# for each cone's open coordinates.

function powboundprim(p::AbstractVector{T}, Δp::AbstractVector{T}) where {T}
    #
    # x₁ > 0, x₂ > 0
    #
    hi = one(T)

    @inbounds begin
        if Δp[1] < 0
            hi = min(hi, -p[1] / Δp[1])
        end

        if Δp[2] < 0
            hi = min(hi, -p[2] / Δp[2])
        end
    end

    return hi
end

function powbounddual(d::AbstractVector{T}, Δd::AbstractVector{T}) where {T}
    #
    # s₁ > 0, s₂ > 0
    #
    hi = one(T)

    @inbounds begin
        if Δd[1] < 0
            hi = min(hi, -d[1] / Δd[1])
        end

        if Δd[2] < 0
            hi = min(hi, -d[2] / Δd[2])
        end
    end

    return hi
end

# compute the primal "shadow" iterate p*, solving
#
#   -f'(p*) = d,
#
# via the scale-free normal form
#
#   1 - u = W₀ (1 + k₁u)ᵃ (1 + k₂u)ᵇ
#
# in the reciprocal variable u = 1/ρ ∈ (0, 1], with
#
#   a = 2α,  b = 2 - 2α,  k₁ = (1-α)/a,  k₂ = α/b
#
# and the constant
#
#   W₀ = ¼ d₃² (a/d₁)ᵃ (b/d₂)ᵇ
#
# computed once. W₀ < 1 is exactly dual-cone membership and
# f(0) = 1 - W₀ is the margin; on [0, u*] the increasing factor
# W stays below 1 and k₁a + k₂b = 1, so
#
#   |f'| ∈ (1, 2)   ⇒   u* ∈ [(1-W₀)/2, 1-W₀]
#
# brackets the root for free — the same slope window and bracket
# structure as the exp normal form. Every per-evaluation quantity
# is O(1) whatever the scales of d (all magnitude sits in W₀),
# with one pow per evaluation at the exact exponent a = 2α; the
# two-pow form with b̃ = fl(2-2α) carries an exponent bias and
# raced 15-30% slower.
#
# PRECONDITION: powdualmargin(d, α) > sqrt(eps(T)).
# Callers must gate on the margin (see tdscale!); the solver stays
# well behaved below the gate, but accuracy degrades as eps/margin
# and the shadow primal is numerically meaningless.
#
function powdualgrad!(sp::AbstractVector{T}, seed::T, d::AbstractVector{T}, α::T) where {T}
    @inbounds d1, d2, d3 = d[1], d[2], d[3]

    a =     2α
    b = 2 - 2α
    #
    # d₃ = 0 gives ρ = 1 and p₃* = 0 exactly: no root-find
    #
    if iszero(d3)
        @inbounds sp[1] = (a + 1 - α) / d1
        @inbounds sp[2] = (b +     α) / d2
        @inbounds sp[3] = zero(T)

        return one(T)
    end

    k1 = (1 - α) / a
    k2 =      α  / b
    t1 = a / d1
    t2 = b / d2
    #
    # compute W₀ once, in the balanced grouping
    #
    #   W₀ = ((d₃/2) t₂)² (t₁/t₂)ᵃ,
    #
    # requiring every intermediate to be NORMAL: a subnormal
    # intermediate (d₃² for |d₃| < ~1e-154, or the ratio at
    # extreme d₁/d₂) quantizes catastrophically yet can launder
    # into a normal-looking W₀. Otherwise fall back to log space,
    # whose Σ|log| eps rounding matters only beyond ~1e±150
    #
    z1  = (d3 / 2) * t2
    z1s = z1 * z1
    z2  = t1 / t2
    z2a = z2^a
    W0  = z1s * z2a

    normal = floatmin(T) <= z1s < T(Inf) &&
             floatmin(T) <= z2  < T(Inf) &&
             floatmin(T) <= z2a < T(Inf) &&
             floatmin(T) <= W0  < T(Inf)

    if !normal
        W0 = exp(2log(abs(d3)) - log(T(4)) +
                 a * (log(a) - log(d1)) + b * (log(b) - log(d2)))
    end
    #
    # bracket for free, flooring the margin f(0) = 1 - W₀ for
    # graceful out-of-contract behavior (no root for W₀ ≥ 1):
    #
    #   u* ∈ [f₀/2, f₀]
    #
    f0 = 1 - W0

    if !(f0 > 0)
        f0 = floatmin(T)
    end

    lo = f0 / 2
    hi = f0
    #
    # map the seed into u-space; outside the bracket, take a free
    # Halley step from u = 0 built from the precomputed constants
    #
    #   f(0)   = f₀
    #   f'(0)  = -(1 + W₀)
    #   f''(0) = -W₀ (1 - k₁²a - k₂²b),
    #
    # which lands in the bracket except at extreme α, where the
    # Halley denominator changes sign and 0.7 f₀ stands in
    #
    u = inv(max(seed, one(T)))

    if !(lo < u < hi)
        fp0  = -(1 + W0)
        fpp0 = -W0 * (1 - (k1^2 * a + k2^2 * b))

        u = -2f0 * fp0 / (2 * fp0^2 - f0 * fpp0)

        if !(lo < u < hi)
            u = T(0.7) * f0
        end
    end
    #
    # solve 1 - u = W₀ (1 + k₁u)ᵃ (1 + k₂u)ᵇ by one-division
    # Halley; f'' reuses W and the slope sum σ, no extra pow
    #
    for _ in 1:60
        r1 = 1 + k1 * u
        r2 = 1 + k2 * u
        W  = W0 * (r2 * r2) * (r1 / r2)^a
        fv = 1 - u - W

        σ  = k1 * a / r1 + k2 * b / r2
        fp = -1 - W * σ

        if fv > 0
            lo = u
        elseif fv < 0
            hi = u
        else
            break
        end

        σp  = -(k1^2 * a / r1^2 + k2^2 * b / r2^2)
        fpp = -W * (σ^2 + σp)
        un  = u - 2fv * fp / (2 * fp^2 - fv * fpp)

        if !(lo < un < hi)
            #
            # a step that rounds onto a bracket endpoint means
            # the root is within one ulp of it: accept and stop
            #
            if un == lo || un == hi
                u = un
                break
            end

            un = (lo + hi) / 2
        end
        #
        # f carries the absolute eps of its leading 1, so u is
        # resolvable only to ~eps absolute — which is the
        # eps/margin conditioning of ρ = 1/u, intrinsic to the
        # problem in any parametrization
        #
        converged = abs(un - u) < 4 * eps(T) * (1 + abs(un))

        u = un

        if converged
            break
        end
    end

    ρ = inv(u)
    #
    # recover p* as
    #
    #   p₁* = (aρ + 1 - α) / d₁
    #   p₂* = (bρ +     α) / d₂
    #   p₃* = -2(ρ - 1) / d₃
    #
    @inbounds sp[1] = (a * ρ + 1 - α) / d1
    @inbounds sp[2] = (b * ρ +     α) / d2
    @inbounds sp[3] = -2(ρ - one(T)) / d3

    return ρ
end

#
# AbstractTDCone Interface
#

function tddet!(cache::PowerConeCache, p::AbstractVector)
    @inbounds p1, p2, p3 = p[1], p[2], p[3]
    pv, ρ, invdet = powdet(p1, p2, p3, cache.cone.α)
    cache.d1[] = pv
    cache.d2[] = ρ
    cache.d3[] = invdet
    return cache
end

function tdfact!(L::AbstractMatrix, p::AbstractVector, cache::PowerConeCache{T}) where {T}
    return powfact!(L, p, cache.cone.α, cache.d1[], cache.d2[], cache.d3[])
end

function tdbarrgrad!(g::AbstractVector, p::AbstractVector, cache::PowerConeCache{T}) where {T}
    return powbarrgrad!(g, p, cache.cone.α, cache.d2[], cache.d3[])
end

function tdbarrthird!(D::AbstractMatrix, p::AbstractVector, u::AbstractVector, cache::PowerConeCache{T}) where {T}
    return powbarrthird!(D, p, u, cache.cone.α, cache.d1[], cache.d2[], cache.d3[])
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

# The primal power-cone jet (g, gp, gpp, noise, c), where
#
#   g = x₁^α x₂^(1-α) - |x₃|
#
# on the closed domain x₁, x₂ ≥ 0 (boundary gives finite g = -|x₃|),
# evaluated in the exact-α form
#
#   t₁ = x₂ (x₁/x₂)^α:
#
# only α — an exact input — ever appears as an exponent (one ^, not
# two). The naive two-pow form exponentiates β̃ = fl(1-α), a
# systematic relative bias
#
#   |β̃ - β| |log x₂| ≤ (ε/4) |log x₂|
#
# — up to ~175ε at extreme scales, outside any fixed noise constant;
# here β enters only linearly, where fl(1-α) is exact by Sterbenz for
# α ≥ 1/2 and a u-level perturbation otherwise. Noise 16 ε (t₁ + t₂);
# adversarial one-sided error ≤ 1.1 ε·base under the normality guard:
# refusal (noise = NaN) when a coordinate, the ratio x₁/x₂, or t₁ is
# subnormal, or the ratio is non-finite — a rounding into the
# subnormal range carries O(1) relative error which (·)^α damps only
# to ~α·O(1) ≈ 10¹²·ε, silently unsound under any ε-scale band. The
# kink at x₃ = 0 needs only a supergradient (sgn = 0 there).
#
#   g″ = -αβ t₁ (Δ₁/x₁ - Δ₂/x₂)².
#
# No overflow (powers of positives are finite). No candidates (K = 0):
# the frozen-coordinate root extractions xᵢ (t₂/t₁)^(1/exponent)
# measured 0.4-0.5 evals of value against two ^ per accept — a
# wall-time loss — and pow has no starved region for candidates to
# serve. Requires a faithfully-rounded ^. Not scale-equivariant in
# floats (2^(kα) is non-dyadic): results drift ulps under power-of-2
# input scaling, within the noise band.
#
function powjetprim(τ::T, p::AbstractVector{T}, Δp::AbstractVector{T}, α::T) where {T}
    @inbounds p1, p2, p3 = p[1], p[2], p[3]
    @inbounds Δ1, Δ2, Δ3 = Δp[1], Δp[2], Δp[3]

    β = one(T) - α

    x1 = fma(τ, Δ1, p1)
    x2 = fma(τ, Δ2, p2)
    x3 = fma(τ, Δ3, p3)

    g, gp, gpp, n = T(-Inf), T(NaN), T(NaN), zero(T)

    if x1 > 0 && x2 > 0
        #
        # exact-α: one ^
        #
        q  = x1 / x2
        t1 = x2 * q^α
        t2 = abs(x3)

        g = t1 - t2

        if iszero(x3)
            sgn = zero(T)
        else
            sgn = copysign(one(T), x3)
        end

        gp = t1 * (α * Δ1 / x1 + β * Δ2 / x2) - sgn * Δ3

        qd  = Δ1 / x1 - Δ2 / x2
        gpp = -α * β * t1 * qd * qd
        #
        # refusal: nothing certifies outside the rounding
        # model's preconditions
        #
        if issubnormal(x1) | issubnormal(x2) | issubnormal(x3) |
           issubnormal(q) | issubnormal(t1) | !isfinite(q)
            n = T(NaN)
        else
            n = 16 * eps(T) * (t1 + t2)
        end
    elseif (iszero(x1) && x2 >= 0) || (x1 >= 0 && iszero(x2))
        g = -abs(x3)

        if iszero(x3)
            sgn = zero(T)
        else
            sgn = copysign(one(T), x3)
        end

        gp  = -sgn * Δ3
        gpp = zero(T)

        if issubnormal(x1) | issubnormal(x2) | issubnormal(x3)
            n = T(NaN)
        end
    end

    return g, gp, gpp, n, T(NaN)
end

# The dual power-cone jet (g, gp, gpp, noise, c), where
#
#   g = (s₁/α)^α (s₂/(1-α))^(1-α) - |s₃|
#
# on closed s₁, s₂ ≥ 0, evaluated in the exact-α form
#
#   t₁ = u₂ (u₁/u₂)^α,   u₁ = s₁/α,   u₂ = s₂/β:
#
# only α exponentiates (one ^); β enters only linearly (see
# powjetprim for the β̃-exponent bias this removes — it also
# supersedes the former max-scaling, whose job of keeping pow-argument
# logs bounded is now done explicitly by the normality guard:
# |log(u₁/u₂)| ≤ ~708 whenever the eval certifies). Isomorphic to the
# primal under diag(1/α, 1/(1-α), 1); same noise 16 ε (t₁ + t₂), same
# guard (coordinates, the ratio u₁/u₂, t₁), same g″, same
# candidate-free K = 0 path. Adversarial one-sided error ≤ 1.6 ε·base
# under the guard.
#
# The guard is not decorative: without it, straddled scales drive the
# ratio subnormal and the 64ε-era certificate is falsifiable — e.g.
# α = 0.005, s = (9.87e251, 1.13e-60, 4.23e-59),
# Δs = (0, -2.60e-61, 4.18e-63): fl certifies g = 1.26e-72 ≥ noise at
# τ = 0 while exact g = -2.91e-74 < 0; measured 31/250 false
# certificates on that family, 0 with the guard (147 explicit
# refusals). Requires a faithfully-rounded ^.
#
function powjetdual(τ::T, d::AbstractVector{T}, Δd::AbstractVector{T}, α::T) where {T}
    @inbounds d1, d2, d3 = d[1], d[2], d[3]
    @inbounds Δ1, Δ2, Δ3 = Δd[1], Δd[2], Δd[3]

    β = one(T) - α

    s1 = fma(τ, Δ1, d1)
    s2 = fma(τ, Δ2, d2)
    s3 = fma(τ, Δ3, d3)

    g, gp, gpp, n = T(-Inf), T(NaN), T(NaN), zero(T)

    if s1 > 0 && s2 > 0
        #
        # exact-α: one ^
        #
        u1 = s1 / α
        u2 = s2 / β

        q  = u1 / u2
        t1 = u2 * q^α
        t2 = abs(s3)

        g = t1 - t2

        if iszero(s3)
            sgn = zero(T)
        else
            sgn = copysign(one(T), s3)
        end

        gp = t1 * (α * Δ1 / s1 + β * Δ2 / s2) - sgn * Δ3

        qd  = Δ1 / s1 - Δ2 / s2
        gpp = -α * β * t1 * qd * qd
        #
        # refusal: nothing certifies outside the rounding
        # model's preconditions
        #
        if issubnormal(s1) | issubnormal(s2) | issubnormal(s3) |
           issubnormal(q) | issubnormal(t1) | !isfinite(q)
            n = T(NaN)
        else
            n = 16 * eps(T) * (t1 + t2)
        end
    elseif (iszero(s1) && s2 >= 0) || (s1 >= 0 && iszero(s2))
        g = -abs(s3)

        if iszero(s3)
            sgn = zero(T)
        else
            sgn = copysign(one(T), s3)
        end

        gp  = -sgn * Δ3
        gpp = zero(T)

        if issubnormal(s1) | issubnormal(s2) | issubnormal(s3)
            n = T(NaN)
        end
    end

    return g, gp, gpp, n, T(NaN)
end

function tdjetprim(τ::T, p::AbstractVector{T}, Δp::AbstractVector{T}, cache::PowerConeCache{T}) where {T}
    return powjetprim(τ, p, Δp, cache.cone.α)
end

function tdjetdual(τ::T, d::AbstractVector{T}, Δd::AbstractVector{T}, cache::PowerConeCache{T}) where {T}
    return powjetdual(τ, d, Δd, cache.cone.α)
end
