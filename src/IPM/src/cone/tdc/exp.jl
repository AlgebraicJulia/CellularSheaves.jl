"""
    ExponentialCone <: AbstractTDCone

The exponential cone, consisting of all triples (x, y, z)
such that x > 0, y > 0, and y log(x/y) ≥ z.
"""
struct ExponentialCone <: AbstractTDCone end

const ExponentialConeCache{T} = AbstractTDConeCache{ExponentialCone, T}

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

# Compute the "determinant"
#
#   det(x) = x₂ log(x₁/x₂) − x₃.
#
# in double precision. This function
# satisfies
#
#   f(x) = -log(det(x)) - log(x₁) - log(x₂).
#
# Returns a quadruple
#
#   (det(x), 1/det(x), log(x₁/x₂), log(x₁/x₂) - 1)
#
function expdet(x1::T, x2::T, x3::T) where {T <: AbstractFloat}
    if T(0.5) * x2 <= x1 <= 2x2
        d = x1 - x2                              # Sterbenz: exact
        qh, ql = twodiv(d, zero(T), x2, zero(T))
        lrh, lrl = twolog1p(qh, ql)
    else
        l1h, l1l = twolog(x1)
        l2h, l2l = twolog(x2)
        rh, re = twosum(l1h, -l2h)
        lrh, lrl = twosum(rh, re + (l1l - l2l))
    end

    mh, me = twosum(lrh, -one(T))
    lrm1 = T(mh + (me + lrl))                    # lr − 1 off the two-word lr
    ph, pe = twoprod(lrh, x2)
    pe = muladd(lrl, x2, pe)                     # x₂·lr, two words
    sh, se = twosum(ph, -x3)
    ψh, ψl = twosum(sh, se + pe)
    ih, il = twodiv(one(T), zero(T), ψh, ψl)
    return T(ψh + ψl), T(ih + il), T(lrh + lrl), lrm1
end

# The gradient
#
#   f'(x) = -ψ'(x)/ψ(x) - (1/x₁, 1/x₂, 0)
#
# of the barrier function, from the compensated determinant chain
# (invψ, lrm1).
function expbarrgrad!(g::AbstractVector, x::AbstractVector, invψ, lrm1)
    x1, x2 = x[1], x[2]

    g[1] = -(x2 / x1) * invψ - inv(x1)
    g[2] = -lrm1 * invψ - inv(x2)
    g[3] = invψ

    return g
end

# Compute the pivoted Cholesky factor L
#
#   P f''(p) Pᵀ = L Lᵀ.
#
# Returns true on success, false if x is too close to the boundary.
function expfact!(L::AbstractMatrix{T}, x::AbstractVector{T}, ψ, invψ, lrm1) where {T}
    x1, x2 = x[1], x[2]

    r1 = 1 + x2 * invψ           # (ψ + x₂) / ψ
    r2 = (ψ + 2x2) / (ψ + x2)

    flag = ψ > 0 && r1 > 0 && r2 > 0 && isfinite(r2)

    if flag
        L[1,1] = invψ
        L[2,1] = -x2 * invψ / x1
        L[3,1] = -lrm1 * invψ
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
function expbarrthird!(D::AbstractMatrix, x::AbstractVector, u::AbstractVector, ψ, invψ, lrm1)
    x1, x2     = x[1], x[2]
    u1, u2, u3 = u[1], u[2], u[3]

    ψg1 = x2 / x1
    ψg2 = lrm1
    ψgu = ψg1 * u1 + ψg2 * u2 - u3

    ψH11 = -x2 / x1^2
    ψH21 =  inv(x1)
    ψH22 = -inv(x2)

    ψHu1 = ψH11 * u1 + ψH21 * u2
    ψHu2 = ψH21 * u1 + ψH22 * u2

    invψ2 = invψ * invψ
    invψ3 = invψ2 * invψ

    α = -2ψgu * invψ3

    D[1,1] =  α * ψg1^2 + (2ψg1 * ψHu1 + ψH11 * ψgu) * invψ2 - (2x2 * u1 / x1^3 - u2 / x1^2) * invψ - 2u1 / x1^3
    D[2,1] =  α * ψg1 * ψg2 + (ψg1 * ψHu2 + ψg2 * ψHu1 + ψH21 * ψgu) * invψ2 + u1 / x1^2 * invψ
    D[3,1] = -α * ψg1 - ψHu1 * invψ2
    D[2,2] =  α * ψg2^2 + (2ψg2 * ψHu2 + ψH22 * ψgu) * invψ2 - u2 / x2^2 * invψ - 2u2 / x2^3
    D[3,2] = -α * ψg2 - ψHu2 * invψ2
    D[3,3] =  α

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

function expboundprim(p::AbstractVector{T}, Δp::AbstractVector{T}) where {T}
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

function expbounddual(d::AbstractVector{T}, Δd::AbstractVector{T}) where {T}
    #
    # z₁ > 0, z₃ < 0
    #
    hi = one(T)

    @inbounds begin
        if Δd[1] < 0
            hi = min(hi, -d[1] / Δd[1])
        end

        if Δd[3] > 0
            hi = min(hi, -d[3] / Δd[3])
        end
    end

    return hi
end

# compute the "shadow" primal, solving
#
#   f'(p*) = -d
#
# via the scale-free normal form
#
#   φ(v) = v + log1p(v) = R
#
# in the reciprocal variable
#
#   v = 1 / (p₂* (-d₃)),
#
# with the constant
#
#   R = h(∞) / (-d₃) = d₂ / (-d₃) - (log(-d₃/d₁) - 1)
#
# computed once. Globally φ' = 1 + 1/(1+v) ∈ (1, 2) and φ is
# concave, and
#
#   0 ≤ log1p(v) ≤ v   ⇒   v* ∈ [R/2, R]
#
# brackets the root for free, so Halley needs no bracketing loop
# and bisection is a one-line safety net (never fired in testing).
# The small-margin regime R → 0 is the well-conditioned one here,
# with v* ≈ R/2; the previous p₂-space iteration degenerated to
# bisection there (20-38 evaluations vs 1-3 for this form).
#
# PRECONDITION: expdualmargin(d) > sqrt(eps(T)).
# Callers must gate on the margin (see tdscale!); the solver stays
# well behaved below the gate, but accuracy degrades as eps/margin
# and the shadow primal is numerically meaningless.
#
function expdualgrad!(sp::AbstractVector{T}, seed::T, d::AbstractVector{T}) where {T}
    @inbounds d1, d2, d3 = d[1], d[2], d[3]

    #
    # d₃ = 0 degenerates h to -1/p₂ + d₂: root p₂ = 1/d₂, no
    # root-find; sp[3] is the d₃ → 0⁻ limit -Inf (the incumbent's
    # inv(d₃) gave the wrong-side +Inf on a positive zero)
    #
    if iszero(d3)
        p2 = inv(d2)
        s  = d2 / d1

        @inbounds sp[1] = s * p2
        @inbounds sp[2] = p2
        @inbounds sp[3] = T(-Inf)

        return p2
    end

    md3 = -d3
    #
    # compute the constant
    #
    #   R = d₂ / (-d₃) - (log(-d₃/d₁) - 1),
    #
    # splitting the log when the ratio leaves the normal range,
    # and flooring R for graceful out-of-contract behavior (no
    # root exists for R ≤ 0)
    #
    rat = md3 / d1

    if floatmin(T) <= rat < T(Inf)
        lr = log(rat)
    else
        lr = log(md3) - log(d1)
    end

    R = d2 / md3 - (lr - one(T))

    if !(R > 0)
        R = floatmin(T)
    end
    #
    # bracket for free:
    #
    #   v* ∈ [R/2, R]
    #
    lo = R / 2
    hi = R
    #
    # map the seed into v-space; outside the bracket, start from
    # the series / asymptotic inverse of φ
    #
    #   v₀ = R/2 + R²/16     (R ≤ 2)
    #   v₀ = R - log R       (R > 2)
    #
    # which lies in the bracket for all finite R except when
    # log R < ulp(R) rounds v₀ onto R, where 0.7 R stands in
    #
    v = inv(max(seed, floatmin(T)) * md3)

    if !(lo < v < hi)
        if R <= 2
            v = R / 2 + R^2 / 16
        else
            v = R - log(R)
        end

        if !(lo < v < hi)
            v = T(0.7) * R
        end
    end
    #
    # solve φ(v) = R by one-division Halley; φ'' = -(φ'-1)² is
    # free from φ'
    #
    for _ in 1:60
        fv = v + log1p(v) - R
        fp = one(T) + inv(one(T) + v)

        if fv < 0
            lo = v
        elseif fv > 0
            hi = v
        else
            break
        end

        fpp = -(fp - one(T))^2
        vn  = v - 2fv * fp / (2 * fp^2 - fv * fpp)

        if !(lo < vn < hi)
            #
            # a step that rounds onto a bracket endpoint means
            # the root is within one ulp of it: accept and stop
            #
            if vn == lo || vn == hi
                v = vn
                break
            end

            vn = (lo + hi) / 2
        end
        #
        # the form is scale-free in v, so the stop is purely
        # relative; accuracy of p₂* is then limited only by the
        # conditioning of R (≈ eps/margin, intrinsic)
        #
        converged = abs(vn - v) < 4 * eps(T) * abs(vn)

        v = vn

        if converged
            break
        end
    end
    #
    # recover p* as
    #
    #   p₂* = 1 / ((-d₃) v)
    #   p₁* / p₂* = (-d₃)(v + 1) / d₁
    #   p₃* = p₂* log(p₁*/p₂*) + 1/d₃
    #
    q  = md3 * v
    p2 = inv(q)

    s = md3 * (v + one(T)) / d1

    if !(zero(T) < s < T(Inf))
        s = exp(lr + log1p(v))
    end

    @inbounds sp[1] = s * p2
    @inbounds sp[2] = p2
    @inbounds sp[3] = p2 * log(s) + inv(d3)

    return p2
end

#
# AbstractTDCone Interface
#

function tddet!(cache::ExponentialConeCache, x::AbstractVector)
    @inbounds x1, x2, x3 = x[1], x[2], x[3]
    ψ, invψ, lr, lrm1 = expdet(x1, x2, x3)
    cache.d1[] = ψ
    cache.d2[] = invψ
    cache.d3[] = lr
    cache.d4[] = lrm1
    return cache
end

function tdfact!(L::AbstractMatrix, p::AbstractVector, cache::ExponentialConeCache)
    return expfact!(L, p, cache.d1[], cache.d2[], cache.d4[])
end

function tdbarrgrad!(g::AbstractVector, p::AbstractVector, cache::ExponentialConeCache)
    return expbarrgrad!(g, p, cache.d2[], cache.d4[])
end

function tdbarrthird!(D::AbstractMatrix, p::AbstractVector, u::AbstractVector, cache::ExponentialConeCache)
    return expbarrthird!(D, p, u, cache.d1[], cache.d2[], cache.d4[])
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

# The primal exp-cone jet (g, gp, gpp, noise, c), where
#
#   g = x₂ log(x₁/x₂) - x₃
#
# on x₁, x₂ > 0, with the closure ray x₂ = 0, x₁ ≥ 0 carrying
#
#   g = -x₃.
#
# The noise needs the input-amplification floor
#
#   16 ε (x₂ (1 + |s|) + |x₃|)
#
# (the |t₁| + |t₂| form certifies exactly infeasible points at
# x₁ ≈ x₂, x₃ ≈ 0 crossings); adversarial one-sided error
# ≤ 1.0 ε·base. When the ratio x₁/x₂ leaves the normal range (rounds
# to 0, subnormal, or Inf across ~616 straddled decades) its one
# rounding carries up to O(1) relative error, so the jet switches to
# the split form
#
#   s = log x₁ - log x₂
#
# — no ratio rounding, no spurious -Inf, and the |s| ≥ 708 there
# absorbs the split's cancellation within the band. Refusal
# (noise = NaN) only on subnormal or non-finite coordinates.
#
#   g″ = -W² / (x₁² x₂),   W = Δ₁x₂ - Δ₂x₁  (constant Wronskian).
#
# Requires a faithfully-rounded log. No overflow (log cannot
# overflow); no candidates (K = 0): the x₁-solve bought 0.8 evals for
# one exp per accept — a time wash inside measurement noise — and the
# free x₃-solve measured inert; both dropped for parsimony.
#
function expjetprim(τ::T, p::AbstractVector{T}, Δp::AbstractVector{T}) where {T}
    @inbounds p1, p2, p3 = p[1], p[2], p[3]
    @inbounds Δ1, Δ2, Δ3 = Δp[1], Δp[2], Δp[3]

    x1 = fma(τ, Δ1, p1)
    x2 = fma(τ, Δ2, p2)
    x3 = fma(τ, Δ3, p3)

    g, gp, gpp, n = T(-Inf), T(NaN), T(NaN), zero(T)

    if x1 > 0 && x2 > 0
        #
        # split the log when the ratio leaves the normal range
        #
        q = x1 / x2

        if floatmin(T) <= q < T(Inf)
            s = log(q)
        else
            s = log(x1) - log(x2)
        end

        t1 = x2 * s

        g  = t1 - x3
        gp = Δ1 * x2 / x1 + Δ2 * (s - one(T)) - Δ3

        W   = Δ1 * x2 - Δ2 * x1
        gpp = -(W / x1) * (W / x1) / x2
        #
        # refusal: nothing certifies outside the rounding
        # model's preconditions
        #
        if issubnormal(x1) | issubnormal(x2) | issubnormal(x3) |
           !isfinite(x1) | !isfinite(x2)
            n = T(NaN)
        else
            n = 16 * eps(T) * (x2 * (one(T) + abs(s)) + abs(x3))
        end
    elseif x1 >= 0 && iszero(x2)
        g   = -x3
        gp  = -Δ3
        gpp = zero(T)

        if issubnormal(x1) | issubnormal(x3)
            n = T(NaN)
        end
    end

    return g, gp, gpp, n, T(NaN)
end

# The dual exp-cone jet (g, gp, gpp, noise, c), where
#
#   g = e z₁ + z₃ e^(z₂/z₃)
#
# on z₁ > 0, z₃ < 0, with the closure ray z₃ = 0, z₁ ≥ 0, z₂ ≥ 0
# carrying
#
#   g = e z₁.
#
# The noise
#
#   16 ε (t₁ + |t₂| (1 + |r|))
#
# needs no floor; adversarial one-sided error ≤ 1.4 ε·base. Refusal
# (noise = NaN) on subnormal coordinates; a subnormal ratio r needs
# none — exp is absolutely insensitive at 0, so r's O(1) relative
# error there is harmless.
#
#   g″ = s V² / z₃,   V = Δ₂ - Δ₃r.
#
# The essential face (z₂(τ) < 0 < z₁(τ)) carries an overflow strip
# where naive fl(g) = -Inf; exact g there may be finite, or even
# positive (dom fl ⊊ dom cl). The jet continues such evals in log
# coordinates,
#
#   g̃ = t₁ h,   h = log t₁ - r - log(-z₃),
#
# one extra log since log t₁ = 1 + log z₁ — which never overflows,
# has the same root, and matches g to first order there (gpp = NaN on
# this path: the surrogate's curvature is not g's, so the upper step
# runs Newton). Strip rejects carry usable derivatives; strip accepts
# are sound under the log path's noise
#
#   16 ε t₁ (1 + |log t₁| + |r| + |log(-z₃)|)
#
# (the wide logs, up to ~1400, coarsen certification there by ~3
# decades — still certifiable, where naive fl certifies nothing).
# Measured on strips containing exact τ*: witness gaps fell from the
# strip width (~10¹² tol units) to ~10¹, and evals from ~50 (bisecting
# a region that cannot certify) to ~9; on ordinary essential-face
# strips, ~4% fewer evals.
#
# One candidate: the Lambert z₃-solve (2 Newton steps on
#
#   log y - y = log a;
#
# trajectory-equivalent to 4 steps within 0.03 evals, and the
# construction runs on most essential-face accepts, where the 4-step
# form measured 2.4× whole-solve time — a candidate only needs to beat
# the midpoint, not converge). It is the strip-jumper: dropping it
# costs 8.6 → 23 evals on strips containing exact τ*. The z₂-solve
# (0.2 evals against a log per accept) and the free z₁-solve (marginal
# value given the Lambert ≤ 0.11 evals — inert) were dropped.
#
function expjetdual(τ::T, d::AbstractVector{T}, Δd::AbstractVector{T}) where {T}
    @inbounds d1, d2, d3 = d[1], d[2], d[3]
    @inbounds Δ1, Δ2, Δ3 = Δd[1], Δd[2], Δd[3]

    e = T(ℯ)

    z1 = fma(τ, Δ1, d1)
    z2 = fma(τ, Δ2, d2)
    z3 = fma(τ, Δ3, d3)

    g, gp, gpp, n = T(-Inf), T(NaN), T(NaN), zero(T)
    c3 = T(NaN)

    if z1 > 0 && z3 < 0
        r = z2 / z3
        #
        # exp overflow gives s = Inf and g = -Inf, routed to the
        # log-coordinate continuation below
        #
        s = exp(r)

        t1 = e * z1
        t2 = z3 * s

        g = t1 + t2

        bad = issubnormal(z1) | issubnormal(z2) | issubnormal(z3)

        if g == T(-Inf)
            #
            # fl overflow inside the domain: continue in log
            # coordinates, where g̃ = t₁ h matches g to first
            # order at the root; the surrogate carries no
            # trustworthy curvature, so gpp stays NaN
            #
            lg  = one(T) + log(z1)
            lt3 = log(-z3)
            h   = lg - r - lt3

            g  = t1 * h
            rp = (Δ2 - Δ3 * r) / z3
            gp = e * Δ1 * h + t1 * (Δ1 / z1 - rp - Δ3 / z3)

            if bad
                n = T(NaN)
            else
                n = 16 * eps(T) * t1 *
                    (one(T) + abs(lg) + abs(r) + abs(lt3))
            end
        else
            gp = e * Δ1 + Δ2 * s + Δ3 * s * (one(T) - r)

            V   = Δ2 - Δ3 * r
            gpp = s * V * V / z3

            if bad
                #
                # refusal: nothing certifies
                #
                n = T(NaN)
            else
                n = 16 * eps(T) * (t1 + abs(t2) * (one(T) + abs(r)))
                #
                # the Lambert z₃-solve, built only on a certified
                # accept with z₂ < 0 and a = -z₂/t₁ < 1/e
                #
                if g >= n
                    if z2 < 0
                        a = -z2 / t1

                        if a < T(0.36787944117144233)
                            la = log(a)

                            if -la > one(T)
                                y = -la + log(-la)
                            else
                                y = T(1.5)
                            end

                            for _ in 1:2
                                y -= (log(y) - y - la) * y / (one(T) - y)
                            end

                            if y > one(T)
                                c3 = (z2 / y - d3) / Δ3
                            end
                        end
                    end
                end
            end
        end
    elseif z1 >= 0 && z2 >= 0 && iszero(z3)
        g   = e * z1
        gp  = e * Δ1
        gpp = zero(T)

        if issubnormal(z1) | issubnormal(z2)
            n = T(NaN)
        end
    end

    return g, gp, gpp, n, c3
end

function tdjetprim(τ::T, p::AbstractVector{T}, Δp::AbstractVector{T}, ::ExponentialConeCache) where {T}
    return expjetprim(τ, p, Δp)
end

function tdjetdual(τ::T, d::AbstractVector{T}, Δd::AbstractVector{T}, ::ExponentialConeCache) where {T}
    return expjetdual(τ, d, Δd)
end

