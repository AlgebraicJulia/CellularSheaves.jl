#
# SVD and eigenvalue computations for small matrices
#
# svd1!/2!/3!/4!: one-sided Jacobi SVD computing only U and σ (not V)
# eigmin1/2/3/4: analytic minimum eigenvalue for symmetric matrices
# Sturm bisection for tridiagonal, inverse iteration, Kato-Temple certified Ritz bounds.
#

# ============================================================================
# One-sided Jacobi SVD for small matrices
# ============================================================================
#
# svd1!/2!/3!/4!(U, σ, A, B) computes A = U Σ Vᵀ but only returns U (n×n) and
# σ (descending singular values). V is never formed. B is an n×n work buffer.
#
# Method: one-sided Jacobi on Aᵀ. Rotations J chosen so that AᵀJ has orthogonal
# columns; then A = J Σ Wᵀ, so U = J and σᵢ = ‖(AᵀJ)[:,i]‖.
#

@propagate_inbounds function svd1!(U::AbstractMatrix{T}, σ::AbstractVector{T}, A::AbstractMatrix{T}) where T
    @boundscheck checkbounds(A, 1, 1)
    @boundscheck checkbounds(U, 1, 1)
    @boundscheck checkbounds(σ, 1)
    @inbounds a = A[1,1]
    @inbounds σ[1] = abs(a)
    @inbounds U[1,1] = a ≥ zero(T) ? one(T) : -one(T)
    return U, σ
end

@propagate_inbounds function svd2!(U::AbstractMatrix{T}, σ::AbstractVector{T}, A::AbstractMatrix{T}) where T
    @boundscheck checkbounds(A, 2, 2)
    @boundscheck checkbounds(U, 2, 2)
    @boundscheck checkbounds(σ, 2)

    @inbounds m = max(abs(A[1,1]), abs(A[1,2]), abs(A[2,1]), abs(A[2,2]))

    if iszero(m)
        σ₁ = σ₂ = zero(T)
        U₁₁ = U₂₂ = one(T)
        U₁₂ = U₂₁ = zero(T)
    else
        @inbounds B₁₁ = A[1,1] / m
        @inbounds B₂₁ = A[1,2] / m
        @inbounds B₁₂ = A[2,1] / m
        @inbounds B₂₂ = A[2,2] / m

        α = muladd(B₁₁, B₁₁, B₂₁*B₂₁)
        β = muladd(B₁₂, B₁₂, B₂₂*B₂₂)
        γ = muladd(B₁₁, B₁₂, B₂₁*B₂₂)

        flip = α < β

        if flip
            B₁₁, B₁₂ = B₁₂, B₁₁
            B₂₁, B₂₂ = B₂₂, B₂₁
            α, β = β, α
        end

        if γ*γ > (T(16) * eps(T)^2) * (α * β)
            ζ = (β - α) / (2γ)
            t = copysign(one(T), ζ) / (abs(ζ) + sqrt(muladd(ζ, ζ, one(T))))
            c = inv(sqrt(muladd(t, t, one(T))))
            s = c * t

            σ₁ = m * sqrt(muladd(-t, γ, α))
            r₁₂ = muladd(s, B₁₁, c*B₁₂)
            r₂₂ = muladd(s, B₂₁, c*B₂₂)
            σ₂ = m * hypot(r₁₂, r₂₂)

            U₁₁ = c; U₂₁ = -s; U₁₂ = s; U₂₂ = c
        else
            σ₁ = m * sqrt(α)
            σ₂ = m * hypot(B₁₂, B₂₂)

            U₁₁ = one(T); U₂₁ = zero(T); U₁₂ = zero(T); U₂₂ = one(T)
        end

        if flip
            U₁₁, U₂₁ = U₂₁, U₁₁
            U₁₂, U₂₂ = U₂₂, U₁₂
        end
    end

    @inbounds σ[1] = σ₁
    @inbounds σ[2] = σ₂
    @inbounds U[1,1] = U₁₁
    @inbounds U[1,2] = U₁₂
    @inbounds U[2,1] = U₂₁
    @inbounds U[2,2] = U₂₂

    return U, σ
end

# Rotate one column pair of B (and U) if their columns aren't orthogonal yet.
# α, β, γ are computed FRESH from B's entries — this is what preserves relative
# accuracy for small singular values; never cache/update these across rotations.
# Pure function on scalars: no `!`, everything stays in registers.
@inline function svd3pair(Bp1::T, Bp2::T, Bp3::T, Bq1::T, Bq2::T, Bq3::T,
                          Up1::T, Up2::T, Up3::T, Uq1::T, Uq2::T, Uq3::T) where T
    α = muladd(Bp1, Bp1, muladd(Bp2, Bp2, Bp3*Bp3))
    β = muladd(Bq1, Bq1, muladd(Bq2, Bq2, Bq3*Bq3))
    γ = muladd(Bp1, Bq1, muladd(Bp2, Bq2, Bp3*Bq3))

    flag = γ*γ > (T(16) * eps(T)^2) * (α * β)

    if flag
        ζ = (β - α) / (2γ)
        t = copysign(one(T), ζ) / (abs(ζ) + sqrt(muladd(ζ, ζ, one(T))))
        c = inv(sqrt(muladd(t, t, one(T))))
        s = c * t

        Bp1, Bp2, Bp3, Bq1, Bq2, Bq3, Up1, Up2, Up3, Uq1, Uq2, Uq3 = (
            muladd(c, Bp1, -s*Bq1), muladd(c, Bp2, -s*Bq2), muladd(c, Bp3, -s*Bq3),
            muladd(s, Bp1, c*Bq1), muladd(s, Bp2, c*Bq2), muladd(s, Bp3, c*Bq3),
            muladd(c, Up1, -s*Uq1), muladd(c, Up2, -s*Uq2), muladd(c, Up3, -s*Uq3),
            muladd(s, Up1, c*Uq1), muladd(s, Up2, c*Uq2), muladd(s, Up3, c*Uq3),
        )
    end

    return Bp1, Bp2, Bp3, Bq1, Bq2, Bq3, Up1, Up2, Up3, Uq1, Uq2, Uq3, flag
end

@propagate_inbounds function svd3!(U::AbstractMatrix{T}, σ::AbstractVector{T}, A::AbstractMatrix{T}; maxsweep::Int = 5) where T
    @boundscheck checkbounds(A, 3, 3)
    @boundscheck checkbounds(U, 3, 3)
    @boundscheck checkbounds(σ, 3)

    @inbounds m = max(abs(A[1,1]), abs(A[1,2]), abs(A[1,3]),
                      abs(A[2,1]), abs(A[2,2]), abs(A[2,3]),
                      abs(A[3,1]), abs(A[3,2]), abs(A[3,3]))

    if iszero(m)
        σ₁ = σ₂ = σ₃ = zero(T)
        U₁₁ = U₂₂ = U₃₃ = one(T)
        U₁₂ = U₁₃ = U₂₁ = U₂₃ = U₃₁ = U₃₂ = zero(T)
    else
        @inbounds B₁₁ = A[1,1]/m
        @inbounds B₂₁ = A[1,2]/m
        @inbounds B₃₁ = A[1,3]/m
        @inbounds B₁₂ = A[2,1]/m
        @inbounds B₂₂ = A[2,2]/m
        @inbounds B₃₂ = A[2,3]/m
        @inbounds B₁₃ = A[3,1]/m
        @inbounds B₂₃ = A[3,2]/m
        @inbounds B₃₃ = A[3,3]/m

        U₁₁ = U₂₂ = U₃₃ = one(T)
        U₁₂ = U₁₃ = U₂₁ = U₂₃ = U₃₁ = U₃₂ = zero(T)

        for _ in 1:maxsweep
            B₁₁, B₂₁, B₃₁, B₁₂, B₂₂, B₃₂, U₁₁, U₂₁, U₃₁, U₁₂, U₂₂, U₃₂, r₁₂ =
                svd3pair(B₁₁, B₂₁, B₃₁, B₁₂, B₂₂, B₃₂, U₁₁, U₂₁, U₃₁, U₁₂, U₂₂, U₃₂)

            B₁₁, B₂₁, B₃₁, B₁₃, B₂₃, B₃₃, U₁₁, U₂₁, U₃₁, U₁₃, U₂₃, U₃₃, r₁₃ =
                svd3pair(B₁₁, B₂₁, B₃₁, B₁₃, B₂₃, B₃₃, U₁₁, U₂₁, U₃₁, U₁₃, U₂₃, U₃₃)

            B₁₂, B₂₂, B₃₂, B₁₃, B₂₃, B₃₃, U₁₂, U₂₂, U₃₂, U₁₃, U₂₃, U₃₃, r₂₃ =
                svd3pair(B₁₂, B₂₂, B₃₂, B₁₃, B₂₃, B₃₃, U₁₂, U₂₂, U₃₂, U₁₃, U₂₃, U₃₃)

            (r₁₂ | r₁₃ | r₂₃) || break
        end

        # σ from B's column norms — the entries hold small columns at full
        # relative precision; cached Gram quantities do not.
        σ₁ = m * sqrt(B₁₁*B₁₁ + B₂₁*B₂₁ + B₃₁*B₃₁)
        σ₂ = m * sqrt(B₁₂*B₁₂ + B₂₂*B₂₂ + B₃₂*B₃₂)
        σ₃ = m * sqrt(B₁₃*B₁₃ + B₂₃*B₂₃ + B₃₃*B₃₃)

        if σ₂ > σ₁
            σ₁, σ₂ = σ₂, σ₁
            U₁₁, U₁₂ = U₁₂, U₁₁
            U₂₁, U₂₂ = U₂₂, U₂₁
            U₃₁, U₃₂ = U₃₂, U₃₁
        end

        if σ₃ > σ₂
            σ₂, σ₃ = σ₃, σ₂
            U₁₂, U₁₃ = U₁₃, U₁₂
            U₂₂, U₂₃ = U₂₃, U₂₂
            U₃₂, U₃₃ = U₃₃, U₃₂
        end

        if σ₂ > σ₁
            σ₁, σ₂ = σ₂, σ₁
            U₁₁, U₁₂ = U₁₂, U₁₁
            U₂₁, U₂₂ = U₂₂, U₂₁
            U₃₁, U₃₂ = U₃₂, U₃₁
        end
    end

    @inbounds σ[1] = σ₁
    @inbounds σ[2] = σ₂
    @inbounds σ[3] = σ₃
    @inbounds U[1,1] = U₁₁
    @inbounds U[1,2] = U₁₂
    @inbounds U[1,3] = U₁₃
    @inbounds U[2,1] = U₂₁
    @inbounds U[2,2] = U₂₂
    @inbounds U[2,3] = U₂₃
    @inbounds U[3,1] = U₃₁
    @inbounds U[3,2] = U₃₂
    @inbounds U[3,3] = U₃₃

    return U, σ
end

@inline function svd4pair(Bp1::T, Bp2::T, Bp3::T, Bp4::T, Bq1::T, Bq2::T, Bq3::T, Bq4::T,
                          Up1::T, Up2::T, Up3::T, Up4::T, Uq1::T, Uq2::T, Uq3::T, Uq4::T) where T
    α = muladd(Bp1, Bp1, muladd(Bp2, Bp2, muladd(Bp3, Bp3, Bp4*Bp4)))
    β = muladd(Bq1, Bq1, muladd(Bq2, Bq2, muladd(Bq3, Bq3, Bq4*Bq4)))
    γ = muladd(Bp1, Bq1, muladd(Bp2, Bq2, muladd(Bp3, Bq3, Bp4*Bq4)))

    flag = γ*γ > (T(16) * eps(T)^2) * (α * β)

    if flag
        ζ = (β - α) / (2γ)
        t = copysign(one(T), ζ) / (abs(ζ) + sqrt(muladd(ζ, ζ, one(T))))
        c = inv(sqrt(muladd(t, t, one(T))))
        s = c * t

        Bp1, Bp2, Bp3, Bp4, Bq1, Bq2, Bq3, Bq4, Up1, Up2, Up3, Up4, Uq1, Uq2, Uq3, Uq4 = (
            muladd(c, Bp1, -s*Bq1), muladd(c, Bp2, -s*Bq2), muladd(c, Bp3, -s*Bq3), muladd(c, Bp4, -s*Bq4),
            muladd(s, Bp1, c*Bq1), muladd(s, Bp2, c*Bq2), muladd(s, Bp3, c*Bq3), muladd(s, Bp4, c*Bq4),
            muladd(c, Up1, -s*Uq1), muladd(c, Up2, -s*Uq2), muladd(c, Up3, -s*Uq3), muladd(c, Up4, -s*Uq4),
            muladd(s, Up1, c*Uq1), muladd(s, Up2, c*Uq2), muladd(s, Up3, c*Uq3), muladd(s, Up4, c*Uq4),
        )
    end

    return Bp1, Bp2, Bp3, Bp4, Bq1, Bq2, Bq3, Bq4, Up1, Up2, Up3, Up4, Uq1, Uq2, Uq3, Uq4, flag
end

@propagate_inbounds function svd4!(U::AbstractMatrix{T}, σ::AbstractVector{T}, A::AbstractMatrix{T}; maxsweep::Int = 8) where T
    @boundscheck checkbounds(A, 4, 4)
    @boundscheck checkbounds(U, 4, 4)
    @boundscheck checkbounds(σ, 4)

    @inbounds m = max(abs(A[1,1]), abs(A[1,2]), abs(A[1,3]), abs(A[1,4]),
                      abs(A[2,1]), abs(A[2,2]), abs(A[2,3]), abs(A[2,4]),
                      abs(A[3,1]), abs(A[3,2]), abs(A[3,3]), abs(A[3,4]),
                      abs(A[4,1]), abs(A[4,2]), abs(A[4,3]), abs(A[4,4]))

    if iszero(m)
        σ₁ = σ₂ = σ₃ = σ₄ = zero(T)
        U₁₁ = U₂₂ = U₃₃ = U₄₄ = one(T)
        U₁₂ = U₁₃ = U₁₄ = U₂₁ = U₂₃ = U₂₄ = U₃₁ = U₃₂ = U₃₄ = U₄₁ = U₄₂ = U₄₃ = zero(T)
    else
        @inbounds B₁₁ = A[1,1]/m
        @inbounds B₂₁ = A[1,2]/m
        @inbounds B₃₁ = A[1,3]/m
        @inbounds B₄₁ = A[1,4]/m
        @inbounds B₁₂ = A[2,1]/m
        @inbounds B₂₂ = A[2,2]/m
        @inbounds B₃₂ = A[2,3]/m
        @inbounds B₄₂ = A[2,4]/m
        @inbounds B₁₃ = A[3,1]/m
        @inbounds B₂₃ = A[3,2]/m
        @inbounds B₃₃ = A[3,3]/m
        @inbounds B₄₃ = A[3,4]/m
        @inbounds B₁₄ = A[4,1]/m
        @inbounds B₂₄ = A[4,2]/m
        @inbounds B₃₄ = A[4,3]/m
        @inbounds B₄₄ = A[4,4]/m

        U₁₁ = U₂₂ = U₃₃ = U₄₄ = one(T)
        U₁₂ = U₁₃ = U₁₄ = U₂₁ = U₂₃ = U₂₄ = U₃₁ = U₃₂ = U₃₄ = U₄₁ = U₄₂ = U₄₃ = zero(T)

        for _ in 1:maxsweep
            B₁₁, B₂₁, B₃₁, B₄₁, B₁₂, B₂₂, B₃₂, B₄₂, U₁₁, U₂₁, U₃₁, U₄₁, U₁₂, U₂₂, U₃₂, U₄₂, r₁₂ =
                svd4pair(B₁₁, B₂₁, B₃₁, B₄₁, B₁₂, B₂₂, B₃₂, B₄₂, U₁₁, U₂₁, U₃₁, U₄₁, U₁₂, U₂₂, U₃₂, U₄₂)

            B₁₁, B₂₁, B₃₁, B₄₁, B₁₃, B₂₃, B₃₃, B₄₃, U₁₁, U₂₁, U₃₁, U₄₁, U₁₃, U₂₃, U₃₃, U₄₃, r₁₃ =
                svd4pair(B₁₁, B₂₁, B₃₁, B₄₁, B₁₃, B₂₃, B₃₃, B₄₃, U₁₁, U₂₁, U₃₁, U₄₁, U₁₃, U₂₃, U₃₃, U₄₃)

            B₁₁, B₂₁, B₃₁, B₄₁, B₁₄, B₂₄, B₃₄, B₄₄, U₁₁, U₂₁, U₃₁, U₄₁, U₁₄, U₂₄, U₃₄, U₄₄, r₁₄ =
                svd4pair(B₁₁, B₂₁, B₃₁, B₄₁, B₁₄, B₂₄, B₃₄, B₄₄, U₁₁, U₂₁, U₃₁, U₄₁, U₁₄, U₂₄, U₃₄, U₄₄)

            B₁₂, B₂₂, B₃₂, B₄₂, B₁₃, B₂₃, B₃₃, B₄₃, U₁₂, U₂₂, U₃₂, U₄₂, U₁₃, U₂₃, U₃₃, U₄₃, r₂₃ =
                svd4pair(B₁₂, B₂₂, B₃₂, B₄₂, B₁₃, B₂₃, B₃₃, B₄₃, U₁₂, U₂₂, U₃₂, U₄₂, U₁₃, U₂₃, U₃₃, U₄₃)

            B₁₂, B₂₂, B₃₂, B₄₂, B₁₄, B₂₄, B₃₄, B₄₄, U₁₂, U₂₂, U₃₂, U₄₂, U₁₄, U₂₄, U₃₄, U₄₄, r₂₄ =
                svd4pair(B₁₂, B₂₂, B₃₂, B₄₂, B₁₄, B₂₄, B₃₄, B₄₄, U₁₂, U₂₂, U₃₂, U₄₂, U₁₄, U₂₄, U₃₄, U₄₄)

            B₁₃, B₂₃, B₃₃, B₄₃, B₁₄, B₂₄, B₃₄, B₄₄, U₁₃, U₂₃, U₃₃, U₄₃, U₁₄, U₂₄, U₃₄, U₄₄, r₃₄ =
                svd4pair(B₁₃, B₂₃, B₃₃, B₄₃, B₁₄, B₂₄, B₃₄, B₄₄, U₁₃, U₂₃, U₃₃, U₄₃, U₁₄, U₂₄, U₃₄, U₄₄)

            (r₁₂ | r₁₃ | r₁₄ | r₂₃ | r₂₄ | r₃₄) || break
        end

        σ₁ = m * sqrt(B₁₁*B₁₁ + B₂₁*B₂₁ + B₃₁*B₃₁ + B₄₁*B₄₁)
        σ₂ = m * sqrt(B₁₂*B₁₂ + B₂₂*B₂₂ + B₃₂*B₃₂ + B₄₂*B₄₂)
        σ₃ = m * sqrt(B₁₃*B₁₃ + B₂₃*B₂₃ + B₃₃*B₃₃ + B₄₃*B₄₃)
        σ₄ = m * sqrt(B₁₄*B₁₄ + B₂₄*B₂₄ + B₃₄*B₃₄ + B₄₄*B₄₄)

        # Sort descending: optimal 5-comparator network (1,2)(3,4)(1,3)(2,4)(2,3)
        if σ₂ > σ₁
            σ₁, σ₂ = σ₂, σ₁
            U₁₁, U₁₂ = U₁₂, U₁₁; U₂₁, U₂₂ = U₂₂, U₂₁; U₃₁, U₃₂ = U₃₂, U₃₁; U₄₁, U₄₂ = U₄₂, U₄₁
        end
        if σ₄ > σ₃
            σ₃, σ₄ = σ₄, σ₃
            U₁₃, U₁₄ = U₁₄, U₁₃; U₂₃, U₂₄ = U₂₄, U₂₃; U₃₃, U₃₄ = U₃₄, U₃₃; U₄₃, U₄₄ = U₄₄, U₄₃
        end
        if σ₃ > σ₁
            σ₁, σ₃ = σ₃, σ₁
            U₁₁, U₁₃ = U₁₃, U₁₁; U₂₁, U₂₃ = U₂₃, U₂₁; U₃₁, U₃₃ = U₃₃, U₃₁; U₄₁, U₄₃ = U₄₃, U₄₁
        end
        if σ₄ > σ₂
            σ₂, σ₄ = σ₄, σ₂
            U₁₂, U₁₄ = U₁₄, U₁₂; U₂₂, U₂₄ = U₂₄, U₂₂; U₃₂, U₃₄ = U₃₄, U₃₂; U₄₂, U₄₄ = U₄₄, U₄₂
        end
        if σ₃ > σ₂
            σ₂, σ₃ = σ₃, σ₂
            U₁₂, U₁₃ = U₁₃, U₁₂; U₂₂, U₂₃ = U₂₃, U₂₂; U₃₂, U₃₃ = U₃₃, U₃₂; U₄₂, U₄₃ = U₄₃, U₄₂
        end
    end

    @inbounds σ[1] = σ₁
    @inbounds σ[2] = σ₂
    @inbounds σ[3] = σ₃
    @inbounds σ[4] = σ₄
    @inbounds U[1,1] = U₁₁; @inbounds U[1,2] = U₁₂; @inbounds U[1,3] = U₁₃; @inbounds U[1,4] = U₁₄
    @inbounds U[2,1] = U₂₁; @inbounds U[2,2] = U₂₂; @inbounds U[2,3] = U₂₃; @inbounds U[2,4] = U₂₄
    @inbounds U[3,1] = U₃₁; @inbounds U[3,2] = U₃₂; @inbounds U[3,3] = U₃₃; @inbounds U[3,4] = U₃₄
    @inbounds U[4,1] = U₄₁; @inbounds U[4,2] = U₄₂; @inbounds U[4,3] = U₄₃; @inbounds U[4,4] = U₄₄

    return U, σ
end

# ============================================================================
# Analytic eigmin for small symmetric matrices
# ============================================================================

@propagate_inbounds function eigmin1(M::AbstractMatrix)
    @boundscheck checkbounds(M, 1, 1)
    @inbounds return M[1,1]
end

@propagate_inbounds function eigmin2(M::AbstractMatrix{T}) where T
    @boundscheck checkbounds(M, 2, 2)
    @inbounds M₁₁ = M[1,1]
    @inbounds M₂₂ = M[2,2]
    @inbounds M₂₁ = M[2,1]
    return (M₁₁ + M₂₂) / 2 - hypot((M₁₁ - M₂₂) / 2, M₂₁)
end

@propagate_inbounds function eigmin3(M::AbstractMatrix{T}) where T
    @boundscheck checkbounds(M, 3, 3)
    @inbounds M₁₁ = M[1,1]
    @inbounds M₂₂ = M[2,2]
    @inbounds M₃₃ = M[3,3]
    @inbounds M₂₁ = M[2,1]
    @inbounds M₃₁ = M[3,1]
    @inbounds M₃₂ = M[3,2]

    m = max(abs(M₁₁), abs(M₂₂), abs(M₃₃), abs(M₂₁), abs(M₃₁), abs(M₃₂))

    if iszero(m)
        return zero(T)
    elseif !isfinite(m)
        return typemin(T)
    else
        B₁₁ = M₁₁ / m; B₂₂ = M₂₂ / m; B₃₃ = M₃₃ / m
        B₂₁ = M₂₁ / m; B₃₁ = M₃₁ / m; B₃₂ = M₃₂ / m

        q = (B₁₁ + B₂₂ + B₃₃) / 3

        C₁₁ = B₁₁ - q
        C₂₂ = B₂₂ - q
        C₃₃ = B₃₃ - q

        p = sqrt((C₁₁^2 + C₂₂^2 + C₃₃^2 + 2(B₂₁^2 + B₃₁^2 + B₃₂^2)) / 6)

        if p ≤ 4eps(T) * max(abs(q), one(T))
            r = 4eps(T) * abs(q)
        else
            D₁₁ = C₁₁ / p; D₂₂ = C₂₂ / p; D₃₃ = C₃₃ / p
            D₂₁ = B₂₁ / p; D₃₁ = B₃₁ / p; D₃₂ = B₃₂ / p

            r = ( D₁₁ * (D₂₂ * D₃₃ - D₃₂ * D₃₂)
                - D₂₁ * (D₂₁ * D₃₃ - D₃₂ * D₃₁)
                + D₃₁ * (D₂₁ * D₃₂ - D₂₂ * D₃₁) ) / 2

            r = clamp(r, -one(T), one(T))

            r = 2p * (-cos((acos(r) + π + π) / 3) + 16eps(T) / sqrt(max(one(T) - r^2, 16eps(T))))
        end

        return m * (q - r)
    end
end

# Smallest eigenvalue of a symmetric 4×4 matrix (lower triangle accessed),
# via the closed-form (Euler) quartic solution + one verified Newton polish.
#
# Method:
#   1. scale by max-abs entry, shift by tr/4  ->  trace-free B, char poly
#      g(λ) = λ⁴ + pλ² + qλ + r  with  p = -tr(B²)/2, q = -tr(B³)/3, r = det(B)
#   2. resolvent cubic (roots zₖ = (λ₁+λⱼ)² ≥ 0) solved in *depressed* form,
#      with P, Q formed directly from p, q, r (avoids the cascaded
#      cancellations of the a₂, a₁ = p²-4r, a₀ route)
#   3. the largest root z₁ = 2ρcosθ + s₀ comes from the trig formula, which
#      never cancels (both terms ≥ 0).  The remaining pair is then:
#        - if subordinate (z₂+z₃ < z₁/2): recovered via Vieta,
#            z₂+z₃ = -2p - z₁,   √(z₂z₃) = |q|/√z₁  (AM-GM clipped),
#            s₂+s₃ = √(z₂+z₃ + 2√(z₂z₃))            (no cancellation)
#          This avoids the trig formula's catastrophic cancellation when
#          z₂, z₃ ≪ s₀ (e.g. spectra with two double eigenvalues).
#        - otherwise: all three roots from trig as before — for clustered
#          cubic roots the trig errors are anticorrelated and cancel in the
#          sums s₁+s₂+s₃, which Vieta deflation would forfeit.
#   4. λmin = -(√z₁+√z₂+√z₃)/2 if q ≥ 0, else min(√zₖ) - (√z₁+√z₂+√z₃)/2
#   5. one Newton step on g, accepted only if it decreases |g|
#
# Accuracy (vs LAPACK eigvals, error relative to ‖M‖₂):
#   simple λmin (incl. scales 1e±300):     ~1e-16 .. 1e-14
#   multiple / near-multiple λmin:         ~1e-8 worst case, ~1e-16 typical
#     (√eps is the intrinsic floor of the characteristic-polynomial
#      representation: an O(eps) coefficient perturbation splits a double
#      root by O(√eps).  Applies to double, triple, and two-double spectra.)
#   graded spectra (λmin simple but g'(λmin) ≈ ∏gaps tiny): limited by
#     eps/|g'|, again intrinsic to the coefficient representation.
@propagate_inbounds function eigmin4(M::AbstractMatrix{T}) where T
    @boundscheck checkbounds(M, 4, 4)
    @inbounds M₁₁ = M[1,1]
    @inbounds M₂₂ = M[2,2]
    @inbounds M₃₃ = M[3,3]
    @inbounds M₄₄ = M[4,4]
    @inbounds M₂₁ = M[2,1]
    @inbounds M₃₁ = M[3,1]
    @inbounds M₄₁ = M[4,1]
    @inbounds M₃₂ = M[3,2]
    @inbounds M₄₂ = M[4,2]
    @inbounds M₄₃ = M[4,3]

    m = max(abs(M₁₁), abs(M₂₂), abs(M₃₃), abs(M₄₄),
            abs(M₂₁), abs(M₃₁), abs(M₄₁), abs(M₃₂), abs(M₄₂), abs(M₄₃))

    if iszero(m)
        return zero(T)
    elseif !isfinite(m)
        return typemin(T)
    else
        N₁₁ = M₁₁ / m; N₂₂ = M₂₂ / m; N₃₃ = M₃₃ / m; N₄₄ = M₄₄ / m
        e = M₂₁ / m; f = M₃₁ / m; g = M₄₁ / m
        h = M₃₂ / m; i = M₄₂ / m; j = M₄₃ / m

        q₀ = (N₁₁ + N₂₂ + N₃₃ + N₄₄) / 4

        # trace-free part B (diagonal a,b,c,d; entries now in [-2,2])
        a = N₁₁ - q₀; b = N₂₂ - q₀; c = N₃₃ - q₀; d = N₄₄ - q₀

        # g(λ) = λ⁴ + pλ² + qλ + r  =  det(λI - B)
        p = -(a^2 + b^2 + c^2 + d^2 + 2(e^2 + f^2 + g^2 + h^2 + i^2 + j^2)) / 2

        q = -( a^3 + b^3 + c^3 + d^3
             + 3((a + b)*e^2 + (a + c)*f^2 + (a + d)*g^2
               + (b + c)*h^2 + (b + d)*i^2 + (c + d)*j^2)
             + 6(e*f*h + e*g*i + f*g*j + h*i*j) ) / 3

        C₁ = b*(c*d - j*j) - h*(h*d - j*i) + i*(h*j - c*i)
        C₂ = e*(c*d - j*j) - h*(f*d - j*g) + i*(f*j - c*g)
        C₃ = e*(h*d - j*i) - b*(f*d - j*g) + i*(f*i - h*g)
        C₄ = e*(h*j - c*i) - b*(f*j - c*g) + h*(f*i - h*g)
        r = a*C₁ - e*C₂ + f*C₃ - g*C₄

        # resolvent cubic, depressed form built directly from p, q, r:
        #   z = y + s₀,  y³ + Py + Q = 0
        q² = q * q
        s₀ = -2p / 3                                # ≥ 0 since p ≤ 0
        P  = -(p*p / 3 + 4r)
        Q  = (8p*r - 2(p*p)*p / 9) / 3 - q²

        ρ² = max(-P / 3, zero(T))
        ρ  = sqrt(ρ²)

        if ρ ≤ 4eps(T) * max(s₀, one(T))
            # numerically triple root: z₁ = z₂ = z₃ = s₀
            s₁ = sqrt(max(s₀, zero(T)))

            if q ≥ zero(T)
                λ = -3s₁ / 2
            else
                λ = -s₁ / 2
            end
        else
            w = clamp(-Q / (2ρ*ρ²), -one(T), one(T))
            st, ct = sincos(acos(w) / 3)

            z₁ = max(2ρ*ct + s₀, zero(T))           # largest root: never cancels
            t  = max(-2p - z₁, zero(T))             # z₂ + z₃ by Vieta (trace)

            if t < z₁ / 2
                # pair subordinate to top root: Vieta reconstruction
                s₁  = sqrt(z₁)
                gm  = min(sqrt(q² / z₁), t / 2)     # √(z₂z₃), AM-GM consistent
                Σ₂₃ = sqrt(t + 2gm)                 # s₂ + s₃: no cancellation
                Δ₂₃ = sqrt(max(t - 2gm, zero(T)))   # s₂ - s₃
                s₂  = (Σ₂₃ + Δ₂₃) / 2

                if iszero(s₂)
                    s₃ = zero(T)
                else
                    s₃ = gm / s₂    # min of pair (and overall)
                end

                if q ≥ zero(T)
                    λ = -(s₁ + Σ₂₃) / 2
                else
                    λ = s₃ - (s₁ + Σ₂₃) / 2
                end
            else
                # comparable/clustered roots: symmetric trig — anticorrelated
                # errors cancel in the sums below
                u = sqrt(T(3)) / 2 * st
                v = ct / 2
                z₂ = max(2ρ*(u - v) + s₀, zero(T))
                z₃ = max(2ρ*(-u - v) + s₀, zero(T))

                s₁ = sqrt(z₁); s₂ = sqrt(z₂); s₃ = sqrt(z₃)
                S = s₁ + s₂ + s₃

                if q ≥ zero(T)
                    λ = -S / 2
                else
                    λ = min(s₁, s₂, s₃) - S / 2
                end
            end
        end

        # one verified Newton step on g(λ)
        g₀ = muladd(muladd(muladd(λ, λ, p), λ, q), λ, r)
        g₁ = muladd(muladd(4λ, λ, 2p), λ, q)

        if !iszero(g₁)
            λ′ = λ - g₀ / g₁
            gₙ = muladd(muladd(muladd(λ′, λ′, p), λ′, q), λ′, r)

            if abs(gₙ) ≤ abs(g₀)
                λ = λ′
            end
        end

        return m * (q₀ + λ)
    end
end

# ============================================================================
# Sturm sequences for eigenvalue counting
# ============================================================================

# Sturm test: are there ≥ k eigenvalues of SymTridiagonal(α[1:j], β[1:j-1]) below μ?
# Early-exits once k sign changes are found.
@propagate_inbounds function sturm_ge(α::AbstractVector{T}, β::AbstractVector{T}, j::Int, μ::T, k::Int; pivmin::T = floatmin(T)) where T
    @assert j ≥ 1
    @boundscheck checkbounds(α, j)
    @boundscheck isone(j) || checkbounds(β, j - 1)

    count = 0

    @inbounds q = α[1] - μ

    if q < 0
        count += 1
        count >= k && return true
    end

    @inbounds for i in 1:j - 1
        if abs(q) < pivmin
            q = -pivmin
        end

        q = α[i + 1] - μ - β[i]^2 / q

        if q < 0
            count += 1
            count >= k && return true
        end
    end

    return false
end

# ============================================================================
# Eigenvalue computation via Sturm bisection
# ============================================================================

# Find k-th smallest eigenvalue of SymTridiagonal(α[1:j], β[1:j-1]) via Sturm bisection.
# Returns (lo, hi) bracket for the eigenvalue.
@propagate_inbounds function eigval_tridiag(α::AbstractVector{T}, β::AbstractVector{T}, j::Int, k::Int; tol::T = T(1e-6), hi0::T = typemax(T)) where T
    @assert j ≥ 1
    @boundscheck checkbounds(α, j)
    @boundscheck isone(j) || checkbounds(β, j - 1)

    if isone(j)
        @inbounds lo = hi = α[1]
    else
        #
        # compute Gershgorin bounds
        #
        @inbounds lo = α[1] - abs(β[1])
        @inbounds hi = α[1] + abs(β[1])

        @inbounds for i in 2:j - 1
            r = abs(β[i - 1]) + abs(β[i])
            lo = min(lo, α[i] - r)
            hi = max(hi, α[i] + r)
        end

        @inbounds r = abs(β[j - 1])
        @inbounds lo = min(lo, α[j] - r)
        @inbounds hi = max(hi, α[j] + r)

        lo -= abs(lo) * tol
        hi += abs(hi) * tol
        #
        # interlacing guarantees that
        #
        #   λⱼ ≤ λᵢ
        #
        # where i = j - 1.
        #
        hi = min(hi, hi0)
        #
        # compute minimum pivot:
        #
        #   ϵ maxᵢ ‖βᵢ‖²
        #
        maxβsq = floatmin(T)

        @inbounds for i in 1:j - 1
            maxβsq = max(maxβsq, β[i]^2)
        end

        pivmin = floatmin(T) * max(one(T), maxβsq)
        #
        # compute λ using bisection
        #
        while hi - lo > tol * max(one(T), abs(lo), abs(hi))
            mid = (lo + hi) / 2

            if sturm_ge(α, β, j, mid, k; pivmin)
                hi = mid
            else
                lo = mid
            end
        end
    end
    #
    # return bracket:
    #
    #   lo ≤ λ ≤ hi
    #
    return (lo, hi)
end

# ============================================================================
# Pivoted shifted-tridiagonal solve (§3.1 of spec)
# ============================================================================

# Solve (T − λI)·y = eⱼ where T = SymTridiagonal(α[1:j], β[1:j-1]).
# Uses partial pivoting per LAPACK dlagtf/dlagts to avoid unbounded element growth.
# Writes unnormalized solution into y[1:j]; caller normalizes.
@propagate_inbounds function tridiag_solve_ej!(y::AbstractVector{T}, α::AbstractVector{T},
                                               β::AbstractVector{T}, j::Int, λ::T,
                                               work::AbstractVector{T}) where T
    @assert j ≥ 1
    @boundscheck checkbounds(y, j)
    @boundscheck checkbounds(α, j)
    @boundscheck isone(j) || checkbounds(β, j - 1)
    @boundscheck isone(j) || checkbounds(work, 5j - 3)

    @inbounds d = α[1] - λ

    if isone(j)
        @inbounds y[1] = inv(signfloor(d, floatmin(T)))
    else
        #
        # compute offsets
        #
        o1 = 0j
        o2 = 1j
        o3 = 2j - 1
        ol = 3j - 3
        op = 4j - 3

        anorm = zero(T)

        @inbounds e = β[1]
        #
        # factorize the tri-diagonal matrix
        #
        #   T - λI = PLU
        #
        @inbounds for k in 1:j - 1
            dnxt = α[k + 1] - λ

            if k < j - 1
                enxt = β[k + 1]
            else
                enxt = zero(T)
            end

            c = β[k]

            anorm = max(anorm, abs(d) + abs(e))

            if abs(d) ≥ abs(c)
                #
                # no interchange
                #
                if iszero(d)
                    m = zero(T)
                else
                    m = c / d
                end

                work[op + k] = zero(T)
                work[o1 + k] = d
                work[o2 + k] = e

                if k < j - 1
                    work[o3 + k] = zero(T)
                end

                d = dnxt - m * e
                e = enxt
            else
                #
                # interchange rows k and k + 1
                #
                m = d / c

                work[op + k] = one(T)
                work[o1 + k] = c
                work[o2 + k] = dnxt

                if k < j - 1
                    work[o3 + k] = enxt
                end

                d = -m * dnxt + e
                e = -m * enxt
            end

            work[ol + k] = m
        end

        @inbounds work[o1 + j] = d

        anorm = max(anorm, abs(d))
        #
        # forward substitution: solve
        #
        #   PLw = eⱼ
        #
        @inbounds for k in 1:j - 1
            y[k] = zero(T)
        end

        @inbounds y[j] = one(T)

        @inbounds for k in 1:j - 1
            if isone(work[op + k])
                y[k], y[k + 1] = y[k + 1], y[k]
            end

            y[k + 1] -= work[ol + k] * y[k]
        end
        #
        # backward substitution: solve
        #
        #   Ux = w
        #
        pm = max(eps(T) * anorm, floatmin(T))

        @inbounds y[j] /= signfloor(work[o1 + j], pm)
        @inbounds y[j - 1] = (y[j - 1] - work[o2 + j - 1] * y[j]) / signfloor(work[o1 + j - 1], pm)

        @inbounds for k in j - 2:-1:1
            y[k] = (y[k] - work[o2 + k] * y[k + 1] - work[o3 + k] * y[k + 2]) / signfloor(work[o1 + k], pm)
        end
    end

    return y
end

# ============================================================================
# Inverse iteration for eigenvector components
# ============================================================================

# One step of inverse iteration on (T - λI) to get |y[j]|
# T is SymTridiagonal(α[1:j], β[1:j-1]), λ is an approximate eigenvalue
# work must have length ≥ 6j-3 (5j-3 for pivoted solve + j for solution)
# Returns |y[j]| where y is the normalized solution of (T - λI)y ≈ 0
@propagate_inbounds function inviter_last(α::AbstractVector{T}, β::AbstractVector{T}, j::Int, λ::T, work::AbstractVector{T}) where T
    @assert j ≥ 1
    @boundscheck checkbounds(α, j)
    @boundscheck isone(j) || checkbounds(β, j - 1)
    @boundscheck checkbounds(work, 6j - 3)

    w = view(work,      1:5j - 3)
    y = view(work, 5j - 2:6j - 3)
    #
    # solve the tri-diagonal system
    #
    #   (T - λI)y = eⱼ
    #
    tridiag_solve_ej!(y, α, β, j, λ, w)
    #
    # compute the normalized absolute value
    #
    #   |yⱼ| / ‖y‖
    #
    nsq = zero(T)

    @inbounds for i in 1:j
        nsq += y[i]^2
    end

    if nsq < floatmin(T)
        return zero(T)
    else
        return abs(y[j]) / sqrt(nsq)
    end
end

# Compute Ritz vector q = Q * y where y is the eigenvector of T for λ
# Q is n×(j+1), α[1:j] and β[1:j-1] define T, λ is the target eigenvalue
# work must have length ≥ 6j-3, q is output (length n)
@propagate_inbounds function ritzvec!(q::AbstractVector{T}, Q::AbstractMatrix{T},
                  α::AbstractVector{T}, β::AbstractVector{T}, j::Int,
                  λ::T, work::AbstractVector{T}) where T
    @assert j ≥ 1
    @boundscheck checkbounds(α, j)
    @boundscheck isone(j) || checkbounds(β, j - 1)
    @boundscheck checkbounds(work, 6j - 3)
    @boundscheck checkbounds(Q, axes(q, 1), 1:j)

    w = view(work,      1:5j - 3)
    y = view(work, 5j - 2:6j - 3)
    #
    # solve the tri-diagonal system
    #
    #   (T - λI)y = eⱼ
    #
    tridiag_solve_ej!(y, α, β, j, λ, w)
    #
    # compute the norm
    #
    #   ‖y‖
    #
    nsq = zero(T)

    @inbounds for i in 1:j
        nsq += y[i]^2
    end
    #
    # compute the Ritz vector
    #
    #   q = Qⱼ (y / ‖y‖)
    #
    fill!(q, zero(T))

    if nsq ≥ eps(T)
        inrm = inv(sqrt(nsq))

        @inbounds for k in 1:j
            yk = y[k] * inrm

            for i in axes(Q, 1)
                q[i] += Q[i, k] * yk
            end
        end
    end

    return q
end

# ============================================================================
# Lanczos iteration for smallest eigenvalue
# ============================================================================

# Lanczos iteration to find smallest eigenvalue of symmetric matrix A.
# Q[:,1] must be initialized and normalized before calling.
# Returns (λlo, λ1, flag, jlst) where:
#   λlo: certified lower bound on smallest eigenvalue
#   λ1: estimate of smallest eigenvalue
#   flag: true if tolerance met
#   jlst: number of iterations performed
#
# nlock: number of locked (deflated) vectors in Q[:, 1:nlock]. The Krylov
# basis lives at Q[:, nlock+1:end]. DGKS sweeps orthogonalize against both
# the locked set and the growing Krylov basis.
@propagate_inbounds function eigminlanczos!(
        A::Symmetric{T},
        Q::AbstractMatrix{T},
        α::AbstractVector{T},
        β::AbstractVector{T},
        work::AbstractVector{T};
        tol::T = T(1e-3),
        nlock::Int = 0
    ) where T
    n    = size(A, 1)
    jmax = size(Q, 2) - 1 - nlock
    jmin = min(4, n)

    @boundscheck n == size(Q, 1) || error()
    @boundscheck checkbounds(α, jmax)
    @boundscheck checkbounds(β, jmax)
    @boundscheck checkbounds(work, 6jmax)

    λ1hi  = typemax(T)
    λ2hi  = typemax(T)
    λlo   = typemin(T)
    λ1    = zero(T)
    njp1  = zero(T)
    nT    = one(T)

    flag = false
    jlst = 0

    @inbounds for j in 1:jmax
        jlst = j

        Qj   = view(Q, :, nlock + j)
        Qjp1 = view(Q, :, nlock + j + 1)

        mul!(Qjp1, A, Qj)

        if !isone(j)
            Qjm1 = view(Q, :, nlock + j - 1)
            axpy!(-njp1, Qjm1, Qjp1)
        end
        #
        # two-pass DGKS re-orthogonalization (includes locked columns)
        #
        αj = zero(T)
        pjp1 = sqrt(normsq(Qjp1))

        for _ in (1, 2)
            for k in 1:nlock + j - 1
                Qk = view(Q, :, k)
                axpy!(-dot(Qk, Qjp1), Qk, Qjp1)
            end

            αj += h = dot(Qj, Qjp1)
            axpy!(-h, Qj, Qjp1)
            njp1 = sqrt(normsq(Qjp1))

            njp1 > pjp1 / 2 && break
            pjp1 = njp1
        end

        α[j] = αj
        nT = max(nT, abs(αj), njp1)

        if j ≥ jmin
            λlo, λ1, flag, λ1hi, λ2hi = ritzbound!(α, β, j, njp1, work;
                                                   tol=tol, jmin=jmin,
                                                   λhi1=λ1hi, λhi2=λ2hi)
            flag && break
        end

        njp1 < eps(T) * nT && break

        β[j] = njp1
        ldiv!(njp1, Qjp1)
    end

    return (λlo, λ1, flag, jlst)
end

# ============================================================================
# Kato-Temple certified Ritz bounds
# ============================================================================

# Compute Kato-Temple certified lower bound for smallest eigenvalue
# α[1:j], β[1:j-1] define the Lanczos tridiagonal
# βnxt is the norm of the next Lanczos vector (residual norm)
# work must have length ≥ 6j
# jmin: minimum iteration before certification (default 4, relaxed to min(4,n) by caller)
# λhi1, λhi2: warm brackets from previous iteration (§3.3), typemax if cold start
# Returns (λlo, λmd1, flag, λhi1, λhi2) — last two for warm bracket threading
@propagate_inbounds function ritzbound!(α::AbstractVector{T}, β::AbstractVector{T}, j::Int, βnxt::T,
                    work::AbstractVector{T}; tol::T = T(1e-3), jmin::Int = 4,
                    λhi1::T = typemax(T), λhi2::T = typemax(T)) where T
    @assert j ≥ 1
    @boundscheck checkbounds(α, j)
    @boundscheck isone(j) || checkbounds(β, j - 1)
    @boundscheck isone(j) || checkbounds(work, 6j - 3)

    if isone(j)
        @inbounds λlo  = α[1] - βnxt
        @inbounds λmd1 = α[1]

        flag = jmin ≤ 1 && βnxt ≤ tol * max(one(T), abs(α[1]))

        λhi1 = typemax(T)
        λhi2 = typemax(T)
    else
        #
        # find the two smallest Ritz values
        #
        @inbounds λlo1, λhi1 = eigval_tridiag(α, β, j, 1; tol=tol/10, hi0=λhi1)
        @inbounds λlo2, λhi2 = eigval_tridiag(α, β, j, 2; tol=tol/10, hi0=λhi2)

        λ1 = λlo1
        λ2 = λlo2

        λmd1 = (λlo1 + λhi1) / 2
        λmd2 = (λlo2 + λhi2) / 2
        #
        # get last eigenvector components via inverse iteration
        #
        @inbounds y1 = inviter_last(α, β, j, λmd1, work)
        @inbounds y2 = inviter_last(α, β, j, λmd2, work)
        #
        # compute residual norms
        #
        r1 = βnxt * y1
        r2 = βnxt * y2
        #
        # compute Kato-Temple bound
        #
        gap = (λ2 - r2) - λ1

        if gap > zero(T)
            bnd = min(r1, r1^2 / gap)
        else
            bnd = r1
        end

        λlo = λ1 - bnd
        flag = j ≥ jmin && bnd ≤ tol * max(one(T), abs(λ1))
    end

    return (λlo, λmd1, flag, λhi1, λhi2)
end

# ============================================================================
# LAPACK eigmin with eigenvector
# ============================================================================

function eigminlapack!(
        A::AbstractMatrix{T},
        z::AbstractVector{T},
        W::AbstractVector{T},
        work::Vector{T},
        iwork::Vector{BlasInt}
    ) where {T <: BlasFloat}
    n = size(A, 1)

    @assert size(A, 2) == n
    @assert length(z) >= n
    @assert length(W) >= n

    require_one_based_indexing(A, z, W, work, iwork)
    chkstride1(A)
    m = Ref{BlasInt}()
    info = Ref{BlasInt}()
    lwork = BlasInt(-1)
    liwork = BlasInt(-1)
    isuppz = Vector{BlasInt}(undef, 2)

    for i in 1:2
        ccall((@blasfunc(dsyevr_), libblastrampoline), Cvoid,
              (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
               Ptr{T}, Ref{BlasInt}, Ref{T}, Ref{T},
               Ref{BlasInt}, Ref{BlasInt}, Ref{T}, Ptr{BlasInt},
               Ptr{T}, Ptr{T}, Ref{BlasInt}, Ptr{BlasInt},
               Ptr{T}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
               Ref{BlasInt}, Clong, Clong, Clong),
              'V', 'I', 'L', n,
              A, max(1, stride(A, 2)), zero(T), zero(T),
              1, 1, -one(T), m,
              W, z, n, isuppz,
              work, lwork, iwork, liwork,
              info, 1, 1, 1)

        chklapackerror(info[])

        if i == 1
            lwork = round(BlasInt, nextfloat(real(work[1])))
            liwork = iwork[1]

            if lwork > length(work)
                resize!(work, lwork)
            end

            if liwork > length(iwork)
                resize!(iwork, liwork)
            end
        end
    end

    return W[1]
end

# ============================================================================
# Negative-curvature seed from failed Cholesky (Gill-Murray)
# ============================================================================

# Extract z such that zᵀBz ≤ 0 from a partial Cholesky factor that failed at
# pivot k. z has guaranteed overlap with the eigenspace below the threshold.
@propagate_inbounds function negcurv!(z::AbstractVector{T}, F::AbstractMatrix{T}, k::Int) where T
    n = size(F, 1)
    @inbounds for i in 1:k - 1
        z[i] = F[k, i]
    end
    ldiv!(LowerTriangular(view(F, 1:k - 1, 1:k - 1))', view(z, 1:k - 1))
    @inbounds for i in 1:k - 1
        z[i] = -z[i]
    end
    @inbounds z[k] = one(T)
    @inbounds for i in k + 1:n
        z[i] = zero(T)
    end
    return z
end

# ============================================================================
# eigmin!: Lanczos with Cholesky gate and deflated recovery
# ============================================================================

const MAXDEFLATE = 2

# Compute τ = 1/max(1, -λmin(M)) with Cholesky-gated Lanczos and deflated
# recovery. Q0 is the warm start (already fuzzed by caller); updated in place
# with the new warm start on success.
#
# Arguments:
#   M: the matrix (will be read but not modified)
#   C: Cholesky workspace (n×n, will be overwritten)
#   Q: Lanczos basis workspace (n × (jmax+1))
#   α, β: tridiagonal coefficients
#   ritzwork: workspace for ritzbound!/ritzvec!
#   Q0: warm start vector (updated in place)
#   work, iwork: LAPACK workspace for fallback
#   tol: convergence tolerance
#
function eigmin!(
        M::Symmetric{T},
        C::AbstractMatrix{T},
        Q::AbstractMatrix{T},
        α::AbstractVector{T},
        β::AbstractVector{T},
        ritzwork::AbstractVector{T},
        Q0::AbstractVector{T},
        work::Vector{T},
        iwork::Vector{BlasInt};
        tol::T = T(1e-3)
    ) where T
    n = size(M, 1)
    m = size(Q, 2)
    P = parent(M)

    @inbounds for nlp1 in 1:min(MAXDEFLATE + 1, m - 4)
        Qk = view(Q, :, nlp1:m  )
        Qr = view(Q, :, nlp1    )
        Qs = view(Q, :, nlp1 + 1)

        nlock = nlp1 - 1

        λlo, λ1, flag, jlst = eigminlanczos!(M, Q, α, β, ritzwork; tol=tol, nlock=nlock)
        flag || break

        τ = inv(max(one(T), -λlo))

        copyto!(C, I)
        axpy!(τ, P, C)

        F = cholesky!(Symmetric(C, :L); check=false)

        if issuccess(F)
            if isone(jlst)
                copyto!(Q0, Qr)
            else
                ritzvec!(Q0, Qk, α, β, jlst, λ1, ritzwork)
            end

            return τ
        end

        if nlp1 ≤ MAXDEFLATE
            #
            # recovery: (a) extract Ritz vector, (b) extract neg-curv seed, (c) lock
            #
            ritzvec!(Q0, Qk, α, β, jlst, λ1, ritzwork)
            negcurv!(Qs, C, F.info)
            copyto!(Qr, Q0)
            #
            # project seed against locked columns and normalize
            #

            for l in 1:nlp1
                Ql = view(Q, :, l)
                axpy!(-dot(Ql, Qs), Ql, Qs)
            end

            nQs = normsq(Qs)
            nQs < eps(T)^2 && break

            ldiv!(sqrt(nQs), Qs)
        end
    end
    #
    # LAPACK fallback (Q is dead, reuse column 1 as eigenvalue workspace)
    #
    copyto!(C, P)
    W = view(Q, :, 1)
    λ = eigminlapack!(C, Q0, W, work, iwork)
    return inv(max(one(T), -λ))
end
