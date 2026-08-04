# CRAIG — Golub–Kahan bidiagonalization for the minimum-norm least-squares problem
#
#   minimize ‖x‖  subject to  B x = g
#
# Vendored from Krylov.jl and specialized to our single use: real T, no left preconditioner (M = I), the
# right preconditioner (F Fᵀ)⁻¹ applied inline by two triangular solves against the augmented-matrix
# Cholesky factor F (so we pass F + its DivisionWorkspace rather than a generic operator), λ = 0, and no
# warm start / timing / logging / callbacks. We own it so a later pass can harvest the bidiagonal
# coefficients αₖ, βₖ (the pricing-model Ritz / spectral readings) without an external internal API.
#
# The Krylov @k* BLAS wrappers become their stdlib equivalents (kaxpy!→axpy!, kdot→dot, kscal!→rmul!,
# kdiv!→rdiv!, kmul!→mul!, kfill!→fill!, kcopy!→copyto!, knorm→norm). norm and rdiv! are not bit-identical
# to Krylov's raw nrm2 / multiply-by-reciprocal — they perturb the recurrence in the last bits — so this
# port is validated by convergence and CRAIG iteration count, not by diff.
#
# The workspace form deconstructs its scratch and forwards to the array form; the solution x (n) and y (m)
# are passed in, not owned (as in solveuzw!): x is accumulated onto (craig adds Σ ξ v — pass the base
# solution to get base + δx), y is overwritten with the min-norm δy. g (m) is in/out: the RHS on entry, on exit
# the residual g − B x = −β ξ u — which CRAIG forms internally (‖g − B x‖ = β·|ξ| is its stopping norm),
# so returning it is one length-m scale, not a matvec. The caller uses it to recover y without re-forming
# B x. btol is the backward-error stop (bkwerr ≤ btol); we pass 0 because the pricing model uses an
# absolute residual tolerance (atol). The knob costs one comparison per iteration (bkwerr is computed
# anyway for the machine-precision stop), so it is kept.

# CRAIG_CONTINUE is the in-progress sentinel — the iteration runs `while status == CRAIG_CONTINUE` and
# writes a terminal status the moment a stop fires. The terminals: CRAIG_SOLVED covers both a met
# residual/backward-error tolerance and the zero-residual start; the rest are the diagnostic exits
# (inconsistent right-hand side, operator too ill-conditioned, dimension-based iteration cap hit).
@enum CraigStatus begin
    CRAIG_CONTINUE
    CRAIG_SOLVED
    CRAIG_INCONSISTENT
    CRAIG_ILL_CONDITIONED
    CRAIG_ITMAX
end

struct CraigWorkspace{T}
    Nv::FVector{T}     # Golub–Kahan v-side, before preconditioning (n)
    v::FVector{T}      # (F Fᵀ)⁻¹ N v (n)
    w::FVector{T}      # y-recurrence direction (m)
end

function CraigWorkspace{T}(m::Integer, n::Integer) where {T}
    return CraigWorkspace{T}(FVector{T}(undef, n), FVector{T}(undef, n), FVector{T}(undef, m))
end

function craig!(wrk::CraigWorkspace{T}, B, F, divwrk, x, y, g; kwargs...) where {T}
    return craig!(wrk.Nv, wrk.v, wrk.w, B, F, divwrk, x, y, g; kwargs...)
end

function craig!(
        Nv::AbstractVector{T},
        v::AbstractVector{T},
        w::AbstractVector{T},
        B::AbstractMatrix{T},
        F::AbstractMatrix{T},
        divwrk,
        x::AbstractVector{T},
        y::AbstractVector{T},
        g::AbstractVector{T};
        atol::T = √eps(T),
        rtol::T = √eps(T),
        btol::T = √eps(T),
        ctol::T = √eps(T),
        itmax::Int = 0,
    ) where {T}
    m, n = size(B)
    Bᵀ = B'

    # x is accumulated onto, not zeroed: craig adds its correction Σ ξ v, so the caller passes the base
    # solution and receives base + δx. y is a fresh output (the min-norm δy).
    fill!(y, zero(T))
    #
    # g carries the Golub–Kahan u-side in place (u ≡ M u ≡ g, M = I): the RHS on entry, the u-recurrence
    # through the loop, the residual on exit. g ← g / ‖g‖; x = 0 is a zero-residual solution when g = 0.
    #
    β₁ = norm(g)
    nr = β₁

    if iszero(β₁)
        return 0, CRAIG_SOLVED
    end

    β₁² = β₁^2
    β = β₁
    θ = β₁
    ξ = -one(T)
    ρprv = one(T)

    rdiv!(g, β₁)
    fill!(Nv, zero(T))
    fill!(w, zero(T))

    nB² = zero(T); nB = zero(T)
    nD² = zero(T); condB = zero(T)
    nx² = zero(T); nx = zero(T)

    iter = 0

    if iszero(itmax)
        itmax = m + n
    end

    εc = atol + rtol * nr             # consistent-system residual tolerance
    bkwerr = one(T)                   # backward error ‖r‖ / √(‖b‖² + ‖B‖²‖x‖²)

    # a zero-residual / already-within-tolerance start converges before the first iteration
    if (one(T) + bkwerr ≤ one(T)) || (bkwerr ≤ btol) || (nr ≤ εc) || (nr ≤ btol + atol * nB * nx / β₁)
        status = CRAIG_SOLVED
    else
        status = CRAIG_CONTINUE
    end

    while status == CRAIG_CONTINUE
        #
        # αₖ₊₁ N vₖ₊₁ = Bᵀ uₖ − βₖ N vₖ,   then precondition  v ← (F Fᵀ)⁻¹ N v
        #
        mul!(Nv, Bᵀ, g, one(T), -β)   # Nv ← Bᵀ g − β Nv
        copyto!(v, Nv)
        ldiv!(divwrk, F, v)
        ldiv!(divwrk, F', v)
        α² = dot(v, Nv)               # (F Fᵀ)⁻¹-elliptic norm, squared
        α = sqrt(α²)

        if iszero(α)
            status = CRAIG_INCONSISTENT
        else
            rdiv!(v, α)
            rdiv!(Nv, α)
            nB² += α²
            #
            # advance the solution recurrences for x and y
            #
            ρ = α
            ξ *= -θ / ρ
            axpy!(ξ, v, x)                     # x ← x + ξ v
            axpby!(one(T), g, -θ / ρprv, w)    # w ← u − (θ/ρprv) w
            axpy!(ξ / ρ, w, y)                 # y ← y + (ξ/ρ) w
            nD² += norm(w)
            #
            # βₖ₊₁ M uₖ₊₁ = B vₖ − αₖ M uₖ
            #
            mul!(g, B, v, one(T), -α)          # u ← B v − α u  (u ≡ g)
            β² = dot(g, g)
            β = sqrt(β²)

            if !iszero(β)
                rdiv!(g, β)
            end

            θ = β

            nB² += β²
            nB = sqrt(nB²)
            condB = nB * sqrt(nD²)
            nx² += ξ * ξ
            nx = sqrt(nx²)
            nr = β * abs(ξ)                    # r = −β ξ u
            iter += 1

            bkwerr = nr / sqrt(β₁² + nB² * nx²)
            ρprv = ρ

            if (one(T) + inv(condB) ≤ one(T)) || (inv(condB) ≤ ctol)
                status = CRAIG_ILL_CONDITIONED
            elseif (one(T) + bkwerr ≤ one(T)) || (bkwerr ≤ btol) || (nr ≤ εc) || (nr ≤ btol + atol * nB * nx / β₁)
                status = CRAIG_SOLVED
            elseif iter ≥ itmax
                status = CRAIG_ITMAX
            end
        end
    end

    #
    # g holds the final u; the residual g − B x = −β ξ u, so scaling in place turns g into the residual
    # it must return. Exact at 0 iterations too (β = β₁, ξ = −1, g = g₀/β₁ ⟹ −β ξ g = g₀); the g = 0
    # early return already left g = 0.
    #
    rmul!(g, -β * ξ)

    return iter, status
end
