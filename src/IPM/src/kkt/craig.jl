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
# are passed in, not owned (as in solveuzw!). g (m) is in/out: the right-hand side on entry, and on exit
# the residual g − B x = −β ξ u — which CRAIG forms internally (‖g − B x‖ = β·|ξ| is its stopping norm),
# so returning it is one length-m scale, not a matvec. The caller uses it to recover y without re-forming
# B x. btol is the backward-error stop (bkwerr ≤ btol); we pass 0 because the pricing model uses an
# absolute residual tolerance (atol). The knob costs one comparison per iteration (bkwerr is computed
# anyway for the machine-precision stop), so it is kept.

# Why the CRAIG iteration halted. CRAIG_SOLVED covers both a met residual/backward-error tolerance and
# the zero-residual start; the rest are the diagnostic exits (inconsistent right-hand side, operator too
# ill-conditioned, dimension-based iteration cap hit).
@enum CraigStatus begin
    CRAIG_SOLVED
    CRAIG_INCONSISTENT
    CRAIG_ILL_CONDITIONED
    CRAIG_ITMAX
end

struct CraigWorkspace{T}
    Nv::FVector{T}     # Golub–Kahan v-side, before preconditioning (n)
    Btu::FVector{T}    # Bᵀ u (n)  [u ≡ g, the in/out buffer]
    v::FVector{T}      # (F Fᵀ)⁻¹ N v (n)
    w::FVector{T}      # y-recurrence direction (m)
    Bv::FVector{T}     # B v (m)
end

function CraigWorkspace{T}(m::Integer, n::Integer) where {T}
    return CraigWorkspace{T}(FVector{T}(undef, n), FVector{T}(undef, n), FVector{T}(undef, n),
                             FVector{T}(undef, m), FVector{T}(undef, m))
end

function craig!(wrk::CraigWorkspace{T}, B, F, divwrk, x, y, g; kwargs...) where {T}
    return craig!(wrk.Nv, wrk.Btu, wrk.v, wrk.w, wrk.Bv, B, F, divwrk, x, y, g; kwargs...)
end

function craig!(
        Nv::AbstractVector{T},
        Btu::AbstractVector{T},
        v::AbstractVector{T},
        w::AbstractVector{T},
        Bv::AbstractVector{T},
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

    fill!(x, zero(T))
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

    solved = (one(T) + bkwerr ≤ one(T)) | (bkwerr ≤ btol) | (nr ≤ εc) |
             (nr ≤ btol + atol * nB * nx / β₁)
    illcond = false
    inconsistent = false

    while !(solved || inconsistent || illcond || iter ≥ itmax)
        #
        # αₖ₊₁ N vₖ₊₁ = Bᵀ uₖ − βₖ N vₖ,   then precondition  v ← (F Fᵀ)⁻¹ N v
        #
        mul!(Btu, Bᵀ, g)
        axpby!(one(T), Btu, -β, Nv)
        copyto!(v, Nv)
        ldiv!(divwrk, F, v)
        ldiv!(divwrk, F', v)
        α = sqrt(dot(v, Nv))          # (F Fᵀ)⁻¹-elliptic norm
        if iszero(α)
            inconsistent = true
        else
            rdiv!(v, α)
            rdiv!(Nv, α)
            nB² += α * α
            #
            # advance the solution recurrences for x and y
            #
            ρ = α
            ξ = -θ / ρ * ξ
            axpy!(ξ, v, x)                     # x ← x + ξ v
            axpby!(one(T), g, -θ / ρprv, w)    # w ← u − (θ/ρprv) w
            axpy!(ξ / ρ, w, y)                 # y ← y + (ξ/ρ) w
            nD² += norm(w)
            #
            # βₖ₊₁ M uₖ₊₁ = B vₖ − αₖ M uₖ
            #
            mul!(Bv, B, v)
            axpby!(one(T), Bv, -α, g)          # u ← B v − α u
            β = norm(g)
            !iszero(β) && rdiv!(g, β)
            θ = β

            nB² += β * β
            nB = sqrt(nB²)
            condB = nB * sqrt(nD²)
            nx² += ξ * ξ
            nx = sqrt(nx²)
            nr = β * abs(ξ)                    # r = −β ξ u
            iter += 1

            bkwerr = nr / sqrt(β₁² + nB² * nx²)
            ρprv = ρ

            solved = (one(T) + bkwerr ≤ one(T)) | (bkwerr ≤ btol) | (nr ≤ εc) |
                     (nr ≤ btol + atol * nB * nx / β₁)
            illcond = (one(T) + inv(condB) ≤ one(T)) | (inv(condB) ≤ ctol)
        end
    end

    #
    # g holds the final u; the residual g − B x = −β ξ u, so scaling in place turns g into the residual
    # it must return. Exact at 0 iterations too (β = β₁, ξ = −1, g = g₀/β₁ ⟹ −β ξ g = g₀); the g = 0
    # early return already left g = 0.
    #
    rmul!(g, -β * ξ)

    if inconsistent
        status = CRAIG_INCONSISTENT
    elseif illcond
        status = CRAIG_ILL_CONDITIONED
    elseif solved
        status = CRAIG_SOLVED
    else
        status = CRAIG_ITMAX
    end
    return iter, status
end
