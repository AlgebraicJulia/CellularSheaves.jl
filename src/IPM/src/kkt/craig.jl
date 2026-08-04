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
# the residual g − B x = −ng·ξ·u — which CRAIG forms internally (‖g − B x‖ = ng·|ξ| is its stopping norm),
# so returning it is one length-m scale, not a matvec. The caller uses it to recover y without re-forming
# B x. atol is the only tolerance: an absolute bound on the true primal residual ‖g − B x‖₂ (the left
# preconditioner is I, so nr = ng·|ξ| is exactly that residual's norm). Terminals — CRAIG_SOLVED: nr ≤ atol
# or a zero-residual start; CRAIG_STAGNATED: nr at the machine floor (1 + bkwerr ≤ 1); CRAIG_ITMAX: the
# dimension cap m + n; CRAIG_NUMERICAL_FAILURE: α-breakdown (nv = 0, v-side Krylov space collapsed). CRAIG
# assumes a consistent system (full-row-rank B); it does not detect inconsistency (that shows as a stalled
# residual, not α = 0), so an inconsistent RHS yields a garbage x under a SOLVED/STAGNATED label, not a
# failure — see the §10.4 finding.

# CRAIG_CONTINUE is the in-progress sentinel — the iteration runs `while status == CRAIG_CONTINUE` and
# writes a terminal status the moment a stop fires.
@enum CraigStatus begin
    CRAIG_CONTINUE
    CRAIG_SOLVED             # ‖r‖ ≤ atol, or a zero-residual start
    CRAIG_NUMERICAL_FAILURE  # α-breakdown (nv = 0): the v-side Krylov space collapsed
    CRAIG_STAGNATED          # ‖r‖ at the machine floor for this precision
    CRAIG_ITMAX
end

struct CraigWorkspace{T}
    Av::FVector{T}     # Golub–Kahan v-side, before preconditioning (n)
    v::FVector{T}      # (F Fᵀ)⁻¹ N v (n)
    w::FVector{T}      # y-recurrence direction (m)
end

function CraigWorkspace{T}(m::Integer, n::Integer) where {T}
    return CraigWorkspace{T}(FVector{T}(undef, n), FVector{T}(undef, n), FVector{T}(undef, m))
end

function craig!(wrk::CraigWorkspace{T}, B, F, divwrk, x, y, g; kwargs...) where {T}
    return craig!(wrk.Av, wrk.v, wrk.w, B, F, divwrk, x, y, g; kwargs...)
end

function craig!(
        Av::AbstractVector{T},
        v::AbstractVector{T},
        w::AbstractVector{T},
        B::AbstractMatrix{T},
        F::AbstractMatrix{T},
        divwrk,
        x::AbstractVector{T},
        y::AbstractVector{T},
        g::AbstractVector{T};
        atol::T  = √eps(T),
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
    ng0 = norm(g)
    nr = ng0

    if iszero(ng0)
        return 0, CRAIG_SOLVED
    end

    ng0² = ng0^2
    ng = ng0
    ξ = -one(T)
    ρprv = one(T)

    rdiv!(g, ng0)
    fill!(Av, zero(T))
    fill!(w, zero(T))

    nB² = zero(T); nB = zero(T)
    nx² = zero(T); nx = zero(T)

    iter = 0

    if iszero(itmax)
        itmax = m + n
    end

    bkwerr = one(T)                   # backward error ‖r‖ / √(‖g‖² + ‖B‖²‖x‖²)

    # a zero-residual / already-within-tolerance start converges before the first iteration
    if nr ≤ atol
        status = CRAIG_SOLVED
    else
        status = CRAIG_CONTINUE
    end

    while status == CRAIG_CONTINUE
        #
        # αₖ₊₁ N vₖ₊₁ = Bᵀ uₖ − βₖ N vₖ,   then precondition  v ← (F Fᵀ)⁻¹ N v
        #
        mul!(Av, Bᵀ, g, one(T), -ng)  # Av ← Bᵀ g − ng·Av
        copyto!(v, Av)
        ldiv!(divwrk, F, v)
        ldiv!(divwrk, F', v)
        nv² = dot(v, Av); nv = sqrt(nv²)  # (F Fᵀ)⁻¹-elliptic norm of the v-side

        if iszero(nv)
            status = CRAIG_NUMERICAL_FAILURE
        else
            rdiv!(v, nv)
            rdiv!(Av, nv)
            #
            # advance the solution recurrences for x and y
            #
            ρ = nv
            ξ *= -ng / ρ

            axpy!(ξ, v, x)                     # x ← x + ξ v
            axpby!(one(T), g, -ng / ρprv, w)   # w ← u − (ng/ρprv) w
            axpy!(ξ / ρ, w, y)                 # y ← y + (ξ/ρ) w
            #
            # βₖ₊₁ M uₖ₊₁ = B vₖ − αₖ M uₖ
            #
            mul!(g, B, v, one(T), -nv)         # u ← B v − nv·u  (u ≡ g)

            ng² = dot(g, g); ng = sqrt(ng²)

            if !iszero(ng)
                rdiv!(g, ng)
            end

            nB² += nv² + ng²; nB = sqrt(nB²)
            nx² += ξ * ξ;     nx = sqrt(nx²)

            nr = ng * abs(ξ)                   # r = −ng ξ u

            iter += 1

            bkwerr = nr / sqrt(ng0² + nB² * nx²)
            ρprv = ρ

            if nr ≤ atol
                status = CRAIG_SOLVED
            elseif one(T) + bkwerr ≤ one(T)
                status = CRAIG_STAGNATED
            elseif iter ≥ itmax
                status = CRAIG_ITMAX
            end
        end
    end

    #
    # g holds the final u; the residual g − B x = −ng·ξ·u, so scaling in place turns g into the residual
    # it must return. Exact at 0 iterations too (ng = ng0, ξ = −1, g = g₀/ng0 ⟹ −ng·ξ·g = g₀); the g = 0
    # early return already left g = 0.
    #
    rmul!(g, -ng * ξ)

    return iter, status
end

# ─── Split-preconditioned variant (experimental; A/B against craig!) ────────────────────────────────
# Instead of applying (F Fᵀ)⁻¹ as an operator (two back-to-back solves feeding the κ(F)²-conditioned
# elliptic norm nv² = Avᵀ(F Fᵀ)⁻¹Av), bidiagonalize the explicitly-factored M = B F⁻ᵀ directly (standard
# Golub–Kahan, split preconditioning). The v-side coefficient is then nv = ‖vs‖₂ — a plain 2-norm,
# κ(F)-conditioned. Change of variable z = Fᵀδx turns min ‖δx‖_{FFᵀ} s.t. Bδx = g into min ‖z‖₂ s.t.
# Mz = g, with δx = F⁻ᵀz. The direction v = F⁻ᵀvs equals craig!'s preconditioned direction, so the x/y
# and residual recurrences are identical — same method, same iterates in exact arithmetic, better
# conditioning when F Fᵀ = βA + BᵀB is near-singular (the ρ-ladder regime).
function craig_split!(wrk::CraigWorkspace{T}, B, F, divwrk, x, y, g; kwargs...) where {T}
    return craig_split!(wrk.Av, wrk.v, wrk.w, B, F, divwrk, x, y, g; kwargs...)
end

function craig_split!(
        vs::AbstractVector{T},   # standard Golub–Kahan v-side of M = B F⁻ᵀ (n)
        v::AbstractVector{T},    # v = F⁻ᵀ vs: primal direction (x += ξ v) + forward operand; also Bᵀg temp (n)
        w::AbstractVector{T},    # y-recurrence direction (m)
        B::AbstractMatrix{T},
        F::AbstractMatrix{T},
        divwrk,
        x::AbstractVector{T},
        y::AbstractVector{T},
        g::AbstractVector{T};
        atol::T  = √eps(T),
        itmax::Int = 0,
    ) where {T}
    m, n = size(B)
    Bᵀ = B'

    fill!(y, zero(T))
    ng0 = norm(g)
    nr = ng0

    if iszero(ng0)
        return 0, CRAIG_SOLVED
    end

    ng0² = ng0^2
    ng = ng0
    ξ = -one(T)
    ρprv = one(T)

    rdiv!(g, ng0)
    fill!(vs, zero(T))
    fill!(w, zero(T))

    nB² = zero(T); nB = zero(T)
    nx² = zero(T); nx = zero(T)

    iter = 0

    if iszero(itmax)
        itmax = m + n
    end

    bkwerr = one(T)

    if nr ≤ atol
        status = CRAIG_SOLVED
    else
        status = CRAIG_CONTINUE
    end

    while status == CRAIG_CONTINUE
        #
        # down-leg:  Mᵀ uₖ = F⁻¹ Bᵀ uₖ,   then  αₖ vsₖ = Mᵀ uₖ − βₖ vsₖ₋₁
        #
        mul!(v, Bᵀ, g)                # v = Bᵀ g   (temp)
        ldiv!(divwrk, F, v)           # v = F⁻¹ Bᵀ g = Mᵀ u
        axpby!(one(T), v, -ng, vs)    # vs = Mᵀ u − ng·vs
        nv = norm(vs)                 # α = ‖vs‖₂   (κ(F), not κ(F)²)

        if iszero(nv)
            status = CRAIG_NUMERICAL_FAILURE
        else
            rdiv!(vs, nv)
            nv² = nv * nv
            #
            # advance the solution recurrences for x and y  (v = F⁻ᵀ vs is the primal direction)
            #
            ρ = nv
            ξ *= -ng / ρ

            copyto!(v, vs)
            ldiv!(divwrk, F', v)          # v = F⁻ᵀ vs  (= craig!'s preconditioned direction)
            axpy!(ξ, v, x)                # x ← x + ξ v
            axpby!(one(T), g, -ng / ρprv, w)
            axpy!(ξ / ρ, w, y)
            #
            # up-leg:  βₖ₊₁ uₖ₊₁ = M vsₖ − αₖ uₖ = B v − nv·u
            #
            mul!(g, B, v, one(T), -nv)    # u ← B v − nv·u  (u ≡ g)

            ng² = dot(g, g); ng = sqrt(ng²)

            if !iszero(ng)
                rdiv!(g, ng)
            end

            nB² += nv² + ng²; nB = sqrt(nB²)
            nx² += ξ * ξ;     nx = sqrt(nx²)

            nr = ng * abs(ξ)

            iter += 1

            bkwerr = nr / sqrt(ng0² + nB² * nx²)
            ρprv = ρ

            if nr ≤ atol
                status = CRAIG_SOLVED
            elseif one(T) + bkwerr ≤ one(T)
                status = CRAIG_STAGNATED
            elseif iter ≥ itmax
                status = CRAIG_ITMAX
            end
        end
    end

    rmul!(g, -ng * ξ)

    return iter, status
end
