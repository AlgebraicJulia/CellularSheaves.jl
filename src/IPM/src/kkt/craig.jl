#
# CRAIG — minimum-norm least-squares by
# Golub–Kahan bidiagonalization
#
#   min ‖δx‖_{FFᵀ}  s.t.  B δx = g,
#   FFᵀ = β A + Bᵀ B
#
# Split-preconditioned form: bidiagonalize
# M = B F⁻ᵀ directly (the v-side coefficient
# is then a plain 2-norm nu = ‖u‖₂) rather
# than applying (FFᵀ)⁻¹ as an operator. Same
# iterates as the operator form up to
# floating-point rearrangement — A/B-
# validated, no conditioning advantage, kept
# for taste. Vendored from Krylov.jl (λ = 0,
# left preconditioner M = I, no warm start /
# timing / logging / callbacks).
#
# x accumulates δx onto the base solution; y
# is overwritten with the min-norm δy; g is
# in/out — the RHS on entry, the residual
# g − B x on exit (one length-m scale, so the
# caller recovers y without a matvec). atol
# bounds ‖g − B x‖₂; itmax = 0 skips the
# iteration entirely (y = 0, g untouched).
# Terminates SOLVED (‖r‖ ≤ atol), STAGNATED
# (machine floor), ITMAX (dimension cap), or
# NUMERICAL_FAILURE (breakdown, nu = 0).
# CRAIG assumes a consistent system (full-
# row-rank B); it does not detect
# inconsistency, so an inconsistent RHS
# yields garbage x under a SOLVED/STAGNATED
# label — see §10.4.
#

# CRAIG_CONTINUE is the in-progress sentinel — the loop runs `while status == CRAIG_CONTINUE` and writes a
# terminal status the moment a stop fires.
@enum CraigStatus begin
    CRAIG_CONTINUE
    CRAIG_SOLVED             # ‖r‖ ≤ atol, or a zero-residual start
    CRAIG_NUMERICAL_FAILURE  # breakdown (nu = 0): the v-side recurrence collapsed
    CRAIG_STAGNATED          # ‖r‖ at the machine floor for this precision
    CRAIG_ITMAX
end

struct CraigWorkspace{T}
    u::FVector{T}      # Golub–Kahan v-side of M = B F⁻ᵀ (n)
    v::FVector{T}      # v = F⁻ᵀ u: primal direction + forward operand; also the Bᵀg temp (n)
    w::FVector{T}      # y-recurrence direction (m)
end

function CraigWorkspace{T}(m::Integer, n::Integer) where {T}
    return CraigWorkspace{T}(FVector{T}(undef, n), FVector{T}(undef, n), FVector{T}(undef, m))
end

function craig!(wrk::CraigWorkspace{T}, B, F, divwrk, x, y, g; kwargs...) where {T}
    return craig!(wrk.u, wrk.v, wrk.w, B, F, divwrk, x, y, g; kwargs...)
end

function craig!(
        u::AbstractVector{T},
        v::AbstractVector{T},
        w::AbstractVector{T},
        B::AbstractMatrix{T},
        F::AbstractMatrix{T},
        divwrk,
        x::AbstractVector{T},
        y::AbstractVector{T},
        g::AbstractVector{T};
        atol::T  = √eps(T),
        itmax::Int = sum(size(B)),
    ) where {T}
    @assert length(u) == length(v) == length(x) == size(B, 2)
    @assert length(w) == length(y) == length(g) == size(B, 1)
    @assert size(F, 1) == size(B, 2)

    fill!(y, zero(T))
    #
    # entry residual norm (the residual at x = base is g itself):
    #
    #   ng0 = ‖ g ‖
    #
    ng0² = dot(g, g); ng0 = sqrt(ng0²)
    iter = 0
    status = CRAIG_SOLVED
    #
    # itmax = 0 skips the iteration entirely — the base solution x is returned with y = 0 and g (the
    # residual) untouched. Used for a bare backsolve. ng0 ≤ atol is likewise already converged.
    #
    if itmax > 0 && ng0 > atol
        ng = ng0
        ξ = -one(T)
        pnu = one(T)
        #
        # normalize the Golub–Kahan u-side (u ≡ g); ng0 > atol > 0 here, so this is safe
        #
        #   g ← g / ng0
        #
        rdiv!(g, ng0)
        fill!(u, zero(T))
        fill!(w, zero(T))

        nB² = zero(T); nB = zero(T)
        nx² = zero(T); nx = zero(T)

        η = one(T)
        status = CRAIG_CONTINUE

        while status == CRAIG_CONTINUE
            #
            # v-side bidiagonalization step (M = B F⁻ᵀ):
            #
            #   v  ← F⁻¹ Bᵀ g
            #   u  ← v - ng u,   nu ← ‖ u ‖
            #
            mul!(v, B', g)
            ldiv!(divwrk, F, v)
            axpby!(one(T), v, -ng, u)
            nu² = dot(u, u); nu = sqrt(nu²)

            if iszero(nu)
                status = CRAIG_NUMERICAL_FAILURE
            else
                #
                # advance the solution recurrences (v = F⁻ᵀ u is the primal direction):
                #
                #   u  ← u / nu
                #   ξ  ← -ng/nu ξ
                #   x  ← x + ξ v
                #   w  ← g - ng/pnu w,   y ← y + ξ/nu w
                #
                rdiv!(u, nu)
                ξ *= -ng / nu

                copyto!(v, u)
                ldiv!(divwrk, F', v)
                axpy!(ξ, v, x)
                axpby!(one(T), g, -ng / pnu, w)
                axpy!(ξ / nu, w, y)
                #
                # u-side bidiagonalization step, and accumulate the norm estimates:
                #
                #   g   ← B v - nu g,   ng ← ‖ g ‖,   g ← g / ng
                #   nB  ← √(nB² + nu² + ng²),   nx ← √(nx² + ξ²)
                #
                mul!(g, B, v, one(T), -nu)
                ng² = dot(g, g);  ng = sqrt(ng²)
                nB² += nu² + ng²; nB = sqrt(nB²)
                nx² += ξ * ξ;     nx = sqrt(nx²)

                if !iszero(ng)
                    rdiv!(g, ng)
                end
                #
                # residual norm and normwise backward error (Rigal–Gaches):
                #
                #   nr ← ng |ξ| = ‖ g - B x ‖
                #   η  ← nr / √(ng0² + nB² nx²)
                #
                nr = ng * abs(ξ)
                iter += 1
                η = nr / sqrt(ng0² + nB² * nx²)
                pnu = nu
                #
                # solved (nr ≤ atol), stagnated (1 + η ≤ 1, i.e. η below 1's ULP), or out of iterations
                #
                if nr ≤ atol
                    status = CRAIG_SOLVED
                elseif one(T) + η ≤ one(T)
                    status = CRAIG_STAGNATED
                elseif iter ≥ itmax
                    status = CRAIG_ITMAX
                end
            end
        end
        #
        # write the residual back into g (u holds the final GK u-side):
        #
        #   g ← -ng ξ g = g - B x
        #
        rmul!(g, -ng * ξ)
    end

    return iter, status
end
