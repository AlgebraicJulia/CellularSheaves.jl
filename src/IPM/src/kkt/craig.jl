# CRAIG — Golub–Kahan bidiagonalization for the minimum-norm least-squares problem
#
#   minimize ‖x‖  subject to  A x = b
#
# vendored from Krylov.jl and specialized to our use: real T, no left preconditioner (M = I), and the
# right preconditioner (F Fᵀ)⁻¹ applied inline by two triangular solves against the augmented-matrix
# Cholesky factor F (so we pass F + its DivisionWorkspace rather than a generic operator), λ = 0, and no
# warm start / timing / logging / callbacks. We own it so a later pass can harvest the bidiagonal
# coefficients αₖ, βₖ (the pricing-model Ritz / spectral readings) without an external internal API.
#
# The Krylov @k* BLAS wrappers become their idiomatic standard-library equivalents: kaxpy!→axpy!,
# kaxpby!→axpby!, kdot→dot, kscal!→rmul!, kdiv!(x,s)→rdiv!(x,s), kmul!→mul!, kfill!→fill!,
# kcopy!→copyto!, knorm→norm. Two of these are NOT bit-identical to Krylov: `norm` of a short BlasFloat
# vector (length < NRM2_CUTOFF = 32, as craig's dual-space u/w often are — row dim of B, e.g. m=21 for
# e04) takes the generic_norm2 branch, which differs from Krylov's raw BLAS.nrm2 by 1 ULP; and rdiv!(x,s)
# divides where kdiv! multiplies by the reciprocal. Either gap perturbs the Golub–Kahan recurrence in the
# last bits, but CRAIG still converges to the same atol — so this port is validated by convergence and
# solver iteration cost, not by byte-identity. (The byte-identical intermediate used nrm2 + rmul!(inv).)
#
# btol is the backward-error stopping tolerance (bkwerr ≤ btol); we pass 0 — the pricing model uses an
# ABSOLUTE residual tolerance (atol), so the backward-error stop is disabled. The knob is kept because it
# costs one comparison per iteration (bkwerr is computed anyway for the machine-precision stop).

# Why the CRAIG iteration halted. CRAIG_SOLVED covers both a met residual/backward-error tolerance and
# the zero-residual start; the rest are the diagnostic exits (inconsistent right-hand side, operator too
# ill-conditioned, dimension-based iteration cap hit).
@enum CraigStatus begin
    CRAIG_SOLVED
    CRAIG_INCONSISTENT
    CRAIG_ILL_CONDITIONED
    CRAIG_ITMAX
end

mutable struct CraigWorkspace{T}
    m::Int
    n::Int
    x::Vector{T}     # solution x (n)  — the min-norm correction
    y::Vector{T}     # solution y (m)
    Nv::Vector{T}    # Golub–Kahan v-side, before preconditioning (n)
    Aᴴu::Vector{T}   # Aᴴ u (n)
    v::Vector{T}     # N·Nv (n)
    w::Vector{T}     # y-recurrence direction (m)
    u::Vector{T}     # Golub–Kahan u-side (m)  [ = Mu, M = I ]
    Av::Vector{T}    # A v (m)
    niter::Int
end

function CraigWorkspace{T}(m::Integer, n::Integer) where {T}
    return CraigWorkspace{T}(m, n, zeros(T, n), zeros(T, m), zeros(T, n), zeros(T, n),
                             zeros(T, n), zeros(T, m), zeros(T, m), zeros(T, m), 0)
end

function craig!(wrk::CraigWorkspace{T}, A, F, divwrk, b;
                atol::T = √eps(T), rtol::T = √eps(T), btol::T = √eps(T),
                conlim::T = 1 / √eps(T), itmax::Int = 0) where {T}
    m, n = size(A)
    Aᴴ = A'
    x, y, Nv, Aᴴu, v, w, u, Av = wrk.x, wrk.y, wrk.Nv, wrk.Aᴴu, wrk.v, wrk.w, wrk.u, wrk.Av

    fill!(x, zero(T))
    fill!(y, zero(T))

    copyto!(u, b)                    # u ← b   (u ≡ Mu, M = I)
    β₁ = norm(u)
    rNorm = β₁
    if β₁ == 0
        wrk.niter = 0
        return 0, CRAIG_SOLVED       # x = 0 is a zero-residual solution
    end
    β₁² = β₁^2
    β  = β₁
    θ  = β₁
    ξ  = -one(T)
    ρ_prev = one(T)

    rdiv!(u, β₁)                     # u ← u / β₁
    fill!(Nv, zero(T))
    fill!(w,  zero(T))

    Anorm² = zero(T); Anorm = zero(T)
    Dnorm² = zero(T); Acond = zero(T)
    xNorm² = zero(T); xNorm = zero(T)

    iter = 0
    itmax == 0 && (itmax = m + n)

    ɛ_c  = atol + rtol * rNorm       # consistent-system residual tolerance
    ctol = conlim > 0 ? 1 / conlim : zero(T)
    bkwerr = one(T)                  # backward error ‖r‖ / √(‖b‖² + ‖A‖²‖x‖²)

    solved = (one(T) + bkwerr ≤ one(T)) | (bkwerr ≤ btol) | (rNorm ≤ ɛ_c) |
             (rNorm ≤ btol + atol * Anorm * xNorm / β₁)
    ill_cond = false
    inconsistent = false
    tired = iter ≥ itmax

    while !(solved || inconsistent || ill_cond || tired)
        # 1.  αₖ₊₁ N vₖ₊₁ = Aᴴ uₖ − βₖ N vₖ
        mul!(Aᴴu, Aᴴ, u)
        axpby!(one(T), Aᴴu, -β, Nv)  # Nv ← Aᴴu − β Nv
        # preconditioner  v ← (F Fᵀ)⁻¹ Nv  (F is the augmented-matrix Cholesky factor)
        copyto!(v, Nv)
        ldiv!(divwrk, F,  v)
        ldiv!(divwrk, F', v)
        α = sqrt(dot(v, Nv))         # (F Fᵀ)⁻¹-elliptic norm
        if α == 0
            inconsistent = true
            continue
        end
        rdiv!(v,  α)
        rdiv!(Nv, α)
        Anorm² += α * α

        ρ = α
        ξ = -θ / ρ * ξ
        axpy!(ξ, v, x)                       # x ← x + ξ v
        axpby!(one(T), u, -θ / ρ_prev, w)    # w ← u − (θ/ρ_prev) w
        axpy!(ξ / ρ, w, y)                   # y ← y + (ξ/ρ) w
        Dnorm² += norm(w)

        # 2.  βₖ₊₁ M uₖ₊₁ = A vₖ − αₖ M uₖ
        mul!(Av, A, v)
        axpby!(one(T), Av, -α, u)            # u ← A v − α u
        β = norm(u)
        β ≠ 0 && rdiv!(u, β)
        θ = β

        Anorm² += β * β
        Anorm  = sqrt(Anorm²)
        Acond  = Anorm * sqrt(Dnorm²)
        xNorm² += ξ * ξ
        xNorm  = sqrt(xNorm²)
        rNorm  = β * abs(ξ)                   # r = −β ξ u
        iter += 1

        bkwerr = rNorm / sqrt(β₁² + Anorm² * xNorm²)
        ρ_prev = ρ

        solved = (one(T) + bkwerr ≤ one(T)) | (bkwerr ≤ btol) | (rNorm ≤ ɛ_c) |
                 (rNorm ≤ btol + atol * Anorm * xNorm / β₁)
        ill_cond = (one(T) + inv(Acond) ≤ one(T)) | (inv(Acond) ≤ ctol)
        tired = iter ≥ itmax
    end

    status = inconsistent ? CRAIG_INCONSISTENT :
             ill_cond     ? CRAIG_ILL_CONDITIONED :
             solved       ? CRAIG_SOLVED :
                            CRAIG_ITMAX
    wrk.niter = iter
    return iter, status
end
