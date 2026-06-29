abstract type AbstractTDCone <: AbstractCone end

function workspacesize(::AbstractTDCone, n::Integer)
    @assert n == 3
    return 15
end

struct AbstractTDConeCache{C <: AbstractTDCone, T} <: AbstractCache{C}
    cone::C
    #
    # The factor of the barrier Hessian
    #
    #   P f''(p) Pᵀ = L Lᵀ
    #
    L::FMatrixView{T}
    #
    # The dual "shadow" iterate
    #
    #   d* = -f'(p)
    #
    sd::FVectorView{T}
    #
    # warm-start for computing the primal
    # "shadow" iterate p*, which solves
    #
    #   -f'(p*) = d
    #
    seed::FScalarView{T}
end

# Compute the coefficient t in the rank-1 term
# tzzᵀ of the Tuncel scaling matrix M:
#
#   M = ⟨p, d⟩⁻¹   d  dᵀ
#     + ⟨δp,δd⟩⁻¹ δd δdᵀ
#     + t          z  zᵀ,
#
# where
#
#   δp = p - μ p*
#   δd = d - μ d*
#
function tdbfgs(
        L::AbstractMatrix{T},
        sp::AbstractVector{T},
        sd::AbstractVector{T},
        z::AbstractVector{T},
        p::AbstractVector{T},
        μv::T,
        μt::T,
        cache::AbstractTDConeCache,
        w::AbstractVector{T},
        Hw::AbstractVector{T},
        Hz::AbstractVector{T},
    ) where {T}
    t = zero(T)
    #
    # compute the gap direction
    #
    #   w = p* - μ* p
    #
    # where is the dual centrality parameter
    #
    #   μ* = ⟨p*, d*⟩ / ν.
    #
    copy3!(w, sp)
    axpy3!(-μt, p, w)
    #
    # compute the norm
    #
    #   ⟨w, f''(p) w⟩ = ‖Rᵀ w‖²
    #
    tdhessmul!(Hw, L, w, cache)
    wHw = dot3(w, Hw)

    if wHw > 0
        #
        # compute the norm
        #
        #   ⟨z, f''(p) z⟩ = ‖Rᵀ z‖²
        #
        tdhessmul!(Hz, L, z, cache)
        fppzz = dot3(z, Hz)
        #
        # compute the dot product
        #
        #   ⟨d*, z⟩
        #
        sdz = dot3(sd, z)
        #
        # compute the dot product
        #
        #   ⟨w, f''(p) z⟩ = ⟨Rᵀ w, Rᵀ z⟩
        #
        wHz = dot3(Hw, z)
        #
        # compute t:
        #
        #   t = μ ⟨z, f''(p) z⟩
        #     - μ ⟨d*,       z⟩² / ν
        #     - μ ⟨w, f''(p) z⟩² / ⟨w, f''(p) w⟩
        #
        t = μv * (fppzz - sdz^2 / 3 - wHz^2 / wHw)
    end

    return t
end

# Assemble the Tuncel scaling matrix
#
#   M = ⟨p, d⟩⁻¹   d  dᵀ
#     + ⟨δp,δd⟩⁻¹ δd δdᵀ
#     + t          z  zᵀ,
#
# where
#
#   δp = p - μ p*
#   δd = d - μ d*
#
function tdscale!(
        H::AbstractMatrix{T},
        L::AbstractMatrix{T},
        sd::AbstractVector{T},
        seed::T,
        p::AbstractVector{T},
        d::AbstractVector{T},
        cache::AbstractTDConeCache,
        wrk::ConeWorkspace{T},
    ) where {T}

    sp = view(wrk.data,  1:3)
    z  = view(wrk.data,  4:6)
    w  = view(wrk.data,  7:9)
    δp = view(wrk.data, 10:12)
    δd = view(wrk.data, 13:15)
    #
    # compute the analytic factor
    #
    #   f''(p) = R Rᵀ
    #
    tdfact!(L, p, cache)
    #
    # compute the "shadow" dual
    #
    #   d* = -f'(p)
    #
    tdbarrgrad!(sd, p, cache)
    lmul3!(-1, sd)
    #
    # compute the "shadow" primal, solving
    #
    #   d = -f'(p*)
    #
    seed_new = tddualgrad!(sp, seed, d, cache)
    #
    # compute the centrality parameters
    #
    #   μ  = ⟨p,  d ⟩ / ν
    #   μ* = ⟨p*, d*⟩ / ν
    #
    μv = dot3(p,  d)  / 3
    μt = dot3(sp, sd) / 3
    #
    # compute the cross-product
    #
    #   z = p × p*
    #
    cross3!(z, p, sp)
    #
    # compute the sine of the angle θ between p and p*:
    #
    #   ‖p × p*‖ / (‖p‖ ‖p*‖) = sin(θ)
    #
    # when this quantity is small, the iterate is close to
    # the central path and the term δd δdᵀ / ⟨δp, δd⟩ term in M
    # becomes innaccurate due to cancellation in the difference
    #
    #   δp = p - μ p*
    #
    # in this case, we fall back to the approximation
    #
    #   M ≈ μ f''(p)
    #
    nz  = norm3(z)
    np  = norm3(p)
    nsp = norm3(sp)

    if nz < eps(T) * (np * nsp + eps(T))
        #
        # approximate M by
        #
        #    M ≈ μ f''(p) = μ R Rᵀ
        #
        tdgram!(H, L, μv, cache)
    else
        #
        # normalize z:
        #
        #   z = z / ‖z‖
        #
        ldiv3!(nz, z)
        #
        # compute the coefficent t in the rank-1 term
        #
        #   tzzᵀ
        #
        t = tdbfgs(L, sp, sd, z, p, μv, μt, cache, w, δp, δd)

        if t ≤ 0 || !isfinite(t)
            #
            # approximate M by
            #
            #    M ≈ μ f''(p) = μ R Rᵀ
            #
            tdgram!(H, L, μv, cache)
        else
            #
            # construct M:
            #
            #   M = ⟨p, d⟩⁻¹   d  dᵀ
            #     + ⟨δp,δd⟩⁻¹ δd δdᵀ
            #     + t          z  zᵀ,
            #
            copy3!(δp, p); axpy3!(-μv, sp, δp)
            copy3!(δd, d); axpy3!(-μv, sd, δd)

            pd_dot = 3μv
            δ_dot  = dot3(δp, δd)

            ger3!(H,  d,  d, inv(pd_dot), 0)
            ger3!(H, δd, δd, inv(δ_dot),  1)
            ger3!(H,  z,  z, t,           1)
        end
    end

    return seed_new
end

# Compute the Mehrotra corrector term
#
#   -d - σμ f'(p) - η,
#
# where η is the third-order correction
#
#   η = -½ f'''(p)[Δp, f''(p)⁻¹ Δd].
#
function tdcorr!(
        r::AbstractVector{T},
        L::AbstractMatrix{T},
        sd::AbstractVector{T},
        p::AbstractVector{T},
        d::AbstractVector{T},
        Δp::AbstractVector{T},
        Δd::AbstractVector{T},
        σμ::Real,
        cache::AbstractTDConeCache,
        wrk::ConeWorkspace{T},
    ) where {T}

    v =         view(wrk.data, 1:3)
    η =         view(wrk.data, 4:6)
    D = reshape(view(wrk.data, 7:15), 3, 3)
    #
    # solve for v in
    #
    #   f''(p) v =  Δd
    #
    # using the analytic factorization
    #
    #   f''(p) = R Rᵀ.
    #
    copy3!(v, Δd)
    tdhessldiv!(L, v, cache)
    #
    # compute the third-order correction
    #
    #   η = -½ f'''(p)[Δp, v]
    #
    tdbarrthird!(D, p, Δp, cache)
    mul3!(η, D, v, -0.5, 0)
    #
    # compute the Mehrotra corrector term
    #
    #   -d - σμ f'(p) - η
    #
    copy3!(r, d)
    axpby3!(σμ, sd, -1, r)
    axpy3!(-1, η, r)

    return r
end

# use bisection to find the largest number 0 < τ ≤ 1
# such that
#
#   p + τ Δp
#
# is in the cone (or its dual).
function tdmaxstep(incone, p::AbstractVector{T}, Δp::AbstractVector{T}, cache::AbstractTDConeCache, wrk::ConeWorkspace{T}) where {T}
    w = view(wrk.data, 1:3)

    τ = binarysearchlast(zero(T), one(T), eps(T), 53) do τ
        copy3!(w, p)
        axpy3!(τ, Δp, w)
        return incone(w, cache)
    end

    return τ
end

function tdmaxsteps(
        p::AbstractVector{T},
        Δp::AbstractVector{T},
        d::AbstractVector{T},
        Δd::AbstractVector{T},
        cache::AbstractTDConeCache,
        wrk::ConeWorkspace{T},
    ) where {T}
    return tdmaxstep(tdincone, p, Δp, cache, wrk), tdmaxstep(tdindual, d, Δd, cache, wrk)
end

#
# AbstractCone Interface
#

function degree(::AbstractTDCone, n::Integer)
    @assert n == 3
    return 3
end

function cachesize(::AbstractTDCone, n::Integer)
    @assert n == 3
    return 13
end

function scale!(H::AbstractMatrix{T}, p::AbstractVector{T}, d::AbstractVector{T}, cache::AbstractTDConeCache{C, T}, wrk::ConeWorkspace{T}) where {C, T}
    cache.seed[] = tdscale!(H, cache.L, cache.sd, cache.seed[], p, d, cache, wrk)
    return H
end

function corr!(r::AbstractVector{T}, p::AbstractVector{T}, d::AbstractVector{T}, Δp::AbstractVector{T}, Δd::AbstractVector{T}, σμ::Real, cache::AbstractTDConeCache{C, T}, wrk::ConeWorkspace{T}) where {C, T}
    return tdcorr!(r, cache.L, cache.sd, p, d, Δp, Δd, σμ, cache, wrk)
end

function maxsteps(p::AbstractVector{T}, Δp::AbstractVector{T}, d::AbstractVector{T}, Δd::AbstractVector{T}, cache::AbstractTDConeCache{C, T}, wrk::ConeWorkspace{T}) where {C, T}
    return tdmaxsteps(p, Δp, d, Δd, cache, wrk)
end

include("exp.jl")
include("pow.jl")
