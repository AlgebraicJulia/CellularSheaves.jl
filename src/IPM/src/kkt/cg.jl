@enum CGStatus begin
    CG_CONTINUE
    CG_SOLVED
    CG_NUMERICAL_FAILURE
    CG_ITMAX
end

struct CGWorkspace{T}
    p::FVector{T}
    w::FVector{T}
    hist::FVector{T}
end

function CGWorkspace{T}(m::Integer, n::Integer; itmax::Integer = 2n) where {T}
    p = FVector{T}(undef, m)
    w = FVector{T}(undef, n)
    hist = FVector{T}(undef, itmax + 1)
    return CGWorkspace{T}(p, w, hist)
end

# Solve the KKT system
#
#   [ A -Bᵀ ] [ δx ] = [ 0 ]
#   [ B  0  ] [ δy ]   [ u ]
#
# where A is a n x n positive-definite
# matrix and B is a m × n matrix and L
# is the n × n Cholesky factor of A. The
# solution (δx, δy) meets the following
# tolerances:
#
#   ‖ A δx - Bᵀ δy ‖ ≤ c u ( ‖A‖  ‖δx‖ + ‖B‖ ‖δy‖ )
#   ‖ B δx -     u ‖ ≤ gtol
#
# where c is a small constant depending on n. Then
# update x, y, and u:
#
#   x ← x +   δx
#   y ←       δy
#   u ← u - B δx
#
function cg!(wrk::CGWorkspace{T}, B, L, divwrk, x, y, u; kwargs...) where {T}
    return cg!(wrk.p, wrk.w, B, L, divwrk, x, y, u, wrk.hist; kwargs...)
end

function cg!(
        p::AbstractVector{T},
        w::AbstractVector{T},
        B::AbstractMatrix{T},
        L::AbstractMatrix{T},
        divwrk::DivisionWorkspace{T},
        x::AbstractVector{T},
        y::AbstractVector{T},
        u::AbstractVector{T},
        hist::AbstractVector{T};
        atol::T  = √eps(T),
        itmax::Int = 2size(B, 2),
    ) where {T}

    @assert size(L, 1) == size(B, 2)
    @assert size(L, 2) == size(B, 2)

    @assert length(w) == size(B, 2)
    @assert length(x) == size(B, 2)

    @assert length(p) == size(B, 1)
    @assert length(y) == size(B, 1)
    @assert length(u) == size(B, 1)

    niter = 0

    nu² = dot(u, u); hist[1] = nu = sqrt(nu²)
    np² = nu²

    fill!(y, zero(T))

    status = CG_SOLVED

    if nu > atol
        status = CG_ITMAX

        if itmax > 0
            copyto!(p, u)

            status = CG_CONTINUE

            while status == CG_CONTINUE
                niter += 1
                # solve for w:
                #
                #   L w = Bᵀ p
                #
                mul!(w, B', p)
                ldiv!(divwrk, L, w)
                #
                # compute the squared elliptic norm of the direction p
                #
                #   ep² ← ‖ p ‖²_M = pᵀ B A⁻¹ Bᵀ p = pᵀ M p
                #
                ep² = dot(w, w)

                if ep² ≤ eps(T) * np²
                    status = CG_NUMERICAL_FAILURE
                else
                    # solve for w:
                    #
                    #   Lᵀ w = w
                    #
                    ldiv!(divwrk, L', w)
                    #
                    # solution update: advance the CG iterate
                    #
                    #   α ← ep⁻² nu²,  y ← y + α p,  x ← x + α w,  u ← u − α B w
                    #
                    α = nu² / ep²
                    axpy!(α, p, y)
                    axpy!(α, w, x)
                    mul!(u, B, w, -α, one(T))
                    #
                    # compute the norm of the residual u:
                    #
                    #   nu ← ‖ u ‖
                    #
                    pnu² = nu²
                    nu² = dot(u, u); hist[niter + 1] = nu = sqrt(nu²)

                    if nu ≤ atol || one(T) + nu ≤ one(T)
                        status = CG_SOLVED
                    elseif niter ≥ itmax
                        status = CG_ITMAX
                    else
                        # direction update: compute the next conjugate direction
                        #
                        #   β ← pnu⁻² nu²,  p ← u + β p,  np² ← nu² + β² np²
                        #
                        β = nu² / pnu²
                        axpby!(one(T), u, β, p)
                        np² = nu² + β^2 * np²
                    end
                end
            end
        end
    end

    return niter, status
end
