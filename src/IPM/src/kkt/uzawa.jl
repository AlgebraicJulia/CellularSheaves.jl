struct UzawaSolver{UPLO, T, I <: Integer} <: KKTSolver{T}
    F::FChordalTriangular{:N, UPLO, T, I}
    L::BlockSparseMatrix{T, I}
    facwrk::FactorizationWorkspace{T, I}
    divwrk::DivisionWorkspace{T, I}
    itrwrk::CGWorkspace{T}
    hist::FVector{T}
    δ::FScalar{T}
    δf::FVector{T}
    δg::FVector{T}
    δx::FVector{T}
    δy::FVector{T}
end

function UzawaSolver(S::ChordalSymbolic{I}, B::BlockSparseMatrix{T, I}; cgmax::Integer = 2size(B, 2), rfmax::Integer = 10) where {T, I <: Integer}
    F = FChordalTriangular{:N, :L, T, I}(S)
    return UzawaSolver(F, B; cgmax, rfmax)
end

function UzawaSolver(F::FChordalTriangular{:N, UPLO, T, I}, B::BlockSparseMatrix{T, I}; cgmax::Integer = 2size(B, 2), rfmax::Integer = 10) where {UPLO, T, I <: Integer}
    m, n = size(B)
    facwrk = FactorizationWorkspace(F)
    divwrk = DivisionWorkspace(F, 1)
    itrwrk = CGWorkspace{T}(m, n; itmax = cgmax)
    hist = FVector{T}(undef, rfmax + 2)
    δ = FScalar{T}(undef)
    δf = FVector{T}(undef, n)
    δg = FVector{T}(undef, m)
    δx = FVector{T}(undef, n)
    δy = FVector{T}(undef, m)
    return UzawaSolver(F, B' * B, facwrk, divwrk, itrwrk, hist, δ, δf, δg, δx, δy)
end

function initkkt!(wrk::UzawaSolver{UPLO, T}, A::BlockSparseMatrix; δ::T, rgmin::T) where {UPLO, T}
    wrk.δ[] = δ
    return initkkt!(wrk.facwrk, wrk.F, wrk.L, A, δ, rgmin)
end

function initkkt!(
        facwrk::FactorizationWorkspace{T},
        F::ChordalTriangular{:N, UPLO, T},
        L::BlockSparseMatrix{T},
        A::BlockSparseMatrix{T},
        δ::T,
        rgmin::T,
    ) where {UPLO, T}
    @assert size(F, 1) == size(L, 1) == size(A, 1)
    #
    # factorize the augmented matrix
    #
    #   F Fᵀ = δ A + Bᵀ B
    #
    copyto!(F, L)
    axpy!(δ, A, F)
    info = cholesky!(facwrk, F; check=false)

    ρ = zero(T)

    if !iszero(info)
        #
        # on failure, factorize the perturbed matrix
        #
        #   F Fᵀ = δ A + Bᵀ B + ρ I
        #
        # with an increasing sequence of perturbations
        #
        #   ρ = rgmin, 8 rgmin, 8² rgmin, ...
        #
        ρ = rgmin

        for _ in 1:10
            #
            # factorize the perturbed matrix
            #
            #   F Fᵀ = δ A + Bᵀ B + ρ I
            #
            copyto!(F, L)
            axpy!(δ, A, F)
            axpy!(ρ, I, F)
            info = cholesky!(facwrk, F; check=false)

            iszero(info) && break
            #
            # on failure, increase ρ
            #
            #   ρ ← 8 ρ
            #
            ρ *= 8
        end
    end

    return iszero(info), ρ
end

#
# Estimate the minimum-cost interval
#
#   [δmin, δmax];
#
# any parameter δ = 1/α in this interval would have
# solved the KKT system
#
#   [ A -Bᵀ ] [ x ] = [ f ]
#   [ B  0  ] [ y ]   [ g ]
#
# to the tolerances
#
#   ‖ A x - Bᵀ y - f ‖ ≤ ftol
#   ‖ B x        - g ‖ ≤ gtol
#
# with the minimum possible number of triangular solves.
# The interval is estimated using the equations
#
#   δmin ≈ δ    fres  / ftol
#   δᵢ   ≈ δ ⁱ√(gtol  / gresᵢ)
#
# where δmin is the smallest number δ such that ‖ B x - g ‖ ≤ gtol,
# and δᵢ is the smallest number δ such that exactly 2i triangular
# solves are performed. δmax defined to be
#
#   δmax = min {δᵢ : δᵢ ≥ δmin and 1 ≤ i ≤ niter + 1}
#
# Beware! This estimate can be very poor.
#
function kktwindow(rhist::AbstractVector{T}, chist::AbstractVector{T}, δ::T, nir::Int, ncg::Int; gtol::T, ftol::T) where {T}
    δmin = δmax = T(NaN)

    if nir ≥ 1
        logδ    = log(δ)
        logfres = log(rhist[2])
        loggtol = log(gtol)
        logftol = log(ftol)

        i = 0
        logδmax = typemin(T)
        logδmin = logδ + logfres - logftol

        while logδmin > logδmax && i < ncg + 1
            i += 1
            loggres = log(chist[i])
            logδmax = logδ + (loggtol - loggres) / i
        end

        δmin = exp(logδmin)
        δmax = exp(logδmax)
    end

    return δmin, δmax
end

#
# Solve the KKT system
#
#   [ A -Bᵀ ] [ x ] = [ f ]
#   [ B  0  ] [ y ]   [ g ]
#
# where A is a n x n positive-definite
# matrix and B is a m × n matrix, δ is a
# positive number, and L is the n × n
# Cholesky factor of the augmented matrix
#
#   δ A + Bᵀ B = L Lᵀ.
#
# The solution (x, y) meets the following
# tolerances:
#
#   ‖ A x - Bᵀ y - f ‖ ≤ ftol
#   ‖ B x        - g ‖ ≤ gtol.
#
function solvekkt!(wrk::UzawaSolver, x, y, A, B, f, g; kw...)
    return solvekkt!(wrk.divwrk, wrk.itrwrk, wrk.hist, wrk.F, wrk.δf, wrk.δg,
            wrk.δx, wrk.δy, wrk.δ[], x, y, A, B, f, g; kw...)
end

function solvekkt!(
        divwrk::DivisionWorkspace{T},
        itrwrk::CGWorkspace{T},
        hist::AbstractVector{T},
        L::AbstractMatrix{T},
        δf::AbstractVector{T},
        δg::AbstractVector{T},
        δx::AbstractVector{T},
        δy::AbstractVector{T},
        δ::T,
        x::AbstractVector{T},
        y::AbstractVector{T},
        A::BlockSparseMatrix{T},
        B::BlockSparseMatrix{T},
        f::AbstractVector{T},
        g::AbstractVector{T};
        warm::Bool = false,
        gtol::T = √eps(T),
        ftol::T = √eps(T),
        stall::T = T(0.5),
        rfmax::Int = 10,
        cgmax::Int = 100,
    ) where {T}
    m, n = size(B)

    @assert length(x) == n
    @assert length(y) == m
    @assert length(f) == n
    @assert length(g) == m
    @assert size(L, 1) == n
    @assert size(A, 1) == n

    @assert δ > 0
    @assert gtol ≥ 0
    @assert ftol ≥ 0
    @assert rfmax ≥ 0
    @assert cgmax ≥ 0

    α = inv(δ)
    nir = 0
    ncg = 0
    #
    # compute the residuals
    #
    #   [ δf ] = [ f ] - [ A  -Bᵀ ] [ x ]
    #   [ δg ]   [ g ]   [ B   0  ] [ y ]
    #
    # and the residual norm
    #
    #   fres = ‖ δf ‖
    #
    if warm
        mulkkt!(δf, δg, A, B, x, y)
        axpby!(one(T), f, -one(T), δf)
        axpby!(one(T), g, -one(T), δg)
    else
        copyto!(δf, f)
        copyto!(δg, g)
        fill!(x, zero(T))
        fill!(y, zero(T))
    end

    fprv = fres = norm(δf); hist[1] = fres

    if fres ≤ ftol
        status = KKT_SOLVED
    elseif nir ≥ rfmax
        status = KKT_RFMAX
    else
        status = KKT_CONTINUE
    end
    #
    # refine the KKT system
    #
    #   [ A -Bᵀ ] [ x ] = [ f ]
    #   [ B  0  ] [ y ]   [ g ]
    #
    # until
    #
    #   ‖ A x - Bᵀ y - f ‖ ≤ ftol
    #
    while status === KKT_CONTINUE
        nir += 1
        #
        # solve for the corrections δx and δy
        #
        #   [ A -Bᵀ ] [ δx ] = [ δf ]
        #   [ B  δ  ] [ δy ]   [ δg ]
        #
        copyto!(δx, δf)
        rmul!(δx, δ)
        mul!(δx, B', δg, one(T), one(T))
        ldiv!(divwrk, L,  δx)
        ldiv!(divwrk, L', δx)
        mul!(δg, B, δx, -one(T), one(T))
        #
        # update x and y:
        #
        #   x ← x + δx
        #   y ← y + δy
        #
        axpy!(one(T), δx, x)
        axpy!(α, δg, y)
        #
        # compute the residuals
        #
        #   [ δf ] = [ f ] - [ A  -Bᵀ ] [ x ]
        #   [ δg ]   [ g ]   [ B   0  ] [ y ]
        #
        # and the residual norm
        #
        #   fres = ‖ δf ‖
        #
        mulkkt!(δf, δg, A, B, x, y)
        axpby!(one(T), f, -one(T), δf)
        axpby!(one(T), g, -one(T), δg)
        fres = norm(δf)
        hist[nir + 1] = fres

        if fres ≤ ftol
            status = KKT_SOLVED
        elseif fres > (one(T) - stall) * fprv
            status = KKT_STAGNATED
        elseif nir ≥ rfmax
            status = KKT_RFMAX
        else
            status = KKT_CONTINUE
        end

        fprv = fres
    end
    #
    # solve for δx and δy
    #
    #   [ δ A + Bᵀ B  -Bᵀ ] [ δx ] = [ 0  ]
    #   [          B   0  ] [ δy ]   [ δg ]
    #
    # to the tolerances
    #
    #   ‖ A δx - Bᵀ δy ‖ ≤ c u ( ( δ ‖A‖ + ‖B‖² ) ‖δx‖ + ‖B‖ ‖δy‖ )
    #   ‖ B δx -    δg ‖ ≤ gtol
    #
    # and update x and y:
    #
    #   x ← x +     δx
    #   y ← y + 1/δ δy - 1/δ B δx
    #
    if status === KKT_SOLVED
        axpy!(-α, δg, y)
        ncg, cgstat = cg!(itrwrk, B, L, divwrk, x, δy, δg; atol = gtol, itmax = cgmax)
        axpy!(α, δy, y)
        axpy!(α, δg, y)

        if cgstat === CG_NUMERICAL_FAILURE
            status = KKT_NUMERICAL_FAILURE
        elseif cgstat === CG_ITMAX
            status = KKT_ITMAX
        end
    end
    #
    # estimate the minimum-cost interval
    #
    #   [δmin, δmax];
    #
    # any parameter δ in this interval would have
    # solved the KKT system
    #
    #   [ A -Bᵀ ] [ x ] = [ f ]
    #   [ B  0  ] [ y ]   [ g ]
    #
    # to the tolerances
    #
    #   ‖ A x - Bᵀ y - f ‖ ≤ ftol
    #   ‖ B x        - g ‖ ≤ gtol
    #
    # with no IR iterations and the minimum number
    # of CG iterations.
    #
    # beware! this estimate can be very poor.
    #
    if status === KKT_NUMERICAL_FAILURE
        δmin = δmax = δ
    else
        δmin, δmax = kktwindow(hist, itrwrk.hist, δ, nir, ncg; gtol, ftol)
    end

    return ncg, nir, status, δmin, δmax
end
