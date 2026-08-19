struct UzawaSolver{UPLO, T, I <: Integer} <: KKTSolver{T}
    F::FChordalTriangular{:N, UPLO, T, I}
    L::BlockSparseMatrix{T, I}
    facwrk::FactorizationWorkspace{T, I}
    divwrk::DivisionWorkspace{T, I}
    itrwrk::CGWorkspace{T}
    r::FVector{T}
    δ::FScalar{T}
    δf::FVector{T}
    δg::FVector{T}
    δx::FVector{T}
    δy::FVector{T}
end

function UzawaSolver(S::ChordalSymbolic{I}, B::BlockSparseMatrix{T, I}) where {T, I <: Integer}
    F = FChordalTriangular{:N, :L, T, I}(S)
    return UzawaSolver(F, B)
end

function UzawaSolver(F::FChordalTriangular{:N, UPLO, T, I}, B::BlockSparseMatrix{T, I}) where {UPLO, T, I <: Integer}
    m, n = size(B)
    facwrk = FactorizationWorkspace(F)
    divwrk = DivisionWorkspace(F, 1)
    itrwrk = CGWorkspace{T}(m, n)
    r = FVector{T}(undef, m)
    δ = FScalar{T}(undef)
    δf = FVector{T}(undef, n)
    δg = FVector{T}(undef, m)
    δx = FVector{T}(undef, n)
    δy = FVector{T}(undef, m)
    return UzawaSolver(F, B' * B, facwrk, divwrk, itrwrk, r, δ, δf, δg, δx, δy)
end

function initkkt!(wrk::UzawaSolver{UPLO, T}, A::BlockSparseMatrix; δ::T, rgmin::T) where {UPLO, T}
    wrk.δ[] = δ
    return inituzw!(wrk.facwrk, wrk.F, wrk.L, A, δ, rgmin)
end

function inituzw!(
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
#   δmax = min {δᵢ : δᵢ ≥ δmin and 1 ≤ i ≤ niter}
#
# Beware! This estimate can be very poor.
#
function kktwindow(hist::AbstractVector{T}, δ::T, fres::T, niter::Int; gtol::T, ftol::T) where {T}
    logδ    = log(δ)
    logfres = log(fres)
    loggtol = log(gtol)
    logftol = log(ftol)

    i = 0
    logδmax = typemin(T)
    logδmin = logδ + logfres - logftol

    while logδmin > logδmax && i < niter
        i += 1
        loggres = log(hist[i])
        logδmax = logδ + (loggtol - loggres) / i
    end

    δmin = exp(logδmin)
    δmax = exp(logδmax)

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
# The solution (x, y) meets the
# following tolerances:
#
#   ‖ A x - Bᵀ y - f ‖ ≤ ftol
#   ‖ B x        - g ‖ ≤ gtol.
#
function solvekkt!(wrk::UzawaSolver, args...; kw...)
    return solvekkt!(wrk.divwrk, wrk.itrwrk, wrk.r, wrk.F, wrk.δf,
            wrk.δg, wrk.δx, wrk.δy, wrk.δ[], args...; kw...)
end

function solvekkt!(
        divwrk::DivisionWorkspace{T},
        itrwrk::CGWorkspace{T},
        r::AbstractVector{T},
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
        xwarm::Bool = false,
        ywarm::Bool = false,
        gtol::T = √eps(T),
        ftol::T = √eps(T),
        stall::T = T(0.5),
        rfmax::Int = 10,
        cgmax::Int = 100,
    ) where {T}

    atol = gtol
    rtol = convert(T, INTERIOR_THETA)

    niter = 0
    npass = 0
    #
    # compute the residuals
    #
    #   [ δf ]   [ f ]   [ A  -Bᵀ ] [ x ]
    #   [ δg ] = [ g ] - [ B   0  ] [ y ]
    #
    copyto!(δf, f)
    copyto!(δg, g)

    if xwarm
        mul!(δf, A, x, -one(T), one(T))
        mul!(δg,           B,      x, -one(T), one(T))
    end

    if ywarm
        mul!(δf, B', y, one(T), one(T))
    end
    #
    # solve for the corrections δx and δy
    #
    #   [ A -Bᵀ ] [ δx ] = [ δf ] (1)
    #   [ B  0  ] [ δy ]   [ δg ]
    #
    n, cgstat = solveuzw!(divwrk, itrwrk, δx, δy, r, L, A, B, δf, δg, δ; atol, itmax = cgmax)
    niter += n

    if cgstat === CG_NUMERICAL_FAILURE
        status = KKT_NUMERICAL_FAILURE
        δmin = δmax = δ
    else
        #
        # update x and y:
        #
        #   x ← x + δx
        #   y ← y + δy
        #
        if xwarm
            axpy!(one(T), δx, x)
        else
            copyto!(x, δx)
        end

        if ywarm
            axpy!(one(T), δy, y)
        else
            copyto!(y, δy)
        end
        #
        # compute the residuals
        #
        #   [ δf ]   [ f ]   [ A  -Bᵀ ] [ x ]
        #   [ δg ] = [ g ] - [ B   0  ] [ y ]
        #
        mulkkt!(δf, δg, A, B, x, y)
        axpby!(one(T), f, -one(T), δf)
        axpby!(one(T), g, -one(T), δg)
        #
        # compute the residual norms
        #
        #   gres = ‖ δg ‖
        #   fres = ‖ δf ‖
        #
        gres = norm(δg)
        fres = norm(δf)

        if gres ≤ gtol && fres ≤ ftol
            status = KKT_SOLVED
        elseif npass == rfmax
            status = KKT_ITMAX
        else
            status = KKT_CONTINUE
        end

        gprv = gres
        fprv = fres
        #
        # estimate the minimum-cost interval
        #
        #   [δmin, δmax];
        #
        # any parameter δ in this interval would have
        # solved (1) to the tolerances
        #
        #   ‖ A δx - Bᵀ δy - δf ‖ ≤ ftol
        #   ‖ B δx         - δg ‖ ≤ gtol
        #
        # with the minimum possible number of triangular solves.
        # beware! this estimate can be very poor.
        #
        δmin, δmax = kktwindow(itrwrk.hist, δ, fres, niter; gtol, ftol)

        while status == KKT_CONTINUE
            #
            # solve for the corrections δx and δy
            #
            #   [ A -Bᵀ ] [ δx ] = [ δf ] (1)
            #   [ B  0  ] [ δy ]   [ δg ]
            #
            n, cgstat = solveuzw!(divwrk, itrwrk, δx, δy, r, L, A, B, δf, δg, δ; atol, rtol, itmax = cgmax)
            niter += n
            npass += 1

            if cgstat === CG_NUMERICAL_FAILURE
                status = KKT_NUMERICAL_FAILURE
            else
                #
                # update x and y:
                #
                #   x ← x + δx
                #   y ← y + δy
                #
                axpy!(one(T), δx, x)
                axpy!(one(T), δy, y)
                #
                # compute the residuals
                #
                #   [ δf ]   [ f ]   [ A  -Bᵀ ] [ x ]
                #   [ δg ] = [ g ] - [ B   0  ] [ y ]
                #
                mulkkt!(δf, δg, A, B, x, y)
                axpby!(one(T), f, -one(T), δf)
                axpby!(one(T), g, -one(T), δg)
                #
                # compute the residual norms
                #
                #   gres = ‖ δg ‖
                #   fres = ‖ δf ‖
                #
                gres = norm(δg)
                fres = norm(δf)

                if gres ≤ gtol && fres ≤ ftol
                    status = KKT_SOLVED
                elseif gres > (1 - stall) * gprv && fres > (1 - stall) * fprv
                    status = KKT_STAGNATED
                elseif npass == rfmax
                    status = KKT_ITMAX
                else
                    status = KKT_CONTINUE
                end

                gprv = gres
                fprv = fres
            end
        end
    end

    return niter, npass, status, δmin, δmax
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
# The solution (x, y) meets the
# following tolerances
#
#   ‖ A x - Bᵀ y - f ‖ ≤ c u ( (‖A‖ + ‖B‖²/δ) ‖x‖ + ‖B‖ ‖y‖ + ‖f‖ ) (1)
#   ‖ B x        - g ‖ ≤ ϵ                                          (2)
#
# where ϵ is a specified tolerance, u is the unit
# roundoff, and c is a small constant depending on
# m and n. The bound in (1) is asymptotically
# linear in 1/δ:
#
#   c u ( (‖A‖ + ‖B‖²/δ) ‖x‖ + ‖B‖ ‖y‖ + ‖f‖ ) ~ O(1/δ).
#
function solveuzw!(
        divwrk::DivisionWorkspace{T},
        itrwrk::CGWorkspace{T},
        x::AbstractVector{T},
        y::AbstractVector{T},
        r::AbstractVector{T},
        L::AbstractMatrix{T},
        A::BlockSparseMatrix{T},
        B::BlockSparseMatrix{T},
        f::AbstractVector{T},
        g::AbstractVector{T},
        δ::T;
        atol::T = √eps(T),
        rtol::T = zero(T),
        itmax::Int = 2size(B, 2),
    ) where {T}
    m, n = size(B)

    @assert length(x) == n
    @assert length(y) == m
    @assert length(r) == m
    @assert length(f) == n
    @assert length(g) == m
    @assert size(L, 1) == n
    @assert size(A, 1) == n

    @assert δ ≥ 0
    @assert atol ≥ 0
    @assert rtol ≥ 0
    @assert itmax ≥ 0

    α = inv(δ)
    #
    # solve for x:
    #
    #   (δ A + Bᵀ B) x =  δ f + Bᵀ g
    #
    copyto!(r, g)
    copyto!(x, f)
    mul!(x, B', r, one(T), δ)
    ldiv!(divwrk, L,  x)
    ldiv!(divwrk, L', x)
    #
    # compute the residual
    #
    #   r ← g - B x
    #
    copyto!(r, g)
    mul!(r, B, x, -one(T), one(T))
    #
    # update tolerance:
    #
    #   atol ← max(atol, rtol ‖r‖)
    #
    if !iszero(rtol)
        atol = max(atol, rtol * norm(r))
    end
    #
    # solve for δx and δy:
    #
    #   [ δ A + Bᵀ B  -Bᵀ ] [ δx ] = [ 0 ]
    #   [          B   0  ] [ δy ]   [ r ]
    #
    # and update
    #
    #   x ← x +   δx
    #   y ←       δy
    #   r ← r - B δx
    #
    niter, status = cg!(itrwrk, B, L, divwrk, x, y, r; atol, itmax)
    #
    # recover y:
    #
    #   y ← α (y + r)
    #
    axpy!(one(T), r, y)
    lmul!(α, y)

    return niter + 1, status
end
