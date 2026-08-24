abstract type KKTSolver{T} end

struct UzawaSolver{UPLO, T, I <: Integer} <: KKTSolver{T}
    F::FChordalTriangular{:N, UPLO, T, I}
    L::BlockSparseMatrix{T, I}
    facwrk::FactorizationWorkspace{T, I}
    divwrk::DivisionWorkspace{T}
    itrwrk::CGWorkspace{T}
    hist::FVector{T}
    δ::FScalar{T}
    δf::FVector{T}
    δg::FVector{T}
    δx::FVector{T}
    δy::FVector{T}
end

struct PivotedUzawaSolver{UPLO, T, I <: Integer} <: KKTSolver{T}
    F::FChordalTriangular{:N, UPLO, T, I}
    L::BlockSparseMatrix{T, I}
    facwrk::FactorizationWorkspace{T, I}
    divwrk::DivisionWorkspace{T}
    itrwrk::CGWorkspace{T}
    hist::FVector{T}
    δ::FScalar{T}
    δf::FVector{T}
    δg::FVector{T}
    δx::FVector{T}
    δy::FVector{T}
    S::ChordalSymbolic{I}
    P::FPermutation{I}
    w::FVector{T}
end

function UzawaSolver(S::ChordalSymbolic{I}, B::BlockSparseMatrix{T, I}; cgmax::Integer = 2size(B, 2), irmax::Integer = 10) where {T, I <: Integer}
    F = FChordalTriangular{:N, :L, T, I}(S)
    m, n = size(B)
    facwrk = FactorizationWorkspace(F)
    divwrk = DivisionWorkspace{T}(F.S, 1)
    itrwrk = CGWorkspace{T}(m, n; itmax = cgmax)
    hist = FVector{T}(undef, irmax + 2)
    δ = FScalar{T}(undef)
    δf = FVector{T}(undef, n)
    δg = FVector{T}(undef, m)
    δx = FVector{T}(undef, n)
    δy = FVector{T}(undef, m)
    return UzawaSolver(F, B' * B, facwrk, divwrk, itrwrk, hist, δ, δf, δg, δx, δy)
end

function PivotedUzawaSolver(S::ChordalSymbolic{I}, B::BlockSparseMatrix{T, I}; cgmax::Integer = 2size(B, 2), irmax::Integer = 10) where {T, I <: Integer}
    F = FChordalTriangular{:N, :L, T, I}(S)
    m, n = size(B)
    facwrk = FactorizationWorkspace(F)
    divwrk = DivisionWorkspace{T}(F.S, 1)
    itrwrk = CGWorkspace{T}(m, n; itmax = cgmax)
    hist = FVector{T}(undef, irmax + 2)
    δ = FScalar{T}(undef)
    δf = FVector{T}(undef, n)
    δg = FVector{T}(undef, m)
    δx = FVector{T}(undef, n)
    δy = FVector{T}(undef, m)
    P = FPermutation{I}(n)
    w = FVector{T}(undef, n)
    return PivotedUzawaSolver(F, B' * B, facwrk, divwrk, itrwrk, hist, δ, δf, δg, δx, δy, copy(F.S), P, w)
end

# ============================================================================
# initkkt!
# ============================================================================

function initkkt!(wrk::UzawaSolver{UPLO, T}, A::BlockSparseMatrix; δ::T) where {UPLO, T}
    wrk.δ[] = δ
    return initkkt!(wrk.facwrk, wrk.F, wrk.L, A, δ)
end

function initkkt!(
        facwrk::FactorizationWorkspace,
        F::ChordalTriangular,
        L::BlockSparseMatrix,
        A::BlockSparseMatrix,
        δ::Real,
    )
    @assert size(F, 1) == size(L, 1) == size(A, 1)
    #
    # factorize the augmented matrix
    #
    #   F Fᵀ = δ A + Bᵀ B
    #
    copyto!(F, L)
    axpy!(δ, A, F)
    info = cholesky!(facwrk, F; check=false)

    return iszero(info)
end

function initkkt!(wrk::PivotedUzawaSolver{UPLO, T}, A::BlockSparseMatrix; δ::T) where {UPLO, T}
    wrk.δ[] = δ
    return initkkt!(wrk.facwrk, wrk.F, wrk.S, wrk.P, wrk.L, A, δ)
end

function initkkt!(
        facwrk::FactorizationWorkspace,
        F::ChordalTriangular,
        S::ChordalSymbolic,
        P::Permutation,
        L::BlockSparseMatrix,
        A::BlockSparseMatrix,
        δ::Real,
    )
    @assert size(F, 1) == size(L, 1) == size(A, 1)

    copyto!(F.S, S)
    copyto!(P.perm, axes(F, 1))
    copyto!(P.invp, axes(F, 1))
    #
    # factorize the augmented matrix
    #
    #   F Fᵀ = P (δ A + Bᵀ B) Pᵀ
    #
    copyto!(F, L)
    axpy!(δ, A, F)
    cholesky!(facwrk, F, P, RowMaximum(); check = false, tol = 0)

    return true
end

# ============================================================================
# kktwindow
# ============================================================================

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
function kktwindow(rhist::AbstractVector{T}, chist::AbstractVector{T}, δ, nir::Int, ncg::Int, ftol, gtol) where {T}
    return kktwindow(rhist, chist, convert(T, δ), nir, ncg, convert(T, ftol), convert(T, gtol))
end

function kktwindow(rhist::AbstractVector{T}, chist::AbstractVector{T}, δ::T, nir::Int, ncg::Int, ftol::T, gtol::T) where {T}
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

# ============================================================================
# solvekkt!
# ============================================================================

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
function solvekkt!(
        wrk::UzawaSolver,
        x::AbstractVector,
        y::AbstractVector,
        A::AbstractMatrix,
        B::AbstractMatrix,
        f::AbstractVector,
        g::AbstractVector;
        kw...,
    )
    F = wrk.F
    dw = wrk.divwrk

    function L!(v)
        return ldiv!(dw, F, v)
    end

    function U!(v)
        return ldiv!(dw, F', v)
    end

    return solvekkt!(L!, U!, wrk.itrwrk, wrk.hist, wrk.δf, wrk.δg,
            wrk.δx, wrk.δy, wrk.δ[], x, y, A, B, f, g; kw...)
end

function solvekkt!(
        wrk::PivotedUzawaSolver,
        x::AbstractVector,
        y::AbstractVector,
        A::AbstractMatrix,
        B::AbstractMatrix,
        f::AbstractVector,
        g::AbstractVector;
        kw...,
    )
    F = wrk.F
    dw = wrk.divwrk
    P = wrk.P
    w = wrk.w

    function L!(v)
        mul!(w, P, v)
        lpdiv!(dw, F, w)
        return copyto!(v, w)
    end

    function U!(v)
        lpdiv!(dw, F', v)
        mul!(w, P', v)
        return copyto!(v, w)
    end

    return solvekkt!(L!, U!, wrk.itrwrk, wrk.hist, wrk.δf, wrk.δg,
            wrk.δx, wrk.δy, wrk.δ[], x, y, A, B, f, g; kw...)
end

function solvekkt!(
        L!::Function,
        U!::Function,
        itrwrk::CGWorkspace,
        hist::AbstractVector,
        δf::AbstractVector,
        δg::AbstractVector,
        δx::AbstractVector,
        δy::AbstractVector,
        δ::Real,
        x::AbstractVector,
        y::AbstractVector,
        A::AbstractMatrix,
        B::AbstractMatrix,
        f::AbstractVector,
        g::AbstractVector;
        warm::Bool,
        gtol::Real,
        ftol::Real,
        stall::Real,
        irmax::Int,
        cgmax::Int,
    )
    m, n = size(B)

    @assert length(x) == n
    @assert length(y) == m
    @assert length(f) == n
    @assert length(g) == m
    @assert size(A, 1) == n

    @assert δ > 0
    @assert gtol ≥ 0
    @assert ftol ≥ 0
    @assert irmax ≥ 0
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
        axpby!(1, f, -1, δf)
        axpby!(1, g, -1, δg)
    else
        copyto!(δf, f)
        copyto!(δg, g)
        fill!(x, 0)
        fill!(y, 0)
    end

    fprv = fres = norm(δf); hist[1] = fres

    if fres ≤ ftol
        irstat = KKT_SOLVED
    elseif nir ≥ irmax
        irstat = KKT_IRMAX
    else
        irstat = KKT_CONTINUE
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
    while irstat === KKT_CONTINUE
        nir += 1
        #
        # solve for the corrections δx and δy
        #
        #   [ A -Bᵀ ] [ δx ] = [ δf ]
        #   [ B  δ  ] [ δy ]   [ δg ]
        #
        copyto!(δx, δf)
        rmul!(δx, δ)
        mul!(δx, B', δg, 1, 1)
        L!(δx)
        U!(δx)
        mul!(δg, B, δx, -1, 1)
        #
        # update x and y:
        #
        #   x ← x + δx
        #   y ← y + δy
        #
        axpy!(1, δx, x)
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
        axpby!(1, f, -1, δf)
        axpby!(1, g, -1, δg)
        fres = norm(δf)
        hist[nir + 1] = fres

        if fres ≤ ftol
            irstat = KKT_SOLVED
        elseif fres > (1 - stall) * fprv
            irstat = KKT_STAGNATED
        elseif nir ≥ irmax
            irstat = KKT_IRMAX
        else
            irstat = KKT_CONTINUE
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
    axpy!(-α, δg, y)
    ncg, cgstat = cg!(L!, U!, itrwrk, B, x, δy, δg; atol = gtol, itmax = cgmax)
    axpy!(α, δy, y)
    axpy!(α, δg, y)

    if cgstat === CG_NUMERICAL_FAILURE
        status = KKT_NUMERICAL_FAILURE
    elseif irstat === KKT_STAGNATED
        status = KKT_STAGNATED
    elseif irstat === KKT_IRMAX || cgstat === CG_ITMAX
        status = KKT_ITMAX
    else
        status = KKT_SOLVED
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
    δmin, δmax = kktwindow(hist, itrwrk.hist, δ, nir, ncg, ftol, gtol)

    return ncg, nir, status, δmin, δmax
end
