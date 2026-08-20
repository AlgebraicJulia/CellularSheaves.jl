struct BorderedSolver{UPLO, T, I} <: KKTSolver{T}
    kkt::UzawaSolver{UPLO, T, I}
    d::FVector{T}
    κ::FScalar{T}
    φ::FScalar{T}
end

function BorderedSolver(S::ChordalSymbolic{I}, B::BlockSparseMatrix{T, I}; cgmax::Integer = 2size(B, 2), irmax::Integer = 10) where {T, I <: Integer}
    F = FChordalTriangular{:N, :L, T, I}(S)
    return BorderedSolver(F, B; cgmax, irmax)
end

function BorderedSolver(F::FChordalTriangular{:N, UPLO, T, I}, B::BlockSparseMatrix{T, I}; cgmax::Integer = 2size(B, 2), irmax::Integer = 10) where {UPLO, T, I <: Integer}
    n = size(B, 2)

    kkt = UzawaSolver(F, B; cgmax, irmax)

    d = FVector{T}(undef, n)
    κ = FScalar{T}(undef)
    φ = FScalar{T}(undef)

    return BorderedSolver{UPLO, T, I}(kkt, d, κ, φ)
end

function initkkt!(
        bw::BorderedSolver,
        x2::AbstractVector,
        y2::AbstractVector,
        A::BlockSparseMatrix,
        B::BlockSparseMatrix,
        c::AbstractVector,
        e::AbstractVector,
        d::AbstractVector,
        φ::Real;
        δ::Real,
    )
    bw.φ[] = φ

    flag, witer, wpass, κ = initkkt!(bw.kkt, x2, y2, A, B, c, e, d, φ; δ)
    copyto!(bw.d, d)
    bw.κ[] = κ

    return flag, witer, wpass
end

function initkkt!(
        kkt::UzawaSolver,
        x2::AbstractVector,
        y2::AbstractVector,
        A::BlockSparseMatrix,
        B::BlockSparseMatrix,
        c::AbstractVector,
        e::AbstractVector,
        d::AbstractVector,
        φ::Real;
        δ::Real,
    )
    m, n = size(B)

    @assert length(x2) == n
    @assert length(y2) == m
    @assert length(c)  == n
    @assert length(e)  == m
    @assert length(d)  == n
    @assert size(A, 1) == n
    @assert δ ≥ 0

    witer = 0
    wpass = 0
    #
    # factor the augmented matrix
    #
    #   F Fᵀ = δ A + Bᵀ B
    #
    flag = initkkt!(kkt, A; δ)

    if flag
        #
        # solve the Woodbury system
        #
        #   [ A -Bᵀ ] [ x₂ ] = [ c ]
        #   [ B  0  ] [ y₂ ]   [ e ]
        #
        ncg, nir, _ = solvekkt!(kkt, x2, y2, A, B, c, e;
                                warm = true, gtol = 0, ftol = 0, stall = 0, irmax = 1, cgmax = 0)
        witer = ncg
        wpass = nir
    end
    #
    # compute the capacitance
    #
    #   κ ← dᵀ x₂ + eᵀ y₂ + φ
    #
    κ = dot(d, x2) + dot(e, y2) + φ

    return flag, witer, wpass, κ
end

#
# Solve the bordered KKT system
#
#   [ A    -Bᵀ  -c ] [ x ]   [ f ]
#   [ B     0   -e ] [ y ] = [ g ]
#   [ dᵀ    eᵀ   φ ] [ ζ ]   [ η ]
#
#
# where A is a n x n positive-definite
# matrix, B is a m × n matrix, δ is a
# positive number, and L is the n × n
# Cholesky factor of the augmented matrix
#
#   δ A + Bᵀ B = L Lᵀ.
#
# The solution (x, y, ζ) meets the tolerances
#
#   ‖ A x - Bᵀ y - c ζ - f ‖ ≤ ftol
#   ‖ B x        - e ζ - g ‖ ≤ gtol
#   | dᵀx + eᵀy  + φ ζ - η | ≤ ηtol.
#
function solvekkt!(bw::BorderedSolver, args...; kw...)
    niter, npass, witer, wpass, status, ζ, δmin, δmax, κ =
        solvekkt!(bw.kkt, bw.d, bw.κ[], bw.φ[], args...; kw...)
    bw.κ[] = κ
    return niter, npass, witer, wpass, status, ζ, δmin, δmax
end

function solvekkt!(
        kkt::UzawaSolver,
        d::AbstractVector,
        κ::Real,
        φ::Real,
        x::AbstractVector,
        y::AbstractVector,
        x2::AbstractVector,
        y2::AbstractVector,
        A::AbstractMatrix,
        B::AbstractMatrix,
        c::AbstractVector,
        e::AbstractVector,
        f::AbstractVector,
        g::AbstractVector,
        η::Real;
        warm::Bool,
        gtol::Real,
        ftol::Real,
        ηtol::Real,
        stall::Real,
        irmax::Int,
        cgmax::Int,
    )
    m, n = size(B)

    @assert length(x)  == n
    @assert length(y)  == m
    @assert length(x2) == n
    @assert length(y2) == m
    @assert length(f)  == n
    @assert length(g)  == m
    @assert length(c)  == n
    @assert length(e)  == m
    @assert length(d)  == n
    @assert size(A, 1) == n

    witer = 0
    wpass = 0
    #
    # solve the KKT system
    #
    #   [ A -Bᵀ ] [ x ] = [ f ]
    #   [ B  0  ] [ y ]   [ g ]
    #
    # to the tolerances
    #
    #   ‖ Ax - Bᵀy - f ‖ ≤ 1/2 ftol
    #   ‖ Bx       - g ‖ ≤ 1/2 gtol
    #
    niter, npass, nstat, δmin, δmax = solvekkt!(kkt, x, y, A, B, f, g;
                                warm, ftol = ftol / 2, gtol = gtol / 2, stall, irmax, cgmax)
    #
    # compute the residual
    #
    #   δη = η - dᵀx - eᵀy
    #
    δη = η - dot(d, x) - dot(e, y)
    #
    # compute ζ
    #
    #   ζ = κ⁻¹ δη
    #
    ζ  = δη / κ

    wstat = KKT_SOLVED
    #
    # refine the Woodbury system
    #
    #   [ A -Bᵀ ] [ x₂ ] = [ c ]
    #   [ B  0  ] [ y₂ ]   [ e ]
    #
    # until
    #
    #   ‖ A x₂ - Bᵀy₂ - c ‖ ≤ 1/(2|ζ|) ftol
    #   ‖ B x₂        - e ‖ ≤ 1/(2|ζ|) gtol
    #
    if !iszero(ζ)
        ncg2, nir2, wstat = solvekkt!(kkt, x2, y2, A, B, c, e; warm = true,
                                ftol = ftol / 2abs(ζ), gtol = gtol / 2abs(ζ), stall, irmax, cgmax)
        witer += ncg2
        wpass += nir2
        #
        # compute the capacitance
        #
        #   κ ← dᵀ x₂ + eᵀ y₂ + φ
        #
        κ = dot(d, x2) + dot(e, y2) + φ
        #
        # re-compute ζ
        #
        #   ζ = κ⁻¹ δη
        #
        ζ = δη / κ
    end
    #
    # lift x and y:
    #
    #   x ← x + ζ x₂
    #   y ← y + ζ y₂
    #
    axpy!(ζ, x2, x)
    axpy!(ζ, y2, y)
    #
    # compute the residual norm
    #
    #   ηres = | dᵀ x + eᵀ y + φ ζ - η |
    #
    ηres = abs(η - dot(d, x) - dot(e, y) - φ * ζ)

    if nstat === KKT_NUMERICAL_FAILURE || wstat === KKT_NUMERICAL_FAILURE
        status = KKT_NUMERICAL_FAILURE
    elseif nstat !== KKT_SOLVED
        status = nstat
    elseif wstat !== KKT_SOLVED
        status = wstat
    elseif ηres > ηtol
        status = KKT_STAGNATED
    else
        status = KKT_SOLVED
    end

    return niter, npass, witer, wpass, status, ζ, δmin, δmax, κ
end
