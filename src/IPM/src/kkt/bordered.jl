struct BorderedSolver{T, UPLO, I} <: KKTSolver{T}
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
    d::FVector{T}
    κ::FScalar{T}
    φ::FScalar{T}
end

function BorderedSolver(S::ChordalSymbolic{I}, B::BlockSparseMatrix{T, I}; cgmax::Integer = 2size(B, 2), rfmax::Integer = 10) where {T, I <: Integer}
    F = FChordalTriangular{:N, :L, T, I}(S)
    return BorderedSolver(F, B; cgmax, rfmax)
end

function BorderedSolver(F::FChordalTriangular{:N, UPLO, T, I}, B::BlockSparseMatrix{T, I}; cgmax::Integer = 2size(B, 2), rfmax::Integer = 10) where {UPLO, T, I <: Integer}
    m, n = size(B)

    L = B' * B

    facwrk = FactorizationWorkspace(F)
    divwrk = DivisionWorkspace(F, 1)
    itrwrk = CGWorkspace{T}(m, n; itmax = cgmax)

    hist = FVector{T}(undef, rfmax + 2)
    d = FVector{T}(undef, n)

    δf = FVector{T}(undef, n)
    δg = FVector{T}(undef, m)
    δx = FVector{T}(undef, n)
    δy = FVector{T}(undef, m)

    δ = FScalar{T}(undef)
    κ = FScalar{T}(undef)
    φ = FScalar{T}(undef)

    return BorderedSolver{T, UPLO, I}(F, L, facwrk, divwrk, itrwrk, hist, δ,
                                      δf, δg, δx, δy,
                                      d, κ, φ)
end

function initkkt!(
        bw::BorderedSolver{T},
        x2::AbstractVector{T},
        y2::AbstractVector{T},
        A::BlockSparseMatrix{T},
        B::BlockSparseMatrix{T},
        c::AbstractVector{T},
        e::AbstractVector{T},
        d::AbstractVector{T},
        φ::T;
        δ::T,
        rgmin::T,
    ) where {T}
    bw.δ[] = δ
    bw.φ[] = φ

    ok, ρ, nwood, κ = initkkt!(bw.divwrk, bw.itrwrk, bw.hist, bw.F, bw.L, bw.facwrk,
                               bw.δf, bw.δg, bw.δx, bw.δy,
                               x2, y2, A, B, c, e, d, φ; δ, rgmin)
    copyto!(bw.d, d)
    bw.κ[] = κ

    return ok, ρ, nwood
end

function initkkt!(
        divwrk::DivisionWorkspace{T},
        itrwrk::CGWorkspace{T},
        hist::AbstractVector{T},
        F::AbstractMatrix{T},
        L::BlockSparseMatrix{T},
        facwrk::FactorizationWorkspace{T},
        δf::AbstractVector{T},
        δg::AbstractVector{T},
        δx::AbstractVector{T},
        δy::AbstractVector{T},
        x2::AbstractVector{T},
        y2::AbstractVector{T},
        A::BlockSparseMatrix{T},
        B::BlockSparseMatrix{T},
        c::AbstractVector{T},
        e::AbstractVector{T},
        d::AbstractVector{T},
        φ::T;
        δ::T,
        rgmin::T,
    ) where {T}
    m, n = size(B)

    @assert length(x2) == n
    @assert length(y2) == m
    @assert length(c)  == n
    @assert length(e)  == m
    @assert length(d)  == n
    @assert size(A, 1) == n
    @assert size(L, 1) == n
    @assert δ     ≥ 0
    @assert rgmin ≥ 0

    nwood = 0
    κ = zero(T)
    #
    # factor the augmented matrix
    #
    #   F Fᵀ = δ A + Bᵀ B + ρ I
    #
    ok, ρ = initkkt!(facwrk, F, L, A, δ, rgmin)

    if ok
        #
        # solve the Woodbury system
        #
        #   [ A -Bᵀ ] [ x₂ ] = [ c ]
        #   [ B  0  ] [ y₂ ]   [ e ]
        #
        ncg, nir, _ = solvekkt!(divwrk, itrwrk, hist, F, δf, δg, δx, δy, δ,
                                x2, y2, A, B, c, e; warm = true, rfmax = 1, cgmax = 0)
        #
        # compute the capacitance
        #
        #   κ ← dᵀ x₂ + eᵀ y₂ + φ
        #
        κ = dot(d, x2) + dot(e, y2) + φ
        nwood = nir + ncg
    end

    return ok, ρ, nwood, κ
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
    niter, nwood, npass, status, ζ, δmin, δmax, κ =
        solvekkt!(bw.divwrk, bw.itrwrk, bw.hist, bw.F, bw.δf, bw.δg, bw.δx, bw.δy,
                  bw.d, bw.κ[], bw.φ[], bw.δ[],
                  args...; kw...)
    bw.κ[] = κ
    return niter, nwood, npass, status, ζ, δmin, δmax
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
        d::AbstractVector{T},
        κ::T,
        φ::T,
        δ::T,
        x::AbstractVector{T},
        y::AbstractVector{T},
        x2::AbstractVector{T},
        y2::AbstractVector{T},
        A::BlockSparseMatrix{T},
        B::BlockSparseMatrix{T},
        c::AbstractVector{T},
        e::AbstractVector{T},
        f::AbstractVector{T},
        g::AbstractVector{T},
        η::T;
        warm::Bool = false,
        gtol::T = √eps(T),
        ftol::T = √eps(T),
        ηtol::T = √eps(T),
        stall::T = T(0.5),
        rfmax::Int = 10,
        cgmax::Int = 100,
    ) where {T}
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
    @assert size(L, 1) == n

    nwood = 0; ζ = T(NaN)
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
    niter, npass, status, δmin, δmax = solvekkt!(divwrk, itrwrk, hist, L, δf, δg, δx, δy, δ,
                                x, y, A, B, f, g; warm, ftol = ftol / 2, gtol = gtol / 2, stall, rfmax, cgmax)

    if status === KKT_SOLVED
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
            ncg2, nir2, status = solvekkt!(divwrk, itrwrk, hist, L, δf, δg, δx, δy, δ,
                                    x2, y2, A, B, c, e; warm = true,
                                    ftol = ftol / 2abs(ζ), gtol = gtol / 2abs(ζ), stall, rfmax, cgmax)
            nwood += nir2 + ncg2
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

        if ηres > ηtol
            status = KKT_STAGNATED
        end
    end

    return niter, nwood, npass, status, ζ, δmin, δmax, κ
end
