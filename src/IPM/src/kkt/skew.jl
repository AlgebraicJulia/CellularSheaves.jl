struct SkewSolver{T, KKT}
    kkt::KKT
    d::FVector{T}
    σ::FScalar{T}
    γ::FScalar{T}
end

function SkewSolver(S::ChordalSymbolic{I}, B::BlockSparseMatrix{T, I}; cgmax::Integer = 2size(B, 2), irmax::Integer = 10, pivot::Bool = false) where {T, I <: Integer}
    n = size(B, 2)

    if pivot
        kkt = PivotedUzawaSolver(S, B; cgmax, irmax)
    else
        kkt = UzawaSolver(S, B; cgmax, irmax)
    end

    d = FVector{T}(undef, n)
    σ = FScalar{T}(undef)
    γ = FScalar{T}(undef)

    return SkewSolver(kkt, d, σ, γ)
end

function initkkt!(
        bw::SkewSolver,
        x2::AbstractVector,
        y2::AbstractVector,
        A::BlockSparseMatrix,
        B::BlockSparseMatrix,
        c::AbstractVector,
        e::AbstractVector,
        d::AbstractVector,
        γ::Real;
        δ::Real,
    )
    bw.γ[] = γ

    flag, witer, wpass, σ = initkkt!(bw.kkt, x2, y2, A, B, c, e, d, γ; δ)
    copyto!(bw.d, d)
    bw.σ[] = σ

    return flag, witer, wpass
end

function initkkt!(
        kkt::KKTSolver,
        x2::AbstractVector,
        y2::AbstractVector,
        A::BlockSparseMatrix,
        B::BlockSparseMatrix,
        c::AbstractVector,
        e::AbstractVector,
        d::AbstractVector,
        γ::Real;
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
        # solve the τ column
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
    #   σ ← dᵀ x₂ + eᵀ y₂ + γ
    #
    σ = dot(d, x2) + dot(e, y2) + γ

    return flag, witer, wpass, σ
end

#
# Solve the bordered KKT system
#
#   [ A    -Bᵀ  -c ] [ x ]   [ h ]
#   [ B     0   -e ] [ y ] = [ i ]
#   [ dᵀ    eᵀ   γ ] [ ζ ]   [ κ ]
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
#   ‖ A x - Bᵀ y - c ζ - h ‖ ≤ htol
#   ‖ B x        - e ζ - i ‖ ≤ itol
#
function solvekkt!(bw::SkewSolver, args...; kw...)
    niter, npass, witer, wpass, status, ζ, δmin, δmax, σ =
        solvekkt!(bw.kkt, bw.d, bw.σ[], bw.γ[], args...; kw...)
    bw.σ[] = σ
    return niter, npass, witer, wpass, status, ζ, δmin, δmax
end

function solvekkt!(
        kkt::KKTSolver,
        d::AbstractVector,
        σ::Real,
        γ::Real,
        x::AbstractVector,
        y::AbstractVector,
        x2::AbstractVector,
        y2::AbstractVector,
        A::AbstractMatrix,
        B::AbstractMatrix,
        c::AbstractVector,
        e::AbstractVector,
        h::AbstractVector,
        i::AbstractVector,
        κ::Real;
        warm::Bool,
        htol::Real,
        itol::Real,
        stall::Real,
        irmax::Int,
        cgmax::Int,
    )
    m, n = size(B)

    @assert length(x)  == n
    @assert length(y)  == m
    @assert length(x2) == n
    @assert length(y2) == m
    @assert length(h)  == n
    @assert length(i)  == m
    @assert length(c)  == n
    @assert length(e)  == m
    @assert length(d)  == n
    @assert size(A, 1) == n

    witer = 0
    wpass = 0
    #
    # solve the KKT system
    #
    #   [ A -Bᵀ ] [ x ] = [ h ]
    #   [ B  0  ] [ y ]   [ i ]
    #
    # to the tolerances
    #
    #   ‖ Ax - Bᵀy - h ‖ ≤ 1/2 htol
    #   ‖ Bx       - i ‖ ≤ 1/2 itol
    #
    niter, npass, nstat, δmin, δmax = solvekkt!(kkt, x, y, A, B, h, i;
                                warm, ftol = htol / 2, gtol = itol / 2, stall, irmax, cgmax)
    #
    # compute the residual
    #
    #   δκ = κ - dᵀx - eᵀy
    #
    δκ = κ - dot(d, x) - dot(e, y)
    #
    # compute ζ
    #
    #   ζ = σ⁻¹ δκ
    #
    ζ  = δκ / σ

    wstat = KKT_SOLVED
    #
    # refine the τ column
    #
    #   [ A -Bᵀ ] [ x₂ ] = [ c ]
    #   [ B  0  ] [ y₂ ]   [ e ]
    #
    # until
    #
    #   ‖ A x₂ - Bᵀy₂ - c ‖ ≤ 1/(2|ζ|) htol
    #   ‖ B x₂        - e ‖ ≤ 1/(2|ζ|) itol
    #
    if !iszero(ζ)
        ncg2, nir2, wstat = solvekkt!(kkt, x2, y2, A, B, c, e; warm = true,
                                ftol = htol / 2abs(ζ), gtol = itol / 2abs(ζ), stall, irmax, cgmax)
        witer += ncg2
        wpass += nir2
        #
        # compute the capacitance
        #
        #   σ ← dᵀ x₂ + eᵀ y₂ + γ
        #
        σ = dot(d, x2) + dot(e, y2) + γ
        #
        # re-compute ζ
        #
        #   ζ = σ⁻¹ δκ
        #
        ζ = δκ / σ
    end
    #
    # lift x and y:
    #
    #   x ← x + ζ x₂
    #   y ← y + ζ y₂
    #
    axpy!(ζ, x2, x)
    axpy!(ζ, y2, y)

    if nstat === KKT_NUMERICAL_FAILURE || wstat === KKT_NUMERICAL_FAILURE
        status = KKT_NUMERICAL_FAILURE
    elseif nstat !== KKT_SOLVED
        status = nstat
    elseif wstat !== KKT_SOLVED
        status = wstat
    else
        status = KKT_SOLVED
    end

    return niter, npass, witer, wpass, status, ζ, δmin, δmax, σ
end
