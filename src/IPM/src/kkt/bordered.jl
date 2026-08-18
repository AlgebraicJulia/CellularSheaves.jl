struct BorderedSolver{T, UPLO, I} <: KKTSolver{T}
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
    d::FVector{T}
    κ::FScalar{T}
    φ::FScalar{T}
end

function BorderedSolver(S::ChordalSymbolic{I}, B::BlockSparseMatrix{T, I}) where {T, I <: Integer}
    F = FChordalTriangular{:N, :L, T, I}(S)
    return BorderedSolver(F, B)
end

function BorderedSolver(F::FChordalTriangular{:N, UPLO, T, I}, B::BlockSparseMatrix{T, I}) where {UPLO, T, I <: Integer}
    m, n = size(B)

    L = B' * B

    facwrk = FactorizationWorkspace(F)
    divwrk = DivisionWorkspace(F, 1)
    itrwrk = CGWorkspace{T}(m, n)

    r = FVector{T}(undef, m)
    d = FVector{T}(undef, n)

    δf = FVector{T}(undef, n)
    δg = FVector{T}(undef, m)
    δx = FVector{T}(undef, n)
    δy = FVector{T}(undef, m)

    δ = FScalar{T}(undef)
    κ = FScalar{T}(undef)
    φ = FScalar{T}(undef)

    return BorderedSolver{T, UPLO, I}(F, L, facwrk, divwrk, itrwrk, r, δ,
                                      δf, δg, δx, δy,
                                      d, κ, φ)
end

function initkkt!(
        bw::BorderedSolver{T},
        x2::AbstractVector{T},
        y2::AbstractVector{T},
        r2::AbstractVector{T},
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

    ok, ρ, nwood, κ = initkkt!(bw.facwrk, bw.F, bw.L, bw.divwrk, bw.δf, bw.δx,
                               x2, y2, r2, A, B, c, e, d, φ; δ, rgmin)
    copyto!(bw.d, d)
    bw.κ[] = κ

    return ok, ρ, nwood
end

function initkkt!(
        facwrk::FactorizationWorkspace{T},
        F::AbstractMatrix{T},
        L::BlockSparseMatrix{T},
        divwrk::DivisionWorkspace{T},
        δf::AbstractVector{T},
        δx::AbstractVector{T},
        x2::AbstractVector{T},
        y2::AbstractVector{T},
        r2::AbstractVector{T},
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
    @assert length(r2) == m
    @assert length(c)  == n
    @assert length(e)  == m
    @assert length(d)  == n
    @assert size(A, 1) == n
    @assert size(L, 1) == n
    @assert δ     ≥ 0
    @assert rgmin ≥ 0

    α = inv(δ)

    nwood = 0
    κ = zero(T)
    #
    # factor the augmented matrix
    #
    #   F Fᵀ = δ A + Bᵀ B + ρ I
    #
    ok, ρ = inituzw!(facwrk, F, L, A, δ, rgmin)

    if ok
        #
        # compute the residual
        #
        #   [ δf ] = [ c ] - [ A  -Bᵀ ] [ x₂ ]
        #   [ r₂ ]   [ e ]   [ B   0  ] [ y₂ ]
        #
        mulkkt!(δf, r2, A, B, x2, y2)
        axpby!(one(T), c, -one(T), δf)
        axpby!(one(T), e, -one(T), r2)
        #
        # solve for the direction δx
        #
        #   (δ A + Bᵀ B) δx = δ δf + Bᵀ r₂
        #
        # and update x₂ and r₂:
        #
        #   x₂ ← x₂ +   δx
        #   r₂ ← r₂ - B δx
        #
        copyto!(δx, δf)
        rmul!(δx, δ)
        mul!(δx, B', r2, one(T), one(T))
        ldiv!(divwrk, F,  δx)
        ldiv!(divwrk, F', δx)

        axpy!(one(T), δx, x2)
        mul!(r2, B, δx, -one(T), one(T))
        #
        # compute the capacitance
        #
        #   κ ← dᵀx₂ + eᵀy₂ + α eᵀr₂ + φ
        #
        κ = dot(d, x2) + dot(e, y2) + α * dot(e, r2) + φ
        nwood = 1
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
# The solution (x, y, ζ) meets the tolerances
#
#   ‖ A x - Bᵀ y - c ζ - f ‖ ≤ ftol
#   ‖ B x        - e ζ - g ‖ ≤ gtol
#   | dᵀx + eᵀy  + φ ζ - η | ≤ ηtol
#
function solvekkt!(bw::BorderedSolver, args...; kw...)
    niter, nwood, npass, status, ζ, δmin, δmax, κ =
        solvekkt!(bw.divwrk, bw.itrwrk, bw.r, bw.F, bw.δf, bw.δg, bw.δx, bw.δy,
                  bw.d, bw.κ[], bw.φ[], bw.δ[],
                  args...; kw...)
    bw.κ[] = κ
    return niter, nwood, npass, status, ζ, δmin, δmax
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
        d::AbstractVector{T},
        κ::T,
        φ::T,
        δ::T,
        x::AbstractVector{T},
        y::AbstractVector{T},
        x2::AbstractVector{T},
        y2::AbstractVector{T},
        r2::AbstractVector{T},
        A::BlockSparseMatrix{T},
        B::BlockSparseMatrix{T},
        c::AbstractVector{T},
        e::AbstractVector{T},
        f::AbstractVector{T},
        g::AbstractVector{T},
        η::T;
        xwarm::Bool = false,
        ywarm::Bool = false,
        gtol::T = √eps(T),
        ftol::T = √eps(T),
        ηtol::T = √eps(T),
        stall::T = T(0.5),
        rfmax::Int = 10,
        cgmax::Int = 100,
    ) where {T}
    m, n = size(B)

    α = inv(δ)

    @assert length(x)  == n
    @assert length(y)  == m
    @assert length(x2) == n
    @assert length(y2) == m
    @assert length(f)  == n
    @assert length(g)  == m
    @assert length(c)  == n
    @assert length(e)  == m
    @assert length(r2) == m
    @assert length(d)  == n
    @assert size(A, 1) == n
    @assert size(L, 1) == n

    niter = 0
    nwood = 0
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
        mul!(δg, B,                x, -one(T), one(T))
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
    n, cgstat = solveuzw!(divwrk, itrwrk, δx, δy, r, L, A, B, δf, δg, δ; atol = gtol / 2, itmax = cgmax)
    niter += n

    if cgstat === CG_NUMERICAL_FAILURE
        status = KKT_NUMERICAL_FAILURE
        δmin = δmax = δ
        ζ = zero(T)
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
        #   [ δf ]   [ f ]   [ A   -Bᵀ ] [ x ]
        #   [ δg ] = [ g ] - [ B    0  ] [ y ]
        #   [ δη ]   [ η ]   [ dᵀ   eᵀ ]
        #
        mulkkt!(δf, δg, A, B, x, y)
        axpby!(one(T), f, -one(T), δf)
        axpby!(one(T), g, -one(T), δg)
        δη = η - dot(d, x) - dot(e, y)
        #
        # compute the residual norms
        #
        #   gres = ‖ δg ‖
        #   fres = ‖ δf ‖
        #
        gres = norm(δg)
        fres = norm(δf)
        #
        # estimate the minimum-cost interval
        #
        #   [δmin, δmax];
        #
        # any parameter δ in this interval would have
        # solved (1) to the tolerances
        #
        #   ‖ A δx + Bᵀ δy - δf ‖ ≤ ftol
        #   ‖ B δx         - δg ‖ ≤ gtol
        #
        # with the minimum possible number of triangular solves.
        # beware! this estimate can be very poor.
        #
        δmin, δmax = kktwindow(itrwrk.hist, δ, fres, niter; gtol = gtol / 2, ftol)
        #
        # compute ζ:
        #
        #   ζ = κ⁻¹ δη
        #
        ζ = δη / κ
        #
        # ensure that
        #
        #   ‖ e - B x₂ ‖ ≤ (gtol - gres) / | ζ |
        #
        # by solving for δx₂ and δy₂
        #
        #   [ δ A + Bᵀ B  -Bᵀ ] [ δx₂ ] = [ 0  ]
        #   [          B   0  ] [ δy₂ ]   [ r₂ ]
        #
        # and updating
        #
        #   x₂ ← x₂ +   δx₂
        #   y₂ ← y₂ + α δy₂
        #   r₂ ← r₂ - B δx₂
        #
        nwood, wstat = cg!(itrwrk, B, L, divwrk, x2, δy, r2; atol = (gtol - gres) / abs(ζ), itmax = cgmax)
        axpy!(α, δy, y2)
        #
        # re-compute the capacitance
        #
        #   κ ← dᵀx₂ + eᵀy₂ + α eᵀr₂ + φ
        #
        κ = dot(d, x2) + dot(e, y2) + α * dot(e, r2) + φ
        #
        # re-compute ζ:
        #
        #   ζ = κ⁻¹ δη
        #
        ζ = δη / κ
        #
        # lift x and y
        #
        #   x ← x + ζ x2
        #   y ← y + ζ y2 + α ζ r2
        #
        axpy!(ζ,     x2, x)
        axpy!(ζ,     y2, y)
        axpy!(ζ * α, r2, y)
        #
        # compute the residuals
        #
        #   [ δf ]   [ f ]   [ A   -Bᵀ  -c ] [ x ]
        #   [ δg ] = [ g ] - [ B    0   -e ] [ y ]
        #   [ δη ]   [ η ]   [ dᵀ   eᵀ   φ ] [ ζ ]
        #
        mulkkt!(δf, δg, A, B, x, y)
        axpy!(-ζ, c, δf)
        axpy!(-ζ, e, δg)
        axpby!(one(T), f, -one(T), δf)
        axpby!(one(T), g, -one(T), δg)
        δη = η - dot(d, x) - dot(e, y) - φ * ζ
        #
        # compute the residual norms
        #
        #   gres = ‖ δg ‖
        #   fres = ‖ δf ‖
        #   ηres = | δη |
        #
        gres = norm(δg)
        fres = norm(δf)
        ηres = abs(δη)

        if wstat === CG_NUMERICAL_FAILURE
            status = KKT_NUMERICAL_FAILURE
        elseif gres ≤ gtol && fres ≤ ftol && ηres ≤ ηtol
            status = KKT_SOLVED
        elseif npass == rfmax
            status = KKT_ITMAX
        else
            status = KKT_CONTINUE
        end

        gprv = gres
        fprv = fres
        ηprv = ηres

        while status == KKT_CONTINUE
            #
            # solve for the corrections δx and δy
            #
            #   [ A  -Bᵀ ] [ δx ]   [ δf ]
            #   [ B   0  ] [ δy ] = [ δg ]
            #
            n, cgstat = solveuzw!(divwrk, itrwrk, δx, δy, r, L, A, B, δf, δg, δ; atol = gtol, rtol = T(INTERIOR_THETA), itmax = cgmax)
            niter += n
            npass += 1

            if cgstat === CG_NUMERICAL_FAILURE
                status = KKT_NUMERICAL_FAILURE
            else
                #
                # compute the correction δζ:
                #
                #   δζ = (δη - dᵀδx - eᵀδy) / κ
                #
                δζ = (δη - dot(d, δx) - dot(e, δy)) / κ
                #
                # lift x, y, and ζ:
                #
                #   x ← x + δx + δζ x₂
                #   y ← y + δy + δζ y₂ + α δζ r₂
                #   ζ ← ζ + δζ
                #
                axpy!(δζ, x2, δx)
                axpy!(δζ, y2, δy)
                axpy!(δζ * α, r2, δy)
                axpy!(one(T), δx, x)
                axpy!(one(T), δy, y)
                ζ += δζ
                #
                # compute the residuals
                #
                #   [ δf ]   [ f ]   [ A   -Bᵀ  -c ] [ x ]
                #   [ δg ] = [ g ] - [ B    0   -e ] [ y ]
                #   [ δη ]   [ η ]   [ dᵀ   eᵀ   φ ] [ ζ ]
                #
                mulkkt!(δf, δg, A, B, x, y)
                axpy!(-ζ, c, δf)
                axpy!(-ζ, e, δg)
                axpby!(one(T), f, -one(T), δf)
                axpby!(one(T), g, -one(T), δg)
                δη = η - dot(d, x) - dot(e, y) - φ * ζ
                #
                # compute the residual norms
                #
                #   gres = ‖ δg ‖
                #   fres = ‖ δf ‖
                #   ηres = | δη |
                #
                gres = norm(δg)
                fres = norm(δf)
                ηres = abs(δη)

                if gres ≤ gtol && fres ≤ ftol && ηres ≤ ηtol
                    status = KKT_SOLVED
                elseif gres > (1 - stall) * gprv && fres > (1 - stall) * fprv && ηres > (1 - stall) * ηprv
                    status = KKT_STAGNATED
                elseif npass == rfmax
                    status = KKT_ITMAX
                else
                    status = KKT_CONTINUE
                end

                gprv = gres
                fprv = fres
                ηprv = ηres
            end
        end
    end

    return niter, nwood, npass, status, ζ, δmin, δmax, κ
end
