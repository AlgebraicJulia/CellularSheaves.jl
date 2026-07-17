struct HSDSolver{T, I, W, C} <: AbstractSolver{T}
    Q::BlockSparseMatrix{T, I}
    H::BlockSparseMatrix{T, I}
    Hc::BlockSparseMatrix{T, I}  # cone-only Hessian (before Q added)
    B::BlockSparseMatrix{T, I}
    c::FVector{T}
    g::FVector{T}
    p::FVector{T}
    d::FVector{T}
    y::FVector{T}
    K::FVector{C}
    scaling::Scaling{T}
    P::FPermutation{I}
    wrk::HSDWorkspace{T}
    caches::Caches{T, I}
    conewrk::ConeWorkspace{T}
    kkt::W
    hist::HSDHistory{T}
    ν::Int
    settings::HSDSettings{T}
    ρ::FScalar{T}
    τ::FScalar{T}
    κ::FScalar{T}
    Δy0::FVector{T}
    nc::FScalar{T}
    ng::FScalar{T}
    nB::FScalar{T}
    timers::TimerOutput
end

function result(s::HSDSolver{T}, status::IPMStatus) where {T}
    p = Vector{T}(undef, length(s.p))
    d = Vector{T}(undef, length(s.d))
    y = Vector{T}(undef, length(s.y))

    τ = s.τ[]
    κ = s.κ[]

    ldiv!(p, s.P, s.p)
    ldiv!(d, s.P, s.d)
    copyto!(y, s.y)

    if status in (OPTIMAL, NEAR_OPTIMAL, STALLED, ITERATION_LIMIT, NUMERICAL_FAILURE)
        ldiv!(τ, p)
        ldiv!(τ, d)
        ldiv!(τ, y)
    elseif status == PRIMAL_INFEASIBLE
        ldiv!(norm(y), y)
    elseif status == DUAL_INFEASIBLE
        np = norm(p)
        ldiv!(np, p)
        ldiv!(np, d)
    end

    unscale!(p, d, y, s.scaling)

    niter = 0
    npred = 0
    ncorr = 0

    for row in s.hist
        niter += 1
        npred += row.npred
        ncorr += row.ncorr
    end

    return HSDResult{T}(p, d, y, status, niter, npred, ncorr, τ, κ, s.hist, s.timers)
end

############################################################################################
# mulhsd!
############################################################################################

#
# compute the matrix-vector product
#
#   [ u ]   [ H   -Bᵀ  -c ] [ x ]
#   [ v ] = [ B    0   -g ] [ y ]
#   [ w ]   [ cᵀ   gᵀ   0 ] [ z ]
#
function mulhsd!(
        u::AbstractVector{T},
        v::AbstractVector{T},
        H::BlockSparseMatrix{T},
        B::BlockSparseMatrix{T},
        c::AbstractVector{T},
        g::AbstractVector{T},
        x::AbstractVector{T},
        y::AbstractVector{T},
        z::T,
    ) where {T}
    #
    # compute the matrix-vector product
    #
    #   [ u ] ← [ H  -Bᵀ ] [ x ]
    #   [ v ]   [ B   0  ] [ y ]
    #
    mulkkt!(u, v, H, B, x, y)
    #
    # subtract cz and gz:
    #
    #   u ← u - cz
    #   v ← v - gz
    #
    axpy!(-z, c, u)
    axpy!(-z, g, v)
    #
    # compute w:
    #
    #   w ← cᵀx + gᵀy
    #
    return dot(c, x) + dot(g, y)
end

############################################################################################
# residuals!
############################################################################################

#
# compute the negated residuals
#
#   [ rd ]   [ d  ]   [ Q   -Bᵀ  -c ] [ p ]
#   [ rp ] = [ 0  ] - [ B    0   -g ] [ y ]
#   [ rτ ]   [ κ  ]   [ cᵀ   gᵀ   0 ] [ τ ]
#
# then correct rτ:
#
#   rτ ← rτ + pᵀQp/τ
#
function residuals!(
        rd::AbstractVector{T},
        rp::AbstractVector{T},
        Qp::AbstractVector{T},
        τ::T,
        κ::T,
        B::BlockSparseMatrix,
        p::AbstractVector,
        d::AbstractVector,
        y::AbstractVector,
        c::AbstractVector,
        g::AbstractVector,
        Q::BlockSparseMatrix,
    ) where {T}
    #
    # compute the matrix-vector product
    #
    #   [ rd ]   [ Q   -Bᵀ  -c ] [ p ]
    #   [ rp ] = [ B    0   -g ] [ y ]
    #   [ rτ ]   [ cᵀ   gᵀ   0 ] [ τ ]
    #
    rτ = mulhsd!(rd, rp, Q, B, c, g, p, y, τ)
    #
    # correct rd and rp
    #
    #   rd ← -rd + d
    #   rp ← -rp  
    #
    axpby!(one(T), d, -one(T), rd)
    lmul!(-one(T), rp)
    #
    # correct rτ:
    #
    #   rτmκ ← pᵀQp/τ - rτ
    #
    return dot(p, Qp) / τ - rτ
end

function residuals!(s::HSDSolver)
    w = s.wrk
    mul!(w.Qp, Symmetric(s.Q, :L), s.p)
    return residuals!(w.rd, w.rp, w.Qp, s.τ[], s.κ[], s.B, s.p, s.d, s.y, s.c, s.g, s.Q)
end

############################################################################################
# woodbury! / capacitance! / newton!
############################################################################################

# solve for the Woodbury auxiliary directions
#
#   [ H  -Bᵀ ] [ Δp2 ]   [ c ]
#   [ B   0  ] [ Δy2 ] = [ g ]
#
function woodbury!(s::HSDSolver{T}; force_tol::T, floor_tol::T, y0 = nothing) where {T}
    atol = max(force_tol, floor_tol)
    return solve_kkt!(s.kkt, s.settings.kkt, s.wrk.Δp2, s.wrk.Δy2, s.H, s.B, s.c, s.g, y0; atol)
end

function capacitance!(
        QΔp2::AbstractVector{T},
        aτ::AbstractVector{T},
        τ::T,
        κ::T,
        Δp2::AbstractVector{T},
        c::AbstractVector{T},
        Qp::AbstractVector{T},
        p::AbstractVector{T},
        W::BlockSparseMatrix{T},
        Q::BlockSparseMatrix{T},
    ) where {T}
    #
    # compute aτ = c - 2Qp/τ
    #
    copyto!(aτ, c)
    axpby!(-2 / τ, Qp, one(T), aτ)
    #
    # compute Q Δp2
    #
    mul!(QΔp2, Symmetric(Q, :L), Δp2)
    #
    # compute the Woodbury capacitance scalar
    #
    #   S = Δp2ᵀ W Δp2 + (Δp2 - p/τ)ᵀ Q (Δp2 - p/τ) + κ/τ
    #
    S = dot(Δp2, Symmetric(W, :L), Δp2) + κ / τ

    @inbounds for i in eachindex(Δp2)
        S += (Δp2[i] - p[i] / τ) * (QΔp2[i] - Qp[i] / τ)
    end

    return S
end

function capacitance!(s::HSDSolver)
    w = s.wrk
    return capacitance!(w.QΔp2, w.aτ, s.τ[], s.κ[], w.Δp2, s.c, w.Qp, s.p, s.Hc, s.Q)
end

############################################################################################
# refinehsd! — R2: Govaerts–Pryce BE+1 refinement on the 3-row bordered system
############################################################################################

function refinehsd!(
    Δp::AbstractVector{T},
    Δy::AbstractVector{T},
    Δτ::T,
    Δd::AbstractVector{T},
    wrk::KKTWorkspace{T},
    set::KKTSettings{T},
    H::BlockSparseMatrix{T},
    B::BlockSparseMatrix{T},
    Q::BlockSparseMatrix{T},
    c::AbstractVector{T},
    g::AbstractVector{T},
    Qp::AbstractVector{T},
    p::AbstractVector{T},
    rd::AbstractVector{T},
    τ::T,
    κ::T,
    rp::AbstractVector{T},  # RHS for row P
    f::AbstractVector{T},   # RHS for row D
    fτ::T,                  # RHS for row T
    Δp2::AbstractVector{T}, # Woodbury direction
    Δy2::AbstractVector{T},
    S::T,                   # Woodbury capacitance scalar
    sp::AbstractVector{T},  # scratch for primal residual
    sd::AbstractVector{T},  # scratch for dual residual
    dp::AbstractVector{T},  # scratch for correction
    dy::AbstractVector{T},
    nc::T,
    ng::T;
    itmax::Int,
    force_tol::T,
    floor_tol::T,
    stall::T,
    nH::T,
    nB::T,
) where {T}
    niter = 0
    status = REACHED_FORCE
    prv = typemax(T)

    nf = norm(f)
    nrp = norm(rp)
    nQp = norm(Qp)
    S0 = dot(p, Qp) / τ^2 + κ / τ

    for i in 1:itmax
        #
        # compute the residuals
        #
        #   [ sd ]   [ f  ]   [ H  -Bᵀ  -c ] [ Δp ]
        #   [ sp ] = [ rp ] - [ B   0   -g ] [ Δy ]
        #   [ sτ ]   [ fτ ]   [ cᵀ  gᵀ   0 ] [ Δτ ]
        #
        sτ = fτ - mulhsd!(sd, sp, H, B, c, g, Δp, Δy, Δτ)
        axpby!(one(T), f,  -one(T), sd)
        axpby!(one(T), rp, -one(T), sp)
        #
        # correct sτ:
        #
        #   sτ ← sτ + 2pᵀQΔp/τ - (pᵀQp/τ² + κ/τ) Δτ
        #
        sτ += (2dot(Qp, Δp) - Δτ * (dot(p, Qp) / τ + κ)) / τ

        pres = norm(sp, Inf) / (1 + ng)
        dres = norm(sd, Inf) / (1 + nc)
        τres = abs(sτ) / (1 + nc + ng)

        res = max(pres, dres, τres)

        if res ≤ force_tol
            status = REACHED_FORCE
            break
        end
        #
        # compute backward-error threshold:
        #
        #   100ε max(σd/(1+‖c‖), σp/(1+‖g‖), στ/(1+‖c‖+‖g‖))
        #
        # where
        #
        #   σd = ‖H‖ ‖Δp‖  + ‖B‖‖Δy‖  + ‖c‖ |Δτ| + ‖f‖
        #   σp = ‖B‖ ‖Δp‖             + ‖g‖ |Δτ| + ‖rp‖
        #   στ = ‖aτ ‖‖Δp‖ + ‖g‖ ‖Δy‖ + S₀  |Δτ| + |fτ|
        #
        np, ny, aΔτ = norm(Δp), norm(Δy), abs(Δτ)
        naτ = nc + 2nQp / τ

        σd =  nH * np + nB * ny + nc * aΔτ + nf
        σp =  nB * np           + ng * aΔτ + nrp
        στ = naτ * np + ng * ny + S0 * aΔτ + abs(fτ)

        dynam_tol = 100eps(T) * max(σd / (1 + nc), σp / (1 + ng), στ / (1 + nc + ng))

        if res ≤ max(dynam_tol, floor_tol)
            status = REACHED_FLOOR
            break
        end

        if res > stall * prv
            status = REFINE_STALLED
            break
        end

        prv = res
        #
        # solve for dp and dy:
        #
        #   [ H -Bᵀ ] [ dp ] = [ sd ]
        #   [ B  0  ] [ dy ]   [ sp ]
        #
        niter += solve_kkt!(wrk, set, dp, dy, H, B, sd, sp; atol=max(force_tol, dynam_tol, floor_tol))
        #
        # apply the Schur lift (aτ = c - 2Qp/τ)
        #
        #   dτ = (sτ - aτᵀdp - gᵀdy) / S
        #
        dτ = (sτ - dot(c, dp) + 2dot(Qp, dp) / τ - dot(g, dy)) / S
        #
        # apply border correction:
        #
        #   dp ← dp + dτ Δp2
        #   dy ← dy + dτ Δy2
        #
        axpy!(dτ, Δp2, dp)
        axpy!(dτ, Δy2, dy)
        #
        # update directions:
        #
        #   Δp ← Δp + dp
        #   Δy ← Δy + dy
        #   Δτ ← Δτ + dτ
        #
        axpy!(one(T), dp, Δp)
        axpy!(one(T), dy, Δy)
        Δτ += dτ
    end
    #
    # re-compute Δd:
    #
    #   Δd ← QΔp - cΔτ - BᵀΔy - rd
    #
    mul!(Δd, Symmetric(Q, :L), Δp)
    axpy!(-Δτ, c, Δd)
    mul!(Δd, B', Δy, -one(T), one(T))
    axpy!(-one(T), rd, Δd)

    return niter, status, Δτ
end

#
# solve for the directions Δp, Δy, and Δτ
#
#   [  H          -Bᵀ   -c  ] [ Δp ]   [ f  ]
#   [  B           0    -g  ] [ Δy ] = [ rp ]
#   [ cᵀ - 2pᵀQ/τ  gᵀ    S₀ ] [ Δτ ]   [ fτ ]
#
# where S₀ = pᵀQp/τ² + κ/τ
#
function newton!(
        Δp::AbstractVector{T},
        Δy::AbstractVector{T},
        Δd::AbstractVector{T},
        τ::T,
        κ::T,
        wrk::KKTWorkspace{T},
        set::KKTSettings{T},
        H::BlockSparseMatrix{T},
        B::BlockSparseMatrix{T},
        Q::BlockSparseMatrix{T},
        c::AbstractVector{T},
        g::AbstractVector{T},
        f::AbstractVector{T},
        rp::AbstractVector{T},
        rd::AbstractVector{T},
        fτ::T,
        Δp2::AbstractVector{T},
        Δy2::AbstractVector{T},
        aτ::AbstractVector{T},
        S::T,
        y0 = nothing;
        atol::T = T(0.1),
    ) where {T}
    #
    # solve for Δp and Δy:
    #
    #   [ H -Bᵀ ] [ Δp ] = [ f  ]
    #   [ B  0  ] [ Δy ]   [ rp ]
    #
    niter = solve_kkt!(wrk, set, Δp, Δy, H, B, f, rp, y0; atol)
    #
    # apply the Schur lift:
    #
    #   Δτ = (fτ - aτᵀΔp - gᵀΔy) / S
    #
    Δτ = (fτ - dot(aτ, Δp) - dot(g, Δy)) / S
    #
    # apply Woodbury update:
    #
    #   Δp ← Δp + Δτ Δp2
    #   Δy ← Δy + Δτ Δy2
    #
    axpy!(Δτ, Δp2, Δp)
    axpy!(Δτ, Δy2, Δy)
    #
    # recover Δd:
    #
    #   Δd ← Q Δp - Bᵀ Δy - c Δτ - rd
    #
    mul!(Δd, Symmetric(Q, :L), Δp)
    axpy!(-Δτ, c, Δd)

    mul!(Δd,           B',     Δy, -1, 1)
    axpy!(-1, rd, Δd)

    return niter, Δτ
end

############################################################################################
# solvepredictor! / solvecorrector!
############################################################################################

#
# solve for the Mehrotra predictor direction
#
#   [  H          -Bᵀ   -c  ] [ Δpa ]   [ rd - d ]
#   [  B           0    -g  ] [ Δya ] = [ rp     ]
#   [ cᵀ - 2pᵀQ/τ  gᵀ    S₀ ] [ Δτa ]   [ rτ - κ ]
#
# where S₀ = pᵀQp/τ² + κ/τ, and recover
#
#   Δκa = -κ (1 + Δτa/τ)
#
function solvepredictor!(s::HSDSolver{T}, rτmκ::T, aτ::AbstractVector{T}, S::T; force_tol::T, floor_tol::T) where {T}
    w = s.wrk
    τ = s.τ[]
    κ = s.κ[]
    atol = max(force_tol, floor_tol)

    axpby!(-1,  s.d, 0, w.f)
    axpby!( 1, w.rd, 1, w.f)
    #
    # solve for the directions Δpa, Δya, and Δτa
    #
    #   [  H          -Bᵀ   -c  ] [ Δpa ]   [ rd - d ]
    #   [  B           0    -g  ] [ Δya ] = [ rp     ]
    #   [ cᵀ - 2pᵀQ/τ  gᵀ    S₀ ] [ Δτa ]   [ rτ - κ ]
    #
    npred, Δτa = newton!(
        w.Δpa, w.Δya, w.Δda, τ, κ,
        s.kkt, s.settings.kkt, s.H, s.B, s.Q, s.c, s.g,
        w.f, w.rp, w.rd, rτmκ, w.Δp2, w.Δy2, aτ, S;
        atol
    )
    #
    # recover Δκa:
    #
    #   Δκa = -κ (1 + 1/τ Δτa)
    #
    Δκa = -κ * (τ + Δτa) / τ
    return npred, Δτa, Δκa
end

#
# solve for the Mehrotra combined direction
#
#   [  H          -Bᵀ   -c  ] [ Δp ]   [ rd*                         ]
#   [  B           0    -g  ] [ Δy ] = [ rp                          ]
#   [ cᵀ - 2pᵀQ/τ  gᵀ    S₀ ] [ Δτ ]   [ rτ - κ + (σμ - Δτa·Δκa) / τ ]
#
# where S₀ = pᵀQp/τ² + κ/τ
#
# where rd* is the corrected dual residual, S is the Woodbury capacitance scalar
#
#   S = Δp2ᵀ W Δp2 + (Δp2 - p/τ)ᵀ Q (Δp2 - p/τ) + κ/τ
#
# and
#
#   Δκ = (σμ - τκ - Δτa·Δκa - κ·Δτ) / τ
#
function solvecorrector!(s::HSDSolver{T}, μ::T, rτmκ::T, Δτa::T, Δκa::T, aτ::AbstractVector{T}, S::T; force_tol::T, floor_tol::T, nH::T, nB::T) where {T}
    w = s.wrk
    τ = s.τ[]
    κ = s.κ[]
    atol = max(force_tol, floor_tol)
    #
    # compute the largest step length αa ∈ (0, 1]
    # such that the perturbed iterates
    #
    #   p + αa Δpa ∈ K
    #   d + αa Δda ∈ K*
    #
    #   τ + αa Δτa ≥ 0
    #   κ + αa Δκa ≥ 0
    #
    # lie within their respective cones
    #
    αa = one(T)

    for v in vtxs(s.B)
        τpv, τdv = maxsteps(s, v, w.Δpa, w.Δda)
        αa = min(αa, τpv, τdv)
    end

    if Δτa < 0
        αa = min(αa, -τ / Δτa)
    end

    if Δκa < 0
        αa = min(αa, -κ / Δκa)
    end
    #
    # compute the centering parameter
    #
    #   σμ ← clamp(μa (μa / μ)², 0, μ)
    #
    # where
    #
    #   μa = (⟨p + αa Δpa, d + αa Δda⟩ + (τ + αa Δτa)(κ + αa Δκa)) / (ν + 1)
    #
    σμ = zero(T)

    for j in cols(s.B)
        σμ += (s.p[j] + αa * w.Δpa[j]) * (s.d[j] + αa * w.Δda[j])
    end

    σμ += (τ + αa * Δτa) * (κ + αa * Δκa)
    σμ /= (s.ν + 1)
    σμ = clamp(σμ * (σμ / μ)^2, zero(T), μ)
    #
    # set f to the Mehrotra corrector term:
    #
    #   f ← -d + (σμ e - Δpa ∘ Δda) / p
    #
    # where e is the Jordan identity element e ∈ K.
    #
    for v in vtxs(s.B)
        initcorrector!(s.K[v], v, w.f, s.caches, s.p, s.d, w.Δpa, w.Δda, σμ, s.B, s.conewrk)
    end

    axpy!(1, w.rd, w.f)
    #
    # compute fκ and fτ:
    #
    #   fκ = σμ - τκ -       Δτa Δκa
    #   fτ = rτ -  κ + (σμ - Δτa Δκa) / τ
    #
    fκ = σμ - τ * κ - Δτa * Δκa
    fτ = rτmκ + (σμ - Δτa * Δκa) / τ

    axpy!(-Δτa, w.Δy2, w.Δya)
    #
    # solve for the directions Δp, Δy, and Δτ
    #
    #   [  H          -Bᵀ   -c  ] [ Δp ]   [ rd*                       ]
    #   [  B           0    -g  ] [ Δy ] = [ rp                        ]
    #   [ cᵀ - 2pᵀQ/τ  gᵀ    S₀ ] [ Δτ ]   [ rτmκ + (σμ - Δτa·Δκa) / τ ]
    #
    ncorr, Δτ = newton!(
        w.Δp, w.Δy, w.Δd, τ, κ,
        s.kkt, s.settings.kkt, s.H, s.B, s.Q, s.c, s.g,
        w.f, w.rp, w.rd, fτ, w.Δp2, w.Δy2, aτ, S, w.Δya;
        atol
    )
    #
    # use iterative refinement to improve
    # the solutions Δp, Δy, and Δτ 
    #
    nrefine, refstat, Δτ = refinehsd!(
        w.Δp, w.Δy, Δτ, w.Δd,
        s.kkt, s.settings.kkt, s.H, s.B, s.Q, s.c, s.g, w.Qp, s.p, w.rd,
        τ, κ, w.rp, w.f, fτ, w.Δp2, w.Δy2, S,
        w.sy, w.sp, w.dp, w.dy, s.nc[], s.ng[];
        itmax=s.settings.refine_itmax, force_tol, floor_tol, stall=s.settings.refine_stall,
        nH, nB
    )
    ncorr += nrefine
    #
    # recover Δκ:
    #
    #   Δκ ← (fκ - κ Δτ) / τ
    #
    Δκ = (fκ - κ * Δτ) / τ
    return ncorr, Δτ, Δκ, refstat
end

############################################################################################
# startingpoint!
############################################################################################

function identitypoint(B::BlockSparseMatrix{T}, K) where {T}
    n = size(B, 2)
    m = size(B, 1)

    p = FVector{T}(undef, n)
    d = FVector{T}(undef, n)
    y = FVector{T}(undef, m)

    for v in vtxs(B)
        r = colrange(B, v)
        identity!(view(p, r), K[v])
        identity!(view(d, r), K[v])
    end

    fill!(y, zero(T))

    return p, d, y
end

############################################################################################
# init
############################################################################################

function CommonSolve.init(prob::IPMProblem{T, I}, settings::HSDSettings{T}) where {T, I}
    n = size(prob.B, 2)
    m = size(prob.B, 1)
    ν = conedegree(prob.K, prob.B)

    scaling = Scaling{T}(n, m)

    if settings.scale_itmax > 0
        B = copy(prob.B)
        Q = copy(prob.Q)
        c = copy(prob.c)
        g = copy(prob.g)

        equilibrate!(scaling, B, Q, c, g; itmax=settings.scale_itmax)
    else
        B = prob.B
        Q = prob.Q
        c = prob.c
        g = prob.g
    end

    R, P, B, kkt = make_kkt(settings.kkt, B)

    c = P * c
    Q = halfselectvtxs(halfselectvtxs(Q, R.perm), R.perm)
    cones = tounion(prob.K, R.perm)

    p, d, y = identitypoint(B, cones)

    caches = Caches(cones, B)

    for v in vtxs(B)
        initcache!(cache(caches, v, cones[v]))
    end

    H = allocblockdiag(B)
    Hc = allocblockdiag(B)
    conewrk = ConeWorkspace{T}(cones, B)
    hsdwrk = HSDWorkspace{T}(m, n)
    hist = HSDHistory{T}()
    ρ = FScalar{T}(undef)
    τ = FScalar{T}(undef)
    κ = FScalar{T}(undef)
    nc = FScalar{T}(undef)
    ng = FScalar{T}(undef)
    nB = FScalar{T}(undef)

    ρ[] = settings.rgmin
    nc[] = norm(c)
    ng[] = norm(g)
    nB[] = norm(B)
    Δy0 = FVector{T}(undef, m); fill!(Δy0, false)

    τ[] = one(T)
    κ[] = one(T)

    solver = HSDSolver(Q, H, Hc, B, c, g, p, d, y, cones,
        scaling, P, hsdwrk, caches, conewrk, kkt,
        hist, ν, settings, ρ, τ, κ, Δy0, nc, ng, nB, TimerOutput()
    )

    return solver
end

############################################################################################
# mu / isoptimal / infeasibility / ill-posed
############################################################################################
#
# compute the centrality parameter
#
#   μ = (pᵀd + τκ) / (ν + 1)
#
function mu(s::HSDSolver{T}) where {T}
    pd = dot(s.p, s.d)
    return (pd + s.τ[] * s.κ[]) / (s.ν + 1)
end

function isoptimal(s::HSDSolver{T}, μ::T, pres::T, dres::T) where {T}
    τ = s.τ[]
    w = s.wrk

    pQp = dot(s.p, w.Qp) / τ^2
    pobj = dot(s.c, s.p) / τ + pQp / 2
    dobj = dot(s.g, s.y) / τ - pQp / 2

    max(pres, dres) < s.settings.feas_tol && (μ < s.settings.gap_tol * τ^2 || pobj - dobj < s.settings.gap_tol * (1 + abs(pobj) + abs(dobj)))
end

function isprimalinfeasible(s::HSDSolver)
    w = s.wrk
    τ = s.τ[]
    κ = s.κ[]

    flag = τ / κ < s.settings.infeas_rel

    if flag
        gy = dot(s.g, s.y)
        ny = norm(s.y)
        flag = gy > s.settings.infeas_abs * ny * (1 + s.ng[])

        if flag
            nQp = norm(w.Qp)
            copyto!(w.f, w.Qp)
            axpy!(-1, s.d, w.f)
            mul!(w.f, s.B', s.y, -1, 1)
            flag = max(nQp, norm(w.f)) < s.settings.infeas_rel * gy * (1 + norm(s.d) / ny)
        end
    end

    return flag
end

function isdualinfeasible(s::HSDSolver)
    w = s.wrk
    τ = s.τ[]
    κ = s.κ[]

    flag = τ / κ < s.settings.infeas_rel

    if flag
        cp = dot(s.c, s.p)
        np = norm(s.p)
        flag = cp < -s.settings.infeas_abs * np * (1 + s.nc[])

        if flag
            mul!(w.sy, s.B, s.p)
            flag = max(norm(w.sy), norm(w.Qp)) < s.settings.infeas_rel * abs(cp)
        end
    end

    return flag
end

function isillposed(s::HSDSolver)
    return isillposed(s.hist; tol=s.settings.illposed_tol)
end

############################################################################################
# step!
############################################################################################

function step!(s::HSDSolver{T}) where {T}
    status = CONTINUE

    npred = 0
    ncorr = 0
    nwood = 0
    refstat = REACHED_FORCE

    α = zero(T)

    w = s.wrk
    τ = s.τ[]
    κ = s.κ[]
    #
    # compute negated residuals
    #
    #   [ rd ]   [ d ]   [  Q          -Bᵀ         c ] [ p ]
    #   [ rp ] = [ 0 ] - [  B           0         -g ] [ y ]
    #   [ rτ ]   [ κ ]   [ -cᵀ - 2pᵀQτ  gᵀ   pᵀQp/τ² ] [ τ ]
    #
    rτmκ = residuals!(s)
    #
    # compute the centrality parameter
    #
    #   μ = (pᵀd + τκ) / (ν + 1)
    #
    μ = mu(s)
    pres = norm(w.rp) / τ / (1 + s.ng[])
    dres = norm(w.rd) / τ / (1 + s.nc[])

    if isoptimal(s, μ, pres, dres)
        status = OPTIMAL
    elseif isprimalinfeasible(s)
        status = PRIMAL_INFEASIBLE
    elseif isdualinfeasible(s)
        status = DUAL_INFEASIBLE
    elseif isillposed(s)
        status = ILL_POSED
    else
        #
        # compute the Hessian
        #
        #   f''(w)
        #
        # of the primal barrier function f
        # at the Nestorov-Todd scaling point w
        #
        # for non-symmetric cones, no such point
        # exists, so the Hessian is replaced
        # by a Tuncel scaling matrix
        #
        @timeit s.timers "scale" flag = scale!(s)
        #
        # cache H and add the quadratic term
        #
        #   H ← H + Q
        #
        copyto!(s.Hc, s.H)

        for v in vtxs(s.B)
            axpy!(true, block(s.Q, v, v, v), block(s.H, v, v, v))
        end

        if !flag
            if s.settings.verbose > 1
                @warn "Scaling failed."
            end

            if isnearoptimal(s)
                status = NEAR_OPTIMAL
            else
                status = NUMERICAL_FAILURE
            end
        else
            nH = norm(Symmetric(s.H, :L))

            if !@timeit s.timers "initkkt" initkkt!(s, nH)
                if s.settings.verbose > 1
                    @warn "Failed to initialize KKT solver."
                end

                if isnearoptimal(s)
                    status = NEAR_OPTIMAL
                else
                    status = NUMERICAL_FAILURE
                end
            else
                #
                # compute tolerances for predictor and corrector solves
                #
                #   force: min(θμ, ceil)
                #   floor: 100ϵ (1 + max(‖rp‖, ‖rd‖))
                #
                force_tol = min(s.settings.forcing_frac * μ, s.settings.forcing_ceil)
                floor_tol = 100eps(T) * (1 + max(norm(w.rp, Inf), norm(w.rd, Inf)))
                #
                # solve for the Woodbury auxiliary directions
                #
                #   [ H  -Bᵀ ] [ Δp2 ]   [ c ]
                #   [ B   0  ] [ Δy2 ] = [ g ]
                #
                nwood = @timeit s.timers "woodbury" woodbury!(s; force_tol, floor_tol, y0 = s.Δy0)
                copyto!(s.Δy0, w.Δy2)
                #
                # compute the Woodbury capacitance scalar
                #
                #   S = Δp2ᵀ W Δp2 + (Δp2 - p/τ)ᵀ Q (Δp2 - p/τ) + κ/τ
                #
                S = capacitance!(s)
                #
                # solve for the Mehrota predictor direction
                #
                #   [  H          -Bᵀ             -c ] [ Δpa ]   [ rd - d ]
                #   [  B           0              -g ] [ Δya ] = [ rp     ]
                #   [ cᵀ - 2pᵀQ/τ  gᵀ  pᵀQp/τ² + κ/τ ] [ Δτa ]   [ rτ - κ ]
                #
                # and recover
                #
                #   Δκa = (-τκ - κ Δτa) / τ
                #
                npred, Δτa, Δκa = @timeit s.timers "predictor" solvepredictor!(s, rτmκ, w.aτ, S; force_tol, floor_tol)
                #
                # solve for the Mehrotra combined direction
                #
                #   [  H          -Bᵀ             -c ] [ Δp ]   [ rd*                         ]
                #   [  B           0              -g ] [ Δy ] = [ rp                          ]
                #   [ cᵀ - 2pᵀQ/τ  gᵀ  pᵀQp/τ² + κ/τ ] [ Δτ ]   [ rτ - κ + (σμ - Δτa·Δκa) / τ ]
                #
                # where rd* is the corrected dual residual, and recover
                #
                #   Δκ = (σμ - τκ - Δτa·Δκa - κ·Δτ) / τ
                #
                ncorr, Δτ, Δκ, refstat = @timeit s.timers "corrector" solvecorrector!(s, μ, rτmκ, Δτa, Δκa, w.aτ, S; force_tol, floor_tol, nH, nB=s.nB[])

                if s.settings.verbose > 1 && refstat != REACHED_FORCE
                    @info "KKT solve above target tolerance" refstat
                end
                #
                # find the largest step size α ∈ (0, 1] such that
                #
                #   p + α Δp ∈ K
                #   d + α Δd ∈ K*
                #
                #   τ + α Δτ ≥ 0
                #   κ + α Δκ ≥ 0
                #
                @timeit s.timers "maxsteps" begin
                    τp, τd = maxsteps(s, w.Δp, w.Δd; step_frac=s.settings.step_frac)
                    α = min(τp, τd)
                end

                if Δτ < 0
                    α = min(α, s.settings.step_frac * (-τ / Δτ))
                end

                if Δκ < 0
                    α = min(α, s.settings.step_frac * (-κ / Δκ))
                end
                #
                # compute the updated iterates
                #
                #   p ← p + α Δp ∈ K
                #   d ← d + α Δd ∈ K*
                #
                #   τ ← α Δτ ≥ 0
                #   κ ← α Δκ ≥ 0
                #
                axpy!(α, w.Δp, s.p)
                axpy!(α, w.Δd, s.d)
                axpy!(α, w.Δy, s.y)

                s.τ[] = τ + α * Δτ
                s.κ[] = κ + α * Δκ

                if isstalled(s)
                    if s.settings.verbose > 1
                        @warn "Stalling detected."
                    end

                    if isnearoptimal(s)
                        status = NEAR_OPTIMAL
                    else
                        status = STALLED
                    end
                elseif isnumfail(s)
                    if s.settings.verbose > 1
                        @warn "Step collapse detected."
                    end

                    if isnearoptimal(s)
                        status = NEAR_OPTIMAL
                    else
                        status = NUMERICAL_FAILURE
                    end
                end
            end
        end
    end

    push!(s.hist, (; μ, step=α, pres, dres, npred, ncorr, nwood, τ=s.τ[], κ=s.κ[], refstat))

    if status == CONTINUE && atfloor(s.hist; patience=s.settings.floor_patience)
        if s.settings.verbose > 1
            @warn "Refinement floor reached $(s.settings.floor_patience) consecutive times"
        end

        if isnearoptimal(s)
            status = NEAR_OPTIMAL
        else
            status = NUMERICAL_FAILURE
        end
    end

    return status
end

function reinit!(solver::HSDSolver{T}, prob::IPMProblem{T}; frac::Real=0.1, rgfrac::Real=0.0) where {T}
    @assert ncols(prob.B) == ncols(solver.B)
    @assert nrows(prob.B) == nrows(solver.B)
    @assert nvtxs(prob.B) == nvtxs(solver.B)
    @assert nouts(prob.B) == nouts(solver.B)

    c = copy(prob.c)
    g = copy(prob.g)

    for j in cols(solver.B)
        c[j] *= solver.scaling.cscl[j]
    end

    for i in rows(solver.B)
        g[i] *= solver.scaling.rscl[i]
    end

    mul!(solver.c, solver.P, c)
    copyto!(solver.g, g)

    solver.nc[] = norm(solver.c)
    solver.ng[] = norm(solver.g)

    τ = solver.τ[]
    κ = solver.κ[]

    if τ > sqrt(eps(T)) * κ
        f = convert(T, frac)
    else
        f = one(T)
    end

    if f < 1
        ldiv!(τ, solver.p)
        ldiv!(τ, solver.d)
        ldiv!(τ, solver.y)
    end

    p, d, y = identitypoint(solver.B, solver.K)
    axpby!(f, p, 1 - f, solver.p)
    axpby!(f, d, 1 - f, solver.d)
    axpby!(f, y, 1 - f, solver.y)

    solver.τ[] = one(T)
    solver.κ[] = dot(solver.p, solver.d) / solver.ν

    empty!(solver.hist)

    for v in vtxs(solver.B)
        initcache!(cache(solver.caches, v, solver.K[v]))
    end

    ρprev = solver.ρ[]
    ρcold = solver.settings.rgmin
    solver.ρ[] = rgfrac * ρprev + (1 - rgfrac) * ρcold

    return solver
end
