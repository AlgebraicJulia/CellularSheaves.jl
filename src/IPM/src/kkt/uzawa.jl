# Augmentation KKT backend (rebuilt — rebuild_spec.md §1). There is no settings object: α is a
# per-call argument (a trajectory-level parameter owned by the IPM, held as solver state, fixed at construction), and the ρ-shift ladder bounds are baked into the workspace at construction. The
# augmented system is F = (1/α)·H + BᵀB. There is no CRAIG iteration budget of our own: CRAIG terminates
# in exact arithmetic by iteration n+m, and Krylov.jl caps at its dimension-based default when itmax is
# unset — the hang guard belongs to the library that owns the iteration, at its termination bound.
#
# A different KKT backend is a different workspace type behind the same α/atol argument interface —
# there is no settings-level polymorphism (the KKTWorkspace abstract type is the only extension point).

# Bundled Gauss-quadrature (Ritz) reading of B on one base solve's rhs subspace — the σ̂² extremes plus the
# top-10 nodes/weights of the rhs spectral measure. Produced TWICE per base solve: at the fixed kmax=10
# budget and at kmax = the actual CRAIG iteration count (`k`), the latter a stand-in for a future inline-
# CRAIG harvest that would reuse the Lanczos vectors instead of re-bidiagonalizing off-books.
#   θ[i] = μ̂ᵢ  (Ritz nodes, DESCENDING, top-10 by μ, NaN-padded)   w[i] = quadrature weight (Σw=1)
#   k     = full GK step count (may exceed 10; θ/w keep only the leading 10)
struct SpectralHarvest{T}
    s2min::T
    s2max::T
    ritz_beta::T
    omm::T
    k::Int
    θ::Vector{T}
    w::Vector{T}
end

SpectralHarvest{T}() where {T} =
    SpectralHarvest{T}(T(NaN), T(NaN), T(NaN), T(NaN), 0, fill(T(NaN), 10), fill(T(NaN), 10))

# build from a ritz_spectral result — θraw/wraw are full-length (descending); retain the leading 10.
function SpectralHarvest(s2min::T, s2max::T, ritz_beta::T, omm::T,
                         θraw::AbstractVector{T}, wraw::AbstractVector{T}) where {T}
    θ = fill(T(NaN), 10); w = fill(T(NaN), 10)
    kk = min(length(θraw), 10)
    @views θ[1:kk] .= θraw[1:kk]
    @views w[1:kk] .= wraw[1:kk]
    return SpectralHarvest{T}(s2min, s2max, ritz_beta, omm, length(θraw), θ, w)
end

# Per-refinement-pass trace (E-series per-pass instrumentation). Fixed capture of the first
# PASSTRACE_LEN refinement passes (observed pass-count median 3, max 10; 6 covers the modelled regime).
# For refinement pass i = 1..PASSTRACE_LEN (pass i = the i-th refinement invocation AFTER the base solve):
#   dres[i]/pres[i]/tres[i] = the dual/primal/τ residual components at the TOP of refinement iteration i,
#       i.e. the ENTRY residual of pass i (the residual left by base + passes 1..i-1, which pass i's inner
#       solve consumes). NaN where the loop never reached iteration i. tres is the HSD 3-row τ-component
#       (NaN on the IPM 2-row system, which has no τ row).
#   nkry[i] = CRAIG iterations spent by pass i's inner solve (solve_kkt! return − 1, since the return
#       counts the direct solve as 1); −1 when pass i did NOT fire (the loop broke at the force/floor/stall
#       test before solving at iteration i).
# The base invocation (k=0) is NOT in here — its entry residual is the r0_* column, its Krylov is
# (pbase−1), and its post-solve residual equals this trace's pass-1 entry (dres[1]/pres[1]).
# Buffers are workspace-resident and filled in place (zero per-call allocation); a caller passes the
# per-role trace object (ptrace/ctrace/wtrace) so predictor/corrector/woodbury never clobber each other.
struct PassTrace{T}
    dres::Vector{T}
    pres::Vector{T}
    tres::Vector{T}
    nkry::Vector{Int}
end

const PASSTRACE_LEN = 6
PassTrace{T}() where {T} = PassTrace{T}(fill(T(NaN), PASSTRACE_LEN), fill(T(NaN), PASSTRACE_LEN),
                                        fill(T(NaN), PASSTRACE_LEN), fill(-1, PASSTRACE_LEN))

# reset to the empty state (no iteration reached, no pass fired) before a refinement loop fills it.
function reset!(t::PassTrace{T}) where {T}
    fill!(t.dres, T(NaN)); fill!(t.pres, T(NaN)); fill!(t.tres, T(NaN)); fill!(t.nkry, -1)
    return t
end

struct UzawaWorkspace{UPLO, T, I <: Integer} <: KKTWorkspace{T}
    F::FChordalTriangular{:N, UPLO, T, I}
    L::BlockSparseMatrix{T, I}
    facwrk::FactorizationWorkspace{T, I}
    divwrk::DivisionWorkspace{T, I}
    itrwrk::CraigWorkspace{T, T, Vector{T}}
    r::Vector{T}
    α::Scalar{T}
    r0::Scalar{T}                                  # pre-CRAIG residual ‖g-Bx‖ of the last base solve
    r1::Scalar{T}                                  # post-CRAIG residual ‖g-Bx‖ of the last base solve
    harvest10::Base.RefValue{SpectralHarvest{T}}   # last base solve's spectral reading at kmax=10
    harvestN::Base.RefValue{SpectralHarvest{T}}    # ... and at kmax = the base solve's actual CRAIG iters
    ptrace::PassTrace{T}                            # per-pass refinement trace, last predictor solve
    ctrace::PassTrace{T}                            # ... last corrector solve
    wtrace::PassTrace{T}                            # ... last woodbury solve (HSD; unused/empty on IPM)
    rgmin::T                                        # baked: ρ-shift ladder lower bound
    rgmax::T                                        # baked: ρ-shift ladder upper bound
end

function UzawaWorkspace(F::FChordalTriangular{:N, UPLO, T, I}, L::BlockSparseMatrix{T, I}, B::BlockSparseMatrix{T, I};
                        rgmin::T, rgmax::T) where {UPLO, T, I <: Integer}
    m, n = size(B)
    facwrk = FactorizationWorkspace(F)
    divwrk = DivisionWorkspace(F, 1)
    itrwrk = CraigWorkspace(m, n, Vector{T})
    r = zeros(T, m)
    α = ones(T)
    r0 = fill(T(NaN))
    r1 = fill(T(NaN))
    harvest10 = Ref(SpectralHarvest{T}())
    harvestN  = Ref(SpectralHarvest{T}())
    ptrace = PassTrace{T}(); ctrace = PassTrace{T}(); wtrace = PassTrace{T}()
    return UzawaWorkspace(F, L, facwrk, divwrk, itrwrk, r, α, r0, r1, harvest10, harvestN,
                          ptrace, ctrace, wtrace, rgmin, rgmax)
end

function make_kkt(B::BlockSparseMatrix{T, I}; elim::EliminationAlgorithm = DEFAULT_ELIMINATION_ALGORITHM,
                  rgmin::T, rgmax::T) where {T, I}
    weights, graph = weightedgraph(B)
    R, P, S = symbolic(weights, graph; alg=elim)
    B = selectvtxs(B, R.perm)
    F = FChordalTriangular{:N, :L, T, I}(S)
    L = B' * B
    wrk = UzawaWorkspace(F, L, B; rgmin, rgmax)
    return R, P, B, wrk
end

# build & factor F = (1/α)·A + BᵀB (BᵀB precomputed as wrk.L), with the ρ-shift ladder on failure.
function initkkt!(wrk::UzawaWorkspace{UPLO, T}, A::BlockSparseMatrix; α::T) where {UPLO, T}
    wrk.α[] = α
    return init_uzw!(wrk.facwrk, wrk.F, wrk.L, A, α, wrk.rgmin, wrk.rgmax)
end

function init_uzw!(
        facwrk::FactorizationWorkspace{T},
        F::ChordalTriangular{:N, UPLO, T},
        L::BlockSparseMatrix{T},
        A::BlockSparseMatrix{T},
        α::T,
        rgmin::T,
        rgmax::T
    ) where {UPLO, T}
    @assert size(F, 1) == size(L, 1) == size(A, 1)

    β = inv(α)
    ρ = rgmin
    #
    # assemble and factor the augmented matrix
    #
    #   F = β A + Bᵀ B
    #
    copyto!(F, L)
    axpy!(β, A, F)
    info = cholesky!(facwrk, F; check=false)
    #
    # ρ_applied is the diagonal shift ACTUALLY added — zero when the unshifted factorization succeeds,
    # otherwise the rung the ρ-ladder stopped on. (rgmin is the ladder's first rung, not "no shift"; the
    # two must be distinguishable, so we report the applied shift rather than the ladder's lower bound.)
    #
    ρ_applied = zero(T)
    #
    # on failure, add a diagonal shift ρ I and retry:
    #
    #   F ← β A + Bᵀ B + ρ I
    #
    if !iszero(info)
        copyto!(F, L)
        axpy!(β, A, F)
        axpy!(ρ, I, F)
        info = cholesky!(facwrk, F; check=false)
        ρ_applied = ρ
    end

    while !iszero(info) && 8ρ ≤ rgmax
        ρ *= 8
        copyto!(F, L)
        axpy!(β, A, F)
        axpy!(ρ, I, F)
        info = cholesky!(facwrk, F; check=false)
        ρ_applied = ρ
    end

    return iszero(info), ρ_applied
end

# Gauss-quadrature representation of the rhs spectral measure from a Golub-Kahan lower-bidiagonal (dv=diag
# αₖ, ev=subdiag βₖ). Ritz μ = singular(Bk)² approximate the preconditioned spectrum B F⁻¹ Bᵀ; recover to
# σ² via μ/(α(1-μ)) (μ∈(0,1) for F=(1/α)A+BᵀB). Returns (s2min, s2max, ritz_beta, ritz_gap, θ, w, omm):
#   θ[i] = svᵢ²           — Ritz nodes (μ-units), DESCENDING, length k (k<kmax on early GK termination)
#   w[i] = U[1,i]²         — normalized quadrature weight of node i w.r.t. the rhs (Σw=1); the u-basis is
#                            seeded by the rhs, so U[1,·] is the seed's coordinate. Absolute mass = ‖seed‖²
#                            = r0² (r0_p column), same point/norm as r0ref.
#   omm  = 1 - sv[1]²      — 1-μ̂max computed DIRECTLY (not reconstructed from the s2max sentinel; nodes
#                            near 1 lose their small complement to roundoff, so it is logged separately).
# s2min/s2max and the σ²max saturation sentinel (μ_max→1 or overshoot μ≥1 → −μ_max) are UNCHANGED.
function ritz_spectral(dv::AbstractVector{T}, ev::AbstractVector{T}, α_aug::T) where {T}
    k = length(dv)
    k ≥ 1 || return (T(NaN), T(NaN), T(NaN), T(NaN), T[], T[], T(NaN))
    # svd(Bk) calls LAPACK, which errors (not NaN-returns) on any non-finite entry — refuse those inputs.
    (all(isfinite, dv) && all(isfinite, @view ev[1:k-1])) ||
        return (T(NaN), T(NaN), T(NaN), T(NaN), T[], T[], T(NaN))
    Bk = k == 1 ? Bidiagonal(T[dv[1]], T[], :L) : Bidiagonal(collect(dv), collect(@view ev[1:k-1]), :L)
    # need S.U (left singular vectors) for the weights. gesdd (divide-and-conquer) can throw
    # LAPACKException on a pathological but finite bidiagonal (non-convergence) — the harvest is
    # diagnostic-only, so degrade to an all-NaN reading rather than aborting the whole solve/sweep.
    S = try
        svd(Matrix(Bk))
    catch err
        err isa LinearAlgebra.LAPACKException || rethrow()
        return (T(NaN), T(NaN), T(NaN), T(NaN), T[], T[], T(NaN))
    end
    θ = S.S .^ 2                                          # Ritz nodes μ, descending
    w = vec(S.U[1, :]) .^ 2                               # normalized weights, Σw=1 (U orthogonal, row 1 unit)
    μmin = θ[end]; μmax = θ[1]
    s2min = (zero(T) < μmin < one(T)) ? μmin / (α_aug * (one(T) - μmin)) : T(NaN)
    s2max = if μmax ≤ zero(T); T(NaN)
            elseif μmax ≥ one(T) - T(1e-12); -μmax
            else μmax / (α_aug * (one(T) - μmax)) end
    ritz_beta = isempty(ev) ? T(NaN) : ev[end]
    ritz_gap  = k ≥ 2 ? (θ[end-1] - θ[end]) : T(NaN)
    omm = one(T) - μmax                                   # 1 - μ̂max, raw (may be ≤0 on numerical overshoot)
    return (s2min, s2max, ritz_beta, ritz_gap, θ, w, omm)
end

# k-step preconditioned Golub-Kahan bidiagonalization of B under N = F⁻¹, seeded by `seed` (copied, not
# mutated), reusing the caller's factor F / divwrk. Internally-allocated work vectors (hack; a real
# inline-CRAIG harvest would reuse the Krylov workspace). Returns ritz_spectral(...).
function gk_spectral(divwrk::DivisionWorkspace{T}, F::ChordalTriangular{:N, UPLO, T},
                     B::BlockSparseMatrix{T}, seed::AbstractVector{T}, α_aug::T; kmax::Int = 10) where {UPLO, T}
    m, n = size(B)
    u = copy(seed)
    Nv = zeros(T, n); v = zeros(T, n); Atu = zeros(T, n); Av = zeros(T, m)
    dv = T[]; ev = T[]
    β = norm(u); (β == zero(T) || !isfinite(β)) && return (T(NaN), T(NaN), T(NaN), T(NaN), T[], T[], T(NaN))
    u ./= β
    for _ in 1:kmax
        mul!(Atu, B', u)
        @. Nv = Atu - β * Nv
        copyto!(v, Nv); ldiv!(divwrk, F, v); ldiv!(divwrk, F', v)
        a = sqrt(max(real(dot(v, Nv)), zero(T)))
        # stop the harvest before a non-finite bidiagonal entry can reach svd (a diverging refinement
        # direction at an extreme α seeds a non-finite rhs; svd(Bk) throws in LAPACK on NaN/Inf). The
        # harvest is diagnostic-only, so truncating to the finite prefix is the correct degradation.
        (a == zero(T) || !isfinite(a)) && break
        push!(dv, a)
        v ./= a; Nv ./= a
        mul!(Av, B, v)
        @. u = Av - a * u
        β = norm(u)
        !isfinite(β) && break
        push!(ev, β)
        β == zero(T) && break
        u ./= β
    end
    return ritz_spectral(dv, ev, α_aug)
end

# base KKT solve to `atol`; returns the base CRAIG iteration count (Krylov's dimension-based cap on hang).
function solve_kkt!(
    wrk::UzawaWorkspace{UPLO, T},
    x::AbstractVector{T},
    y::AbstractVector{T},
    A::BlockSparseMatrix{T},
    B::BlockSparseMatrix{T},
    f::AbstractVector{T},
    g::AbstractVector{T},
    y0 = nothing;
    atol::T,
) where {UPLO, T}
    return solve_uzw!(wrk.divwrk, wrk.itrwrk, x, y, wrk.r, wrk.F, A, B,
                      f, g, wrk.α[], atol, y0, wrk.r0, wrk.r1, wrk.harvest10, wrk.harvestN)
end

#
# Solve the KKT system
#
#   [ A -Bᵀ ] [ x ] = [ f ]
#   [ B  0  ] [ y ]   [ g ]
#
function solve_uzw!(
        divwrk::DivisionWorkspace{T},
        itrwrk::CraigWorkspace{T},
        x::AbstractVector{T},
        y::AbstractVector{T},
        r::AbstractVector{T},
        F::ChordalTriangular{:N, UPLO, T},
        A::BlockSparseMatrix{T},
        B::BlockSparseMatrix{T},
        f::AbstractVector{T},
        g::AbstractVector{T},
        α::T,
        atol::T,
        y0 = nothing,
        r0ref = nothing,
        r1ref = nothing,
        harvest10ref = nothing,
        harvestNref = nothing,
    ) where {UPLO, T}
    m, n = size(B)

    @assert length(x) == n
    @assert length(y) == m
    @assert length(r) == m
    @assert length(f) == n
    @assert length(g) == m
    @assert size(F, 1) == n
    @assert size(A, 1) == n
    niter = 1; β = inv(α)
    #
    # solve for x:
    #
    #   (β A + Bᵀ B) x =  β f + Bᵀ (g + β y₀)
    #
    copyto!(r, g)

    if !isnothing(y0)
        axpy!(β, y0, r)
    end

    copyto!(x, f)
    mul!(x, B', r, one(T), β)
    ldiv!(divwrk, F,  x)
    ldiv!(divwrk, F', x)
    #
    # compute the residual
    #
    #   r = g - B x
    #
    copyto!(r, g)
    mul!(r, B, x, -one(T), one(T))
    #
    # solve for δx and δy:
    #
    #   [ β A + Bᵀ B  -Bᵀ ] [ δx ] = [ 0 ]
    #   [ B            0  ] [ δy ]   [ r ]
    #
    function prec!(u, v)
        #
        # solve for u:
        #
        #   (β A + Bᵀ B) u = v
        #
        copyto!(u, v)
        ldiv!(divwrk, F,  u)
        ldiv!(divwrk, F', u)
    end

    N = LinearOperator(T, n, n, true, true, prec!)
    r0ref === nothing || (r0ref[] = norm(r))    # pre-CRAIG residual of the base solve
    #
    # Snapshot the EXACT rhs craig! is about to consume (the post-update recompute below overwrites r),
    # so the spectral harvest runs on craig!'s own rhs subspace.
    #
    harvesting = harvest10ref !== nothing || harvestNref !== nothing
    rseed = harvesting ? copy(r) : r
    craig!(itrwrk, B, r; ldiv = false, btol = zero(T), N, atol)
    ncraig = itrwrk.stats.niter
    niter += ncraig
    #
    # HACK (stand-in for a future inline-CRAIG harvest): AFTER craig! — now the true CRAIG iteration count
    # is known — run preconditioned Golub-Kahan on the snapshotted rhs, reusing the same F/divwrk, at BOTH
    # the fixed kmax=10 budget and kmax=ncraig (what an instrumented CRAIG would have produced in-flight).
    #
    if harvest10ref !== nothing
        s2min, s2max, rβ, _, θ, w, omm = gk_spectral(divwrk, F, B, rseed, α; kmax = 10)
        harvest10ref[] = SpectralHarvest(s2min, s2max, rβ, omm, θ, w)
    end
    if harvestNref !== nothing
        s2min, s2max, rβ, _, θ, w, omm = gk_spectral(divwrk, F, B, rseed, α; kmax = ncraig)
        harvestNref[] = SpectralHarvest(s2min, s2max, rβ, omm, θ, w)
    end
    #
    # update x:
    #
    #   x = x + δx
    #
    axpy!(one(T), itrwrk.x, x)
    #
    # recover y:
    #
    #   y = y₀ + 1/β (δy + (g - B x))
    #
    copyto!(r, g)
    mul!(r, B, x, -one(T), one(T))
    r1ref === nothing || (r1ref[] = norm(r))    # post-CRAIG residual of the base solve
    copyto!(y, itrwrk.y)
    axpy!(one(T), r, y)
    lmul!(α, y)

    if !isnothing(y0)
        axpy!(one(T), y0, y)
    end

    return niter
end
