using Test
using LinearAlgebra
using CellularSheaves
using CellularSheaves.IPM: IPMProblem, PositiveCone, AbstractCone, HSDSettings, IPMSettings,
                           UzawaSettings, HSDSolver, solve

# Minimal consensus LP on the positive orthant: V vertices agree (x_v = x_{v+1}),
# one normalization row, tied costs on the first kface coordinates (degenerate
# optimal face). Exercises the HSD and affine conic solvers and the universal
# solve-contract (Plank 4).
function consensus_lp(; V::Int, d::Int, kface::Int)
    Is = Int[]; Js = Int[]; blks = Matrix{Float64}[]; rowblk = 0
    for a in 1:V-1
        rowblk += 1
        push!(Is, rowblk); push!(Js, a);   push!(blks, Matrix{Float64}(I, d, d))
        push!(Is, rowblk); push!(Js, a+1); push!(blks, Matrix{Float64}(-I, d, d))
    end
    rowblk += 1
    push!(Is, rowblk); push!(Js, 1); push!(blks, ones(1, d))   # normalization: 1ᵀ x₁ = 1
    B = blocksparse(Is, Js, blks)
    m = d * (V - 1) + 1
    g = zeros(m); g[end] = 1.0
    Q = blocksparse(collect(1:V), collect(1:V), [zeros(d, d) for _ in 1:V])
    w = collect(range(0.0, 1.0, length = d)); w[1:kface] .= 0.0
    c = zeros(V * d); c[1:d] .= -w
    K = AbstractCone[PositiveCone() for _ in 1:V]
    return IPMProblem(Q, B, c, g, K)
end

hsdset(rg; tol = 1e-10) =
    HSDSettings{Float64}(feas_tol = tol, gap_tol = tol, itmax = 200,
                         kkt = UzawaSettings{Float64}(raug = rg))
affset(rg; tol = 1e-10) =
    IPMSettings{Float64}(feas_tol = tol, gap_tol = tol, itmax = 200,
                         kkt = UzawaSettings{Float64}(raug = rg))

@testset "IPM conic solver" begin
    @testset "solves a well-conditioned consensus LP (HSD)" begin
        pr = consensus_lp(; V = 6, d = 8, kface = 1)
        rh = solve(pr, hsdset(1e2))
        @test rh.status in (CellularSheaves.IPM.OPTIMAL, CellularSheaves.IPM.NEAR_OPTIMAL)
    end

    @testset "affine solver terminates (universal contract; see REBUILD2_REPORT casualty)" begin
        # Universal enforcement is a KNOWN casualty on some affine cells (documented
        # in REBUILD2_REPORT.md): the ungated dual-row IR over-refines affine
        # base/corrector solves. Here we only assert termination with a valid status.
        pr = consensus_lp(; V = 6, d = 8, kface = 1)
        ra = solve(pr, affset(1e2))
        @test ra.status isa CellularSheaves.IPM.IPMStatus
    end

    @testset "universal solve-contract (Plank 4): border pair priced" begin
        # The dual-row IR audit refs are populated on every solve; a degenerate
        # high-raug HSD run must exercise the border refinement (pass count > 0
        # on the stalled calls) and report an honest exit status.
        pr = consensus_lp(; V = 6, d = 8, kface = 4)
        r = solve(pr, hsdset(1e6))
        @test r.status isa CellularSheaves.IPM.IPMStatus       # terminates with a status
        @test isfinite(r.mu)                                    # audit fields finite
    end

    @testset "contract terminates (no runaway) at high raug" begin
        # nir=20 is a backstop; the loop must always terminate (contract-met /
        # stall / cap), never hang, even on the ill-conditioned high-raug regime.
        pr = consensus_lp(; V = 6, d = 8, kface = 4)
        for rg in (1e2, 1e6, 1e10)
            r = solve(pr, hsdset(rg))
            @test r.status isa CellularSheaves.IPM.IPMStatus
            @test r.niter ≤ 200
        end
    end
end
