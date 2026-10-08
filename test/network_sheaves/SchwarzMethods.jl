using Test
using CellularSheaves
using LinearAlgebra
using SparseArrays
using Graphs
using BlockArrays

# 5-point finite-difference Laplacian for -Δu = f on an m×m interior grid of
# the unit square with homogeneous Dirichlet boundary conditions.
function poisson2d(m)
    h = 1 / (m + 1)
    T = spdiagm(-1 => fill(-1.0, m - 1), 0 => fill(2.0, m), 1 => fill(-1.0, m - 1))
    Id = sparse(1.0I, m, m)
    return (kron(Id, T) + kron(T, Id)) / h^2
end

# Label each grid dof by which of the p×p boxes it falls in.
function box_partition(m, p)
    parts = Vector{Int}(undef, m * m)
    for jy in 1:m, jx in 1:m
        bx = cld(jx * p, m)
        by = cld(jy * p, m)
        parts[(jy - 1) * m + jx] = (by - 1) * p + bx
    end
    return parts
end

# A linear system with an overlapping cover, its owner partition, and the
# local matrices used by the textbook reference iterations below.
struct ReferenceSystem
    A::SparseMatrixCSC{Float64,Int}
    f::Vector{Float64}
    doms::Vector{Vector{Int}}
    parts::Vector{Int}
    locals::Vector{Matrix{Float64}}
end

ReferenceSystem(A, f, doms, parts) = ReferenceSystem(A, f, doms, parts, [Matrix(A[d, d]) for d in doms])

# Robin local matrix built independently of the package: drop the coupling to
# outside dofs from the diagonal (algebraic Neumann) and add p.
function robin_local_matrix(A, d, p)
    inside = falses(size(A, 1))
    inside[d] .= true
    Ai = Matrix(A[d, d])
    for (ℓ, k) in enumerate(d)
        outside = [r for r in findall(!iszero, A[:, k]) if !inside[r]]
        isempty(outside) || (Ai[ℓ, ℓ] += p - sum(abs, A[outside, k]))
    end
    return Ai
end

with_robin(sys::ReferenceSystem, p) =
    ReferenceSystem(sys.A, sys.f, sys.doms, sys.parts, [robin_local_matrix(sys.A, d, p) for d in sys.doms])

# Textbook multiplicative Schwarz on a single global vector.
function reference_multiplicative(sys::ReferenceSystem, u, sweeps)
    u = copy(u)
    for _ in 1:sweeps, (d, Ai) in zip(sys.doms, sys.locals)
        u[d] += Ai \ (sys.f - sys.A * u)[d]
    end
    return u
end

# Textbook (optimized) restricted additive Schwarz.
function reference_ras(sys::ReferenceSystem, u, sweeps)
    u = copy(u)
    for _ in 1:sweeps
        r = sys.f - sys.A * u
        du = zeros(length(u))
        for (i, (d, Ai)) in enumerate(zip(sys.doms, sys.locals))
            c = Ai \ r[d]
            owned = sys.parts[d] .== i
            du[d[owned]] = c[owned]
        end
        u += du
    end
    return u
end

stationary(dd, f; kw...) = solve(SchwarzProblem(dd, f), SchwarzIteration(; kw...))
krylov(dd, f; kw...) = solve(SchwarzProblem(dd, f), SchwarzCG(; kw...))

# p×p boxes grouped into (p/2)×(p/2) blocks of 2×2 boxes.
function box_aggregation(p)
    q = p ÷ 2
    return GraphHomomorphism([(cld(by, 2) - 1) * q + cld(bx, 2) for by in 1:p for bx in 1:p])
end

@testset "SchwarzMethods" begin
    @testset "overlap_sheaf of a cover" begin
        doms = [[1, 2, 3], [3, 4, 5], [5, 6, 1], [7]]
        s = overlap_sheaf(doms)
        @test vertex_stalks(s) == [3, 3, 3, 1]
        @test ne(underlying_graph(s)) == 3
        @test !has_edge(underlying_graph(s), 1, 4)
        @test edge_stalks(s)[UnorderedPair(1, 2)] == 1
        @test get_restriction_map(s, 1, 2) == [0.0 0.0 1.0]
        @test get_restriction_map(s, 3, 1) == [1.0 0.0 0.0]
        # H⁰ ≅ functions on the union of the cover: one dimension per dof.
        @test size(nullspace_ldlt(s), 2) == 7

        cover = OverlapCover(doms)
        @test cover.graph == underlying_graph(s)
        @test cover.overlaps[1 => 3] == [1] && cover.overlaps[3 => 1] == [1]
        @test cover.members[1] == [(1, 1), (3, 1)]
        @test overlap_sheaf(cover) == s
    end

    @testset "overlapping_subdomains" begin
        A = poisson2d(6)
        parts = box_partition(6, 2)
        @test overlapping_subdomains(A, parts; overlap=0) == [findall(==(i), parts) for i in 1:4]
        doms = overlapping_subdomains(A, parts; overlap=1)
        @test all(issubset(findall(==(i), parts), d) for (i, d) in enumerate(doms))
        @test length(doms[1]) == 15
        @test_throws ArgumentError overlapping_subdomains(A, parts; overlap=-1)
    end

    m = 12
    A = poisson2d(m)
    n = m * m
    f = [sin(3x) * cos(2y) + 1 for y in range(0, 1; length=m) for x in range(0, 1; length=m)]
    u_exact = A \ f
    parts = box_partition(m, 3)
    doms = overlapping_subdomains(A, parts; overlap=2)
    dd = SchwarzDecomposition(A, doms; owner=parts)
    prob = SchwarzProblem(dd, f)
    sys = ReferenceSystem(A, f, doms, parts)

    @testset "decomposition structure" begin
        @test length(dd.locals) == 9
        @test dd.cover.subdomains == doms
        @test dd.cover.graph == underlying_graph(overlap_sheaf(dd))
        @test has_edge(dd.cover.graph, 1, 5)          # diagonal boxes overlap
        @test !has_edge(dd.cover.graph, 1, 3)
        @test dd.ownership.owner == parts
        @test all(lp.coupling == A[lp.dofs, lp.boundary] for lp in dd.locals)
        @test occursin("9 subdomains", sprint(show, dd))
        # Default owners hold each dof together with its whole stencil.
        dd_default = SchwarzDecomposition(A, doms)
        @test all(1:n) do k
            issubset(findall(!iszero, A[:, k]), doms[dd_default.ownership.owner[k]])
        end
        @test stationary(dd_default, f; sweep=ParallelSweep(), tol=1e-10).u ≈ u_exact rtol = 1e-8
    end

    @testset "localize / glue / disagreement" begin
        x = localize(dd, u_exact)
        @test x isa BlockVector
        @test blocklengths(axes(x, 1)) == length.(doms)
        @test glue(dd, x) ≈ u_exact
        @test overlap_disagreement(dd, x) ≈ 0 atol = 1e-10
        y = copy(x)
        y[Block(1)] .+= 1.0
        d = coboundary_map(overlap_sheaf(dd))
        @test overlap_disagreement(dd, y) ≈ norm(d * Vector(y))
        @test overlap_disagreement(dd, Vector(y)) ≈ overlap_disagreement(dd, y)
    end

    @testset "multiplicative Schwarz" begin
        r = stationary(dd, f; tol=1e-10)
        @test r.converged
        @test r.u ≈ u_exact rtol = 1e-8
        # Pushing along restriction maps keeps every iterate a global section.
        @test maximum(r.disagreements) < 1e-8 * norm(u_exact)
        @test overlap_disagreement(dd, r.x) < 1e-8 * norm(u_exact)

        x = localize(dd, zeros(n))
        for _ in 1:3
            schwarz_step!(x, prob, MultiplicativeSweep())
        end
        @test glue(dd, x) ≈ reference_multiplicative(sys, zeros(n), 3)

        warm = solve(SchwarzProblem(dd, f; u0=u_exact), SchwarzIteration())
        @test warm.iterations == 0 && warm.converged
    end

    @testset "parallel (Lions) Schwarz = RAS" begin
        u0 = randn(n)
        x = localize(dd, u0)
        for _ in 1:4
            schwarz_step!(x, prob, ParallelSweep())
        end
        @test glue(dd, x) ≈ reference_ras(sys, u0, 4)
        @test overlap_disagreement(dd, x) > 1e-6      # copies disagree mid-iteration

        r = stationary(dd, f; sweep=ParallelSweep(), tol=1e-10, maxiter=2000)
        @test r.converged
        @test r.u ≈ u_exact rtol = 1e-8
        @test last(r.disagreements) < 1e-6 * norm(u_exact)
        @test r.disagreements[1] ≈ 0 atol = 1e-12
    end

    @testset "more overlap converges faster" begin
        its = map(1:3) do δ
            dδ = SchwarzDecomposition(A, overlapping_subdomains(A, parts; overlap=δ); owner=parts)
            stationary(dδ, f).iterations
        end
        @test its[1] > its[2] > its[3]
    end

    @testset "input validation" begin
        @test_throws ArgumentError SchwarzDecomposition(A, [doms[1]])            # not a cover
        @test_throws ArgumentError SchwarzDecomposition(A, doms; owner=fill(1, n))
        nonoverlap = overlapping_subdomains(A, parts; overlap=0)
        @test_throws ArgumentError SchwarzDecomposition(A, nonoverlap; owner=parts)
        @test_throws ArgumentError SchwarzProblem(dd, f[1:end-1])
        @test_throws ArgumentError SchwarzDecomposition(sparse([1.0 2.0; 0.0 1.0]), [[1, 2]])
    end

    @testset "multicolor Schwarz" begin
        @test sort(reduce(vcat, dd.colors)) == 1:length(doms)
        @test length(dd.colors) < length(doms)
        for class in dd.colors, i in class, j in class
            i < j || continue
            @test isempty(intersect(doms[i], doms[j]))
            @test iszero(A[doms[i], doms[j]])
        end

        x = localize(dd, zeros(n))
        y = localize(dd, zeros(n))
        for _ in 1:3
            schwarz_step!(x, prob, MulticolorSweep())
            schwarz_step!(y, prob, MultiplicativeSweep(reduce(vcat, dd.colors)))
        end
        @test Vector(x) ≈ Vector(y)

        r = stationary(dd, f; sweep=MulticolorSweep(), tol=1e-10)
        @test r.converged
        @test r.u ≈ u_exact rtol = 1e-8
    end

    @testset "truncated pushforward coarse space" begin
        c = TruncatedPushforwardCoarseSpace(dd)
        @test coarse_dimension(c) == length(doms)
        @test vec(sum(c.basis; dims=2)) ≈ ones(n)                 # partition of unity
        @test c.matrix ≈ c.basis' * A * c.basis
        for (h, d) in enumerate(doms)
            @test issubset(findall(!iszero, c.basis[:, h]), d)    # supported on Ω_h
        end
        @test occursin("dimension 9", sprint(show, c))

        r = stationary(dd, f; sweep=MulticolorSweep(), coarse=c, tol=1e-10)
        @test r.converged
        @test r.u ≈ u_exact rtol = 1e-8
        @test r.iterations < stationary(dd, f; sweep=MulticolorSweep(), tol=1e-10).iterations

        rp = stationary(dd, f; sweep=ParallelSweep(), coarse=c, tol=1e-10)
        @test rp.converged
        @test rp.u ≈ u_exact rtol = 1e-8

        x = localize(dd, randn(n))
        schwarz_step!(x, prob, ParallelSweep())
        before = overlap_disagreement(dd, x)
        u_before = glue(dd, x)
        coarse_correct!(x, prob, c)
        @test overlap_disagreement(dd, x) ≈ before
        @test c.basis' * (f - A * glue(dd, x)) ≈ zeros(coarse_dimension(c)) atol = 1e-8 * norm(f)
        @test glue(dd, x) != u_before

        two_modes = TruncatedPushforwardCoarseSpace(dd; modes=[ones(n) repeat(1:m, m)])
        @test coarse_dimension(two_modes) == 2 * length(doms)
        @test_throws ArgumentError TruncatedPushforwardCoarseSpace(dd, GraphHomomorphism([1, 1, 1, 1, 1, 1, 1, 1, 3]))
    end

    @testset "coarse spaces and the pushforward sheaf" begin
        A8 = poisson2d(8)
        parts8 = box_partition(8, 4)
        dd8 = SchwarzDecomposition(A8, overlapping_subdomains(A8, parts8; overlap=1); owner=parts8)
        hom = box_aggregation(4)
        pf = pushforward_sheaf(hom, overlap_sheaf(dd8))
        exact = ExactPushforwardCoarseSpace(dd8, hom)
        @test vertex_stalks(pf) == length.(exact.decomposition.cover.subdomains)
        @test coarse_dimension(exact) == sum(vertex_stalks(pf))
        @test ne(underlying_graph(pf)) == ne(exact.decomposition.cover.graph)

        tc = TruncatedPushforwardCoarseSpace(dd8, hom)
        @test coarse_dimension(tc) == hom.n_target
        bases = all_fiber_bases(hom, overlap_sheaf(dd8))
        for h in 1:hom.n_target
            fiber = fiber_vertices(hom, h)
            local_section = reduce(vcat, [Vector(tc.basis[d, h]) for d in dd8.cover.subdomains[fiber]])
            B = bases[h]
            @test norm(B * (B \ local_section) - local_section) < 1e-8   # lies in (φ_*F)(h)
        end
    end

    @testset "exact pushforward coarse space" begin
        point = GraphHomomorphism(ones(Int, length(doms)))
        c1 = ExactPushforwardCoarseSpace(dd, point)
        @test coarse_dimension(c1) == n
        r1 = stationary(dd, f; sweep=MulticolorSweep(), coarse=c1, tol=1e-10)
        @test r1.iterations == 1                 # pushforward to a point = direct solve
        @test r1.u ≈ u_exact rtol = 1e-10

        A16 = poisson2d(16)
        parts16 = box_partition(16, 4)
        dd16 = SchwarzDecomposition(A16, overlapping_subdomains(A16, parts16; overlap=1); owner=parts16)
        f16 = ones(16 * 16)
        exact = ExactPushforwardCoarseSpace(dd16, box_aggregation(4))
        tc = TruncatedPushforwardCoarseSpace(dd16, box_aggregation(4))
        its = Dict(name => stationary(dd16, f16; sweep=MulticolorSweep(), coarse=c).iterations
                   for (name, c) in (:none => nothing, :exact => exact, :truncated => tc))
        @test its[:exact] < its[:truncated] < its[:none]
        @test coarse_dimension(exact) > coarse_dimension(tc)
        @test stationary(dd16, f16; coarse=exact, tol=1e-10).u ≈ A16 \ f16 rtol = 1e-8
    end

    @testset "two-level scalability" begin
        its = map((2, 4, 6)) do p
            mp = 6p
            Ap = poisson2d(mp)
            pp = box_partition(mp, p)
            ddp = SchwarzDecomposition(Ap, overlapping_subdomains(Ap, pp; overlap=1); owner=pp)
            fp = ones(mp * mp)
            one_level = stationary(ddp, fp; sweep=MulticolorSweep()).iterations
            two_level = stationary(ddp, fp; sweep=MulticolorSweep(),
                coarse=TruncatedPushforwardCoarseSpace(ddp)).iterations
            (one_level, two_level)
        end
        one, two = first.(its), last.(its)
        @test all(two .< one)
        @test two[end] - two[1] < one[end] - one[1]
    end

    @testset "Robin transmission" begin
        h = 1 / (m + 1)
        p = optimized_robin_parameter(5h) / h           # overlap 2 ⇒ L = 5h
        ddr = SchwarzDecomposition(A, doms; owner=parts, transmission=RobinTransmission(p))
        probr = SchwarzProblem(ddr, f)
        @test all(lp -> isempty(lp.interface), dd.locals)
        @test all(lp -> !isempty(lp.interface), ddr.locals)

        for sweep in (ParallelSweep(), MultiplicativeSweep(), MulticolorSweep())
            x = localize(ddr, u_exact)
            schwarz_step!(x, probr, sweep)
            @test glue(ddr, x) ≈ u_exact                 # exact solution is still the fixed point
        end

        robin_sys = with_robin(sys, p)
        u0 = randn(n)
        x = localize(ddr, u0)
        for _ in 1:3
            schwarz_step!(x, probr, ParallelSweep())
        end
        @test glue(ddr, x) ≈ reference_ras(robin_sys, u0, 3)              # optimized RAS

        x = localize(ddr, zeros(n))
        for _ in 1:2
            schwarz_step!(x, probr, MultiplicativeSweep())
        end
        @test glue(ddr, x) ≈ reference_multiplicative(robin_sys, zeros(n), 2)   # optimized multiplicative

        per_edge = SchwarzDecomposition(A, doms; owner=parts, transmission=RobinTransmission((i, j) -> p))
        @test [lp.shift for lp in per_edge.locals] == [lp.shift for lp in ddr.locals]
        @test_throws ArgumentError RobinTransmission(0.0)
        @test_throws ArgumentError SchwarzDecomposition(A, doms; owner=parts,
            transmission=RobinTransmission((i, j) -> i == 1 ? -1.0 : p))

        @test optimized_robin_parameter(0.1) ≈ cbrt(π^2) / cbrt(0.2)
        @test optimized_robin_parameter(0.1; kmin=2.0, η=1.0) ≈ cbrt(5.0) / cbrt(0.2)
        @test_throws ArgumentError optimized_robin_parameter(0.0)

        # Strips have no cross points; there the optimized parameter is far
        # better than Dirichlet transmission.
        ms = 24
        As = poisson2d(ms)
        strips = [cld(jx * 4, ms) for jy in 1:ms for jx in 1:ms]
        sdoms = overlapping_subdomains(As, strips; overlap=1)
        fs = ones(ms * ms)
        hs = 1 / (ms + 1)
        robin = RobinTransmission(optimized_robin_parameter(3hs) / hs)
        dirichlet = stationary(SchwarzDecomposition(As, sdoms; owner=strips), fs; sweep=ParallelSweep())
        optimized = stationary(SchwarzDecomposition(As, sdoms; owner=strips, transmission=robin), fs;
            sweep=ParallelSweep())
        @test dirichlet.converged && optimized.converged
        @test 3 * optimized.iterations < dirichlet.iterations
        @test optimized.u ≈ As \ fs rtol = 1e-6
    end

    @testset "Schwarz-preconditioned CG" begin
        A8 = poisson2d(8)
        parts8 = box_partition(8, 2)
        doms8 = overlapping_subdomains(A8, parts8; overlap=1)
        for transmission in (DirichletTransmission(), RobinTransmission(50.0))
            d8 = SchwarzDecomposition(A8, doms8; owner=parts8, transmission)
            for coarse in (nothing, TruncatedPushforwardCoarseSpace(d8),
                           ExactPushforwardCoarseSpace(d8, GraphHomomorphism([1, 1, 2, 2])))
                P = SchwarzPreconditioner(d8, coarse)
                Pm = reduce(hcat, [P * e for e in eachcol(Matrix(1.0I, 64, 64))])
                @test Pm ≈ Pm'
                @test eigmin(Symmetric(Pm)) > 0
                @test size(P) == (64, 64) && eltype(P) == Float64
            end
        end

        for coarse in (nothing, TruncatedPushforwardCoarseSpace(dd), ExactPushforwardCoarseSpace(dd, GraphHomomorphism(fill(1, 9))))
            r = krylov(dd, f; coarse, tol=1e-10)
            @test r.converged
            @test r.u ≈ u_exact rtol = 1e-8
            @test r.residuals[1] ≈ 1
            @test length(r.residuals) == r.iterations + 1
            @test glue(dd, r.x) ≈ r.u
        end
        warm = solve(SchwarzProblem(dd, f; u0=u_exact + 1e-6 * randn(n)), SchwarzCG(tol=1e-10))
        @test warm.converged && warm.residuals[1] < 0.1
        @test solve(SchwarzProblem(dd, zeros(n)), SchwarzCG()).iterations == 0

        ddr = SchwarzDecomposition(A, doms; owner=parts,
            transmission=RobinTransmission(optimized_robin_parameter(5 / (m + 1)) * (m + 1)))
        @test krylov(ddr, f; tol=1e-10).u ≈ u_exact rtol = 1e-8

        A128 = poisson2d(128)
        parts128 = box_partition(128, 16)
        d128 = SchwarzDecomposition(A128, overlapping_subdomains(A128, parts128; overlap=1); owner=parts128)
        f128 = ones(128^2)
        @test krylov(d128, f128; coarse=TruncatedPushforwardCoarseSpace(d128)).iterations <
              krylov(d128, f128).iterations
    end
end
