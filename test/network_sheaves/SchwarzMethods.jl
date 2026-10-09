using Test
using CellularSheaves
using LinearAlgebra
using SparseArrays
using Graphs
using BlockArrays
using Random

# 5-point finite-difference Laplacian for -Δu = f on an m×m interior grid of
# the unit square with homogeneous Dirichlet boundary conditions.
function poisson2d(m)
    h = 1 / (m + 1)
    T = spdiagm(-1 => fill(-1.0, m - 1), 0 => fill(2.0, m), 1 => fill(-1.0, m - 1))
    Id = sparse(1.0I, m, m)
    return (kron(Id, T) + kron(T, Id)) / h^2
end

# Label each grid dof by which of the p×p boxes it falls in.
function index_boxes(m, p)
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

# Robin local matrix built independently of the package: replace the coupling
# to each outside dof by a Robin term p (algebraic Neumann plus Robin per face).
function robin_local_matrix(A, d, p)
    inside = falses(size(A, 1))
    inside[d] .= true
    Ai = Matrix(A[d, d])
    for (ℓ, k) in enumerate(d)
        outside = [r for r in findall(!iszero, A[:, k]) if !inside[r]]
        isempty(outside) || (Ai[ℓ, ℓ] += sum(r -> p - abs(A[r, k]), outside))   # one Robin term per face
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
        @test get_restriction_map(s, 1, 2) isa SelectionRestriction{Float64}
        @test Matrix(get_restriction_map(s, 1, 2)) == [0.0 0.0 1.0]
        @test Matrix(get_restriction_map(s, 3, 1)) == [1.0 0.0 0.0]
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
        parts = index_boxes(6, 2)
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
    parts = index_boxes(m, 3)
    doms = overlapping_subdomains(A, parts; overlap=2)
    dd = SchwarzDecomposition(A, doms; owner=parts)
    prob = SchwarzProblem(dd, f)
    sys = ReferenceSystem(A, f, doms, parts)

    @testset "decomposition structure" begin
        @test length(dd.locals) == 9
        @test dd.cover.subdomains == doms
        @test dd.cover.graph == underlying_graph(overlap_sheaf(dd))
        @test has_edge(dd.cover.graph, 1, 5)          # diagonal boxes overlap
        @test isempty(intersect(doms[1], doms[3]))      # interiors of boxes 1 and 3 are disjoint
        @test has_edge(dd.cover.graph, 1, 3) == !isempty(intersect(dd.cover.stalks[1], dd.cover.stalks[3]))
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
        @test blocklengths(axes(x, 1)) == length.(dd.cover.stalks)
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
        @test_throws ArgumentError SchwarzDecomposition(A, nonoverlap; owner=parts,
            transmission=RobinTransmission(1.0))                             # Robin needs overlap
        @test_throws ArgumentError SchwarzProblem(dd, f[1:end-1])
        @test_throws ArgumentError SchwarzDecomposition(sparse([1.0 2.0; 0.0 1.0]), [[1, 2]];
            transmission=RobinTransmission(1.0))                             # Robin needs symmetric A
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
        parts8 = index_boxes(8, 4)
        dd8 = SchwarzDecomposition(A8, overlapping_subdomains(A8, parts8; overlap=1); owner=parts8)
        hom = box_aggregation(4)
        pf = pushforward_sheaf(hom, overlap_sheaf(dd8))
        exact = ExactPushforwardCoarseSpace(dd8, hom)
        @test vertex_stalks(pf) == length.(exact.decomposition.cover.stalks)
        @test coarse_dimension(exact) == sum(length, exact.decomposition.cover.subdomains)
        @test ne(underlying_graph(pf)) == ne(exact.decomposition.cover.graph)

        tc = TruncatedPushforwardCoarseSpace(dd8, hom)
        @test coarse_dimension(tc) == hom.n_target
        bases = all_fiber_bases(hom, overlap_sheaf(dd8))
        for h in 1:hom.n_target
            fiber = fiber_vertices(hom, h)
            local_section = reduce(vcat, [Vector(tc.basis[d, h]) for d in dd8.cover.stalks[fiber]])
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
        parts16 = index_boxes(16, 4)
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
            pp = index_boxes(mp, p)
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
        @test all(lp -> isempty(lp.faces), dd.locals)
        @test all(lp -> !isempty(lp.faces), ddr.locals)
        @test ddr.transmission isa RobinTransmission
        # A corner dof of box 1 has two faces, towards boxes 2 and 4.
        corner = ddr.locals[1]
        @test any(length(unique(face.source for face in corner.faces if face.dof == ℓ)) == 2
                  for ℓ in unique(face.dof for face in corner.faces))
        for (i, lp) in enumerate(ddr.locals), face in lp.faces
            @test has_edge(ddr.cover.graph, i, face.source)
            @test ddr.cover.stalks[face.source][face.source_dof] == lp.dofs[face.dof]
            @test ddr.cover.interior[face.source][face.source_dof]
        end

        for sweep in (ParallelSweep(), MultiplicativeSweep(), MulticolorSweep())
            x = localize(ddr, u_exact)
            schwarz_step!(x, probr, sweep)
            @test glue(ddr, x) ≈ u_exact                 # exact solution is still the fixed point
        end

        # From a global section, one parallel step is one step of optimized RAS.
        robin_sys = with_robin(sys, p)
        u0 = randn(n)
        x = localize(ddr, u0)
        schwarz_step!(x, probr, ParallelSweep())
        @test glue(ddr, x) ≈ reference_ras(robin_sys, u0, 1)

        # No pushes under Robin transmission: copies stay distinct, and the
        # colored sweep is the alternating sweep in color order.
        x = localize(ddr, zeros(n))
        y = localize(ddr, zeros(n))
        for _ in 1:2
            schwarz_step!(x, probr, MulticolorSweep())
            schwarz_step!(y, probr, MultiplicativeSweep(reduce(vcat, ddr.colors)))
        end
        @test Vector(x) ≈ Vector(y)
        @test overlap_disagreement(ddr, x) > 1e-6

        per_edge = SchwarzDecomposition(A, doms; owner=parts, transmission=RobinTransmission((i, j) -> p))
        @test [lp.faces for lp in per_edge.locals] == [lp.faces for lp in ddr.locals]
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

    @testset "Robin transmission at cross points" begin
        # 4×4 boxes meet at nine cross points. Small Robin parameters used to
        # make the stationary iteration diverge there.
        mc = 32
        hc = 1 / (mc + 1)
        Ac = poisson2d(mc)
        pc = index_boxes(mc, 4)
        dc = overlapping_subdomains(Ac, pc; overlap=1)
        fc = ones(mc * mc)
        pstar = optimized_robin_parameter(3hc) / hc
        dirichlet = stationary(SchwarzDecomposition(Ac, dc; owner=pc), fc; sweep=MulticolorSweep())
        for scale in (0.125, 0.5, 1.0, 4.0)
            ddc = SchwarzDecomposition(Ac, dc; owner=pc, transmission=RobinTransmission(scale * pstar))
            for sweep in (ParallelSweep(), MultiplicativeSweep(), MulticolorSweep())
                r = stationary(ddc, fc; sweep, maxiter=2000)
                @test r.converged
                @test r.u ≈ Ac \ fc rtol = 1e-6
            end
        end
        optimized = stationary(SchwarzDecomposition(Ac, dc; owner=pc, transmission=RobinTransmission(pstar)), fc;
            sweep=MulticolorSweep())
        @test 2 * optimized.iterations < dirichlet.iterations
    end

    @testset "notched rectangle" begin
        dom = notched_rectangle(15)
        @test dom.h ≈ 1 / 16
        @test size(dom.inside) == (31, 15)
        @test length(dom.points) == 31 * 15 - 5 * 8     # the slot removes 5 × 8 grid points
        @test !any(abs(x - 1) <= 0.125 && y >= 0.5 for (x, y) in dom.points)
        An = poisson_matrix(dom)
        @test issymmetric(An)
        @test all(>(0), diag(An)) && all(<=(0), An - Diagonal(diag(An)))
        @test isposdef(Matrix(An))
        V = grid_values(dom, collect(1.0:length(dom.points)))
        @test count(isnan, V) == 40 && V[1, 1] == 1.0

        wide = notched_rectangle(15; notch_width=0.5)
        @test sort(unique(box_partition(wide, 8, 2))) == 1:14   # the two boxes inside the slot are dropped
        parts_n = box_partition(dom, 4, 2)
        @test box_partition(unit_square(4), 2, 2) == [1, 1, 2, 2, 1, 1, 2, 2, 3, 3, 4, 4, 3, 3, 4, 4]

        fn = ones(length(dom.points))
        un = An \ fn
        hn = dom.h
        ddn = SchwarzDecomposition(An, overlapping_subdomains(An, parts_n; overlap=2); owner=parts_n)
        robin = RobinTransmission(optimized_robin_parameter(5hn) / hn)
        ddr = SchwarzDecomposition(An, overlapping_subdomains(An, parts_n; overlap=2); owner=parts_n,
            transmission=robin)
        for d in (ddn, ddr), sweep in (MultiplicativeSweep(), MulticolorSweep(), ParallelSweep())
            r = stationary(d, fn; sweep, tol=1e-10, maxiter=3000)
            @test r.converged
            @test r.u ≈ un rtol = 1e-8
        end
        c = TruncatedPushforwardCoarseSpace(ddn)
        @test stationary(ddn, fn; sweep=MulticolorSweep(), coarse=c, tol=1e-10).u ≈ un rtol = 1e-8
        @test krylov(ddn, fn; coarse=c, tol=1e-10).u ≈ un rtol = 1e-8
    end

    @testset "Schwarz-preconditioned CG" begin
        A8 = poisson2d(8)
        parts8 = index_boxes(8, 2)
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
        parts128 = index_boxes(128, 16)
        d128 = SchwarzDecomposition(A128, overlapping_subdomains(A128, parts128; overlap=1); owner=parts128)
        f128 = ones(128^2)
        @test krylov(d128, f128; coarse=TruncatedPushforwardCoarseSpace(d128)).iterations <
              krylov(d128, f128).iterations
    end

    @testset "coordination problems and the pushforward tower" begin
        # A coordination sheaf on a 6×6 grid of agents with 2-dimensional stalks
        # and random (invertible, symmetric) restriction maps, pinned at the four
        # corner agents (targets).
        rng = Random.MersenneTwister(7)
        g = Graphs.grid([6, 6])
        s = sheaf_from_graph(g, 2, d -> randn(rng, d, d); symmetric_edges=true)
        pins = Dict(v => randn(rng, 2) for v in (1, 6, 31, 36))
        x_p, _ = harmonic_extension(s, pins)
        boundary = sort(collect(keys(pins)))
        interior = setdiff(1:36, boundary)
        H, B = restricted_laplacian_blocks(s, interior, boundary)
        p = reduce(vcat, [pins[v] for v in boundary])
        q_star = reduce(vcat, [x_p[Block(v)] for v in interior])
        @test H * q_star ≈ -B * p

        # Teams of agents (quadrants of the grid) are the subdomains; each team
        # solves a pinned harmonic extension with its neighbours as targets.
        team(v) = 1 + (mod1(v, 6) > 3) + 2 * (cld(v, 6) > 3)
        dof_team = reduce(vcat, [fill(team(v), 2) for v in interior])
        doms_c = overlapping_subdomains(H, dof_team; overlap=2)
        ddc = SchwarzDecomposition(H, doms_c; owner=dof_team)
        probc = SchwarzProblem(ddc, -B * p)
        for alg in (SchwarzIteration(sweep=MulticolorSweep(), tol=1e-10, maxiter=5000),
                    SchwarzIteration(sweep=MultiplicativeSweep(), coarse=TruncatedPushforwardCoarseSpace(ddc), tol=1e-10),
                    SchwarzCG(coarse=TruncatedPushforwardCoarseSpace(ddc), tol=1e-10))
            r = solve(probc, alg)
            @test r.converged
            @test r.u ≈ q_star rtol = 1e-7
        end

        # The pushforward of a sheaf along a graph homomorphism has Laplacian
        # Bᵀ L B, where B lifts fiber sections: the hierarchical solve of the
        # nested tower is a Galerkin coarse problem with prolongation B.
        hom = GraphHomomorphism([team(v) for v in 1:36])
        pf = pushforward_sheaf(hom, s)
        bases = all_fiber_bases(hom, s)
        offsets = [0; cumsum(vertex_stalks(s))]
        Blift = zeros(72, sum(size.(bases, 2)))
        col = 0
        for h in 1:4
            k = size(bases[h], 2)
            fiber_rows = reduce(vcat, [collect(offsets[v]+1:offsets[v+1]) for v in fiber_vertices(hom, h)])
            Blift[fiber_rows, col+1:col+k] = bases[h]
            col += k
        end
        @test size(Blift, 2) == 8                       # each team: a 2-dimensional space of sections
        L = Matrix(sheaf_laplacian_matrix(s))
        @test Matrix(sheaf_laplacian_matrix(pf)) ≈ Blift' * L * Blift atol = 1e-8 * norm(L)
    end

    @testset "ghost-layer covers" begin
        cover = dd.cover
        for (i, lp) in enumerate(dd.locals)
            @test cover.stalks[i] == sort(union(doms[i], lp.boundary))
            @test cover.stalks[i][lp.interior] == doms[i]
            @test cover.stalks[i][lp.ghosts] == lp.boundary
            # Every ghost value comes from its owner through a shared edge stalk.
            for k in lp.boundary
                j = dd.ownership.owner[k]
                @test has_edge(cover.graph, i, j)
                @test k in cover.stalks[i][cover.overlaps[i => j]]
            end
        end
        s = overlap_sheaf(dd)
        @test s isa EuclideanSheaf{Float64,SelectionRestriction{Float64}}
        @test coboundary_map(s) isa SparseMatrixCSC
        @test ghost_layer_cover(A, doms).stalks == cover.stalks
        @test OverlapCover(doms).stalks == doms                  # no ghosts without a matrix

        # Ghost layers carry the interface, so non-overlapping subdomains work:
        # the sweeps become block Gauss–Seidel and block Jacobi.
        blocks0 = overlapping_subdomains(A, parts; overlap=0)
        dd0 = SchwarzDecomposition(A, blocks0; owner=parts)
        @test all(lp -> !isempty(lp.ghosts), dd0.locals)
        for sweep in (MultiplicativeSweep(), MulticolorSweep(), ParallelSweep())
            r = stationary(dd0, f; sweep, tol=1e-10, maxiter=5000)
            @test r.converged
            @test r.u ≈ u_exact rtol = 1e-8
        end
    end

    @testset "sheaf ADMM" begin
        objectives = local_objectives(prob)
        K = spzeros(n, n)
        b = zeros(n)
        for (obj, s) in zip(objectives, dd.cover.stalks)
            K[s, s] += obj.matrix
            b[s] += obj.rhs
            @test isposdef(Matrix(obj.matrix) + 1e-12I)      # every local objective is convex
        end
        @test K ≈ A                                          # the energy splits exactly
        @test b ≈ f

        h = 1 / (m + 1)
        ρ = optimized_robin_parameter(5h) / h
        for penalty in (:stalk, :shared)
            r = solve(prob, SheafADMM(rho=ρ, penalty=penalty, tol=1e-9, maxiter=5000))
            @test r.converged
            @test r.u ≈ u_exact rtol = 1e-7
        end
        @test_throws ArgumentError solve(prob, SheafADMM(rho=0.0))
        @test_throws ArgumentError solve(prob, SheafADMM(rho=1.0, penalty=:edges))

        # On strips, sheaf diffusion can replace the exact projection; ADMM
        # converges but optimized Schwarz needs far fewer iterations.
        ms = 24
        As = poisson2d(ms)
        strips = [cld(jx * 4, ms) for jy in 1:ms for jx in 1:ms]
        sdoms = overlapping_subdomains(As, strips; overlap=1)
        fs = ones(ms * ms)
        hs = 1 / (ms + 1)
        ps = optimized_robin_parameter(3hs) / hs
        dds = SchwarzDecomposition(As, sdoms; owner=strips)
        admm = solve(SchwarzProblem(dds, fs), SheafADMM(rho=ps, maxiter=5000))
        diffused = solve(SchwarzProblem(dds, fs), SheafADMM(rho=ps, penalty=:shared, projection_steps=1, maxiter=5000))
        robin = stationary(SchwarzDecomposition(As, sdoms; owner=strips, transmission=RobinTransmission(ps)), fs;
            sweep=ParallelSweep())
        @test admm.converged && diffused.converged && robin.converged
        @test admm.u ≈ As \ fs rtol = 1e-6
        @test 5 * robin.iterations < admm.iterations
    end

    @testset "nonsymmetric M-matrices (upwind transport)" begin
        # ρu + b·∇u = 1 on an m×m grid with a rotating velocity field, upwinded:
        # a nonsingular, nonsymmetric M-matrix.
        m = 24
        h = 1 / (m + 1)
        idx(i, j) = (j - 1) * m + i
        I_, J_, V_ = Int[], Int[], Float64[]
        for j in 1:m, i in 1:m
            x, y = i * h - 0.5, j * h - 0.5
            b = (-y, x)
            diag = 0.5
            for (k, (di, dj)) in enumerate(((1, 0), (0, 1)))
                a = abs(b[k]) / h
                s = b[k] > 0 ? 1 : -1
                diag += a
                ii, jj = i + s * di, j + s * dj
                (1 <= ii <= m && 1 <= jj <= m) && (push!(I_, idx(i, j)); push!(J_, idx(ii, jj)); push!(V_, -a))
            end
            push!(I_, idx(i, j)); push!(J_, idx(i, j)); push!(V_, diag)
        end
        A = sparse(I_, J_, V_, m * m, m * m)
        f = ones(m * m)
        @test !issymmetric(A)
        parts = index_boxes(m, 3)
        doms = overlapping_subdomains(A, parts; overlap=1)
        dd = SchwarzDecomposition(A, doms; owner=parts)
        u = A \ f
        for sweep in (MultiplicativeSweep(), MulticolorSweep(), ParallelSweep())
            r = solve(SchwarzProblem(dd, f), SchwarzIteration(sweep=sweep, tol=1e-10, maxiter=5000))
            @test r.converged
            @test r.u ≈ u rtol = 1e-8
        end
        g = solve(SchwarzProblem(dd, f), SchwarzGMRES(sweep=ParallelSweep(), tol=1e-10, maxiter=200))
        @test g.converged
        @test g.u ≈ u rtol = 1e-8
        @test_throws ArgumentError solve(SchwarzProblem(dd, f), SchwarzCG())
        @test_throws ArgumentError solve(SchwarzProblem(dd, f), SheafADMM(rho=1.0))
        @test_throws ArgumentError TruncatedPushforwardCoarseSpace(dd)
        @test_throws ArgumentError SchwarzDecomposition(A, doms; owner=parts, transmission=RobinTransmission(10.0))

        # refactor: a new matrix on the same cover.
        A2 = A + 0.25I
        dd2 = refactor(dd, A2)
        @test dd2.cover === dd.cover && dd2.colors === dd.colors
        r2 = solve(SchwarzProblem(dd2, f), SchwarzIteration(sweep=MulticolorSweep(), tol=1e-10, maxiter=5000))
        @test r2.u ≈ A2 \ f rtol = 1e-8
        far = copy(A2); far[1, m * m] = -1e-3                          # outside the pattern
        @test_throws ArgumentError refactor(dd, far)
        # A cover built from a wider structure accepts any matrix inside it.
        wide = A + sparse(transpose(A)) + I
        ddw = SchwarzDecomposition(A, doms; owner=parts, structure=wide)
        @test solve(SchwarzProblem(refactor(ddw, sparse(transpose(A))), f),
            SchwarzIteration(sweep=MulticolorSweep(), tol=1e-10, maxiter=5000)).u ≈ sparse(transpose(A)) \ f rtol = 1e-8
        @test_throws ArgumentError SchwarzDecomposition(A, doms; owner=parts, structure=sparse(1.0I, m * m, m * m))

        # refactor!: in place, on a fixed stored pattern with explicit zeros.
        function embed(M)
            Z = copy(wide)
            nonzeros(Z) .= 0
            for (i, j, v) in zip(findnz(M)...)
                Z[i, j] = v
            end
            return Z
        end
        ddi = SchwarzDecomposition(embed(A), doms; owner=parts, dropzeros=false)
        @test nnz(ddi.A) == nnz(wide)
        for M in (A2, sparse(transpose(A)), A, A2)
            Z = embed(M)
            @test refactor!(ddi, Z) === ddi
            @test ddi.At == sparse(transpose(Z))
            for sweep in (MulticolorSweep(), ParallelSweep())
                @test solve(SchwarzProblem(ddi, f), SchwarzIteration(sweep=sweep, tol=1e-10, maxiter=5000)).u ≈
                    M \ f rtol = 1e-8
            end
        end
        Z = embed(sparse(transpose(A)))
        refactor!(ddi, Z; At=copy(transpose(Z)))
        @test solve(SchwarzProblem(ddi, f), SchwarzGMRES(sweep=ParallelSweep(), tol=1e-10, maxiter=200)).u ≈
            sparse(transpose(A)) \ f rtol = 1e-8
        @test_throws ArgumentError refactor!(ddi, A2)                   # pattern differs
        @test_throws ArgumentError refactor!(ddi, Z; At=A)

        # Inexact local solves: symmetric Gauss–Seidel inside each subdomain.
        M = SymmetricGaussSeidel(A)
        D, Lo, Up = Diagonal(A), tril(A, -1), triu(A, 1)
        r = collect(range(-1.0, 1.0; length=m * m))
        @test ldiv!(similar(r), M, r) ≈ (Matrix(D + Up) \ (D * (Matrix(D + Lo) \ r))) rtol = 1e-12
        M3 = SymmetricGaussSeidel(A; sweeps=3)
        y = zeros(m * m)
        for _ in 1:3
            y += ldiv!(similar(r), M, r - A * y)
        end
        @test ldiv!(similar(r), M3, r) ≈ y rtol = 1e-12
        @test_throws ArgumentError SymmetricGaussSeidel(A - Diagonal(A))
        @test_throws ArgumentError SymmetricGaussSeidelLocalSolve(sweeps=0)
        for sweeps in (1, 2)
            dds = SchwarzDecomposition(A, doms; owner=parts,
                local_solver=SymmetricGaussSeidelLocalSolve(; sweeps))
            @test all(lp -> lp.factor isa SymmetricGaussSeidel, dds.locals)
            for sweep in (MultiplicativeSweep(), MulticolorSweep(), ParallelSweep())
                rs = solve(SchwarzProblem(dds, f), SchwarzIteration(sweep=sweep, tol=1e-10, maxiter=20_000))
                @test rs.converged
                @test rs.u ≈ u rtol = 1e-8
            end
            gs = solve(SchwarzProblem(dds, f), SchwarzGMRES(sweep=ParallelSweep(), tol=1e-10, maxiter=500))
            @test gs.converged
            @test gs.u ≈ u rtol = 1e-8
        end
        # Inexact local solves take more iterations than exact ones.
        exact = solve(SchwarzProblem(dd, f), SchwarzIteration(sweep=MulticolorSweep(), tol=1e-10, maxiter=5000))
        inexact = solve(SchwarzProblem(SchwarzDecomposition(A, doms; owner=parts,
                local_solver=SymmetricGaussSeidelLocalSolve()), f),
            SchwarzIteration(sweep=MulticolorSweep(), tol=1e-10, maxiter=20_000))
        @test inexact.iterations > exact.iterations
        # refactor! only gathers the new values into the smoothers.
        ddg = SchwarzDecomposition(embed(A), doms; owner=parts, dropzeros=false,
            local_solver=SymmetricGaussSeidelLocalSolve())
        for M in (A2, sparse(transpose(A)), A2)
            refactor!(ddg, embed(M))
            @test all(lp -> lp.factor isa SymmetricGaussSeidel, ddg.locals)
            @test solve(SchwarzProblem(ddg, f), SchwarzGMRES(sweep=ParallelSweep(), tol=1e-10, maxiter=500)).u ≈
                M \ f rtol = 1e-8
        end
        @test refactor(ddg, A2).locals[1].factor isa SymmetricGaussSeidel
        # The direct RAS application equals one parallel sweep from the zero cochain.
        for d in (dd, ddg)
            P = CellularSheaves.NetworkSheaves.SchwarzMethods.SchwarzSweepPreconditioner(d, ParallelSweep())
            sweep_prob = SchwarzProblem(d, r)
            xz = localize(d, zeros(m * m))
            schwarz_step!(xz, sweep_prob, ParallelSweep())
            @test P * r ≈ glue(d, xz) rtol = 1e-13
        end
    end

    @testset "refactor! with ChordalLDLt local factors" begin
        m = 12
        L = poisson2d(m)
        parts = index_boxes(m, 3)
        doms = overlapping_subdomains(L, parts; overlap=1)
        dd = SchwarzDecomposition(L, doms; owner=parts)
        f = ones(m * m)
        refactor!(dd, 2L)
        refactor!(dd, L + 3I)
        @test solve(SchwarzProblem(dd, f), SchwarzIteration(sweep=MulticolorSweep(), tol=1e-10, maxiter=5000)).u ≈
            (L + 3I) \ f rtol = 1e-8
        robin = SchwarzDecomposition(L, doms; owner=parts, transmission=RobinTransmission(10.0))
        @test_throws ArgumentError refactor!(robin, 2L)
    end
end
