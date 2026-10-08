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

# Textbook reference implementations on a single global vector.
function reference_multiplicative(A, f, doms, u, sweeps)
    u = copy(u)
    for _ in 1:sweeps, d in doms
        u[d] += Matrix(A[d, d]) \ (f - A * u)[d]
    end
    return u
end

function reference_ras(A, f, doms, parts, u, sweeps)
    u = copy(u)
    for _ in 1:sweeps
        r = f - A * u
        du = zeros(length(u))
        for (i, d) in enumerate(doms)
            c = Matrix(A[d, d]) \ r[d]
            owned = parts[d] .== i
            du[d[owned]] = c[owned]
        end
        u += du
    end
    return u
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

    @testset "decomposition structure" begin
        @test length(dd.subdomains) == 9
        @test dd.graph == underlying_graph(overlap_sheaf(dd))
        @test has_edge(dd.graph, 1, 5)          # diagonal boxes overlap
        @test !has_edge(dd.graph, 1, 3)
        @test occursin("9 subdomains", sprint(show, dd))
        # Default owners hold each dof together with its whole stencil.
        dd_default = SchwarzDecomposition(A, doms)
        @test all(1:n) do k
            issubset(findall(!iszero, A[:, k]), dd_default.subdomains[dd_default.owner[k]])
        end
        @test schwarz_solve(dd_default, f; method=:parallel, tol=1e-10).u ≈ u_exact rtol = 1e-8
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
        r = schwarz_solve(dd, f; method=:multiplicative, tol=1e-10)
        @test r.converged
        @test r.u ≈ u_exact rtol = 1e-8
        # Pushing along restriction maps keeps every iterate a global section.
        @test maximum(r.disagreements) < 1e-8 * norm(u_exact)
        @test overlap_disagreement(dd, r.x) < 1e-8 * norm(u_exact)

        x = localize(dd, zeros(n))
        for _ in 1:3
            schwarz_step!(x, dd, f; method=:multiplicative)
        end
        @test glue(dd, x) ≈ reference_multiplicative(A, f, doms, zeros(n), 3)
    end

    @testset "parallel (Lions) Schwarz = RAS" begin
        u0 = randn(n)
        x = localize(dd, u0)
        for _ in 1:4
            schwarz_step!(x, dd, f; method=:parallel)
        end
        @test glue(dd, x) ≈ reference_ras(A, f, doms, parts, u0, 4)
        @test overlap_disagreement(dd, x) > 1e-6      # copies disagree mid-iteration

        r = schwarz_solve(dd, f; method=:parallel, tol=1e-10, maxiter=2000)
        @test r.converged
        @test r.u ≈ u_exact rtol = 1e-8
        @test last(r.disagreements) < 1e-6 * norm(u_exact)
        @test r.disagreements[1] ≈ 0 atol = 1e-12
    end

    @testset "more overlap converges faster" begin
        its = map(1:3) do δ
            dδ = SchwarzDecomposition(A, overlapping_subdomains(A, parts; overlap=δ); owner=parts)
            schwarz_solve(dδ, f; method=:multiplicative, tol=1e-8).iterations
        end
        @test its[1] > its[2] > its[3]
    end

    @testset "input validation" begin
        @test_throws ArgumentError SchwarzDecomposition(A, [doms[1]])            # not a cover
        @test_throws ArgumentError SchwarzDecomposition(A, [d for d in doms]; owner=fill(1, n))
        nonoverlap = overlapping_subdomains(A, parts; overlap=0)
        @test_throws ArgumentError SchwarzDecomposition(A, nonoverlap; owner=parts)
        @test_throws ArgumentError schwarz_solve(dd, f; method=:additive)
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
            schwarz_step!(x, dd, f; method=:multicolor)
            schwarz_step!(y, dd, f; method=:multiplicative, order=reduce(vcat, dd.colors))
        end
        @test Vector(x) ≈ Vector(y)

        r = schwarz_solve(dd, f; method=:multicolor, tol=1e-10)
        @test r.converged
        @test r.u ≈ u_exact rtol = 1e-8
    end

    # p×p boxes grouped into (p/2)×(p/2) blocks of 2×2 boxes.
    function box_aggregation(p)
        q = p ÷ 2
        return GraphHomomorphism([(cld(by, 2) - 1) * q + cld(bx, 2) for by in 1:p for bx in 1:p])
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

        r = schwarz_solve(dd, f; method=:multicolor, coarse=c, tol=1e-10)
        @test r.converged
        @test r.u ≈ u_exact rtol = 1e-8
        @test r.iterations < schwarz_solve(dd, f; method=:multicolor, tol=1e-10).iterations

        rp = schwarz_solve(dd, f; method=:parallel, coarse=c, tol=1e-10)
        @test rp.converged
        @test rp.u ≈ u_exact rtol = 1e-8

        x = localize(dd, randn(n))
        schwarz_step!(x, dd, f; method=:parallel)
        before = overlap_disagreement(dd, x)
        u_before = glue(dd, x)
        coarse_correct!(x, dd, c, f)
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
        @test vertex_stalks(pf) == length.(exact.decomposition.subdomains)
        @test coarse_dimension(exact) == sum(vertex_stalks(pf))
        @test ne(underlying_graph(pf)) == ne(exact.decomposition.graph)

        tc = TruncatedPushforwardCoarseSpace(dd8, hom)
        @test coarse_dimension(tc) == hom.n_target
        bases = all_fiber_bases(hom, overlap_sheaf(dd8))
        for h in 1:hom.n_target
            fiber = fiber_vertices(hom, h)
            local_section = reduce(vcat, [Vector(tc.basis[d, h]) for d in dd8.subdomains[fiber]])
            B = bases[h]
            @test norm(B * (B \ local_section) - local_section) < 1e-8   # lies in (φ_*F)(h)
        end
    end

    @testset "exact pushforward coarse space" begin
        point = GraphHomomorphism(ones(Int, length(doms)))
        c1 = ExactPushforwardCoarseSpace(dd, point)
        @test coarse_dimension(c1) == n
        r1 = schwarz_solve(dd, f; method=:multicolor, coarse=c1, tol=1e-10)
        @test r1.iterations == 1                 # pushforward to a point = direct solve
        @test r1.u ≈ u_exact rtol = 1e-10

        A16 = poisson2d(16)
        parts16 = box_partition(16, 4)
        dd16 = SchwarzDecomposition(A16, overlapping_subdomains(A16, parts16; overlap=1); owner=parts16)
        f16 = ones(16 * 16)
        exact = ExactPushforwardCoarseSpace(dd16, box_aggregation(4))
        tc = TruncatedPushforwardCoarseSpace(dd16, box_aggregation(4))
        its = Dict(name => schwarz_solve(dd16, f16; method=:multicolor, coarse=c, tol=1e-8).iterations
                   for (name, c) in (:none => nothing, :exact => exact, :truncated => tc))
        @test its[:exact] < its[:truncated] < its[:none]
        @test coarse_dimension(exact) > coarse_dimension(tc)
        @test schwarz_solve(dd16, f16; coarse=exact, tol=1e-10).u ≈ A16 \ f16 rtol = 1e-8
    end

    @testset "two-level scalability" begin
        its = map((2, 4, 6)) do p
            mp = 6p
            Ap = poisson2d(mp)
            pp = box_partition(mp, p)
            ddp = SchwarzDecomposition(Ap, overlapping_subdomains(Ap, pp; overlap=1); owner=pp)
            fp = ones(mp * mp)
            one_level = schwarz_solve(ddp, fp; method=:multicolor, tol=1e-8).iterations
            two_level = schwarz_solve(ddp, fp; method=:multicolor, tol=1e-8,
                coarse=TruncatedPushforwardCoarseSpace(ddp)).iterations
            (one_level, two_level)
        end
        one, two = first.(its), last.(its)
        @test all(two .< one)
        @test two[end] - two[1] < one[end] - one[1]
    end

    # Robin local matrix built independently of the package: drop the coupling
    # to outside dofs from the diagonal (algebraic Neumann) and add p.
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

    @testset "Robin transmission" begin
        h = 1 / (m + 1)
        p = optimized_robin_parameter(5h) / h           # overlap 2 ⇒ L = 5h
        ddr = SchwarzDecomposition(A, doms; owner=parts, robin=p)
        @test all(isempty, dd.interface)
        @test all(!isempty, ddr.interface)

        for meth in (:parallel, :multiplicative, :multicolor)
            x = localize(ddr, u_exact)
            schwarz_step!(x, ddr, f; method=meth)
            @test glue(ddr, x) ≈ u_exact                 # exact solution is still the fixed point
        end

        locals = [robin_local_matrix(A, d, p) for d in doms]
        u0 = randn(n)
        x = localize(ddr, u0)
        u = copy(u0)
        for _ in 1:3
            schwarz_step!(x, ddr, f; method=:parallel)
            r = f - A * u
            du = zeros(n)
            for (i, d) in enumerate(doms)
                c = locals[i] \ r[d]
                owned = parts[d] .== i
                du[d[owned]] = c[owned]
            end
            u += du
        end
        @test glue(ddr, x) ≈ u                           # optimized RAS

        x = localize(ddr, zeros(n))
        u = zeros(n)
        for _ in 1:2
            schwarz_step!(x, ddr, f; method=:multiplicative)
            for (i, d) in enumerate(doms)
                u[d] += locals[i] \ (f - A * u)[d]
            end
        end
        @test glue(ddr, x) ≈ u                           # optimized multiplicative Schwarz

        per_edge = SchwarzDecomposition(A, doms; owner=parts, robin=(i, j) -> p)
        @test per_edge.robin_shift == ddr.robin_shift
        @test_throws ArgumentError SchwarzDecomposition(A, doms; owner=parts, robin=0.0)
        @test_throws ArgumentError SchwarzDecomposition(A, doms; owner=parts, robin=(i, j) -> i == 1 ? -1.0 : p)

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
        dirichlet = schwarz_solve(SchwarzDecomposition(As, sdoms; owner=strips), fs; method=:parallel)
        optimized = schwarz_solve(SchwarzDecomposition(As, sdoms; owner=strips,
            robin=optimized_robin_parameter(3hs) / hs), fs; method=:parallel)
        @test dirichlet.converged && optimized.converged
        @test 3 * optimized.iterations < dirichlet.iterations
        @test optimized.u ≈ As \ fs rtol = 1e-6
    end

    @testset "Schwarz-preconditioned CG" begin
        A8 = poisson2d(8)
        parts8 = box_partition(8, 2)
        doms8 = overlapping_subdomains(A8, parts8; overlap=1)
        for robin in (nothing, 50.0)
            d8 = SchwarzDecomposition(A8, doms8; owner=parts8, robin)
            for coarse in (nothing, TruncatedPushforwardCoarseSpace(d8),
                           ExactPushforwardCoarseSpace(d8, GraphHomomorphism([1, 1, 2, 2])))
                P = schwarz_preconditioner(d8; coarse)
                Pm = reduce(hcat, [P * e for e in eachcol(Matrix(1.0I, 64, 64))])
                @test Pm ≈ Pm'
                @test eigmin(Symmetric(Pm)) > 0
                @test size(P) == (64, 64) && eltype(P) == Float64
            end
        end

        for coarse in (nothing, TruncatedPushforwardCoarseSpace(dd), ExactPushforwardCoarseSpace(dd, GraphHomomorphism(fill(1, 9))))
            r = schwarz_cg(dd, f; coarse, tol=1e-10)
            @test r.converged
            @test r.u ≈ u_exact rtol = 1e-8
            @test r.residuals[1] ≈ 1
            @test length(r.residuals) == r.iterations + 1
        end
        ddr = SchwarzDecomposition(A, doms; owner=parts, robin=optimized_robin_parameter(5 / (m + 1)) * (m + 1))
        @test schwarz_cg(ddr, f; tol=1e-10).u ≈ u_exact rtol = 1e-8
        A128 = poisson2d(128)
        parts128 = box_partition(128, 16)
        d128 = SchwarzDecomposition(A128, overlapping_subdomains(A128, parts128; overlap=1); owner=parts128)
        f128 = ones(128^2)
        @test schwarz_cg(d128, f128; coarse=TruncatedPushforwardCoarseSpace(d128)).iterations <
              schwarz_cg(d128, f128).iterations
        z = schwarz_cg(dd, zeros(n))
        @test z.converged && z.iterations == 0 && iszero(z.u)
    end
end
