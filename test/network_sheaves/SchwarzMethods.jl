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
end
