using Test
using CellularSheaves
using LinearAlgebra
using Graphs

@testset "pushout_sheaf" begin

    n = 3
    d = 2

    # ------------------------------------------------------------------
    # Test 1: Pushout of identity span  F --id--> F <--id-- F
    # Pushout should be isomorphic to F (same stalk dims, same maps).
    # ------------------------------------------------------------------
    @testset "Identity span" begin
        F = constant_sheaf(cycle_graph(n), d)
        f = id(SheafMorphism, F)
        g = id(SheafMorphism, F)

        span = SheafSpan(f, g, F, F, F)
        P, iF, iG = pushout_sheaf(span)

        # Stalk dimensions should match F
        @test vertex_stalks(P) == vertex_stalks(F)
        @test edge_stalk_dimensions(P) == edge_stalk_dimensions(F)

        # Injections should be valid morphisms
        @test is_morphism(ComplexMorphism(F, P, iF), F, P)
        @test is_morphism(ComplexMorphism(F, P, iG), F, P)

        # Universal property: iF ∘ f ≈ iG ∘ g
        cm_f  = ComplexMorphism(F, F, f)
        cm_g  = ComplexMorphism(F, F, g)
        cm_iF = ComplexMorphism(F, P, iF)
        cm_iG = ComplexMorphism(F, P, iG)
        lhs = compose(cm_f, cm_iF)
        rhs = compose(cm_g, cm_iG)
        @test norm(lhs.V - rhs.V) < 1e-10
        @test norm(lhs.E - rhs.E) < 1e-10
    end

    # ------------------------------------------------------------------
    # Test 2: Pushout over the zero sheaf  F (+)_0 G = F (+) G
    # K is the zero sheaf (all stalks dimension 0).
    # ------------------------------------------------------------------
    @testset "Zero span (direct sum)" begin
        g_cycle = cycle_graph(n)
        F = constant_sheaf(g_cycle, d)
        G = constant_sheaf(g_cycle, d)
        K = zero_sheaf(g_cycle)

        # Zero morphisms K -> F and K -> G
        f = SheafMorphism(
            collect(1:n), collect(1:n),
            [zeros(Float64, d, 0) for _ in 1:n],
            [zeros(Float64, d, 0) for _ in 1:n],
        )
        g = SheafMorphism(
            collect(1:n), collect(1:n),
            [zeros(Float64, d, 0) for _ in 1:n],
            [zeros(Float64, d, 0) for _ in 1:n],
        )

        span = SheafSpan(f, g, K, F, G)
        P, iF, iG = pushout_sheaf(span)

        # Each vertex stalk of P should be F(v) (+) G(v) = R^(2d)
        @test all(vertex_stalks(P) .== 2d)
        @test all(edge_stalk_dimensions(P) .== 2d)

        # Both injections are valid morphisms
        @test is_morphism(ComplexMorphism(F, P, iF), F, P)
        @test is_morphism(ComplexMorphism(G, P, iG), G, P)
    end

    # ------------------------------------------------------------------
    # Test 3: Round-trip / universal property on a non-trivial span
    # K has d-dim stalks; f and g are non-trivial morphisms.
    # ------------------------------------------------------------------
    @testset "Universal property (non-trivial span)" begin
        g_cycle = cycle_graph(n)
        K = constant_sheaf(g_cycle, d)
        F = constant_sheaf(g_cycle, d)
        G = constant_sheaf(g_cycle, d)

        # f : K -> F is the identity
        f = id(SheafMorphism, K)

        # g : K -> G with rotation by pi/4 at each stalk
        theta = pi / 4
        R = [cos(theta) -sin(theta); sin(theta) cos(theta)]
        g = SheafMorphism(
            collect(1:n), collect(1:n),
            [R for _ in 1:n],
            [R for _ in 1:n],
        )

        span = SheafSpan(f, g, K, F, G)
        P, iF, iG = pushout_sheaf(span)

        # Both injections are valid morphisms F -> P and G -> P
        @test is_morphism(ComplexMorphism(F, P, iF), F, P)
        @test is_morphism(ComplexMorphism(G, P, iG), G, P)

        # Universal property: iF ∘ f ≈ iG ∘ g
        cm_f  = ComplexMorphism(K, F, f)
        cm_g  = ComplexMorphism(K, G, g)
        cm_iF = ComplexMorphism(F, P, iF)
        cm_iG = ComplexMorphism(G, P, iG)
        lhs = compose(cm_f, cm_iF)
        rhs = compose(cm_g, cm_iG)
        @test norm(lhs.V - rhs.V) < 1e-10
        @test norm(lhs.E - rhs.E) < 1e-10
    end

    # ------------------------------------------------------------------
    # Test 4: Dimension check  dim P(v) = dim F(v) + dim G(v) - rank(A_v)
    # ------------------------------------------------------------------
    @testset "Stalk dimension formula" begin
        g_cycle = cycle_graph(n)
        K = constant_sheaf(g_cycle, d)
        F = constant_sheaf(g_cycle, d)
        G = constant_sheaf(g_cycle, d)

        f = id(SheafMorphism, K)
        theta = pi / 3
        R = [cos(theta) -sin(theta); sin(theta) cos(theta)]
        g = SheafMorphism(
            collect(1:n), collect(1:n),
            [R for _ in 1:n],
            [R for _ in 1:n],
        )

        P, _, _ = pushout_sheaf(f, g, K, F, G)

        for v in 1:n
            A_v = vcat(f.Vmaps[v], -g.Vmaps[v])
            expected_dim = vertex_stalks(F)[v] + vertex_stalks(G)[v] - rank(A_v)
            @test vertex_stalks(P)[v] == expected_dim
        end

        for (k, _) in enumerate(edges(g_cycle))
            A_e = vcat(f.Emaps[k], -g.Emaps[k])
            expected_dim = edge_stalk_dimensions(F)[k] + edge_stalk_dimensions(G)[k] - rank(A_e)
            @test edge_stalk_dimensions(P)[k] == expected_dim
        end
    end

    # ------------------------------------------------------------------
    # Test 5: cycle_sheaf constructor with custom rm_pair_fn
    # ------------------------------------------------------------------
    @testset "cycle_sheaf constructor" begin
        # Asymmetric sheaf: src restriction is identity, dst is scaled by 2
        s = cycle_sheaf(n, d, (n, d, k) -> (Matrix{Float64}(I, d, d), 2.0 * Matrix{Float64}(I, d, d)))
        @test vertex_stalks(s) == fill(d, n)
        @test edge_stalk_dimensions(s) == fill(d, n)
    end

end
