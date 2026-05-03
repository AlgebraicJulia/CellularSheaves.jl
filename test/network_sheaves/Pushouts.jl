using Test
using CellularSheaves
using LinearAlgebra
using Graphs

@testset "pushout_sheaf" begin

    # ------------------------------------------------------------------
    # Helper: build a 3-cycle sheaf with given stalk dim and rm function
    # ------------------------------------------------------------------
    function cycle_sheaf(n, d, rm_fn)
        s = EuclideanSheaf{Float64}(Int[])
        for _ in 1:n
            add_vertex_stalk!(s, d)
        end
        for i in 1:n
            v1 = i
            v2 = (i % n) + 1
            rm = rm_fn(d)
            add_sheaf_edge!(s, v1, v2, rm, rm)
        end
        return s
    end

    n = 3
    d = 2

    # ------------------------------------------------------------------
    # Test 1: Pushout of identity span  F --id--> F <--id-- F
    # Pushout should be isomorphic to F (same stalk dims, same maps).
    # ------------------------------------------------------------------
    @testset "Identity span" begin
        F = cycle_sheaf(n, d, d -> Matrix{Float64}(I, d, d))
        f = id(SheafMorphism, F)
        g = id(SheafMorphism, F)

        P, iF, iG = pushout_sheaf(f, g, F, F, F)

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
    # Test 2: Pushout over the zero sheaf  F ⊕_0 G = F ⊕ G
    # K is the zero sheaf (all stalks dimension 0).
    # ------------------------------------------------------------------
    @testset "Zero span (direct sum)" begin
        F = cycle_sheaf(n, d, d -> Matrix{Float64}(I, d, d))
        G = cycle_sheaf(n, d, d -> Matrix{Float64}(I, d, d))

        # Build the zero sheaf on the same graph
        K = EuclideanSheaf{Float64}(Int[])
        for _ in 1:n
            add_vertex_stalk!(K, 0)
        end
        g_graph = underlying_graph(F)
        for e in edges(g_graph)
            rm = zeros(Float64, 0, 0)
            add_sheaf_edge!(K, src(e), dst(e), rm, rm)
        end

        # Zero morphisms K → F and K → G
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

        P, iF, iG = pushout_sheaf(f, g, K, F, G)

        # Each vertex stalk of P should be F(v) ⊕ G(v) = R^(2d)
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
        # Build K, F, G all with d-dim vertex stalks and identity restriction maps
        K = cycle_sheaf(n, d, d -> Matrix{Float64}(I, d, d))
        F = cycle_sheaf(n, d, d -> Matrix{Float64}(I, d, d))
        G = cycle_sheaf(n, d, d -> Matrix{Float64}(I, d, d))

        # f : K → F is the identity
        f = id(SheafMorphism, K)

        # g : K → G — use a non-trivial vertex map (e.g. first coordinate embedding)
        # For each vertex map K(v) → G(v), use the rotation by π/4
        θ = π / 4
        R = [cos(θ) -sin(θ); sin(θ) cos(θ)]
        g = SheafMorphism(
            collect(1:n), collect(1:n),
            [R for _ in 1:n],
            [R for _ in 1:n],
        )

        P, iF, iG = pushout_sheaf(f, g, K, F, G)

        # Both injections are valid morphisms F → P and G → P
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
        K = cycle_sheaf(n, d, d -> Matrix{Float64}(I, d, d))
        F = cycle_sheaf(n, d, d -> Matrix{Float64}(I, d, d))
        G = cycle_sheaf(n, d, d -> Matrix{Float64}(I, d, d))

        f = id(SheafMorphism, K)
        θ = π / 3
        R = [cos(θ) -sin(θ); sin(θ) cos(θ)]
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

        g_graph = underlying_graph(K)
        for (k, _) in enumerate(edges(g_graph))
            A_e = vcat(f.Emaps[k], -g.Emaps[k])
            expected_dim = edge_stalk_dimensions(F)[k] + edge_stalk_dimensions(G)[k] - rank(A_e)
            @test edge_stalk_dimensions(P)[k] == expected_dim
        end
    end

end
