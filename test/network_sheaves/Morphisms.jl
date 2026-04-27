using Test
using CellularSheaves
using LinearAlgebra
using Graphs

@testset "SheafMorphism / ComplexMorphism" begin

    # Build a small 3-cycle; each vertex stalk is R^2, edge stalk is R^2
    # (identity restriction maps).
    n = 3
    d = 2
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

    F = cycle_sheaf(n, d, d -> Matrix{Float64}(I, d, d))
    G = cycle_sheaf(n, d, d -> Matrix{Float64}(I, d, d))

    # Identity morphism F -> F
    Vmap = collect(1:n)
    Emap = collect(1:n)
    Vmaps = [Matrix{Float64}(I, d, d) for _ in 1:n]
    Emaps = [Matrix{Float64}(I, d, d) for _ in 1:n]

    spec = SheafMorphism(Vmap, Emap, Vmaps, Emaps)
    @test length(spec.Vmap) == n
    @test length(spec.Emap) == n

    cm = ComplexMorphism(F, G, spec)
    @test is_morphism(cm, F, G)

    # Convenience helper constructor round-trips
    spec2 = make_sheaf_morphism_spec(Vmap, Emap, Vmaps, Emaps)
    cm2 = ComplexMorphism(F, G, spec2)
    @test is_morphism(cm2, F, G)

    # Explicit (Vmap, Emap, Vmaps, Emaps) overload
    cm3 = ComplexMorphism(F, G, Vmap, Emap, Vmaps, Emaps)
    @test is_morphism(cm3, F, G)

    # Non-morphism: perturb the vertex map
    bad_V = rand(size(cm.V)...)
    bad_cm = ComplexMorphism(bad_V, Array(cm.E))
    @test !is_morphism(bad_cm, F, G)

    # Composition: id ∘ id == id
    comp_spec = compose(spec, spec2)
    comp_cm   = compose(cm, cm2)
    @test is_morphism(comp_cm, F, G)

    # SheafMorphism -> ComplexMorphism composition matches direct compose
    direct_cm = ComplexMorphism(F, G, comp_spec)
    @test norm(direct_cm.V - comp_cm.V) < 1e-12
    @test norm(direct_cm.E - comp_cm.E) < 1e-12

    # -------------------------------------------------------------------
    # Two-step chain: F -> G -> H where H has 1-dimensional edge stalks
    # -------------------------------------------------------------------
    H = EuclideanSheaf{Float64}(Int[])
    for _ in 1:n
        add_vertex_stalk!(H, d)
    end
    for i in 1:n
        v1 = i
        v2 = (i % n) + 1
        rmH = zeros(1, d)
        rmH[1, 1] = 1.0
        add_sheaf_edge!(H, v1, v2, rmH, rmH)
    end

    edimsG = [get_edge_stalk(G, src(e), dst(e)) for e in edges(underlying_graph(G))]
    Emaps_GH = [zeros(1, edimsG[k]) for k in 1:length(edimsG)]
    for k in 1:length(edimsG)
        Emaps_GH[k][1, 1] = 1.0
    end
    specGH = make_sheaf_morphism_spec(collect(1:n), collect(1:n),
                                      [Matrix{Float64}(I, d, d) for _ in 1:n],
                                      Emaps_GH)
    cmGH = ComplexMorphism(G, H, specGH)
    @test is_morphism(cmGH, G, H)

    # Composite F -> G -> H via SheafMorphism compose
    comp_spec_FH = compose(spec, specGH)
    comp_cm_FH   = ComplexMorphism(F, H, comp_spec_FH)
    @test is_morphism(comp_cm_FH, F, H)

    # -------------------------------------------------------------------
    # Identity morphisms
    # -------------------------------------------------------------------

    id_spec = id(SheafMorphism, F)
    @test id_spec.Vmap == collect(1:n)
    @test id_spec.Emap == collect(1:n)
    @test all(id_spec.Vmaps[i] == Matrix{Float64}(I, d, d) for i in 1:n)
    @test all(id_spec.Emaps[k] == Matrix{Float64}(I, d, d) for k in 1:n)
    @test is_morphism(ComplexMorphism(F, F, id_spec), F, F)

    id_cm = id(ComplexMorphism, F)
    total_v = sum(vertex_stalks(F))
    total_e = sum(collect(values(edge_stalks(F))))
    @test id_cm.V == Matrix{Float64}(I, total_v, total_v)
    @test id_cm.E == Matrix{Float64}(I, total_e, total_e)
    @test is_morphism(id_cm, F, F)

    # id is a left and right unit for compose
    @test norm(compose(id_cm, cm).V - cm.V) < 1e-12
    @test norm(compose(cm, id(ComplexMorphism, G)).V - cm.V) < 1e-12
end
