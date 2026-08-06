using Test
using LinearAlgebra
using CellularSheaves
using CellularSheaves.ControlSheaves.NestedSystems
using CellularSheaves.ControlSheaves.NestedDSL
using CellularSheaves.ControlSheaves.AgentControllers

# Helper fragments defined at top level so the composability tests can call them the way a user
# would: an ordinary Julia function taking ordinary Julia values and returning a fragment.
escort_ring(name::Symbol, m::Int, r::Real, tgt::Symbol) = @nested_system begin
    @team $(name) = ring(m; radius=r)
    @target $(tgt)
    @observe $(name) => $(tgt)
end

@testset "NestedDSL" begin

    @testset "fragments and terms" begin
        f = @nested_system begin
            @team a = ring(3; radius=1.0)
            @team b = path(2, 0.5)
            @target t1
            @observe a => t1
        end
        @test f isa SystemFragment
        @test length(f) == 4
        ts = terms(f)
        @test ts[1] isa TeamTerm && ts[1].name == :a && ts[1].kind == :ring
        @test ts[1].n_agents == 3 && ts[1].radius == 1.0 && ts[1].observers == [1]
        @test ts[2] isa TeamTerm && ts[2].kind == :path && ts[2].radius == 0.5
        @test ts[3] isa TargetTerm && ts[3].name == :t1
        @test ts[4] isa ObserveTerm && ts[4].target == :t1
        @test occursin("SystemFragment", sprint(show, f))

        @testset "merge concatenates in order" begin
            g = @nested_system begin
                @team c = star(4)
            end
            h = merge(f, g)
            @test length(h) == 5
            @test terms(h)[5] isa TeamTerm && terms(h)[5].name == :c
            @test length(merge(f)) == length(f)
            @test isempty(@nested_system begin end)
        end
    end

    @testset "@team surface forms" begin
        rad = 2.5
        n = 6
        kindvar = :clique
        f = @nested_system begin
            @team a = ring(n; radius=rad, observers=[1, 3])
            @team b = clique(3, 0.5)
            @team c = star(4; radius=rad)
            @team d = team(kindvar, 3; radius=1.0)
            @team $(Symbol(:e, 1)) = path(2)
        end
        ts = terms(f)
        @test ts[1].n_agents == 6 && ts[1].radius == 2.5 && ts[1].observers == [1, 3]
        @test ts[2].kind == :clique && ts[2].radius == 0.5
        @test ts[3].kind == :star && ts[3].radius == 2.5
        @test ts[4].kind == :clique
        @test ts[5].name == :e1 && ts[5].radius == 1.0   # default radius

        @test_throws NestedDSLError team_term(:x, :hexagon, 3, 1.0)
        @test_throws NestedDSLError team_term(:x, :ring, 0, 1.0)
        @test_throws NestedDSLError team_term(:x, :ring, 3, 1.0; observers=[4])
    end

    @testset "endpoint designators" begin
        M = zeros(4, 12)
        M[1, 1] = 1.0
        spec_expr = NestedSystems.translation_pin(3, 4, 2)
        f = @nested_system begin
            @team a = ring(3)
            @system s begin
                @team inner = ring(3)
                @team inner2 = ring(2)
            end
            @link a => s
            @link centroid(a) => centroid(s)
            @link project(a, 2) => project(s, :inner2)
            @link raw(a, M) => via(s, spec_expr)
        end
        ls = [t for t in terms(f) if t isa LinkTerm]
        @test ls[1].src.map isa NestedSystems.ProjectMember && ls[1].src.map.member == 1
        @test ls[2].src.map isa NestedSystems.Centroid && ls[2].dst.map isa NestedSystems.Centroid
        @test ls[3].src.map.member == 2 && ls[3].dst.map.member == :inner2
        @test ls[4].src.map isa NestedSystems.RawRestriction
        @test ls[4].src.map.M == M
        @test ls[4].dst.map === spec_expr

        @testset "dotted paths in @observe" begin
            g = @nested_system begin
                @team a = ring(3)
                @system s begin
                    @team inner = ring(3)
                end
                @target t1
                @observe centroid(s.inner) => t1
            end
            obs = only(t for t in terms(g) if t isa ObserveTerm)
            @test obs.system.path == [:s, :inner]
        end
    end

    @testset "lowering: structure, ranges, and paths" begin
        c = compile_nested_system(@nested_system begin
            @dim 4
            @affine true
            @team ringA = ring(4; radius=1.0)
            @system mid begin
                @team ringB = ring(5; radius=1.0)
                @team pod = path(2; radius=0.3)
                @link centroid(ringB) => pod
            end
            @target t1 t2
            @observe ringA => t1
            @observe centroid(mid.ringB) => t2
        end)

        @test c.spec.D == 4
        @test c.spec.affine
        @test length(c.spec.root.children) == 2
        @test c.spec.root.children[1] isa LeafTeam
        @test c.spec.root.children[2] isa RefinedSystem
        @test n_agents(c.spec) == 11

        @testset "agent ranges follow depth-first declaration order" begin
            @test agent_range(c, :ringA) == 1:4
            @test agent_range(c, :mid, :ringB) == 5:9
            @test agent_range(c, "mid.pod") == 10:11
            @test agent_range(c, :mid) == 5:11
            @test agent_range(c) == 1:11
            @test agent_vertices(c, :ringA) == c.tower.agent_vertices[1:4]
            @test_throws NestedDSLError agent_range(c, :nope)
        end

        @testset "observation paths resolve to indices" begin
            @test c.spec.observations[1].system_path == [1]
            @test c.spec.observations[1].target_index == 1
            @test c.spec.observations[2].system_path == [2, 1]
            @test c.spec.observations[2].target_index == 2
            @test c.spec.observations[2].system_map isa NestedSystems.Centroid
        end

        @testset "internal edges land on the declaring node" begin
            @test isempty(c.spec.root.internal_edges)
            mid = c.spec.root.children[2]
            e = only(mid.internal_edges)
            @test (e.src, e.dst) == (1, 2)
            @test e.src_map isa NestedSystems.Centroid
        end

        @testset "target name lookups" begin
            @test c.targets == [:t1, :t2]
            @test target_index(c, :t2) == 2
            @test_throws NestedDSLError target_index(c, :nope)
            v = target_vector(c, Dict(:t1 => [1.0, 0, 0, 1], :t2 => [2.0, 0, 0, 1]))
            @test v == [[1.0, 0, 0, 1], [2.0, 0, 0, 1]]
            @test_throws NestedDSLError target_vector(c, Dict(:t1 => [1.0, 0, 0, 1]))
            @test_throws NestedDSLError target_vector(c, Dict(:t1 => [1.0, 0, 0, 1],
                                                              :t2 => [1.0, 0, 0, 1],
                                                              :t9 => [1.0, 0, 0, 1]))
            # String keys normalize to the same Symbol lookups as keyword-style Symbol keys.
            v_str = target_vector(c, Dict("t1" => [1.0, 0, 0, 1], "t2" => [2.0, 0, 0, 1]))
            @test v_str == v
            @test_throws NestedDSLError target_vector(c, Dict("t1" => [1.0, 0, 0, 1]))
        end
    end

    @testset "deep nesting" begin
        c = compile_nested_system(@nested_system begin
            @system wing begin
                @system flight begin
                    @system element begin
                        @team alpha = ring(3; radius=0.5)
                        @team bravo = ring(3; radius=0.5)
                        @link alpha => bravo
                    end
                    @team charlie = ring(3; radius=0.5)
                    @link element => charlie
                end
                @team delta = ring(3; radius=0.5)
                @link flight => delta
            end
            @target t1
            @observe wing.flight.element.alpha => t1
        end)

        # root → wing → flight → element → alpha is five levels of tree, and the tower needs one
        # pushforward per level of collapse plus the finest level itself.
        @test c.tower.depth == 5
        @test n_agents(c.spec) == 12
        @test agent_range(c, "wing.flight.element.alpha") == 1:3
        @test agent_range(c, "wing.flight.element") == 1:6
        @test agent_range(c, "wing.flight") == 1:9
        @test agent_range(c, :wing) == 1:12
        @test only(c.spec.observations).system_path == [1, 1, 1, 1]

        # `project`-wired edges lower all the way to identity restrictions on `H_N`, so with a
        # single target the whole four-level tree has an exact zero-energy solution: everything
        # translated together. Both halves of that are checked — the pinned agent lands exactly
        # on its target, and moving the target moves every one of the twelve agents by the same
        # vector, through four levels of pushforward and lifting.
        tgt = [3.0, 0.0, 1.5, 1.0]
        q0 = solve_hierarchical(c.tower, [tgt])[end]
        @test q0[agent_vertices(c, "wing.flight.element.alpha")[1]] ≈ tgt atol = 1e-8
        shift = [1.0, -2.0, 0.0, 0.0]
        q1 = solve_hierarchical(c.tower, [tgt .+ shift])[end]
        for v in c.tower.agent_vertices
            @test q1[v] ≈ q0[v] .+ shift atol = 1e-8
        end
        @test sheaf_energy(c.tower.levels[end], reduce(vcat, q0)) < 1e-16

        @testset "four levels match a hand-built spec, centroid edges and all" begin
            alpha = LeafTeam(:alpha, :ring, 3, 0.5)
            bravo = LeafTeam(:bravo, :ring, 3, 0.5)
            element = RefinedSystem(:element, AbstractSystemNode[alpha, bravo],
                                    [SystemEdge(1, 2; src_map=centroid(), dst_map=centroid())])
            charlie = LeafTeam(:charlie, :ring, 3, 0.5)
            flight = RefinedSystem(:flight, AbstractSystemNode[element, charlie],
                                   [SystemEdge(1, 2; src_map=centroid(), dst_map=centroid())])
            delta = LeafTeam(:delta, :ring, 3, 0.5)
            wing = RefinedSystem(:wing, AbstractSystemNode[flight, delta],
                                 [SystemEdge(1, 2; src_map=centroid(), dst_map=centroid())])
            hand = build_sheaf_tower(NestedSystemSpec(
                RefinedSystem(:root, AbstractSystemNode[wing]), [TargetSpec(:t1)],
                [Observation([1, 1, 1, 1], 1)], 4, true))

            dsl = compile_nested_system(@nested_system begin
                @system wing begin
                    @system flight begin
                        @system element begin
                            @team alpha = ring(3; radius=0.5)
                            @team bravo = ring(3; radius=0.5)
                            @link centroid(alpha) => centroid(bravo)
                        end
                        @team charlie = ring(3; radius=0.5)
                        @link centroid(element) => centroid(charlie)
                    end
                    @team delta = ring(3; radius=0.5)
                    @link centroid(flight) => centroid(delta)
                end
                @target t1
                @observe wing.flight.element.alpha => t1
            end)
            @test dsl.tower.depth == hand.depth == 5
            @test solve_hierarchical(dsl.tower, [tgt])[end] ≈ solve_hierarchical(hand, [tgt])[end]
        end
    end

    @testset "composability" begin
        @testset "a fragment is position-independent" begin
            shallow = compile_nested_system(merge(escort_ring(:r1, 4, 1.0, :t1),
                                                  escort_ring(:r2, 4, 1.0, :t2)))
            deep = compile_nested_system(@nested_system begin
                @system outer begin
                    @system inner begin
                        @include escort_ring(:r1, 4, 1.0, :t1)
                        @include escort_ring(:r2, 4, 1.0, :t2)
                    end
                end
            end)
            @test [o.system_path for o in shallow.spec.observations] == [[1], [2]]
            @test [o.target_index for o in shallow.spec.observations] == [1, 2]
            @test deep.spec.observations[1].system_path == [1, 1, 1]
            @test deep.spec.observations[2].system_path == [1, 1, 2]
            @test deep.targets == shallow.targets == [:t1, :t2]
            @test agent_range(deep, "outer.inner.r2") == agent_range(shallow, :r2)
        end

        @testset "@system name = fragment" begin
            c = compile_nested_system(@nested_system begin
                @system alpha = escort_ring(:r, 3, 1.0, :t1)
                @system beta = escort_ring(:r, 3, 1.0, :t2)
                @link centroid(alpha) => centroid(beta)
            end)
            @test length(c.spec.root.children) == 2
            @test c.spec.root.children[1].name == :alpha
            @test agent_range(c, "beta.r") == 4:6
        end

        @testset "a Julia loop replaces DSL iteration" begin
            n, m = 4, 3
            looped = @nested_system begin
                for i in 1:n
                    @team $(Symbol(:ring, i)) = ring(m; radius=0.5)
                    @target $(Symbol(:t, i))
                    @observe $(Symbol(:ring, i)) => $(Symbol(:t, i))
                end
                for i in 1:n
                    @link centroid($(Symbol(:ring, i))) => centroid($(Symbol(:ring, mod1(i + 1, n))))
                end
            end
            written = @nested_system begin
                @team ring1 = ring(3; radius=0.5)
                @team ring2 = ring(3; radius=0.5)
                @team ring3 = ring(3; radius=0.5)
                @team ring4 = ring(3; radius=0.5)
                @target t1 t2 t3 t4
                @observe ring1 => t1
                @observe ring2 => t2
                @observe ring3 => t3
                @observe ring4 => t4
                @link centroid(ring1) => centroid(ring2)
                @link centroid(ring2) => centroid(ring3)
                @link centroid(ring3) => centroid(ring4)
                @link centroid(ring4) => centroid(ring1)
            end
            cl, cw = compile_nested_system(looped), compile_nested_system(written)
            tv = [[Float64(i), 0.0, 1.5, 1.0] for i in 1:4]
            @test solve_hierarchical(cl.tower, tv)[end] ≈ solve_hierarchical(cw.tower, tv)[end]
        end

        @testset "conditionals are Julia's too" begin
            for bridged in (true, false)
                c = compile_nested_system(@nested_system begin
                    @team a = ring(3)
                    @team b = ring(3)
                    @target t1
                    @observe a => t1
                    if bridged
                        @link centroid(a) => centroid(b)
                    end
                end)
                @test length(c.spec.root.internal_edges) == (bridged ? 1 : 0)
            end
        end
    end

    @testset "bindings" begin
        dyn = QuadrotorDynamics()
        K = zeros(3, 10)
        K_soft = ones(3, 10)
        p0 = [1.0, 2.0, 3.0]

        c = compile_nested_system(@nested_system begin
            @team a = ring(2; radius=1.0)
            @system mid begin
                @team b = ring(2; radius=1.0)
                @bind b K_lqr=K_soft
            end
            @target t1
            @observe a => t1
            @bind dynamics=dyn K_lqr=K
            @bind a[2] initial_position=p0
        end)

        resolved = resolve_dynamics(c.spec, c.bindings)
        @test length(resolved) == 4
        @test all(r -> r.dynamics === dyn, resolved)
        @test resolved[1].K_lqr == K            # a[1]: root default
        @test resolved[2].initial_position == p0
        @test resolved[3].K_lqr == K_soft       # mid.b: overridden locally
        @test resolved[4].K_lqr == K_soft

        @testset "a descendant path binds the same node as a local one" begin
            from_root = nested_bindings(@nested_system begin
                @team a = ring(2)
                @system mid begin
                    @team b = ring(2)
                end
                @target t1
                @observe a => t1
                @bind mid.b K_lqr=K_soft
            end)
            @test from_root.children[:mid].children[:b].K_lqr == K_soft
        end

        @testset "later declarations win field by field" begin
            b = nested_bindings(@nested_system begin
                @team a = ring(2)
                @target t1
                @observe a => t1
                @bind dynamics=dyn K_lqr=K
                @bind K_lqr=K_soft
            end)
            @test b.dynamics === dyn && b.K_lqr == K_soft
        end

        @testset "empty bindings are not stored" begin
            b = nested_bindings(@nested_system begin
                @team a = ring(2)
                @target t1
                @observe a => t1
                @bind a[1] initial_position=p0
            end)
            @test collect(keys(b.children)) == [:a]
            @test b.dynamics === nothing
            @test collect(keys(b.children[:a].agents)) == [1]
        end
    end

    @testset "validation errors" begin
        no_targets = @nested_system begin
            @team a = ring(3)
        end
        @test_throws NestedDSLError compile_nested_system(no_targets)

        dup_child = @nested_system begin
            @team a = ring(3)
            @team a = ring(4)
            @target t1
            @observe a => t1
        end
        @test_throws NestedDSLError validate_fragment(dup_child)

        dup_target = @nested_system begin
            @team a = ring(3)
            @target t1 t1
            @observe a => t1
        end
        @test_throws NestedDSLError validate_fragment(dup_target)

        unknown_target = @nested_system begin
            @team a = ring(3)
            @target t1
            @observe a => t9
        end
        @test_throws NestedDSLError validate_fragment(unknown_target)

        unknown_system = @nested_system begin
            @team a = ring(3)
            @target t1
            @observe nope => t1
        end
        @test_throws NestedDSLError validate_fragment(unknown_system)

        deep_link = @nested_system begin
            @team a = ring(3)
            @system s begin
                @team b = ring(3)
            end
            @target t1
            @observe a => t1
            @link a => s.b
        end
        @test_throws NestedDSLError validate_fragment(deep_link)

        self_link = @nested_system begin
            @team a = ring(3)
            @team b = ring(3)
            @target t1
            @observe a => t1
            @link a => a
        end
        @test_throws NestedDSLError validate_fragment(self_link)

        project_range = @nested_system begin
            @team a = ring(3)
            @team b = ring(3)
            @target t1
            @observe project(a, 9) => t1
        end
        @test_throws NestedDSLError validate_fragment(project_range)

        project_name_on_team = @nested_system begin
            @team a = ring(3)
            @target t1
            @observe project(a, :inner) => t1
        end
        @test_throws NestedDSLError validate_fragment(project_name_on_team)

        bad_agent = @nested_system begin
            @team a = ring(3)
            @target t1
            @observe a => t1
            @bind a[9] initial_position=[0.0]
        end
        @test_throws NestedDSLError validate_fragment(bad_agent)

        agent_on_system = @nested_system begin
            @system s begin
                @team a = ring(3)
            end
            @target t1
            @observe s.a => t1
            @bind s[1] initial_position=[0.0]
        end
        @test_throws NestedDSLError validate_fragment(agent_on_system)

        stray_position = @nested_system begin
            @team a = ring(3)
            @target t1
            @observe a => t1
            @bind a initial_position=[0.0]
        end
        @test_throws NestedDSLError validate_fragment(stray_position)

        conflicting_dim = @nested_system begin
            @dim 4
            @team a = ring(3)
            @target t1
            @observe a => t1
            @system s begin
                @dim 3
                @team b = ring(2)
            end
        end
        @test_throws NestedDSLError validate_fragment(conflicting_dim)

        empty_system = @nested_system begin
            @team a = ring(3)
            @system s = SystemFragment()
            @target t1
            @observe a => t1
        end
        @test_throws NestedDSLError validate_fragment(empty_system)

        @test_throws NestedDSLError bind_term(Symbol[], nothing)
    end

    @testset "defaults and settings" begin
        f = @nested_system begin
            @team a = ring(3)
            @target t1
            @observe a => t1
        end
        @test fragment_dim(f) == 4
        @test fragment_affine(f)
        s = nested_spec(f)
        @test s.D == 4 && s.affine

        ## A non-affine (linear) stalk cannot carry a translation offset, so `build_escort_topology`
        ## requires radius 0 whenever `affine=false` — hence the degenerate radius here.
        g = @nested_system begin
            @dim 3
            @affine false
            @team a = ring(3; radius=0.0)
            @target t1
            @observe a => t1
        end
        @test nested_spec(g).D == 3
        @test !nested_spec(g).affine
        @test nested_tower(g) isa SheafTower
    end

    @testset "matches a hand-written spec" begin
        # The `centroid_formation_tracking` topology, built both ways. The DSL must be a pure
        # convenience over the type API: same tree, same edges, same numbers.
        D = 4
        ringA = LeafTeam(:ringA, :ring, 4, 1.0)
        ringB = LeafTeam(:ringB, :ring, 5, 1.0)
        mid = RefinedSystem(:mid, AbstractSystemNode[ringA, ringB],
                            [SystemEdge(1, 2; src_map=centroid(), dst_map=centroid())])
        root = RefinedSystem(:root, AbstractSystemNode[mid])
        observations = vcat(
            [Observation([1, 1], 1; system_map=redundant_pin(4, D, k)) for k in 1:2:4],
            [Observation([1, 2], 2; system_map=redundant_pin(5, D, k)) for k in 1:2:5])
        hand = build_sheaf_tower(NestedSystemSpec(root, [TargetSpec(:t1), TargetSpec(:t2)],
                                                  observations, D, true))

        dsl = compile_nested_system(@nested_system begin
            @dim 4
            @system mid begin
                @team ringA = ring(4; radius=1.0)
                @team ringB = ring(5; radius=1.0)
                @link centroid(ringA) => centroid(ringB)
            end
            @target t1 t2
            for k in 1:2:4
                @observe via(mid.ringA, redundant_pin(4, 4, k)) => t1
            end
            for k in 1:2:5
                @observe via(mid.ringB, redundant_pin(5, 4, k)) => t2
            end
        end)

        tv = [[1.0, 0.5, 1.5, 1.0], [-1.0, -0.5, 1.5, 1.0]]
        @test dsl.tower.depth == hand.depth
        @test length(dsl.tower.agent_vertices) == length(hand.agent_vertices)
        @test solve_hierarchical(dsl.tower, tv)[end] ≈ solve_hierarchical(hand, tv)[end]
        @test solve_direct(dsl.tower, tv) ≈ solve_direct(hand, tv)
        gap_dsl = approximation_gap(dsl.tower, tv)
        gap_hand = approximation_gap(hand, tv)
        @test gap_dsl.hierarchical ≈ gap_hand.hierarchical
        @test gap_dsl.gap >= -1e-8
    end

    @testset "escort_problem" begin
        dyn = QuadrotorDynamics()
        Ad, Bd = AgentControllers.discrete_matrices(dyn, 0.05)
        K = AgentControllers.solve_dare(Ad, Bd,
                                        Matrix(Diagonal([500.0, 500.0, 500.0, 150.0, 150.0,
                                                         100.0, 100.0, 100.0, 5.0, 5.0])),
                                        Matrix(Diagonal([0.005, 0.005, 0.005])))
        c = compile_nested_system(@nested_system begin
            @team ring = ring(4; radius=1.0)
            @target t1
            @observe ring => t1
            @bind dynamics=dyn K_lqr=K
        end)
        prob = escort_problem(c, Dict(:t1 => t -> [cos(t), sin(t), 1.5, 1.0]); dt=0.05, steps=5)
        @test prob isa NestedEscortProblem
        @test prob.steps == 5
        res = run_nested_escort_simulation(prob)
        @test length(res.sim_data) == 5
        @test length(res.qstar_history[1]) == 4

        ff = escort_problem(c, Dict(:t1 => t -> [cos(t), sin(t), 1.5, 1.0]);
                            velocities=Dict(:t1 => t -> [-sin(t), cos(t), 0.0, 0.0]),
                            dt=0.05, steps=3)
        @test ff.target_velocities !== nothing
        @test length(run_nested_escort_simulation(ff; use_feedforward=true).sim_data) == 3
    end
end
