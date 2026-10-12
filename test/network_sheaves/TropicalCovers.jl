using Test
using CellularSheaves
using Graphs
using SparseArrays
using Random

# Agents on a cycle of n cells: each step moves to a neighbouring cell or waits, at
# cost 1, except waiting at the goal (cost 0). Two agents collide when they end in
# the same cell or swap cells. Summary: blocks of `b` consecutive cells; two
# summaries are in contact when their blocks are at most one block apart, which is
# conservative for b ≥ 2 (out of contact, agents are at least b + 1 ≥ 3 cells apart).
function cycle_graph_costs(n, goal)
    rows, cols, vals = Int[], Int[], Float64[]
    for s in 1:n, s2 in (s, mod1(s + 1, n), mod1(s - 1, n))
        push!(rows, s); push!(cols, s2); push!(vals, s == s2 == goal ? 0.0 : 1.0)
    end
    return sparse(rows, cols, vals, n, n, min)
end
cycle_blocks(n, b) = GraphHomomorphism([cld(s, b) for s in 1:n], cld(n, b))
cycle_dist(i, j, m) = min(mod(i - j, m), mod(j - i, m))
cycle_collides(u, s, s2, v, t, t2) = s2 == t2 || (s2 == t && t2 == s)

# The joint graph of the agents (non-conflicting tuples, colliding moves removed).
function joint_graph(Ws, collides)
    K = length(Ws)
    n = [size(W, 1) for W in Ws]
    states = [x for x in Iterators.product((1:m for m in n)...)
              if all(!collides(a, x[a], x[a], b, x[b], x[b]) for a in 1:K for b in (a + 1):K)]
    index = Dict(x => i for (i, x) in enumerate(states))
    Wts = [sparse(transpose(W)) for W in Ws]               # out-moves by column, zero costs included
    arcs = Tuple{Int,Int,Float64}[]
    for (i, x) in enumerate(states)
        moves = [[(rowvals(Wt)[k], nonzeros(Wt)[k]) for k in nzrange(Wt, x[a])] for (a, Wt) in enumerate(Wts)]
        for combo in Iterators.product(moves...)
            any(collides(a, x[a], combo[a][1], b, x[b], combo[b][1]) for a in 1:K for b in (a + 1):K) && continue
            push!(arcs, (i, index[ntuple(a -> combo[a][1], K)], sum(c[2] for c in combo)))
        end
    end
    return states, index, arcs
end

@testset "TropicalCovers" begin
    @testset "simplicial complexes" begin
        K4 = clique_complex(complete_graph(4))
        @test length(simplices(K4, 0)) == 4
        @test length(simplices(K4, 1)) == 6
        @test length(simplices(K4, 2)) == 4
        @test isempty(simplices(K4, 3))
        @test length(simplices(clique_complex(complete_graph(4); dim = 3), 3)) == 1
        @test faces([1, 2, 3]) == [[2, 3], [1, 3], [1, 2]]
        K = SimplicialComplex(5, [[3, 1, 2]])
        @test simplices(K, 1) == [[1, 2], [1, 3], [2, 3]]
        @test simplices(K, 0) == [[v] for v in 1:5]
        P = path_graph(4)
        @test isempty(simplices(clique_complex(P), 2))
        @test simplices(clique_complex(P), 1) == [[1, 2], [2, 3], [3, 4]]
        @test_throws ArgumentError SimplicialComplex(3, [[1, 4]])
    end

    n, b = 10, 2
    goals = [7, 2, 4]
    Ws = [cycle_graph_costs(n, g) for g in goals]
    π = cycle_blocks(n, b)
    contact(u, k, v, l) = cycle_dist(k, l, cld(n, b)) <= 1

    @testset "the stalk is the pushforward of the joint graph ($K agents)" for K in (2, 3)
        P = ProductStalk(Tuple(1:K), Ws[1:K], fill(π, K), goals[1:K]; contact, collides = cycle_collides)
        G = stalk_graph(P)
        states, index, arcs = joint_graph(Ws[1:K], cycle_collides)
        μ = [coarse_state(P, x) for x in states]
        @test Set(μ) == Set(G.states)
        pushed = Dict{Tuple{Int,Int},Float64}()
        for (i, j, c) in arcs
            key = (G.index[μ[i]], G.index[μ[j]])
            pushed[key] = min(get(pushed, key, Inf), c)
        end
        stalk = Dict((i, j) => c for (i, j, c) in zip(findnz(G.W)...))
        @test stalk == pushed
        # Symmetry: swapping two agents with the same goal permutes the coarse states.
        Q = ProductStalk((1, 2), [Ws[1], Ws[1]], [π, π], (goals[1], goals[1]); contact, collides = cycle_collides)
        @test all(reverse(coarse_state(Q, (s, t))) == coarse_state(Q, (t, s)) for s in 1:n for t in 1:n if s != t)
        # The cost-to-go is a lower bound on the joint one, and exact when nothing is summarized.
        N = length(states)
        Wj = sparse(first.(arcs), getindex.(arcs, 2), last.(arcs), N, N, min)
        bj = fill(Inf, N); bj[index[Tuple(goals[1:K])]] = 0.0
        V = minplus_closure(Wj, bj)
        @test all(G.Δ[G.index[μ[i]]] <= V[i] + 1e-9 for i in 1:N)
        E = ProductStalk(Tuple(1:K), Ws[1:K], fill(π, K), goals[1:K]; contact = (u, k, v, l) -> true,
            collides = cycle_collides)
        GE = stalk_graph(E)
        @test all(GE.Δ[GE.index[coarse_state(E, x)]] ≈ V[i] for (i, x) in enumerate(states) if isfinite(V[i]))
    end

    @testset "restriction maps" begin
        P3 = ProductStalk((1, 2, 3), Ws, fill(π, 3), goals; contact, collides = cycle_collides)
        G3 = stalk_graph(P3)
        rng = MersenneTwister(1)
        for τ in ([1, 2], [1, 3], [2, 3])
            P2 = ProductStalk(Tuple(τ), Ws[τ], fill(π, 2), goals[τ]; contact, collides = cycle_collides)
            G2 = stalk_graph(P2)
            # The square commutes: restrict ∘ μ₃ = μ₂ ∘ forget.
            for _ in 1:200
                x = Tuple(randperm(rng, n)[1:3])
                @test restrict(P3, P2, coarse_state(P3, x)) == coarse_state(P2, Tuple(x[v] for v in τ))
            end
            # Every arc restricts to an arc that costs no more, so the face's cost-to-go
            # pulled back is a lower bound.
            arcs2 = Dict((i, j) => c for (i, j, c) in zip(findnz(G2.W)...))
            @test all(get(arcs2, (G2.index[restrict(P3, P2, G3.states[i])], G2.index[restrict(P3, P2, G3.states[j])]), Inf) <= c
                      for (i, j, c) in zip(findnz(G3.W)...))
            @test all(G2.Δ[G2.index[restrict(P3, P2, x)]] <= G3.Δ[i] + 1e-9 for (i, x) in enumerate(G3.states))
        end
        P1 = ProductStalk((1,), Ws[1:1], [π], goals[1:1]; contact, collides = cycle_collides)
        @test_throws ArgumentError restrict(P1, P3, coarse_state(P1, (1,)))
    end

    @testset "on-demand stalks and the cover" begin
        K = clique_complex(complete_graph(3))
        C1 = CoverSheaf(K, Ws, fill(π, 3), goals; contact, collides = cycle_collides, explicit_dim = 1)
        C2 = CoverSheaf(K, Ws, fill(π, 3), goals; contact, collides = cycle_collides, explicit_dim = 2)
        G3 = C2.explicit[[1, 2, 3]]
        rng = MersenneTwister(2)
        for _ in 1:20
            x = randperm(rng, n)[1:3]
            v, expanded, exact = section_value(C1, [1, 2, 3], x)
            @test exact
            @test v ≈ first(section_value(C2, [1, 2, 3], x))
            # Without the face heuristic the search agrees and expands at least as much.
            v0, e0, _ = stalk_value(C1.stalks[[1, 2, 3]], coarse_state(C1.stalks[[1, 2, 3]], Tuple(x)))
            @test v0 ≈ v
            @test e0 >= expanded
            # A capped search returns a lower bound.
            vc, _, ex = section_value(C1, [1, 2, 3], x; max_expanded = 2)
            @test vc <= v + 1e-9
        end
        # The cover bound is a lower bound on the joint cost, and triangles only tighten it.
        states, index, arcs = joint_graph(Ws, cycle_collides)
        N = length(states)
        Wj = sparse(first.(arcs), getindex.(arcs, 2), last.(arcs), N, N, min)
        bj = fill(Inf, N); bj[index[Tuple(goals)]] = 0.0
        V = minplus_closure(Wj, bj)
        Ke = clique_complex(complete_graph(3); dim = 1)
        Ce = CoverSheaf(Ke, Ws, fill(π, 3), goals; contact, collides = cycle_collides)
        for (i, x) in enumerate(states)
            isfinite(V[i]) || continue
            full, _ = cover_bound(C1, collect(x))
            edges_only, _ = cover_bound(Ce, collect(x))
            @test edges_only <= full <= V[i] + 1e-9
        end
    end

    @testset "contact growing with the dimension" begin
        wide(u, k, v, l) = cycle_dist(k, l, cld(n, b)) <= 2
        K = clique_complex(complete_graph(3))
        Cn = CoverSheaf(K, Ws, fill(π, 3), goals; contact, collides = cycle_collides, explicit_dim = 2)
        Cw = CoverSheaf(K, Ws, fill(π, 3), goals; contact = k -> k >= 2 ? wide : contact, collides = cycle_collides,
            explicit_dim = 2)
        Gn, Gw = Cn.explicit[[1, 2, 3]], Cw.explicit[[1, 2, 3]]
        @test length(Gw.states) > length(Gn.states)
        @test Cw.explicit[[1, 2]].states == Cn.explicit[[1, 2]].states        # the edges are unchanged
        # Restriction is defined and still gives lower bounds.
        P3 = Cw.stalks[[1, 2, 3]]
        for τ in ([1, 2], [1, 3], [2, 3])
            G2 = Cw.explicit[τ]
            @test all(G2.Δ[G2.index[restrict(P3, Cw.stalks[τ], x)]] <= Gw.Δ[i] + 1e-9 for (i, x) in enumerate(Gw.states))
        end
        # The wider triangle resolves at least as much: its values dominate on every joint state.
        for x in Iterators.product(1:n, 1:n, 1:n)
            allunique(x) || continue
            @test Gw.Δ[Gw.index[coarse_state(P3, x)]] >= Gn.Δ[Gn.index[coarse_state(Cn.stalks[[1, 2, 3]], x)]] - 1e-9
        end
        # A face that resolves earlier than its simplex has no restriction.
        Pn = Cn.stalks[[1, 2, 3]]
        Pw2 = ProductStalk((1, 2), Ws[1:2], fill(π, 2), goals[1:2]; contact = wide, collides = cycle_collides)
        far_apart = coarse_state(Pn, (1, 5, 8))
        @test_throws ArgumentError restrict(Pn, Pw2, far_apart)
    end

    @testset "stable summaries" begin
        # Two routes merged by the summary: 1 → 3 → 5 and 2 → 4 → 6 with {1, 2} and {3, 4}
        # summarized together. The quotient lets a walk start on one route and finish on
        # the other; the stable refinement separates the routes.
        W = sparse([1, 2, 3, 4, 5, 6], [3, 4, 5, 6, 5, 6], zeros(6), 6, 6)
        ψ = GraphHomomorphism([1, 1, 2, 2, 3, 4], 4)
        @test !is_stable(W, ψ)
        ρ = stable_refinement(W, ψ)
        @test is_stable(W, ρ)
        @test ρ.n_target == 6
        # A stable summary is its own refinement (up to renumbering).
        @test is_stable(W, ρ) && stable_refinement(W, ρ).n_target == ρ.n_target
        Wc = cycle_graph_costs(8, 1)
        @test stable_refinement(Wc, GraphHomomorphism(ones(Int, 8), 1); arcs = (s, s2, c) -> s != s2).n_target == 1
        # Only the selected arcs matter.
        Wm = sparse([1, 2, 1], [3, 4, 4], [0.0, 0.0, 1.0], 4, 4)
        ψm = GraphHomomorphism([1, 1, 2, 3], 3)
        @test is_stable(Wm, ψm; arcs = (s, s2, c) -> c == 0) == false
        @test stable_refinement(Wm, ψm; arcs = (s, s2, c) -> c > 0).n_target == 4
        # Random graphs: the refinement is stable and refines ψ.
        rng = MersenneTwister(3)
        for _ in 1:30
            m = 12
            Wr = sprand(rng, m, m, 0.25)
            ψr = GraphHomomorphism(rand(rng, 1:3, m), 3)
            ρr = stable_refinement(Wr, ψr)
            @test is_stable(Wr, ρr)
            @test all(ψr.vertex_map[s] == ψr.vertex_map[t] for s in 1:m for t in 1:m if ρr.vertex_map[s] == ρr.vertex_map[t])
        end
        # Futures seen through a coarser map: the two routes end in different fibres of
        # φ, so their starts are separated; states with the same future are not.
        φ = GraphHomomorphism([1, 1, 1, 1, 2, 3], 3)
        ϕ = future_refinement(W, ψ, φ)
        @test ϕ.vertex_map[1] != ϕ.vertex_map[2] && ϕ.vertex_map[3] != ϕ.vertex_map[4]
        @test ϕ.n_target == 6
        φ1 = GraphHomomorphism([1, 1, 1, 1, 2, 2], 2)
        @test future_refinement(W, ψ, φ1).n_target == 4
        # Random graphs: it refines ψ, and the states of a fibre reach the same fibres of φ.
        function reachable(Wr, φr, s)
            seen, stack = Set([s]), [s]
            while !isempty(stack)
                x = pop!(stack)
                for y in findall(!iszero, Wr[x, :])
                    y in seen || (push!(seen, y); push!(stack, y))
                end
            end
            return Set(φr.vertex_map[collect(seen)])
        end
        for _ in 1:30
            m = 12
            Wr = sprand(rng, m, m, 0.2)
            ψr = GraphHomomorphism(rand(rng, 1:3, m), 3)
            φr = GraphHomomorphism(rand(rng, 1:4, m), 4)
            ϕr = future_refinement(Wr, ψr, φr)
            same = [(s, t) for s in 1:m for t in 1:m if ϕr.vertex_map[s] == ϕr.vertex_map[t]]
            @test all(ψr.vertex_map[s] == ψr.vertex_map[t] for (s, t) in same)
            @test all(reachable(Wr, φr, s) == reachable(Wr, φr, t) for (s, t) in same)
            # Coarsest: states of one fibre of ψ with the same future share a fibre.
            @test all(ϕr.vertex_map[s] == ϕr.vertex_map[t] for s in 1:m for t in 1:m
                      if ψr.vertex_map[s] == ψr.vertex_map[t] && reachable(Wr, φr, s) == reachable(Wr, φr, t))
        end
        # Refining the summaries refines the stalks: their values can only grow.
        ψb = π
        ρb = [stable_refinement(W, ψb; arcs = (s, s2, c) -> s != s2) for W in Ws[1:2]]
        P = ProductStalk((1, 2), Ws[1:2], fill(ψb, 2), goals[1:2]; contact, collides = cycle_collides)
        blocks_of(ρ) = [ψb.vertex_map[findfirst(==(k), ρ.vertex_map)] for k in 1:ρ.n_target]
        bl = blocks_of.(ρb)
        Q = ProductStalk((1, 2), Ws[1:2], ρb, goals[1:2];
            contact = (u, k, v, l) -> contact(u, bl[u][k], v, bl[v][l]), collides = cycle_collides)
        GP, GQ = stalk_graph(P), stalk_graph(Q)
        for s in 1:n, t in 1:n
            s == t && continue
            @test GQ.Δ[GQ.index[coarse_state(Q, (s, t))]] >= GP.Δ[GP.index[coarse_state(P, (s, t))]] - 1e-9
        end
    end

    @testset "simplicial cover" begin
        tri = [([1, 2], 1.0), ([1, 3], 1.0), ([2, 3], 1.0)]
        @test simplicial_cover(3, tri) == 2.0
        @test simplicial_cover(3, [tri; ([1, 2, 3], 3.0)]) == 3.0
        @test simplicial_cover(3, [([1, 3], 2.5)]) == 3.0
        @test simplicial_cover(4, Tuple{Vector{Int},Float64}[]) == 0.0
        @test simplicial_cover(3, [([1, 2], Inf), ([2, 3], 1.0)]) == 1.0
        @test_throws ArgumentError simplicial_cover(17, tri)
    end

    @testset "argument checks" begin
        bad = sparse([1, 2], [2, 1], [-1.0, 1.0], 2, 2)
        @test_throws ArgumentError ProductStalk((1,), [bad], [GraphHomomorphism([1, 1], 1)], (1,);
            contact = (u, k, v, l) -> false, collides = cycle_collides)
        @test_throws ArgumentError ProductStalk((1, 2), Ws[1:2], [π], goals[1:2]; contact, collides = cycle_collides)
    end
end
