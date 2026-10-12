# The interaction sheaf as a sheaf on a 2-complex, built by the package
# (`TropicalCovers`): the clique complex of the interaction graph, with a symmetric
# product stalk on every edge and triangle.
#
#   - Every agent is summarized by its own map, (Reeb arc of its cost-to-go, block of
#     the map), unless it is in contact with another agent of the simplex (blocks at
#     most one apart), and then it keeps its exact state. No agent of a simplex is
#     privileged; the edge stalks are symmetric too.
#   - Edge stalks are built and solved once. Triangle stalks are searched on demand
#     from the start, with the edge values pulled back along the restriction maps as
#     the heuristic.
#   - The cover bound is the simplicial cover of these values.
#
# Compared with exact pair and triple delays (joint search of two and three agents),
# the asymmetric edge stalks of `tropical_pair_sheaf.jl`, and the optimum where
# independence detection finds it.
#
#   julia --project=docs docs/scripts/tropical_cover_triangles.jl [n] [instances]

include(joinpath(@__DIR__, "tropical_multiagent_cover.jl"))
using Graphs

const TRI_CAP = 2_000_000                             # expansions per triangle search

grid_collides(u, s, s2, v, t, t2) = s2 == t2 || (s2 == t && t2 == s)

function cover_sheaf(agents, pairs)
    N = length(agents)
    g = SimpleGraph(N)
    for (i, j) in pairs
        add_edge!(g, i, j)
    end
    summaries = [far_summary(MAP, A, BLOCK) for A in agents]
    blocks = [last(s) for s in summaries]
    contact(u, k, v, l) = cheb(blocks[u][k] .- blocks[v][l]) <= 1
    return CoverSheaf(clique_complex(g), [A.delays for A in agents], [first(s) for s in summaries],
        [A.goal for A in agents]; contact, collides = grid_collides)
end

# The symmetric edge stalk equals the pushforward of the pair graph.
function check_edge(C, agents, e)
    P, G = C.stalks[e], C.explicit[e]
    Ai, Aj = agents[e[1]], agents[e[2]]
    Si = length(Ai.pos)
    W = pair_graph(Ai, Aj)
    pushed = Dict{Tuple{Int,Int},Float64}()
    for (x, y, c) in zip(findnz(W)...)
        s, t = (x - 1) % Si + 1, (x - 1) ÷ Si + 1
        s2, t2 = (y - 1) % Si + 1, (y - 1) ÷ Si + 1
        key = (G.index[coarse_state(P, (s, t))], G.index[coarse_state(P, (s2, t2))])
        pushed[key] = min(get(pushed, key, Inf), c)
    end
    stalk = Dict((i, j) => c for (i, j, c) in zip(findnz(G.W)...))
    @assert stalk == pushed "the symmetric edge stalk is the pushforward of the pair graph"
    println("check: the symmetric edge stalk $e is the pushforward of the pair graph ($(length(G.states)) coarse states, ",
        "pair graph $(count(isfinite, Ai.d) * count(isfinite, Aj.d)))")
end

println("grid $(n)², wall with two 1-cell corridors: $(NF) free cells; block $BLOCK")
rows = NamedTuple[]
rng = MersenneTwister(11)                             # the instances of tropical_multiagent_cover.jl
for N in (4, 6, 8), inst in 1:NI
    starts, goals = instance(rng, N)
    agents = [grid_agent(g) for g in goals]
    upper = prioritized_cost(agents, starts, rng)     # keeps the instance stream aligned
    IS = InteractionSheaf(MAP, agents, starts; block = BLOCK)
    decoupled = sum(agents[i].d[starts[i]] for i in 1:N)
    t_build = @elapsed C = cover_sheaf(agents, IS.pairs)
    N == 4 && inst == 1 && check_edge(C, agents, first(simplices(C.complex, 1)))
    edges, tris = simplices(C.complex, 1), simplices(C.complex, 2)
    # Edges: symmetric stalk, asymmetric stalks, exact.
    sym = Dict(e => first(section_value(C, e, starts)) for e in edges)
    asym = Dict(e => pair_bound(MAP, IS.stalks[(e[1], e[2])], IS.stalks[(e[2], e[1])], agents[e[1]], agents[e[2]],
        starts[e[1]], starts[e[2]]) for e in edges)
    t_pairs = @elapsed exact2 = Dict(e => subset_delay(agents, starts, e) for e in edges)
    # Triangles: on-demand stalk search, exact.
    tri_stats = Dict{Vector{Int},Tuple{Float64,Int,Bool,Float64}}()
    for σ in tris
        t = @elapsed (v, ex, ok) = section_value(C, σ, starts; max_expanded = TRI_CAP)
        tri_stats[σ] = (v, ex, ok, t)
    end
    t_tri = @elapsed exact3 = Dict(σ => subset_delay(agents, starts, σ) for σ in tris)
    face(σ, f) = maximum(f[τ] for τ in faces(σ))
    beyond = [σ for σ in tris if exact3[σ] > face(σ, exact2) + 1e-9]
    caught = [σ for σ in beyond if tri_stats[σ][1] > face(σ, sym) + 1e-9]
    # Covers.
    cover_sym_edges = decoupled + simplicial_cover(N, [(e, sym[e]) for e in edges])
    cover_sheaf_bound = decoupled + simplicial_cover(N, [[(e, sym[e]) for e in edges]; [(σ, tri_stats[σ][1]) for σ in tris]])
    cover_exact = decoupled + simplicial_cover(N, [[(e, exact2[e]) for e in edges]; [(σ, exact3[σ]) for σ in tris]])
    # The optimum, where independence detection finds it.
    pd, _, _ = N <= 6 ? independence_detection(agents, starts, members -> (S = InteractionSheaf(MAP, agents[members],
        starts[members]; block = BLOCK); c -> joint_bound(S, c))) : (nothing, 0, 0)
    opt = pd === nothing ? NaN : sum(plan_cost(agents[i], pd[i]) for i in 1:N)
    @assert all(sym[e] <= exact2[e] + 1e-9 for e in edges)
    @assert all(tri_stats[σ][1] <= exact3[σ] + 1e-9 for σ in tris)
    @assert decoupled <= cover_sym_edges <= cover_sheaf_bound <= cover_exact + 1e-9
    @assert isnan(opt) || cover_exact <= opt + 1e-9
    sizes = [length(C.explicit[e].states) for e in edges]
    push!(rows, (N = N, instance = inst, edges = length(edges), triangles = length(tris),
        edges_exact_sym = count(e -> sym[e] >= exact2[e] - 1e-9, edges),
        edges_exact_asym = count(e -> asym[e] >= exact2[e] - 1e-9, edges),
        triangles_beyond = length(beyond), triangles_caught = length(caught),
        triangles_exact = count(σ -> tri_stats[σ][1] >= exact3[σ] - 1e-9, tris),
        triangles_capped = count(σ -> !tri_stats[σ][3], tris),
        expanded_mean = isempty(tris) ? 0.0 : mean(tri_stats[σ][2] for σ in tris),
        expanded_max = isempty(tris) ? 0 : maximum(tri_stats[σ][2] for σ in tris),
        seconds_triangles_sheaf = sum((tri_stats[σ][4] for σ in tris); init = 0.0), seconds_triangles_exact = t_tri,
        seconds_build = t_build, edge_states_mean = mean(sizes),
        decoupled = decoupled, cover_sym_edges = cover_sym_edges, cover_sheaf = cover_sheaf_bound,
        cover_exact = cover_exact, optimal = opt, upper = upper))
    @printf("N=%d #%d  edges %2d (symmetric stalk exact on %2d, asymmetric on %2d; %.0f states each, built in %.0f s)\n",
        N, inst, length(edges), rows[end].edges_exact_sym, rows[end].edges_exact_asym, mean(sizes), t_build)
    @printf("        triangles %2d: %d beyond their edges, %d of them seen by the stalk; stalk exact on %d; expansions mean %.0f max %d%s (%.0f s; exact triples %.0f s)\n",
        length(tris), length(beyond), length(caught), rows[end].triangles_exact, rows[end].expanded_mean,
        rows[end].expanded_max, rows[end].triangles_capped > 0 ? ", $(rows[end].triangles_capped) capped" : "",
        rows[end].seconds_triangles_sheaf, t_tri)
    for σ in beyond
        println("          $σ: exact $(exact3[σ]) (edges $( [exact2[τ] for τ in faces(σ)] )), stalk $(tri_stats[σ][1]) (edges $([sym[τ] for τ in faces(σ)]))")
    end
    @printf("        bound: decoupled %5.1f  cover of edges %5.1f  cover of edges + triangle stalks %5.1f  with exact triangles %5.1f  optimum %5.1f  prioritized %5.1f\n",
        decoupled, cover_sym_edges, cover_sheaf_bound, cover_exact, opt, upper)
end

println("\nsummary (means over instances)")
@printf("%3s %22s %22s %26s %24s %30s\n", "N", "edges exact: sym/asym", "triangles beyond/seen", "triangle expansions (max)",
    "gap closed: edges/sheaf", "ratio prioritized/bound: sheaf")
for N in (4, 6, 8)
    rs = filter(r -> r.N == N, rows)
    solved = filter(r -> !isnan(r.optimal) && r.optimal > r.decoupled, rs)
    gap(f) = isempty(solved) ? "–" : @sprintf("%.2f", mean((f(r) - r.decoupled) / (r.optimal - r.decoupled) for r in solved))
    @printf("%3d %22s %22s %26s %24s %30.3f\n", N,
        @sprintf("%d/%d of %d", sum(r.edges_exact_sym for r in rs), sum(r.edges_exact_asym for r in rs), sum(r.edges for r in rs)),
        @sprintf("%d/%d of %d", sum(r.triangles_beyond for r in rs), sum(r.triangles_caught for r in rs), sum(r.triangles for r in rs)),
        @sprintf("%.0f (%d)", mean(r.expanded_mean for r in rs), maximum(r.expanded_max for r in rs)),
        gap(r -> r.cover_sym_edges) * " / " * gap(r -> r.cover_sheaf),
        mean(r.upper / r.cover_sheaf for r in rs))
end

mkpath(OUT)
open(joinpath(OUT, "tropical_cover_triangles.csv"), "w") do io
    println(io, join(keys(first(rows)), ","))
    for r in rows
        println(io, join(values(r), ","))
    end
end
