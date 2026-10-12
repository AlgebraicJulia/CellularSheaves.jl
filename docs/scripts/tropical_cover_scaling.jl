# The interaction sheaf in parallel, on larger maps.
#
# The cover's work is a set of independent tasks, one per simplex: an edge stalk is
# built from its two agents' graphs and summaries alone, and a triangle is searched
# from its own coarse state, reading only the finished edge stalks. The agents'
# summaries (vertex data) are independent too. Only the gluing, a small integer
# program on the agents, is sequential. So the construction should scale with the
# number of threads up to the size of the largest single task.
#
# For each instance this reports the stalk sizes against the pair graph, the time to
# build the cover and to evaluate the bound on `Threads.nthreads()` threads, and the
# certificate against prioritized planning. With `profile`, it also times every task
# alone, giving the total work and the longest task, hence the ideal speedup.
#
#   julia -t T --project=docs docs/scripts/tropical_cover_scaling.jl n agents instances [profile]
#
# Environment: COVER_BLOCK (block size, default 2), COVER_TRI_RADIUS (triangle contact
# radius in blocks, default 1). Summaries are (Reeb arc, block) refined by the Reeb arcs
# ahead (`future_refinement`); edges keep agents exact within one block.

include(joinpath(@__DIR__, "tropical_multiagent_cover.jl"))
using Graphs: SimpleGraph, add_edge!, complete_graph

const NA = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 8
const NINST = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 2
const PROFILE = length(ARGS) >= 4 && ARGS[4] == "profile"
const TRI_RADIUS = parse(Int, get(ENV, "COVER_TRI_RADIUS", "1"))
const TRI_CAP = 2_000_000

grid_collides(u, s, s2, v, t, t2) = s2 == t2 || (s2 == t && t2 == s)

# An agent's summary: (Reeb arc, block), refined by the Reeb arcs ahead.
function agent_summary(A)
    π, B = far_summary(MAP, A, BLOCK)
    ρ = future_refinement(A.delays, π, A.arcs; arcs = (s, s2, c) -> c == 0)
    parent = zeros(Int, ρ.n_target)
    parent[ρ.vertex_map] .= π.vertex_map
    return ρ, B[parent]
end

function build_cover(agents, summ, pairs)
    g = SimpleGraph(length(agents))
    for (i, j) in pairs
        add_edge!(g, i, j)
    end
    blocks = [last(s) for s in summ]
    within(r) = (u, k, v, l) -> cheb(blocks[u][k] .- blocks[v][l]) <= r
    return CoverSheaf(clique_complex(g), [A.delays for A in agents], [first(s) for s in summ], [A.goal for A in agents];
        contact = k -> within(k >= 2 ? TRI_RADIUS : 1), collides = grid_collides, keep_arcs = false)
end

T = Threads.nthreads()
println("grid $(n)², wall with two 1-cell corridors: $(NF) free cells; $NA agents; block $BLOCK; ",
    "triangle contact radius $TRI_RADIUS; $T threads")
rng = MersenneTwister(2026 + n)
rows = NamedTuple[]
for inst in 1:NINST
    starts, goals = instance(rng, NA)
    agents = [grid_agent(g) for g in goals]
    if inst == 1                                        # compile everything on a three-agent cover first
        s3 = [agent_summary(A) for A in agents[1:3]]
        C3 = build_cover(agents[1:3], s3, [(1, 2), (1, 3), (2, 3)])
        cover_bound(C3, starts[1:3]; max_expanded = 1000)
    end
    upper = prioritized_cost(agents, starts, rng)
    decoupled = sum(agents[i].d[starts[i]] for i in 1:NA)
    # Vertex data: each agent's summary, in parallel.
    summ = Vector{Tuple{GraphHomomorphism,Vector{NTuple{2,Int}}}}(undef, NA)
    t_summ = @elapsed Threads.@threads :dynamic for i in 1:NA
        summ[i] = agent_summary(agents[i])
    end
    tubes = [tube(agents[i], starts[i], 4.0) for i in 1:NA]
    pairs = [(i, j) for i in 1:NA for j in (i + 1):NA if tubes_meet(MAP, tubes[i], tubes[j])]
    GC.gc()
    t_build = @elapsed C = build_cover(agents, summ, pairs)
    t_bound = @elapsed (bound, vals) = cover_bound(C, starts; max_expanded = TRI_CAP)
    edges, tris = simplices(C.complex, 1), simplices(C.complex, 2)
    edge_states = [length(C.explicit[e].states) for e in edges]
    finite = [count(isfinite, A.d) for A in agents]
    pair_states = [finite[e[1]] * finite[e[2]] for e in edges]
    capped = count(σ -> !vals[σ][3], tris)
    work_edges, max_edge, work_tris, max_tri = NaN, NaN, NaN, NaN
    if PROFILE                                          # every task alone, one after another
        te = [@elapsed(stalk_graph(C.stalks[e])) for e in edges]
        tt = [@elapsed(section_value(C, σ, starts; max_expanded = TRI_CAP)) for σ in tris]
        work_edges, max_edge = sum(te), maximum(te)
        work_tris, max_tri = sum(tt; init = 0.0), maximum(tt; init = 0.0)
    end
    push!(rows, (n = n, block = BLOCK, agents = NA, instance = inst, threads = T, edges = length(edges),
        triangles = length(tris), summaries_mean = mean(s[1].n_target for s in summ), cells = NF,
        edge_states_mean = mean(edge_states), pair_states_mean = mean(pair_states),
        seconds_summaries = t_summ, seconds_build = t_build, seconds_bound = t_bound,
        work_edges = work_edges, max_edge = max_edge, work_triangles = work_tris, max_triangle = max_tri,
        triangles_capped = capped, triangle_expansions_mean = isempty(tris) ? 0.0 : mean(vals[σ][2] for σ in tris),
        decoupled = decoupled, bound = decoupled + bound, upper = upper))
    @printf("#%d  %d edges, %d triangles; summaries %.0f per agent (%d cells); edge stalk %.0f states vs pair graph %.0f (%.0f%%)\n",
        inst, length(edges), length(tris), rows[end].summaries_mean, NF, mean(edge_states), mean(pair_states),
        100 * mean(edge_states ./ pair_states))
    @printf("    %d threads: summaries %.1f s, build %.1f s, bound %.1f s (triangle expansions mean %.0f%s)\n",
        T, t_summ, t_build, t_bound, rows[end].triangle_expansions_mean, capped > 0 ? ", $capped capped" : "")
    if PROFILE
        @printf("    tasks alone: edges %.1f s total, longest %.1f s; triangles %.1f s total, longest %.2f s\n",
            work_edges, max_edge, work_tris, max_tri)
    end
    @printf("    decoupled %.0f  sheaf bound %.0f  prioritized %.0f  certified ratio %.3f (decoupled %.3f)\n",
        decoupled, decoupled + bound, upper, upper / (decoupled + bound), upper / decoupled)
end

mkpath(OUT)
open(joinpath(OUT, "tropical_cover_scaling_n$(n)_b$(BLOCK)_t$(T).csv"), "w") do io
    println(io, join(keys(first(rows)), ","))
    for r in rows
        println(io, join(values(r), ","))
    end
end
