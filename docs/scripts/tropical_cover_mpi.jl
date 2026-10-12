# The interaction sheaf over MPI: the cover's tasks spread over processes and nodes.
#
# On one node the explicit edge stalks are bound by memory bandwidth (each streams
# tens of millions of arcs), so more threads stop helping at about 4x. The tasks
# are independent and what they return is small, so they spread over nodes:
#
#   1. Every rank builds the cheap shared structure itself, deterministically: the
#      instance, the agents' summaries (vertex data), the interaction complex and the
#      product stalk of every simplex (no arcs).
#   2. Edges are dealt to ranks, longest first by the size of their box of summaries;
#      each rank builds its edges over its threads (`stalk_cost_to_go`).
#   3. Each edge's cost-to-go is broadcast from its rank; every rank enumerates the
#      edge's states and index again (`StalkGraph(P, Δ)`), which is cheap.
#   4. Triangles are dealt to ranks and searched with the pulled-back edge values;
#      only their values are reduced.
#   5. Rank 0 glues them: the simplicial cover.
#
#   srun --mpi=pmix in-sandbox julia -t T --project=docs docs/scripts/tropical_cover_mpi.jl n agents instances
#
# Environment: COVER_BLOCK (block size, default 2), COVER_TRI_RADIUS (triangle contact
# radius, default 1). The instances are those of `tropical_cover_scaling.jl`.

using MPI
MPI.Init()
const COMM = MPI.COMM_WORLD
const RANK = MPI.Comm_rank(COMM)
const NRANKS = MPI.Comm_size(COMM)

include(joinpath(@__DIR__, "tropical_multiagent_cover.jl"))
using Graphs: SimpleGraph, add_edge!

const NA = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 8
const NINST = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 2
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

function shared_structure(agents, starts)
    N = length(agents)
    summ = Vector{Tuple{GraphHomomorphism,Vector{NTuple{2,Int}}}}(undef, N)
    Threads.@threads :dynamic for i in 1:N
        summ[i] = agent_summary(agents[i])
    end
    tubes = [tube(agents[i], starts[i], 4.0) for i in 1:N]
    g = SimpleGraph(N)
    for i in 1:N, j in (i + 1):N
        tubes_meet(MAP, tubes[i], tubes[j]) && add_edge!(g, i, j)
    end
    blocks = [last(s) for s in summ]
    within(r) = (u, k, v, l) -> cheb(blocks[u][k] .- blocks[v][l]) <= r
    complex = clique_complex(g)
    stalks = cover_stalks(complex, [A.delays for A in agents], [first(s) for s in summ], [A.goal for A in agents];
        contact = k -> within(k >= 2 ? TRI_RADIUS : 1), collides = grid_collides)
    return complex, stalks, summ
end

# Longest processing time first: tasks by decreasing cost, each to the least loaded rank.
function deal(costs, nranks)
    load = zeros(nranks)
    owner = zeros(Int, length(costs))
    for i in sortperm(costs; rev = true)
        r = argmin(load)
        owner[i] = r - 1
        load[r] += costs[i]
    end
    return owner
end

# The cost-to-go of the given stalks, over this rank's threads (one arc scratch per thread).
function edge_tables(stalks, mine)
    out = Vector{Vector{Float64}}(undef, length(mine))
    pool = Channel{ArcScratch}(Threads.nthreads())
    foreach(_ -> put!(pool, ArcScratch()), 1:Threads.nthreads())
    Threads.@threads :dynamic for k in eachindex(mine)
        sc = take!(pool)
        try
            out[k] = stalk_cost_to_go(stalks[mine[k]]; scratch = sc)
        finally
            put!(pool, sc)
        end
    end
    return Dict(zip(mine, out))
end

# The whole pipeline from the shared structure to the bound; times are maxima over ranks.
function distributed_bound(agents, starts)
    MPI.Barrier(COMM)
    t0 = MPI.Wtime()
    complex, stalks, summ = shared_structure(agents, starts)
    edges, tris = simplices(complex, 1), simplices(complex, 2)
    owner = deal([Float64(summ[e[1]][1].n_target) * summ[e[2]][1].n_target for e in edges], NRANKS)
    mine = [e for (e, r) in zip(edges, owner) if r == RANK]
    t1 = MPI.Wtime()
    local_tables = edge_tables(stalks, mine)
    t2 = MPI.Wtime()
    MPI.Barrier(COMM)                                   # everyone's edges done: the slowest rank sets the pace
    t3 = MPI.Wtime()
    tables = Dict{Vector{Int},Vector{Float64}}()
    for (e, r) in zip(edges, owner)
        len = [r == RANK ? length(local_tables[e]) : 0]
        MPI.Bcast!(len, COMM; root = r)
        buf = r == RANK ? local_tables[e] : Vector{Float64}(undef, len[1])
        MPI.Bcast!(buf, COMM; root = r)
        tables[e] = buf
    end
    t4 = MPI.Wtime()
    graphs = Vector{StalkGraph}(undef, length(edges))
    Threads.@threads :dynamic for k in eachindex(edges)
        graphs[k] = StalkGraph(stalks[edges[k]], tables[edges[k]])
    end
    C = CoverSheaf(complex, stalks, Dict(zip(edges, graphs)))
    t5 = MPI.Wtime()
    tri_owner = deal(ones(length(tris)), NRANKS)
    mytris = findall(==(RANK), tri_owner)
    vals = fill(-Inf, length(tris))
    capped = zeros(Int, length(tris))
    Threads.@threads :dynamic for k in mytris
        v, _, ok = section_value(C, tris[k], starts; max_expanded = TRI_CAP)
        vals[k] = v
        capped[k] = ok ? 0 : 1
    end
    MPI.Allreduce!(vals, max, COMM)
    MPI.Allreduce!(capped, +, COMM)
    t6 = MPI.Wtime()
    bound = simplicial_cover(length(agents), [[(e, C.explicit[e].Δ[C.explicit[e].index[coarse_state(stalks[e], Tuple(starts[e]))]])
        for e in edges]; [(σ, v) for (σ, v) in zip(tris, vals)]])
    t7 = MPI.Wtime()
    worst(x) = MPI.Allreduce(x, max, COMM)
    best(x) = MPI.Allreduce(x, min, COMM)
    times = (shared = worst(t1 - t0), edges_max = worst(t2 - t1), edges_min = best(t2 - t1), wait = worst(t3 - t2),
        exchange = worst(t4 - t3), assemble = worst(t5 - t4), triangles = worst(t6 - t5), glue = worst(t7 - t6),
        total = worst(t7 - t0))
    states = [length(C.explicit[e].states) for e in edges]
    return bound, times, length(edges), length(tris), sum(capped), states, sum(8 * s for s in states)
end

T = Threads.nthreads()
nodes = MPI.Allreduce(Int(MPI.Comm_rank(MPI.Comm_split_type(COMM, MPI.COMM_TYPE_SHARED, RANK)) == 0), +, COMM)   # one leader per node
RANK == 0 && println("grid $(n)², wall with two 1-cell corridors: $(NF) free cells; $NA agents; block $BLOCK; ",
    "triangle contact radius $TRI_RADIUS; $NRANKS ranks on $nodes nodes × $T threads")
rng = MersenneTwister(2026 + n)
rows = NamedTuple[]
for inst in 1:NINST
    starts, goals = instance(rng, NA)
    agents = [grid_agent(g) for g in goals]
    upper = prioritized_cost(agents, starts, rng)       # every rank, to keep the instance stream in step
    decoupled = sum(agents[i].d[starts[i]] for i in 1:NA)
    if inst == 1                                        # compile everything on a three-agent cover first
        distributed_bound(agents[1:3], starts[1:3])
    end
    GC.gc()
    bound, t, ne, nt, capped, states, bytes = distributed_bound(agents, starts)
    if RANK == 0
        push!(rows, (n = n, block = BLOCK, agents = NA, instance = inst, ranks = NRANKS, nodes = nodes, threads = T,
            edges = ne, triangles = nt, edge_states_mean = sum(states) / ne, exchanged_mib = bytes / 2^20,
            seconds_shared = t.shared, seconds_edges_max = t.edges_max, seconds_edges_min = t.edges_min,
            seconds_wait = t.wait, seconds_exchange = t.exchange, seconds_assemble = t.assemble,
            seconds_triangles = t.triangles, seconds_glue = t.glue, seconds_total = t.total, triangles_capped = capped,
            decoupled = decoupled, bound = decoupled + bound, upper = upper))
        @printf("#%d  %d edges (%.0f states each), %d triangles; %d ranks × %d threads\n", inst, ne, sum(states) / ne, nt, NRANKS, T)
        @printf("    total %.1f s: shared %.1f, edges %.1f (fastest rank %.1f), exchange %.1f (%.0f MiB), assemble %.1f, triangles %.1f, glue %.2f\n",
            t.total, t.shared, t.edges_max, t.edges_min, t.exchange, bytes / 2^20, t.assemble, t.triangles, t.glue)
        @printf("    decoupled %.0f  sheaf bound %.0f  prioritized %.0f  certified ratio %.3f%s\n",
            decoupled, decoupled + bound, upper, upper / (decoupled + bound), capped > 0 ? " ($capped triangles capped)" : "")
    end
end

if RANK == 0
    mkpath(OUT)
    open(joinpath(OUT, "tropical_cover_mpi_n$(n)_b$(BLOCK)_r$(NRANKS)_t$(T).csv"), "w") do io
        println(io, join(keys(first(rows)), ","))
        for r in rows
            println(io, join(values(r), ","))
        end
    end
end
MPI.Finalize()
