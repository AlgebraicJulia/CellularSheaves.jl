# The simplicial cover bound at map sizes where no stalk can be tabulated: the
# Warcraft III maps of MovingAI (512², 70–90k free cells), with MAPF scenarios made by
# make_mapf_scenarios.jl.
#
# tropical_cover_triangles.jl tabulates every edge stalk (a pushforward of the pair
# graph, about |cells|² states). The cover bound only needs each stalk's value at one
# section, the instance's start. Here that value is computed on demand, exactly:
#   - edges: the pairs whose tubes meet (delay slack 4, as in InteractionSheaf); the
#     value is the pair's least total delay, by joint A* with the decoupled heuristic
#     Σ dᵢ (ties broken toward deeper states, which matters at this size: every joint
#     state on the agents' own shortest paths ties). If the agents' own shortest paths
#     do not collide, the value is 0 and no search runs.
#   - triangles: the triangles of the interaction graph, solved by joint A* the same way.
#     The obstruction o(σ) = Δ_σ − c_E(σ) is the triangle's value beyond the least
#     cover of its three edges.
#   - the cover: the simplicial cover of the positive values, solved exactly on each
#     connected component (the integer program separates over them). A component with
#     more than 16 agents falls back to a greedy packing of disjoint simplices, which
#     is still a lower bound. A capped search contributes 0, also still a lower bound.
#   - the upper bound: prioritized planning in space-time, with a horizon that scales
#     with the paths. Where it meets the cover bound, the plan is certified optimal.
# No joint optimum is computed: at these sizes independence detection does not finish.
#
#   MAPF_MAP=~/data/movingai/wc3/losttemple.map MAPF_SCEN_DIR=~/data/movingai/wc3/scen-random \
#   LAZY_N=25,50,100 LAZY_INSTANCES=4 julia -t 16 --project=docs docs/scripts/tropical_cover_lazy.jl

isempty(get(ENV, "MAPF_MAP", "")) && error("set MAPF_MAP to a MovingAI map")
include(joinpath(@__DIR__, "tropical_multiagent_cover.jl"))

const SIZES = parse.(Int, split(get(ENV, "LAZY_N", "25,50,100"), ","))
const INSTANCES = parse(Int, get(ENV, "LAZY_INSTANCES", "4"))
const PAIR_CAP = parse(Int, get(ENV, "LAZY_PAIR_CAP", "2000000"))
const TRI_CAP = parse(Int, get(ENV, "LAZY_TRI_CAP", "2000000"))
const MAX_TRIANGLES = parse(Int, get(ENV, "LAZY_MAX_TRIANGLES", "5000"))
const SLACK = 4.0
const ORDERS = parse(Int, get(ENV, "LAZY_ORDERS", "3"))       # priority orders tried (5 replays tropical_cover_triangles.jl)

# The least total delay of `agents` from `starts`, ignoring all others: (value, expanded, exact).
function joint_delay(agents, starts; cap)
    m = length(agents)
    solo = [solo_plan(agents[k], starts[k]) for k in 1:m]
    first_conflict(solo, collect(1:m)) === nothing && return 0.0, 0, true
    h(c) = sum(agents[k].d[c[k]] for k in 1:m)
    goal = key(Tuple(A.goal for A in agents))
    s0 = key(Tuple(starts))
    g = Dict(s0 => 0.0)
    heap = Tuple{Float64,Float64,Int}[]
    TA._push!(heap, (h(Tuple(starts)), 0.0, s0))
    expanded = 0
    while !isempty(heap)
        f, ng, x = TA._pop!(heap)
        -ng > g[x] + 1e-9 && continue                               # stale entry
        x == goal && return g[x] - h(Tuple(starts)), expanded, true
        (expanded += 1) > cap && return 0.0, expanded, false
        c = unkey(x, m)
        for c2 in Iterators.product((NBRS[ci] for ci in c)...)
            any(collides(c, c2, k, l) for k in 1:m for l in (k + 1):m) && continue
            y = key(c2)
            gy = g[x] + sum(step_cost(agents[k].goal, c[k], c2[k]) for k in 1:m)
            if gy < get(g, y, Inf) - 1e-9
                g[y] = gy
                TA._push!(heap, (gy + h(c2), -gy, y))
            end
        end
    end
    return Inf, expanded, true                                       # no joint plan
end

# Each agent's tube and the cells within one king move of it, as bit masks over the map.
# Pair (i, j) interacts when j's tube meets i's neighbourhood: the rule of tubes_meet.
function tube_masks(agents, starts)
    tubes = [falses(NF) for _ in agents]
    near = [falses(NF) for _ in agents]
    Threads.@threads for i in eachindex(agents)
        for c in tube(agents[i], starts[i], SLACK)
            tubes[i][c] = true
            for dy in -1:1, dx in -1:1
                c2 = cell_at(MAP, c, (dx, dy))
                c2 == 0 || (near[i][c2] = true)
            end
        end
    end
    return tubes, near
end

# The simplicial cover, exact per connected component of the positive constraints.
function component_cover(N, cons)
    cons = [(σ, Δ) for (σ, Δ) in cons if isfinite(Δ) && Δ > 1e-9]
    parent = collect(1:N)
    root(v) = parent[v] == v ? v : (parent[v] = root(parent[v]))
    for (σ, _) in cons, v in σ[2:end]
        parent[root(v)] = root(σ[1])
    end
    groups = Dict{Int,Vector{Int}}()
    for (k, (σ, _)) in enumerate(cons)
        push!(get!(groups, root(σ[1]), Int[]), k)
    end
    total, largest, packed = 0.0, 0, 0
    for ks in values(groups)
        verts = sort!(unique(v for k in ks for v in cons[k][1]))
        largest = max(largest, length(verts))
        if length(verts) <= 16
            idx = Dict(v => i for (i, v) in enumerate(verts))
            total += simplicial_cover(length(verts), [([idx[v] for v in cons[k][1]], cons[k][2]) for k in ks])
        else                                                         # greedy packing: a lower bound
            packed += 1
            used = falses(N)
            for k in sort(ks; by = k -> -cons[k][2])
                σ, Δ = cons[k]
                any(used[σ]) && continue
                used[σ] .= true
                total += ceil(Δ - 1e-9)
            end
        end
    end
    return total, largest, packed
end

# Prioritized planning in space-time (as in tropical_multiagent_cover.jl) with horizon T,
# ties broken toward later times.
function prioritized_scaled(agents, starts, order, T)
    vertex = Set{Tuple{Int,Int}}()
    edge = Set{Tuple{Int,Int,Int}}()
    held = Dict{Int,Int}()
    last_use = Dict{Int,Int}()
    plans = Vector{Vector{Int}}(undef, length(agents))
    for i in order
        A = agents[i]
        free(c, t) = !((c, t) in vertex) && t < get(held, c, typemax(Int))
        g = Dict((starts[i], 1) => 0.0)
        parent = Dict{Tuple{Int,Int},Tuple{Int,Int}}()
        heap = Tuple{Float64,Int,Int,Int}[]
        TA._push!(heap, (A.d[starts[i]], -1, starts[i], 1))
        found = nothing
        while !isempty(heap)
            f, _, c, t = TA._pop!(heap)
            f > g[(c, t)] + A.d[c] + 1e-9 && continue
            if c == A.goal && t > get(last_use, c, 0)
                found = (c, t)
                break
            end
            t >= T && continue
            for c2 in NBRS[c]
                (free(c2, t + 1) && !((c2, c, t) in edge)) || continue
                gy = g[(c, t)] + (c == c2 == A.goal ? 0.0 : 1.0)
                if gy < get(g, (c2, t + 1), Inf)
                    g[(c2, t + 1)] = gy
                    parent[(c2, t + 1)] = (c, t)
                    TA._push!(heap, (gy + A.d[c2], -(t + 1), c2, t + 1))
                end
            end
        end
        found === nothing && return nothing
        path = [found]
        while haskey(parent, path[end])
            push!(path, parent[path[end]])
        end
        plan = [c for (c, _) in reverse(path)]
        for (t, c) in enumerate(plan)
            push!(vertex, (c, t))
            last_use[c] = max(get(last_use, c, 0), t)
            t < length(plan) && push!(edge, (c, plan[t + 1], t))
        end
        held[plan[end]] = length(plan)
        plans[i] = plan
    end
    return plans
end

function upper_bound_scaled(agents, starts, rng)
    N = length(agents)
    T = 2 * round(Int, maximum(agents[i].d[starts[i]] for i in 1:N)) + 2N + 50
    return minimum(1:ORDERS) do r
        p = prioritized_scaled(agents, starts, r == 1 ? collect(1:N) : randperm(rng, N), T)
        p === nothing ? Inf : sum(plan_cost(agents[i], p[i]) for i in 1:N)
    end
end

println("MovingAI map $(basename(MAPF_MAP)): $(NF) free cells; sizes $(SIZES), $(INSTANCES) instances each; $(Threads.nthreads()) threads")
rows = NamedTuple[]
rng = MersenneTwister(11)
for N in SIZES, inst in 1:INSTANCES
    starts, goals = instance(rng, N)
    length(starts) == N || error("the scenario has fewer than $N usable agents")
    t_agents = @elapsed begin
        agents = Vector{AgentStalk}(undef, N)
        Threads.@threads for i in 1:N
            agents[i] = grid_agent(goals[i])
        end
    end
    decoupled = sum(agents[i].d[starts[i]] for i in 1:N)
    t_upper = @elapsed upper = upper_bound_scaled(agents, starts, rng)
    tubes, near = tube_masks(agents, starts)
    pairs = [(i, j) for i in 1:N for j in (i + 1):N if any(near[i] .& tubes[j])]
    # edges
    t_pairs = @elapsed begin
        ev = Vector{Tuple{Float64,Int,Bool}}(undef, length(pairs))
        Threads.@threads :dynamic for k in eachindex(pairs)
            i, j = pairs[k]
            ev[k] = joint_delay(agents[[i, j]], starts[[i, j]]; cap = PAIR_CAP)
        end
    end
    edgeval = Dict(pairs[k] => ev[k][1] for k in eachindex(pairs))
    # triangles of the interaction graph
    adj = [Int[] for _ in 1:N]
    for (i, j) in pairs
        push!(adj[i], j)
    end
    tris = [[i, j, k] for (i, j) in pairs for k in adj[j] if k in adj[i]]
    ntris = length(tris)
    # positive edges first: those are where a triangle can matter most
    sort!(tris; by = σ -> -sum(edgeval[(σ[a], σ[b])] for (a, b) in ((1, 2), (1, 3), (2, 3))))
    tris = tris[1:min(end, MAX_TRIANGLES)]
    t_tris = @elapsed begin
        tv = Vector{Tuple{Float64,Int,Bool}}(undef, length(tris))
        Threads.@threads :dynamic for k in eachindex(tris)
            tv[k] = joint_delay(agents[tris[k]], starts[tris[k]]; cap = TRI_CAP)
        end
    end
    edge_cons = [([i, j], edgeval[(i, j)]) for (i, j) in pairs]
    tri_cons = [(tris[k], tv[k][1]) for k in eachindex(tris)]
    obstruction = [tv[k][1] - simplicial_cover(3, [([1, 2], edgeval[(σ[1], σ[2])]), ([1, 3], edgeval[(σ[1], σ[3])]),
                                                   ([2, 3], edgeval[(σ[2], σ[3])])]) for (k, σ) in enumerate(tris)]
    cover_e, comp_e, packed_e = component_cover(N, edge_cons)
    cover_t, comp_t, packed_t = component_cover(N, [edge_cons; tri_cons])
    @assert cover_e <= cover_t + 1e-9
    @assert decoupled + cover_t <= upper + 1e-9 "the cover bound must not exceed a feasible plan"
    push!(rows, (map = splitext(basename(MAPF_MAP))[1], N = N, instance = inst, decoupled = decoupled,
        pairs = length(pairs), pairs_positive = count(v -> v[1] > 0, ev), pairs_capped = count(v -> !v[3], ev),
        pair_expanded_max = isempty(ev) ? 0 : maximum(v[2] for v in ev),
        triangles = ntris, triangles_evaluated = length(tris), triangles_positive = count(v -> v[1] > 0, tv),
        triangles_obstructed = count(>(1e-9), obstruction), triangles_capped = count(v -> !v[3], tv),
        tri_expanded_max = isempty(tv) ? 0 : maximum(v[2] for v in tv),
        largest_component = comp_t, packed_components = packed_t,
        bound_edges = decoupled + cover_e, bound_sheaf = decoupled + cover_t, upper = upper,
        seconds_agents = t_agents, seconds_pairs = t_pairs, seconds_triangles = t_tris, seconds_upper = t_upper))
    r = rows[end]
    @printf("N=%3d #%d  pairs %4d (%d with delay, %d capped; max %d expanded, %.0f s)  triangles %d (%d evaluated, %d with delay, %d obstructed, %d capped; %.0f s)\n",
        N, inst, r.pairs, r.pairs_positive, r.pairs_capped, r.pair_expanded_max, t_pairs, ntris, length(tris),
        r.triangles_positive, r.triangles_obstructed, r.triangles_capped, t_tris)
    @printf("        decoupled %.0f  cover of edges %.0f  + triangles %.0f  prioritized %.0f  (ratio %.4f%s; largest component %d%s)\n",
        decoupled, r.bound_edges, r.bound_sheaf, upper, upper / r.bound_sheaf, upper == r.bound_sheaf ? ", certified optimal" : "",
        comp_t, packed_t > 0 ? ", $(packed_t) packed" : "")
    flush(stdout)
end

println("\nsummary (sums over instances)")
@printf("%4s %8s %10s %10s %12s %14s %18s %12s\n", "N", "pairs", "w/ delay", "triangles", "obstructed",
    "plan − decoupled", "covered: e / e+t", "certified")
for N in SIZES
    rs = filter(r -> r.N == N, rows)
    gap = sum(r.upper - r.decoupled for r in rs)
    cov(f) = gap == 0 ? "–" : @sprintf("%.0f%%", 100 * sum(f(r) - r.decoupled for r in rs) / gap)
    @printf("%4d %8d %10d %10d %12d %14.0f %18s %9d of %d\n", N, sum(r.pairs for r in rs), sum(r.pairs_positive for r in rs),
        sum(r.triangles for r in rs), sum(r.triangles_obstructed for r in rs), gap,
        cov(r -> r.bound_edges) * " / " * cov(r -> r.bound_sheaf), count(r -> r.upper == r.bound_sheaf, rs), length(rs))
end
mkpath(OUT)
open(joinpath(OUT, "tropical_cover_lazy_" * splitext(basename(MAPF_MAP))[1] * ".csv"), "w") do io
    println(io, join(keys(first(rows)), ","))
    for r in rows
        println(io, join(values(r), ","))
    end
end
