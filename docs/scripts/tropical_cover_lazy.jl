# The simplicial cover bound at map sizes where no stalk can be tabulated, run with
# the standard MAPF benchmark protocol.
#
# tropical_cover_triangles.jl tabulates every edge stalk (a pushforward of the pair
# graph, about |cells|² states). The cover bound only needs each stalk's value at one
# section, the instance's start. Here that value is computed on demand, exactly:
#   - edges: the pairs whose tubes meet (delay slack 4, the rule of InteractionSheaf);
#     the value is the pair's least total delay, by joint A* with the decoupled
#     heuristic Σ dᵢ (ties broken toward deeper states, which matters at this size:
#     every joint state on the agents' own shortest paths ties). If the agents' own
#     shortest paths do not collide, the value is 0 and no search runs.
#   - triangles: the triangles of the interaction graph with at least one edge of
#     positive delay (at most LAZY_MAX_TRIANGLES per scenario), by joint A* the same
#     way. The obstruction o(σ) = Δ_σ − c_E(σ) is the triangle's value beyond the least
#     cover of its three edges. Leaving a triangle out only weakens the bound.
#   - the cover: the simplicial cover of the positive values, solved exactly on each
#     connected component (the integer program separates over them). A component with
#     more than 16 agents falls back to a greedy packing of disjoint simplices, still a
#     lower bound; the bound with triangles is taken as at least the cover of the edges.
#     A capped search contributes 0, also still a lower bound.
#   - the upper bound: prioritized planning in space-time, in the scenario's order.
#     Where it meets the cover bound, the plan is certified optimal.
#
# Protocol (Stern et al. 2019, §4.3): for each scenario file, the instance with N
# agents is the file's first N agents, and N grows until the method does not finish
# within the time limit. N grows in steps of LAZY_STEP rather than 1. The instances
# are nested, so everything carries over from one N to the next: the first N agents'
# prioritized plans do not depend on later agents, and pair and triangle values are
# reused. The time of an instance is the cumulative time to reach it, i.e. the time to
# solve it from scratch. Single-agent preprocessing (each agent's cost-to-go and tube)
# is included.
#
#   MAPF_MAP=~/data/movingai/wc3/losttemple.map MAPF_SCEN_DIR=~/data/movingai/wc3/scen-random \
#   LAZY_SCENARIOS=1:25 LAZY_STEP=50 LAZY_MAX_N=1000 LAZY_TIME_LIMIT=300 \
#   julia -t 16 --project=docs docs/scripts/tropical_cover_lazy.jl
# With LAZY_OPT=1 (small maps), the optimum is also found by independence detection,
# where it finishes, and the plans are checked collision-free.

isempty(get(ENV, "MAPF_MAP", "")) && error("set MAPF_MAP to a MovingAI map")
include(joinpath(@__DIR__, "tropical_multiagent_cover.jl"))

const SCENARIOS = eval(Meta.parse(get(ENV, "LAZY_SCENARIOS", "1:25")))::AbstractVector{Int}
const STEP = parse(Int, get(ENV, "LAZY_STEP", "50"))
const FIRST = parse(Int, get(ENV, "LAZY_FIRST", string(STEP)))
const MAX_N = parse(Int, get(ENV, "LAZY_MAX_N", "1000"))
const TIME_LIMIT = parse(Float64, get(ENV, "LAZY_TIME_LIMIT", "300"))
const PAIR_CAP = parse(Int, get(ENV, "LAZY_PAIR_CAP", "2000000"))
const TRI_CAP = parse(Int, get(ENV, "LAZY_TRI_CAP", "2000000"))
const MAX_TRIANGLES = parse(Int, get(ENV, "LAZY_MAX_TRIANGLES", "20000"))
const WITH_OPT = get(ENV, "LAZY_OPT", "0") == "1"
const SLACK = 4.0

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


# ---------------------------------------------------------------------------
# Prioritized planning in space-time, one agent at a time in the scenario's order (as
# in tropical_multiagent_cover.jl, ties broken toward later times). The reservations
# of the first agents do not depend on later ones, so the planner grows with N.
# ---------------------------------------------------------------------------

struct Reservations
    vertex::Set{Tuple{Int,Int}}                                      # (cell, time)
    edge::Set{Tuple{Int,Int,Int}}                                    # (from, to, time) of the move t → t + 1
    held::Dict{Int,Int}                                              # cell -> time from which it is held forever
    last_use::Dict{Int,Int}                                          # cell -> last time an earlier agent uses it
    T::Int
end
Reservations(T) = Reservations(Set{Tuple{Int,Int}}(), Set{Tuple{Int,Int,Int}}(), Dict{Int,Int}(), Dict{Int,Int}(), T)

function plan_next!(R::Reservations, A, start)
    free(c, t) = !((c, t) in R.vertex) && t < get(R.held, c, typemax(Int))
    g = Dict((start, 1) => 0.0)
    parent = Dict{Tuple{Int,Int},Tuple{Int,Int}}()
    heap = Tuple{Float64,Int,Int,Int}[]
    TA._push!(heap, (A.d[start], -1, start, 1))
    found = nothing
    while !isempty(heap)
        f, _, c, t = TA._pop!(heap)
        f > g[(c, t)] + A.d[c] + 1e-9 && continue
        if c == A.goal && t > get(R.last_use, c, 0)
            found = (c, t)
            break
        end
        t >= R.T && continue
        for c2 in NBRS[c]
            (free(c2, t + 1) && !((c2, c, t) in R.edge)) || continue
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
        push!(R.vertex, (c, t))
        R.last_use[c] = max(get(R.last_use, c, 0), t)
        t < length(plan) && push!(R.edge, (c, plan[t + 1], t))
    end
    R.held[plan[end]] = length(plan)
    return plan
end

# The first `n` agents of scenario file k, in the file's order (cells 0-based there).
function scenario_agents(k, n)
    base = splitext(basename(MAPF_MAP))[1]
    file = joinpath(get(ENV, "MAPF_SCEN_DIR", joinpath(dirname(MAPF_MAP), "scen-random")), "$(base)-random-$(k).scen")
    starts, goals = Int[], Int[]
    for l in readlines(file)[2:end]
        length(starts) == n && break
        f = split(l)
        s, g = (parse(Int, f[5]) + 1, parse(Int, f[6]) + 1), (parse(Int, f[7]) + 1, parse(Int, f[8]) + 1)
        (haskey(MAP.index, s) && haskey(MAP.index, g)) || error("$(basename(file)): agent outside the largest region")
        push!(starts, MAP.index[s]); push!(goals, MAP.index[g])
    end
    return starts, goals
end

const W_MAP = maximum(c[1] for c in MAP.cells)
const H_MAP = maximum(c[2] for c in MAP.cells)

println("MovingAI map $(basename(MAPF_MAP)): $(NF) free cells; scenarios $(SCENARIOS); N = $(FIRST):$(STEP):$(MAX_N) until $(TIME_LIMIT) s; $(Threads.nthreads()) threads")
rows = NamedTuple[]
for k in SCENARIOS
    starts_all, goals_all = scenario_agents(k, MAX_N)
    nmax = length(starts_all)
    agents = AgentStalk[]
    tubes, near = BitVector[], BitVector[]
    R = Reservations(4 * (W_MAP + H_MAP) + 2 * nmax)
    plans = Vector{Int}[]
    upper_sum, failed = 0.0, false
    edgeval = Dict{Tuple{Int,Int},Tuple{Float64,Int,Bool}}()
    trival = Dict{Vector{Int},Tuple{Float64,Int,Bool}}()
    adj = Dict{Int,Set{Int}}()
    tri_seen, tri_skipped = 0, 0
    elapsed = 0.0
    n = 0
    for N in FIRST:STEP:nmax
        step_time = @elapsed begin
            new = (n + 1):N
            # single agents: cost-to-go, tube
            built = Vector{AgentStalk}(undef, length(new))
            Threads.@threads for q in eachindex(new)
                built[q] = grid_agent(goals_all[new[q]])
            end
            append!(agents, built)
            tb, nb = tube_masks(built, starts_all[new])
            append!(tubes, tb); append!(near, nb)
            # prioritized planning, the new agents in order
            for i in new
                failed && break
                p = plan_next!(R, agents[i], starts_all[i])
                if p === nothing
                    failed = true
                else
                    push!(plans, p)
                    upper_sum += plan_cost(agents[i], p)
                end
            end
            # new pairs (i, j), j new
            newpairs = [(i, j) for j in new for i in 1:(j - 1) if any(near[i] .& tubes[j])]
            ev = Vector{Tuple{Float64,Int,Bool}}(undef, length(newpairs))
            Threads.@threads :dynamic for q in eachindex(newpairs)
                i, j = newpairs[q]
                ev[q] = joint_delay(agents[[i, j]], starts_all[[i, j]]; cap = PAIR_CAP)
            end
            for (q, (i, j)) in enumerate(newpairs)
                edgeval[(i, j)] = ev[q]
                push!(get!(adj, i, Set{Int}()), j)
                push!(get!(adj, j, Set{Int}()), i)
            end
            # new triangles (i, j, l), l new, with an edge of positive delay
            pos(a, b) = edgeval[minmax(a, b)][1] > 0
            newtris = Vector{Int}[]
            for l in new, i in get(adj, l, Set{Int}()), j in get(adj, l, Set{Int}())
                (i < j < l && j in get(adj, i, Set{Int}())) || continue
                tri_seen += 1
                (pos(i, j) || pos(i, l) || pos(j, l)) || continue
                if length(trival) + length(newtris) >= MAX_TRIANGLES
                    tri_skipped += 1
                    continue
                end
                push!(newtris, [i, j, l])
            end
            tv = Vector{Tuple{Float64,Int,Bool}}(undef, length(newtris))
            Threads.@threads :dynamic for q in eachindex(newtris)
                tv[q] = joint_delay(agents[newtris[q]], starts_all[newtris[q]]; cap = TRI_CAP)
            end
            for (q, σ) in enumerate(newtris)
                trival[σ] = tv[q]
            end
            # the covers
            edge_cons = [([i, j], v[1]) for ((i, j), v) in edgeval]
            tri_cons = [(σ, v[1]) for (σ, v) in trival]
            cover_e, _, _ = component_cover(N, edge_cons)
            cover_t, comp_t, packed_t = component_cover(N, [edge_cons; tri_cons])
            # Triangles can merge components past 16 agents, where the packing fallback is
            # weaker than the exact cover of the edges alone. Both are lower bounds: keep the larger.
            cover_t = max(cover_t, cover_e)
        end
        elapsed += step_time
        n = N
        decoupled = sum(agents[i].d[starts_all[i]] for i in 1:N)
        upper = failed ? Inf : upper_sum
        obstructed = count(trival) do (σ, v)
            e(a, b) = edgeval[(σ[a], σ[b])][1]
            v[1] - simplicial_cover(3, [([1, 2], e(1, 2)), ([1, 3], e(1, 3)), ([2, 3], e(2, 3))]) > 1e-9
        end
        @assert cover_e <= cover_t + 1e-9                                 # by construction
        @assert decoupled + cover_t <= upper + 1e-9 "the cover bound must not exceed a feasible plan"
        opt = NaN
        if WITH_OPT
            @assert failed || first_conflict(plans, collect(1:N)) === nothing "prioritized plans must be collision-free"
            pd, _, _ = independence_detection(agents, starts_all[1:N], m -> (c -> sum(agents[m[q]].d[c[q]] for q in eachindex(m))))
            opt = pd === nothing ? NaN : sum(plan_cost(agents[i], pd[i]) for i in 1:N)
            @assert isnan(opt) || decoupled + cover_t <= opt + 1e-9 <= upper + 2e-9
        end
        r = (map = splitext(basename(MAPF_MAP))[1], scenario = k, N = N, seconds = elapsed,
            pairs = length(edgeval), pairs_positive = count(v -> v[1] > 0, values(edgeval)),
            pairs_capped = count(v -> !v[3], values(edgeval)),
            pair_expanded_max = isempty(edgeval) ? 0 : maximum(v[2] for v in values(edgeval)),
            triangles = tri_seen, triangles_evaluated = length(trival), triangles_skipped = tri_skipped,
            triangles_positive = count(v -> v[1] > 0, values(trival)), triangles_obstructed = obstructed,
            triangles_capped = count(v -> !v[3], values(trival)),
            largest_component = comp_t, packed_components = packed_t,
            decoupled = decoupled, bound_edges = decoupled + cover_e, bound_sheaf = decoupled + cover_t,
            upper = upper, optimal = opt, certified = upper == decoupled + cover_t)
        push!(rows, r)
        @printf("scen %2d N=%4d %7.1f s  pairs %6d (%4d with delay, %d capped)  triangles %6d (%d evaluated, %d with delay, %d obstructed%s)  decoupled %.0f  edges %.0f  +tri %.0f  plan %.0f%s%s\n",
            k, N, elapsed, r.pairs, r.pairs_positive, r.pairs_capped, tri_seen, r.triangles_evaluated, r.triangles_positive,
            obstructed, tri_skipped > 0 ? ", $(tri_skipped) skipped" : "", decoupled, r.bound_edges, r.bound_sheaf, upper,
            r.certified ? " certified" : @sprintf(" (gap %.0f)", upper - r.bound_sheaf), isnan(opt) ? "" : @sprintf("  optimum %.0f", opt))
        flush(stdout)
        (elapsed > TIME_LIMIT || failed) && break
    end
end

println("\nsummary over scenarios (largest N finished within the time limit, and certification at each N)")
for N in sort(unique(r.N for r in rows))
    rs = filter(r -> r.N == N && r.seconds <= TIME_LIMIT, rows)
    isempty(rs) && continue
    gap = sum(r.upper - r.decoupled for r in rs if isfinite(r.upper); init = 0.0)
    cov(f) = gap == 0 ? "–" : @sprintf("%.0f%%", 100 * sum(f(r) - r.decoupled for r in rs if isfinite(r.upper); init = 0.0) / gap)
    @printf("N=%4d  scenarios %2d  certified %2d  plan − decoupled %6.0f  covered: edges %s, + triangles %s  mean time %.1f s\n",
        N, length(rs), count(r -> r.certified, rs), gap, cov(r -> r.bound_edges), cov(r -> r.bound_sheaf), mean(r.seconds for r in rs))
end
reached = [maximum((r.N for r in rows if r.scenario == k && r.seconds <= TIME_LIMIT); init = 0) for k in SCENARIOS]
println("largest N within $(TIME_LIMIT) s per scenario: ", reached)
mkpath(OUT)
tag = length(SCENARIOS) == 25 ? "" : "_s$(first(SCENARIOS))-$(last(SCENARIOS))"
open(joinpath(OUT, "tropical_cover_lazy_" * splitext(basename(MAPF_MAP))[1] * tag * ".csv"), "w") do io
    println(io, join(keys(first(rows)), ","))
    for r in rows
        println(io, join(values(r), ","))
    end
end
