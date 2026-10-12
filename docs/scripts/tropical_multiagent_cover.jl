# N-agent planning on the cover of the interaction sheaf (`tropical_pair_sheaf.jl`), at
# sizes where the joint graph cannot be built. Three pieces:
#
#   - Lower bounds (oplax side), from the edge stalks over the interaction graph of agents
#     whose tubes meet, all added to Σᵢ dᵢ:
#       * matching: the heaviest matching of pair bounds;
#       * vertex cover: the weighted vertex cover of pair bounds (Li, Felner et al. 2019);
#       * cover gauge (`tropical_cover_gauge.jl`): each agent's delay split state by state
#         among its edges, the split optimized from the matching by supergradient ascent;
#       * simplicial cover: the cover over edges and triangles of the interaction graph,
#         each triangle solved exactly (the 2-cells the edge bounds cannot see).
#   - Upper bound (lax side): a feasible plan by prioritized planning (each agent plans
#     in space-time around the reservations of the agents before it). Upper / lower
#     certifies the plan's suboptimality without the optimum.
#   - Search on the cover: independence detection (Standley 2010). Agents follow their
#     own optimal paths; when two groups' plans collide the groups merge and are solved
#     jointly by A* on the fly, with a bound of the group as heuristic: the decoupled
#     Σᵢ dᵢ, the matching, or the larger of the matching and the cover gauge. Merging on
#     conflict is optimal, so the optimum is known and the sandwich can be checked.
#
#   julia --project=docs docs/scripts/tropical_multiagent_cover.jl [n] [instances]

include(joinpath(@__DIR__, "tropical_pair_sheaf.jl"))
using Printf
using Random
using Statistics
import CellularSheaves.NetworkSheaves.TropicalAbstraction as TA
include(joinpath(@__DIR__, "tropical_cover_gauge.jl"))

const n = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 20
const NI = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 4
const BLOCK = 2
const CAP = 300_000                                   # expansions per joint search
const MAXGROUP = 4                                    # joint searches move all agents at once: 5^m successors
const ITERS = 60                                      # supergradient steps of the cover gauge
const OUT = joinpath(@__DIR__, "..", "figures")

const WALL = (n ÷ 2):(n ÷ 2 + 1)
wall(x, y) = x in WALL && !(y in (n ÷ 4, 3n ÷ 4))
const MAP = GridMap([(x, y) for y in 1:n for x in 1:n if !wall(x, y)])
const NF = length(MAP.cells)
const NBRS = [[k; [MAP.index[(x + dx, y + dy)] for (dx, dy) in ((1, 0), (-1, 0), (0, 1), (0, -1))
                   if haskey(MAP.index, (x + dx, y + dy))]] for (k, (x, y)) in enumerate(MAP.cells)]
grid_agent(goal) = AgentStalk(collect(1:NF), ones(Int, NF), NBRS, goal)

collides(c, c2, k, l) = c2[k] == c2[l] || (c2[k] == c[l] && c2[l] == c[k])

# ---------------------------------------------------------------------------
# Plans: one cell sequence per agent, held at the goal after its end.
# ---------------------------------------------------------------------------

at(plan, t) = plan[min(t, length(plan))]
plan_cost(A, plan) = sum(step_cost(A.goal, at(plan, t), at(plan, t + 1)) for t in 1:(length(plan) - 1); init = 0.0)

function solo_plan(A, start)
    plan = [start]
    while plan[end] != A.goal
        c = plan[end]
        push!(plan, argmin(c2 -> step_cost(A.goal, c, c2) + A.d[c2], NBRS[c]))
    end
    return plan
end

# The first time two agents of different groups collide, or nothing.
function first_conflict(plans, group_of)
    T = maximum(length, plans)
    N = length(plans)
    for t in 1:T, k in 1:N, l in (k + 1):N
        group_of[k] == group_of[l] && continue
        c, c2 = (at(plans[k], t), at(plans[l], t)), (at(plans[k], t + 1), at(plans[l], t + 1))
        (c2[1] == c2[2] || (c2[1] == c[2] && c2[2] == c[1])) && return (k, l)
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Joint A* on the fly for a group of agents.
# ---------------------------------------------------------------------------

key(c) = 1 + sum((c[k] - 1) * NF^(k - 1) for k in eachindex(c))
unkey(x, m) = ntuple(k -> ((x - 1) ÷ NF^(k - 1)) % NF + 1, m)

function joint_astar(agents, starts, h)
    m = length(agents)
    goal = key(Tuple(A.goal for A in agents))
    s0 = key(Tuple(starts))
    g = Dict(s0 => 0.0)
    parent = Dict(s0 => 0)
    heap = Tuple{Float64,Int}[]
    TA._push!(heap, (h(Tuple(starts)), s0))
    expanded = 0
    while !isempty(heap)
        f, x = TA._pop!(heap)
        c = unkey(x, m)
        f > g[x] + h(c) + 1e-9 && continue                         # stale entry
        x == goal && break
        (expanded += 1) > CAP && return nothing, expanded
        for c2 in Iterators.product((NBRS[ci] for ci in c)...)
            any(collides(c, c2, k, l) for k in 1:m for l in (k + 1):m) && continue
            y = key(c2)
            gy = g[x] + sum(step_cost(agents[k].goal, c[k], c2[k]) for k in 1:m)
            if gy < get(g, y, Inf)
                g[y] = gy
                parent[y] = x
                TA._push!(heap, (gy + h(c2), y))
            end
        end
    end
    haskey(g, goal) || return nothing, expanded
    path = [goal]
    while parent[path[end]] != 0
        push!(path, parent[path[end]])
    end
    return [unkey(x, m) for x in reverse(path)], expanded
end

# Independence detection with heuristic `heuristic(group) -> (states -> bound)`.
function independence_detection(agents, starts, heuristic)
    N = length(agents)
    plans = [solo_plan(agents[i], starts[i]) for i in 1:N]
    group_of = collect(1:N)
    expanded, largest = 0, 1
    while (conflict = first_conflict(plans, group_of)) !== nothing
        a, b = group_of[conflict[1]], group_of[conflict[2]]
        members = [i for i in 1:N if group_of[i] in (a, b)]
        group_of[members] .= a
        largest = max(largest, length(members))
        length(members) > MAXGROUP && return nothing, expanded, largest
        path, e = joint_astar(agents[members], starts[members], heuristic(members))
        expanded += e
        path === nothing && return nothing, expanded, largest
        for (k, i) in enumerate(members)
            plans[i] = [c[k] for c in path]
        end
    end
    return plans, expanded, largest
end

# ---------------------------------------------------------------------------
# Prioritized planning: space-time A* around the reservations of earlier agents.
# ---------------------------------------------------------------------------

function prioritized(agents, starts, order)
    T = 6n
    vertex = Set{Tuple{Int,Int}}()                                  # (cell, time)
    edge = Set{Tuple{Int,Int,Int}}()                                # (from, to, time) of the move t → t + 1
    held = Dict{Int,Int}()                                          # cell -> time from which it is held forever
    last_use = Dict{Int,Int}()                                      # cell -> last time an earlier agent uses it
    plans = Vector{Vector{Int}}(undef, length(agents))
    for i in order
        A = agents[i]
        free(c, t) = !((c, t) in vertex) && t < get(held, c, typemax(Int))
        g = Dict((starts[i], 1) => 0.0)
        parent = Dict{Tuple{Int,Int},Tuple{Int,Int}}()
        heap = Tuple{Float64,Int,Int}[]
        TA._push!(heap, (A.d[starts[i]], starts[i], 1))
        found = nothing
        while !isempty(heap)
            f, c, t = TA._pop!(heap)
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
                    TA._push!(heap, (gy + A.d[c2], c2, t + 1))
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

# ---------------------------------------------------------------------------

# The best total cost of prioritized planning over a few priority orders.
function prioritized_cost(agents, starts, rng; orders = 5)
    N = length(agents)
    return minimum(1:orders) do r
        p = prioritized(agents, starts, r == 1 ? collect(1:N) : randperm(rng, N))
        p === nothing ? Inf : sum(plan_cost(agents[i], p[i]) for i in 1:N)
    end
end

# The weighted edge stalks with all weights 1 are the edge stalks: a check of GaugeFactor.
function check_factor(IS)
    (i, j) = first(IS.pairs)
    Ai, Aj = IS.agents[i], IS.agents[j]
    F = GaugeFactor(IS.map, Ai, Aj, i, j, 1, IS.stalks[(i, j)].b)
    Δ = weighted_stalk(F, ones(length(Ai.pos)), ones(length(Aj.pos))).Δ
    Δ0 = IS.stalks[(i, j)].Δ
    @assert length(Δ) == length(Δ0) && all(isapprox(a, b; atol = 1e-9) || (isinf(a) && isinf(b)) for (a, b) in zip(Δ, Δ0))
    println("check: the weighted edge stalk with unit weights is the edge stalk ($(count(isfinite, Δ)) coarse states)")
end

# The least total delay of the agents `m` planning together, ignoring the others.
function subset_delay(agents, starts, m)
    S = InteractionSheaf(MAP, agents[m], starts[m]; block = BLOCK)
    path, _ = joint_astar(agents[m], starts[m], c -> joint_bound(S, c))
    path === nothing && return NaN
    return sum(plan_cost(agents[m[k]], [c[k] for c in path]) - agents[m[k]].d[starts[m[k]]] for k in eachindex(m))
end

function instance(rng, N)
    left = [k for (k, (x, y)) in enumerate(MAP.cells) if x < first(WALL) - 2]
    right = [k for (k, (x, y)) in enumerate(MAP.cells) if x > last(WALL) + 2]
    sl, gl = left[randperm(rng, length(left))[1:N]], left[randperm(rng, length(left))[1:N]]
    sr, gr = right[randperm(rng, length(right))[1:N]], right[randperm(rng, length(right))[1:N]]
    half = N ÷ 2
    starts = [sl[1:half]; sr[1:(N - half)]]                         # half cross each way
    goals = [gr[1:half]; gl[1:(N - half)]]
    return starts, goals
end

# The study itself, when run as a script (other studies include the definitions above).
if abspath(PROGRAM_FILE) == @__FILE__
    println("grid $(n)², wall with two 1-cell corridors: $(NF) free cells")
    rows = NamedTuple[]
    rng = MersenneTwister(11)
    for N in (4, 6, 8), inst in 1:NI
        starts, goals = instance(rng, N)
        agents = [grid_agent(g) for g in goals]
        t_sheaf = @elapsed IS = InteractionSheaf(MAP, agents, starts; block = BLOCK)
        decoupled = sum(agents[i].d[starts[i]] for i in 1:N)
        lower = joint_bound(IS, starts)
        N == 4 && inst == 1 && check_factor(IS)
        # Upper bound: the best of a few priority orders.
        t_upper = @elapsed upper = prioritized_cost(agents, starts, rng)
        # Vertex cover and cover gauge.
        pw = zeros(N, N)
        for (i, j) in IS.pairs
            pw[i, j] = pw[j, i] = pair_bound(MAP, IS.stalks[(i, j)], IS.stalks[(j, i)], agents[i], agents[j], starts[i], starts[j])
        end
        cover = decoupled + min_vertex_cover(pw)
        # Triangles of the interaction graph, each solved exactly by joint A* (the identity
        # pushforward of its triple graph), and the cover over edges and triangles.
        E = Set(IS.pairs)
        tris = [[i, j, k] for i in 1:N for j in (i + 1):N for k in (j + 1):N if (i, j) in E && (i, k) in E && (j, k) in E]
        t_tri = @elapsed tri_delays = [subset_delay(agents, starts, m) for m in tris]
        simplices = [[([i, j], pw[i, j]) for (i, j) in IS.pairs]; [(m, Δ) for (m, Δ) in zip(tris, tri_delays) if !isnan(Δ)]]
        simplicial = decoupled + simplicial_cover(N, simplices)
        higher = count(((m, Δ),) -> Δ > maximum(pw[a, b] for a in m for b in m) + 1e-9, zip(tris, tri_delays))
        t_gauge = @elapsed CG = CoverGauge(IS, starts; target = upper - decoupled, iters = ITERS)
        gauge = decoupled + CG.bound
        @assert isapprox(decoupled + CG.matching, lower) "the ascent starts at the matching"
        @assert abs(gauge - (decoupled + gauge_bound(CG, starts))) < 1e-6 "the weighted stalks reproduce the bound"
        # The optimum by independence detection, with the decoupled and the sheaf heuristic.
        sub(members) = InteractionSheaf(MAP, agents[members], starts[members]; block = BLOCK)
        t_dec = @elapsed (pd, ed, gd) = independence_detection(agents, starts,
            members -> (c -> sum(agents[members[k]].d[c[k]] for k in eachindex(members))))
        t_she = @elapsed (ps, es, gs) = independence_detection(agents, starts, members -> (S = sub(members); c -> joint_bound(S, c)))
        gauged(members) = begin
            S = sub(members)
            up = prioritized_cost(agents[members], starts[members], MersenneTwister(1))   # leaves the instance stream alone
            dec = sum(agents[members[k]].d[starts[members[k]]] for k in eachindex(members))
            C = CoverGauge(S, starts[members]; target = up - dec, iters = ITERS)
            c -> sum(S.agents[k].d[c[k]] for k in eachindex(c)) + max(delay_bound(S, c), gauge_bound(C, c))
        end
        t_gau = @elapsed (pg, eg, gg) = independence_detection(agents, starts, gauged)
        cost(p) = p === nothing ? NaN : sum(plan_cost(agents[i], p[i]) for i in 1:N)
        opts = filter(!isnan, [cost(pd), cost(ps), cost(pg)])
        @assert isempty(opts) || all(==(first(opts)), opts) "every search is optimal"
        opt = isempty(opts) ? NaN : first(opts)
        @assert decoupled <= lower + 1e-9 && lower <= gauge + 1e-9
        @assert isnan(opt) || (max(gauge, cover, simplicial) <= opt + 1e-9 && opt <= upper + 1e-9)
        push!(rows, (N = N, instance = inst, pairs = length(IS.pairs), all_pairs = N * (N - 1) ÷ 2,
            decoupled = decoupled, lower = lower, cover = cover, gauge = gauge, simplicial = simplicial,
            triangles = length(tris), higher_triangles = higher, optimal = opt, upper = upper,
            expanded_decoupled = isnothing(pd) ? -ed : ed, expanded_sheaf = isnothing(ps) ? -es : es,
            expanded_gauge = isnothing(pg) ? -eg : eg,
            largest_group = gs, seconds_decoupled = t_dec, seconds_sheaf = t_she, seconds_gauge_search = t_gau,
            seconds_sheaf_build = t_sheaf, seconds_gauge = t_gauge, seconds_triangles = t_tri, seconds_upper = t_upper))
        @printf("N=%d #%d  pairs %2d/%2d  decoupled %5.1f  matching %5.1f  vertex cover %5.1f  cover gauge %6.2f (%.1f s)  simplicial %5.1f (%d triangles, %d beyond their edges; %.0f s)  optimum %5.1f  prioritized %5.1f\n",
            N, inst, length(IS.pairs), N * (N - 1) ÷ 2, decoupled, lower, cover, gauge, t_gauge, simplicial, length(tris), higher, t_tri, opt, upper)
        @printf("        independence detection: decoupled heuristic %9d expanded (%.1f s), matching %9d (%.1f s), gauge %9d (%.1f s), largest group %d%s\n",
            ed, t_dec, es, t_she, eg, t_gau, gs, isnan(opt) ? " (joint search out of reach)" : "")
    end

    avg(xs) = isempty(xs) ? NaN : mean(xs)
    println("\nsummary (means over instances)")
    println("gap closed: (bound - decoupled) / (optimum - decoupled) over solved instances with a gap; ratio: prioritized / bound")
    @printf("%3s %12s %30s %34s %12s %24s %14s\n", "N", "pairs", "gap: match/cover/gauge/simpl", "ratio: match/cover/gauge/simpl",
        "true ratio", "expanded: dec/match/gauge", "largest group")
    for N in (4, 6, 8)
        rs = filter(r -> r.N == N, rows)
        solved = filter(r -> !isnan(r.optimal), rs)
        gap(f) = (g = [(f(r) - r.decoupled) / (r.optimal - r.decoupled) for r in solved if r.optimal > r.decoupled];
                  isempty(g) ? "–" : @sprintf("%.2f", mean(g)))
        @printf("%3d %6.1f/%-5d %30s %34s %12.3f %24s %14.1f\n", N, mean(r.pairs for r in rs), rs[1].all_pairs,
            join((gap(r -> r.lower), gap(r -> r.cover), gap(r -> r.gauge), gap(r -> r.simplicial)), " / "),
            join((@sprintf("%.3f", mean(r.upper / f(r) for r in rs)) for f in (r -> r.lower, r -> r.cover, r -> r.gauge, r -> r.simplicial)), " / "),
            avg([r.upper / r.optimal for r in solved]),
            join((@sprintf("%.0f", avg([f(r) for r in solved])) for f in (r -> r.expanded_decoupled, r -> r.expanded_sheaf, r -> r.expanded_gauge)), " / "),
            mean(r.largest_group for r in rs))
        @printf("    solved by joint search: %d of %d\n", length(solved), length(rs))
    end

    mkpath(OUT)
    open(joinpath(OUT, "tropical_multiagent_cover.csv"), "w") do io
        println(io, join(keys(first(rows)), ","))
        for r in rows
            println(io, join(values(r), ","))
        end
    end
end
