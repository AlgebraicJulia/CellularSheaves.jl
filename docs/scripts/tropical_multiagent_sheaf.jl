# N-agent cooperative planning from the interaction sheaf (`tropical_pair_sheaf.jl`),
# checked where the joint problem is still solvable: three agents on a small grid with a
# wall and two corridors. Compared heuristics on the joint graph:
#
#   - decoupled: Σᵢ dᵢ, the sum of single-agent distances;
#   - the interaction sheaf: Σᵢ dᵢ + the best matching of pair bounds Δᵢⱼ, each the edge stalk
#     built from the two agents' own graphs (no pair graph), checked equal to the pushforward
#     of the pair graph along the same map;
#   - the same with exact pairwise delays (the best any pairwise gluing can do);
#   - the exact joint cost-to-go.
#
#   julia --project=docs docs/scripts/tropical_multiagent_sheaf.jl [n] [queries]

include(joinpath(@__DIR__, "tropical_pair_sheaf.jl"))
using Printf
using Random
using Statistics

const n = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 10
const NQ = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 8
const NA = 3
const BLOCK = 2                                  # block size of the far summary
const OUT = joinpath(@__DIR__, "..", "figures")

const WALL = (n ÷ 2):(n ÷ 2 + 1)
wall(x, y) = x in WALL && !(y in (3, n - 2))
const MAP = GridMap([(x, y) for y in 1:n for x in 1:n if !wall(x, y)])
const NF = length(MAP.cells)
const NBRS = [[k; [MAP.index[(x + dx, y + dy)] for (dx, dy) in ((1, 0), (-1, 0), (0, 1), (0, -1))
                   if haskey(MAP.index, (x + dx, y + dy))]] for (k, (x, y)) in enumerate(MAP.cells)]
grid_agent(goal) = AgentStalk(collect(1:NF), ones(Int, NF), NBRS, goal)

# The joint problem, for checking only.
joint(cells) = 1 + sum((cells[i] - 1) * NF^(i - 1) for i in 1:NA)
cells_of(x) = ntuple(i -> ((x - 1) ÷ NF^(i - 1)) % NF + 1, NA)
conflict(c, c2) = any(c2[i] == c2[j] || (c2[i] == c[j] && c2[j] == c[i]) for i in 1:NA for j in (i + 1):NA)

function joint_graph(agents)
    rows, cols, vals = Int[], Int[], Float64[]
    for x in 1:NF^NA
        c = cells_of(x)
        allunique(c) || continue
        for c2 in Iterators.product((NBRS[ci] for ci in c)...)
            conflict(c, c2) && continue
            push!(rows, x); push!(cols, joint(c2)); push!(vals, sum(step_cost(agents[i].goal, c[i], c2[i]) for i in 1:NA))
        end
    end
    return sparse(rows, cols, vals, NF^NA, NF^NA)
end

# With three agents a matching is a single pair.
const PAIRS = [(i, j) for i in 1:NA for j in (i + 1):NA]

function glued(agents, pair_value)
    h = fill(Inf, NF^NA)
    for x in 1:NF^NA
        c = cells_of(x)
        allunique(c) || continue
        h[x] = sum(agents[i].d[c[i]] for i in 1:NA) + maximum(pair_value(i, j, c[i], c[j]) for (i, j) in PAIRS)
    end
    return h
end

function queries(rng)
    left = [k for (k, (x, y)) in enumerate(MAP.cells) if x < first(WALL) - 1]
    right = [k for (k, (x, y)) in enumerate(MAP.cells) if x > last(WALL) + 1]
    cross = map(1:NQ) do _
        s12, g3 = left[randperm(rng, length(left))[1:2]], left[rand(rng, 1:length(left))]
        g12, s3 = right[randperm(rng, length(right))[1:2]], right[rand(rng, 1:length(right))]
        starts, goals = (s12..., s3), (g12..., g3)
        (allunique(starts) && allunique(goals)) ? ("cross", starts, goals) : nothing
    end
    randoms = [("random", Tuple(p[1:NA]), Tuple(p[(NA + 1):(2NA)])) for p in (randperm(rng, NF) for _ in 1:NQ)]
    return filter(!isnothing, [cross; randoms])
end

println("grid $(n)², wall with two corridors: $(NF) free cells, $(NA) agents, $(NF^NA) joint states")
rows = NamedTuple[]
for (kind, starts, goals) in queries(MersenneTwister(3))
    agents = [grid_agent(g) for g in goals]
    t_sheaf = @elapsed P = Dict((i, j) => PairAbstraction(MAP, agents[i], agents[j], BLOCK) for i in 1:NA for j in 1:NA if i != j)
    t_pushed = @elapsed pushed = Dict((i, j) => pushed_pair(MAP, P[(i, j)], agents[i], agents[j]) for i in 1:NA for j in 1:NA if i != j)
    t_exact = @elapsed exact = Dict((i, j) => cost_to_go(pair_graph(agents[i], agents[j]), [pair_goal(agents[i], agents[j])])
                                    for (i, j) in PAIRS)
    sheaf_states = sum(coarse_states(p) for p in values(P))
    W = joint_graph(agents)
    goal, start = joint(goals), joint(starts)
    V = cost_to_go(W, [goal])
    isfinite(V[start]) || continue
    # The edge stalks are the exact pushforwards of the pair graphs, and bound the pair delays.
    for (i, j) in PAIRS, b in 1:NF, a in 1:NF
        a == b && continue
        e = exact[(i, j)][a + NF * (b - 1)]
        h = pair_bound(MAP, P[(i, j)], P[(j, i)], agents[i], agents[j], a, b)
        @assert h <= e + 1e-9
        @assert isinf(e) || h ≈ max(pushed[(i, j)][a + NF * (b - 1)], pushed[(j, i)][b + NF * (a - 1)])
    end
    dec = glued(agents, (i, j, a, b) -> 0.0)
    sheaf = glued(agents, (i, j, a, b) -> pair_bound(MAP, P[(i, j)], P[(j, i)], agents[i], agents[j], a, b))
    # The same through the interaction sheaf: every pair (its matching reduces to the best
    # pair for three agents), and the sparse interaction graph of meeting tubes.
    full = InteractionSheaf(MAP, agents, starts; block = BLOCK, slack = Inf)
    t_sparse = @elapsed sparse_sheaf = InteractionSheaf(MAP, agents, starts; block = BLOCK)
    sparse_h = fill(Inf, NF^NA)
    for x in 1:NF^NA
        c = cells_of(x)
        allunique(c) || continue
        @assert joint_bound(full, c) ≈ sheaf[x]
        sparse_h[x] = joint_bound(sparse_sheaf, c)
    end
    pairwise = glued(agents, (i, j, a, b) -> exact[(i, j)][a + NF * (b - 1)])
    reach = isfinite.(V)
    out = sparse(transpose(W))
    @printf("%-6s  optimal %3.0f  decoupled %3.0f  (coordination cost %.0f), interacting pairs %d of %d\n", kind,
        V[start], dec[start], V[start] - dec[start], length(sparse_sheaf.pairs), length(PAIRS))
    for (name, h, build, states) in (("decoupled (sum of distances)", dec, 0.0, 0),
                                     ("sheaf: edge stalks, block $BLOCK", sheaf, t_sheaf, sheaf_states),
                                     ("sheaf: interacting pairs only", sparse_h, t_sparse,
                                      sum(coarse_states(p) for p in values(sparse_sheaf.stalks); init = 0)),
                                     ("exact pairwise delays", pairwise, t_exact, length(PAIRS) * NF^2),
                                     ("exact cost-to-go", V, 0.0, NF^NA))
        @assert all(h[reach] .<= V[reach] .+ 1e-9) "$name is not admissible"
        r = shortest_path(W, start, [goal]; heuristic = h, out_edges = out)
        @assert r.cost == V[start]
        push!(rows, (kind = kind, heuristic = name, states = states, optimal = V[start], decoupled = dec[start],
            bound = h[start], expanded = r.expanded, build_seconds = build))
        @printf("    %-44s bound %5.1f  expanded %7d  build %6.3f s\n", name, h[start], r.expanded, build)
    end
end

println("\nsummary (means over queries)")
@printf("%-7s %-44s %9s %10s %11s %9s %9s\n", "queries", "heuristic", "states", "bound/opt", "gap closed", "expanded", "build s")
for kind in ("cross", "random"), name in unique(r.heuristic for r in rows)
    rs = filter(r -> r.kind == kind && r.heuristic == name, rows)
    isempty(rs) && continue
    gap = [(r.bound - r.decoupled) / (r.optimal - r.decoupled) for r in rs if r.optimal > r.decoupled]
    @printf("%-7s %-44s %9d %10.3f %11s %9.0f %9.3f\n", kind, name, round(Int, mean(r.states for r in rs)),
        mean(r.bound / r.optimal for r in rs), isempty(gap) ? "–" : @sprintf("%.2f", mean(gap)),
        mean(r.expanded for r in rs), mean(r.build_seconds for r in rs))
end

mkpath(OUT)
open(joinpath(OUT, "tropical_multiagent_sheaf.csv"), "w") do io
    println(io, join(keys(first(rows)), ","))
    for r in rows
        println(io, join((occursin(',', string(v)) ? "\"$v\"" : v for v in values(r)), ","))
    end
end
