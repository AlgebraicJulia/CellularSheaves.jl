# Cooperative two-agent planning with interaction-aware abstractions. Each agent's
# coarse problem keeps its own cell exactly and a low-dimensional summary of the other
# agent: its relative offset when near (the encounter, kept exactly) and its block on a
# coarse grid when far. Pushed forward along that summary, the joint graph gives an
# admissible lower bound on the joint cost; how much of the coordination cost (waiting,
# detours) does it see, against the decoupled sum of single-agent distances and the
# exact joint cost-to-go?
#
#   julia --project=docs docs/scripts/tropical_multiagent_study.jl [n] [queries]

using CellularSheaves
using SparseArrays
using Printf
using Random
using Statistics

const n = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 32
const NQ = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 8
const OUT = joinpath(@__DIR__, "..", "figures")

# A thick wall (columns 14–19) with two one-cell corridors (rows 16 and 28): agents
# that swap sides must queue through a corridor or detour to the other one.
wall(x, y) = (n ÷ 2 - 2) <= x <= (n ÷ 2 + 3) && !(y in (n ÷ 2, n - 4))
const CELLS = [(x, y) for y in 1:n for x in 1:n if !wall(x, y)]
const INDEX = Dict(c => k for (k, c) in enumerate(CELLS))
const NF = length(CELLS)
const NBRS = [[k; [INDEX[(x + dx, y + dy)] for (dx, dy) in ((1, 0), (-1, 0), (0, 1), (0, -1))
                   if haskey(INDEX, (x + dx, y + dy))]] for (k, (x, y)) in enumerate(CELLS)]

# Joint states (a, b) with a ≠ b, indexed a + NF (b - 1) (states with a = b are unused and
# isolated). A joint move moves both agents one step (or keeps them), with no vertex
# conflict (same cell) and no swap conflict (exchanging cells).
joint(a, b) = a + NF * (b - 1)
agents(x) = ((x - 1) % NF + 1, (x - 1) ÷ NF + 1)

function joint_structure()
    rows, cols = Int[], Int[]
    for b in 1:NF, a in 1:NF
        a == b && continue
        for a2 in NBRS[a], b2 in NBRS[b]
            (a2 == b2 || (a2 == b && b2 == a)) && continue
            push!(rows, joint(a, b)); push!(cols, joint(a2, b2))
        end
    end
    return sparse(rows, cols, ones(length(rows)), NF^2, NF^2)
end

# Each agent pays 1 per step, waiting included, except while it stays at its own goal.
# Returns the joint costs and each agent's own share, on the same stored pattern.
function costs(S, g1, g2)
    v1, v2 = similar(nonzeros(S)), similar(nonzeros(S))
    rv = rowvals(S)
    for y in 1:size(S, 2), p in nzrange(S, y)
        (a, b), (a2, b2) = agents(rv[p]), agents(y)
        v1[p] = a == a2 == g1 ? 0.0 : 1.0
        v2[p] = b == b2 == g2 ? 0.0 : 1.0
    end
    mk(v) = SparseMatrixCSC(size(S)..., S.colptr, rv, v)
    return mk(v1 .+ v2), mk(v1), mk(v2)
end

function bfs(g)
    d = fill(Inf, NF)
    d[g] = 0
    queue = [g]
    while !isempty(queue)
        u = popfirst!(queue)
        for v in NBRS[u]
            isinf(d[v]) && (d[v] = d[u] + 1; push!(queue, v))
        end
    end
    return d
end

# Summaries of the other agent, as seen by the agent whose cell is `own`: a mode in 1:K.
cheb(o) = max(abs(o[1]), abs(o[2]))
offset(own, other) = CELLS[other] .- CELLS[own]
near_index(o, r) = (o[1] + r) + (2r + 1) * (o[2] + r) + 1
block_index(c, B) = (CELLS[c][1] - 1) ÷ cld(n, B) + B * ((CELLS[c][2] - 1) ÷ cld(n, B)) + 1

struct Summary{F}
    name::String
    K::Int
    mode::F                 # (own, other) -> 1:K
end

solo() = Summary("solo (1 mode)", 1, (own, other) -> 1)
near(r) = Summary("near r=$r", (2r + 1)^2 + 1, function (own, other)
    o = offset(own, other)
    cheb(o) <= r ? near_index(o, r) : (2r + 1)^2 + 1
end)
blocks(B) = Summary("blocks $(B)×$(B)", B^2, (own, other) -> block_index(other, B))
blocks_near(B, r) = Summary("blocks $(B)×$(B) + near r=$r", (2r + 1)^2 + B^2, function (own, other)
    o = offset(own, other)
    cheb(o) <= r ? near_index(o, r) : (2r + 1)^2 + block_index(other, B)
end)

# The mode map of agent `i` (1 or 2) on the joint states: (own cell, summary of the other).
function mode_map(s::Summary, i)
    vmap = Vector{Int}(undef, NF^2)
    for x in 1:NF^2
        a, b = agents(x)
        own, other = i == 1 ? (a, b) : (b, a)
        vmap[x] = own + NF * (s.mode(own, other) - 1)
    end
    return GraphHomomorphism(vmap, NF * s.K)
end

# The lower bound of agent i's abstraction along summary s, for costs Wc.
function bound(Wc, s, i, goal)
    A = Abstraction(Wc, mode_map(s, i))
    return lower_bound(A, [goal]), size(A.quotient, 1)
end

# Heuristics: name => (W, W1, W2, goal, decoupled) -> (heuristic over the joint states, coarse states)
const HEURISTICS = [
    "decoupled (sum of distances)" => (W, W1, W2, goal, dec) -> (dec, 0),
    "own costs, near r=2, summed" => function (W, W1, W2, goal, dec)
        (h1, c1), (h2, c2) = bound(W1, near(2), 1, goal), bound(W2, near(2), 2, goal)
        return h1 .+ h2, c1 + c2
    end,
    "joint costs, near r=2" => function (W, W1, W2, goal, dec)
        (h1, c1), (h2, c2) = bound(W, near(2), 1, goal), bound(W, near(2), 2, goal)
        return max.(dec, h1, h2), c1 + c2
    end,
    "joint costs, blocks 4×4" => function (W, W1, W2, goal, dec)
        (h1, c1), (h2, c2) = bound(W, blocks(4), 1, goal), bound(W, blocks(4), 2, goal)
        return max.(dec, h1, h2), c1 + c2
    end,
    "joint costs, blocks 4×4 + near r=2" => function (W, W1, W2, goal, dec)
        (h1, c1), (h2, c2) = bound(W, blocks_near(4, 2), 1, goal), bound(W, blocks_near(4, 2), 2, goal)
        return max.(dec, h1, h2), c1 + c2
    end,
    "joint costs, blocks 8×8 + near r=2" => function (W, W1, W2, goal, dec)
        (h1, c1), (h2, c2) = bound(W, blocks_near(8, 2), 1, goal), bound(W, blocks_near(8, 2), 2, goal)
        return max.(dec, h1, h2), c1 + c2
    end,
]

# Queries: agents swapping sides near the corridors (they must coordinate), and random
# start/goal pairs (they mostly need not).
function queries(rng)
    left = [k for (k, (x, y)) in enumerate(CELLS) if n ÷ 2 - 9 <= x <= n ÷ 2 - 4 && abs(y - n ÷ 2) <= 4]
    right = [k for (k, (x, y)) in enumerate(CELLS) if n ÷ 2 + 5 <= x <= n ÷ 2 + 10 && abs(y - n ÷ 2) <= 4]
    swaps = map(1:NQ) do _
        a1, g2 = randperm(rng, length(left))[1:2]
        a2, g1 = randperm(rng, length(right))[1:2]
        ("swap", left[a1], right[a2], right[g1], left[g2])
    end
    randoms = map(1:NQ) do _
        a1, a2, g1, g2 = randperm(rng, NF)[1:4]
        ("random", a1, a2, g1, g2)
    end
    return [swaps; randoms]
end

println("grid $(n)², wall with two 1-cell corridors: $(NF) free cells, $(NF^2) joint states")
t_structure = @elapsed S = joint_structure()
@printf("joint graph: %d moves, built in %.1f s\n\n", nnz(S), t_structure)

rows = NamedTuple[]
rng = MersenneTwister(7)
for (kind, a1, a2, g1, g2) in queries(rng)
    W, W1, W2 = costs(S, g1, g2)
    goal = joint(g1, g2)
    start = joint(a1, a2)
    V = cost_to_go(W, [goal])
    d1, d2 = bfs(g1), bfs(g2)
    dec = [d1[a] + d2[b] for b in 1:NF for a in 1:NF]
    out = sparse(transpose(W))
    exact = shortest_path(W, start, [goal]; heuristic = V, out_edges = out)
    @assert exact.cost == V[start]
    @printf("%-6s  optimal %3.0f  decoupled %3.0f  (coordination cost %.0f)\n", kind, V[start], dec[start], V[start] - dec[start])
    for (name, heuristic) in HEURISTICS
        build = @elapsed h, coarse = heuristic(W, W1, W2, goal, dec)
        reachable = isfinite.(V)
        @assert all(h[reachable] .<= V[reachable] .+ 1e-9) "$name is not admissible"
        search = @elapsed r = shortest_path(W, start, [goal]; heuristic = h, out_edges = out)
        @assert r.cost == V[start]
        push!(rows, (kind = kind, heuristic = name, coarse_states = coarse, optimal = V[start],
            decoupled = dec[start], bound = h[start], expanded = r.expanded, exact_expanded = exact.expanded,
            build_seconds = build, search_seconds = search))
        @printf("    %-36s bound %5.1f  expanded %8d  build %5.2f s  search %6.3f s\n",
            name, h[start], r.expanded, build, search)
    end
end

println("\nsummary (means over queries)")
@printf("%-8s %-36s %13s %11s %12s %11s %9s %9s\n", "queries", "heuristic", "coarse states", "bound/opt",
    "gap closed", "expanded", "build s", "search s")
for kind in ("swap", "random"), (name, _) in HEURISTICS
    rs = filter(r -> r.kind == kind && r.heuristic == name, rows)
    gap = [(r.bound - r.decoupled) / (r.optimal - r.decoupled) for r in rs if r.optimal > r.decoupled]
    @printf("%-8s %-36s %13d %11.3f %12s %11.0f %9.2f %9.3f\n", kind, name, round(Int, mean(r.coarse_states for r in rs)),
        mean(r.bound / r.optimal for r in rs), isempty(gap) ? "–" : @sprintf("%.2f", mean(gap)),
        mean(r.expanded for r in rs), mean(r.build_seconds for r in rs), mean(r.search_seconds for r in rs))
end
for kind in ("swap", "random")
    rs = filter(r -> r.kind == kind && r.heuristic == first(HEURISTICS)[1], rows)
    @printf("%-8s %-36s %13s %11.3f %12s %11.0f\n", kind, "exact cost-to-go (reference)", "–", 1.0, "1.00",
        mean(r.exact_expanded for r in rs))
end

mkpath(OUT)
open(joinpath(OUT, "tropical_multiagent.csv"), "w") do io
    println(io, join(keys(first(rows)), ","))
    for r in rows
        println(io, join((occursin(',', string(v)) ? "\"$v\"" : v for v in values(r)), ","))
    end
end
