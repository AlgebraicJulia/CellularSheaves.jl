# Pairs of double integrators: the edge stalk of the interaction sheaf
# (`tropical_pair_sheaf.jl`) where an agent's state includes its velocity, so a pair's
# state space is the product of two 4-D lattices. The map is small enough that the pair
# graph can still be solved exactly, to check that the edge stalk built from the two
# agents' own graphs equals the pushforward of the pair graph along its map, and to
# compare it with the exact pair delay, for velocity modes of the near field from exact
# down to rest/moving (the near states are still enumerated exactly, so each is the exact
# pushforward along its own map). Searches
# run on the pair graph in the delay gauge (Dijkstra there is A* with the decoupled
# heuristic Σᵢ dᵢ).
#
#   julia --project=docs docs/scripts/tropical_multiagent_double_integrator.jl [n] [queries]

include(joinpath(@__DIR__, "tropical_pair_sheaf.jl"))
using Printf
using Random
using Statistics

const n = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 16
const NQ = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 8
const BLOCK = 2                                  # block size of the far summary
const OUT = joinpath(@__DIR__, "..", "figures")

const WALL = (n ÷ 2):(n ÷ 2 + 1)
wall(x, y) = x in WALL && !(y in (4, n - 3))
const MAP = GridMap([(x, y) for y in 1:n for x in 1:n if !wall(x, y)])
const NF = length(MAP.cells)

# The double-integrator lattice: state (cell, velocity), velocity in {-1, 0, 1}²; a move
# changes each velocity component by at most one and then moves by the new velocity.
const VELS = [(vx, vy) for vy in -1:1 for vx in -1:1]
const STATES = [(c, w) for w in 1:9 for c in 1:NF]
const STATE = Dict(s => k for (k, s) in enumerate(STATES))
const SUCC = map(STATES) do (c, w)
    out = Int[]
    for (w2, v2) in enumerate(VELS)
        maximum(abs.(v2 .- VELS[w])) <= 1 || continue
        c2 = cell_at(MAP, c, v2)
        c2 == 0 || push!(out, STATE[(c2, w2)])
    end
    out
end
const REST = findfirst(==((0, 0)), VELS)
# Velocity modes of the near field: exact, rest + 4 headings, rest/moving.
heading(v) = v == (0, 0) ? 1 : 2 + mod(round(Int, atan(v[2], v[1]) / (π / 2)), 4)
const VELMODES = [("exact velocity", collect(1:9)), ("rest + 4 headings", heading.(VELS)),
                  ("rest / moving", [v == (0, 0) ? 1 : 2 for v in VELS])]
di_agent(goal_cell) = AgentStalk(first.(STATES), last.(STATES), SUCC, STATE[(goal_cell, REST)])

function queries(rng)
    left = [k for (k, (x, y)) in enumerate(MAP.cells) if x < first(WALL) - 1 && abs(y - 4) <= 3]
    right = [k for (k, (x, y)) in enumerate(MAP.cells) if x > last(WALL) + 1 && abs(y - 4) <= 3]
    swaps = map(1:NQ) do _
        a1, g2 = left[randperm(rng, length(left))[1:2]]
        a2, g1 = right[randperm(rng, length(right))[1:2]]
        ("swap", (a1, a2), (g1, g2))
    end
    randoms = [("random", (p[1], p[2]), (p[3], p[4])) for p in (randperm(rng, NF) for _ in 1:NQ)]
    return [swaps; randoms]
end

S = length(STATES)
println("grid $(n)², wall with two corridors, velocities in {-1, 0, 1}²: $(S) states per agent, $(S^2) pair states")
rows = NamedTuple[]
for (kind, starts, goals) in queries(MersenneTwister(5))
    A1, A2 = di_agent(goals[1]), di_agent(goals[2])
    s0 = STATE[(starts[1], REST)] + S * (STATE[(starts[2], REST)] - 1)
    t_pair = @elapsed W = pair_graph(A1, A2)
    goal = pair_goal(A1, A2)
    t_exact = @elapsed D = cost_to_go(W, [goal])                # exact pair delay
    isfinite(D[s0]) || continue
    base = A1.d[STATE[(starts[1], REST)]] + A2.d[STATE[(starts[2], REST)]]   # Σᵢ dᵢ at the start
    W21 = pair_graph(A2, A1)
    reach = isfinite.(D)
    stalks = map(VELMODES) do (label, vm)
        tb = @elapsed P12, P21 = PairAbstraction(MAP, A1, A2, BLOCK; velmode = vm), PairAbstraction(MAP, A2, A1, BLOCK; velmode = vm)
        h = [pair_bound(MAP, P12, P21, A1, A2, s, t) for t in 1:S for s in 1:S]
        Q12, Q21 = pushed_pair(MAP, P12, A1, A2, W), pushed_pair(MAP, P21, A2, A1, W21)
        pushforward = [max(Q12[s + S * (t - 1)], Q21[t + S * (s - 1)]) for t in 1:S for s in 1:S]
        @assert all(h[reach] .<= D[reach] .+ 1e-9)
        @assert all(h[reach] .≈ pushforward[reach])              # the edge stalk is the exact pushforward
        ("sheaf, near field: $label", h, tb, coarse_states(P12) + coarse_states(P21))
    end
    out = sparse(transpose(W))
    @printf("%-6s  optimal %3.0f  decoupled %3.0f  (coordination cost %.0f)\n", kind, base + D[s0], base, D[s0])
    for (name, h, build, states) in (("decoupled (Dijkstra on delays)", nothing, 0.0, 0), stalks...,
                                     ("exact pair delay", D, t_pair + t_exact, count(reach)))
        r = shortest_path(W, s0, [goal]; heuristic = h, out_edges = out)
        @assert r.cost == D[s0]
        bound = isnothing(h) ? 0.0 : h[s0]
        push!(rows, (kind = kind, heuristic = name, states = states, optimal = base + D[s0], decoupled = base,
            bound = base + bound, expanded = r.expanded, build_seconds = build))
        @printf("    %-40s bound %5.1f  expanded %8d  build %6.2f s\n", name, base + bound, r.expanded, build)
    end
end

println("\nsummary (means over queries)")
@printf("%-7s %-40s %10s %10s %11s %9s %9s\n", "queries", "heuristic", "states", "bound/opt", "gap closed", "expanded", "build s")
for kind in ("swap", "random"), name in unique(r.heuristic for r in rows)
    rs = filter(r -> r.kind == kind && r.heuristic == name, rows)
    isempty(rs) && continue
    gap = [(r.bound - r.decoupled) / (r.optimal - r.decoupled) for r in rs if r.optimal > r.decoupled]
    @printf("%-7s %-40s %10d %10.3f %11s %9.0f %9.2f\n", kind, name, round(Int, mean(r.states for r in rs)),
        mean(r.bound / r.optimal for r in rs), isempty(gap) ? "–" : @sprintf("%.2f", mean(gap)),
        mean(r.expanded for r in rs), mean(r.build_seconds for r in rs))
end

mkpath(OUT)
open(joinpath(OUT, "tropical_multiagent_double_integrator.csv"), "w") do io
    println(io, join(keys(first(rows)), ","))
    for r in rows
        println(io, join((occursin(',', string(v)) ? "\"$v\"" : v for v in values(r)), ","))
    end
end
