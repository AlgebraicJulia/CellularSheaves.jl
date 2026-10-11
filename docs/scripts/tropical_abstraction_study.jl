# Coarse grids in min-plus: the minimum-time double integrator on a lattice,
# solved exactly by A* guided by the quotient graph of a position aggregation, against
# Dijkstra's algorithm. For each depth of the hierarchy (blocks of 2^ℓ × 2^ℓ
# positions) it reports the size of the coarse problem, the time to solve it, and
# the states the exact fine search expands.
#
#   julia --project=docs docs/scripts/tropical_abstraction_study.jl [n] [vmax]

using CellularSheaves
using SparseArrays
using Printf
using Random

function double_integrator_lattice(n, vmax, obstacle)
    nv = 2vmax + 1
    points = (n, n, nv, nv)
    L = LinearIndices(points)
    rows, cols = Int[], Int[]
    for I in CartesianIndices(points)
        qx, qy, ix, iy = Tuple(I)
        vx, vy = ix - vmax - 1, iy - vmax - 1
        obstacle(qx, qy) && continue
        for ax in -1:1, ay in -1:1
            wx, wy = vx + ax, vy + ay
            (abs(wx) <= vmax && abs(wy) <= vmax) || continue
            px, py = qx + wx, qy + wy
            (1 <= px <= n && 1 <= py <= n) || continue
            (obstacle(px, py) || obstacle(qx + wx ÷ 2, qy + wy ÷ 2)) && continue
            push!(rows, L[I]); push!(cols, L[px, py, wx + vmax + 1, wy + vmax + 1])
        end
    end
    N = prod(points)
    return sparse(rows, cols, ones(length(rows)), N, N), points
end

n = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 64
vmax = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 4
# Two walls with offset gaps, so that the time-optimal path winds.
m = n ÷ 3
walls(qx, qy) = (qx in m:(m + 1) && !(qy in 2:(n ÷ 6))) || (qx in (2m):(2m + 1) && !(qy in (n - n ÷ 6):(n - 1)))
W, points = double_integrator_lattice(n, vmax, walls)
L = LinearIndices(points)
stop = vmax + 1
goal = [L[qx, qy, stop, stop] for qx in (n - 3):(n - 1) for qy in 2:4]
println("double integrator lattice: $(n)² positions × $(2vmax + 1)² velocities = $(size(W, 1)) states, $(nnz(W)) edges")

cost_to_go(W, goal; source = goal[1])                                   # compile
t = @elapsed full = cost_to_go(W, goal)
@printf("full value function (Dijkstra): %.3f s, %d states reachable\n", t, full.expanded)

rng = Random.MersenneTwister(1)
candidates = [L[qx, qy, stop, stop] for qx in 1:m-1 for qy in 1:n if isfinite(full.values[L[qx, qy, stop, stop]])]
sources = candidates[randperm(rng, length(candidates))[1:10]]

function measure(f)
    f()
    return @elapsed f()
end

rows = []
dijkstra_expanded = Int[]
dijkstra_time = 0.0
for s in sources
    r = cost_to_go(W, goal; source = s)
    push!(dijkstra_expanded, r.expanded)
    global dijkstra_time += measure(() -> cost_to_go(W, goal; source = s))
end
push!(rows, (level = 0, blocks = size(W, 1), setup_seconds = 0.0, coarse_seconds = 0.0, expanded = sum(dijkstra_expanded) / length(sources),
    fine_seconds = dijkstra_time / length(sources), bound = 1.0))

ψ = GraphHomomorphism(collect(1:size(W, 1)), size(W, 1))
sizes = points
level = 0
while sizes[1] > 2
    global level += 1
    step = aggregation_homomorphism(sizes, (2, 2, 1, 1))
    global ψ = compose(ψ, step)
    global sizes = coarse_points(sizes, (2, 2, 1, 1))
    setup = measure(() -> sparse(transpose(quotient_graph(W, ψ))))   # once per abstraction, any number of queries
    Wr = sparse(transpose(quotient_graph(W, ψ)))
    # How much of the cost the coarse bound sees: the coarse cost-to-go at the source.
    Vc = cost_to_go(sparse(transpose(Wr)), unique(ψ.vertex_map[goal])).values
    expanded, coarse_t, fine_t, ratio = 0, 0.0, 0.0, 0.0
    for s in sources
        h = abstraction_heuristic(W, ψ, s; reversed_quotient = Wr)
        coarse_t += measure(() -> abstraction_heuristic(W, ψ, s; reversed_quotient = Wr))
        r = cost_to_go(W, goal; source = s, heuristic = h)
        r.cost == full.values[s] || error("A* is not exact")
        expanded += r.expanded
        fine_t += measure(() -> cost_to_go(W, goal; source = s, heuristic = h))
        ratio += Vc[ψ.vertex_map[s]] / full.values[s]
    end
    k = length(sources)
    push!(rows, (level = level, blocks = ψ.n_target, setup_seconds = setup, coarse_seconds = coarse_t / k, expanded = expanded / k,
        fine_seconds = fine_t / k, bound = ratio / k))
end

println("\nper query, averaged over $(length(sources)) sources (exact costs $(Int.(full.values[sources])))")
println("level  blocks of positions   coarse states  quotient (once)  coarse solve  bound/cost  fine expansions  fine search  per query")
for r in rows
    b = 2^r.level
    @printf("%5d  %-20s %13d %14.4f s %11.4f s %11.2f %16.0f %10.4f s %8.4f s\n", r.level,
        r.level == 0 ? "none (Dijkstra)" : "$(b) × $(b)", r.blocks, r.setup_seconds, r.coarse_seconds, r.bound, r.expanded,
        r.fine_seconds, r.coarse_seconds + r.fine_seconds)
end

out = joinpath(@__DIR__, "..", "figures", "tropical_abstraction.csv")
open(out, "w") do io
    println(io, "level,block_side,coarse_states,quotient_seconds,coarse_seconds,bound_over_cost,fine_expansions,fine_seconds")
    for r in rows
        println(io, join((r.level, 2^r.level, r.blocks, r.setup_seconds, r.coarse_seconds, r.bound, r.expanded, r.fine_seconds), ","))
    end
end
println("wrote $out")
