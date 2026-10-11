# Layered control on the minimum-time double integrator lattice: a planner on the 2-D
# position grid (the pushforward along (q, v) ↦ q) and a low-level exact search on the
# 4-D states above a corridor around the planned route. What does committing to the
# route cost, against the exact 4-D optimum, as the corridor narrows?
#
#   julia --project=docs docs/scripts/tropical_layered_study.jl [n] [vmax] [queries]

using CellularSheaves
using SparseArrays
using Printf
using Random
using Statistics

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

timed(f) = (f(); t = @elapsed r = f(); (r, t))

n = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 128
vmax = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 6
nq = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 30
m = n ÷ 3
walls(qx, qy) = (qx in m:(m + 1) && !(qy in 2:(n ÷ 6))) || (qx in (2m):(2m + 1) && !(qy in (n - n ÷ 6):(n - 1)))
W, points = double_integrator_lattice(n, vmax, walls)
L = LinearIndices(points)
stop = vmax + 1
goal = [L[qx, qy, stop, stop] for qx in (n - 3):(n - 1) for qy in 2:4]
V = cost_to_go(W, goal)
rng = Random.MersenneTwister(4)
# Sources anywhere in the free space, at rest or moving.
candidates = [x for x in 1:size(W, 1) if isfinite(V[x]) && V[x] > 0]
sources = candidates[randperm(rng, length(candidates))[1:nq]]
println("lattice $(n)² positions × $(2vmax + 1)² velocities = $(size(W, 1)) states; $(nq) queries from random states")
@printf("exact optimal costs: mean %.1f, range %d–%d\n", mean(V[sources]), minimum(V[sources]), maximum(V[sources]))

ψ = projection_homomorphism(points, (1, 2))
P, tbuild = timed(() -> LayeredPlanner(W, ψ, (n, n)))
h, tbound = timed(() -> lower_bound(P.abstraction, goal))
@printf("planner: %d cells, %d quotient moves, built in %.2f s; 2-D bound for the goal in %.3f s (bound/cost %.2f)\n",
    ψ.n_target, nnz(P.abstraction.quotient), tbuild, tbound, mean(h[sources] ./ V[sources]))
troute = mean(timed(() -> plan_route(P, s, goal))[2] for s in sources)
@printf("route planning: %.4f s per query (2-D)\n", troute)

Wout = sparse(transpose(W))
function summarize(name, rs, ts)
    costs = [r.cost for r in rs]
    ratio = costs ./ V[sources]
    first = [hasproperty(r, :first_try) ? r.first_try : true for r in rs]
    states = [hasproperty(r, :corridor_states) ? r.corridor_states : count(isfinite, V) for r in rs]
    @printf("%-26s %6.0f%% %10.3f %8.3f %8.0f%% %12.0f %12.0f %9.4f s\n", name, 100 * mean(first), mean(ratio), maximum(ratio),
        100 * mean(ratio .== 1), mean(states), mean(r.expanded for r in rs), mean(ts))
    return (name = name, first_try = mean(first), mean_ratio = mean(ratio), max_ratio = maximum(ratio),
        optimal = mean(ratio .== 1), states = mean(states), expanded = mean(r.expanded for r in rs), seconds = mean(ts))
end

println("\nmethod                     first try  cost/opt  worst  optimal  search states  expansions   time/query")
rows = []
q(f) = (rs = [f(s) for s in sources]; ts = [timed(() -> f(s))[2] for s in sources]; (rs, ts))
push!(rows, summarize("4-D Dijkstra (exact)", q(s -> shortest_path(W, s, goal; out_edges = Wout))...))
push!(rows, summarize("4-D A*, 2-D bound (exact)", q(s -> shortest_path(W, s, goal; heuristic = h, out_edges = Wout))...))
for width in (0, 1, 2, 4, 8, 16)
    push!(rows, summarize("layered, width $(width)", q(s -> layered_path(P, s, goal; width, heuristic = h))...))
end

open(joinpath(@__DIR__, "..", "figures", "tropical_layered.csv"), "w") do io
    println(io, "method,first_try,mean_cost_over_optimal,max_cost_over_optimal,fraction_optimal,search_states,expansions,seconds_per_query")
    for r in rows
        println(io, join((r.name, r.first_try, r.mean_ratio, r.max_ratio, r.optimal, r.states, r.expanded, r.seconds), ","))
    end
end
