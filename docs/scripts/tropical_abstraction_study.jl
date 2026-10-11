# Coarse grids in idempotent semirings, on the minimum-time double integrator lattice.
#
#   julia --project=docs docs/scripts/tropical_abstraction_study.jl [n] [vmax]
#
# 1. Exact A* queries guided by abstractions: positions only (blocks 2^ℓ, one time
#    step) against space–time (velocity blocks k, position blocks k², k steps), with
#    Dijkstra's algorithm as the baseline.
# 2. The coarse closure for many targets: Dijkstra on the coarse graph per target set
#    against PathPack's chordal min-plus LU factorized once.
# 3. Reach-avoid in the bottleneck semiring: the safest route to the goal, checked
#    against min-plus reachability.

using CellularSheaves
using PATHPACK
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

function timed(f)
    f()
    t = @elapsed r = f()
    return r, t
end

n = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 128
vmax = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 6
m = n ÷ 3
walls(qx, qy) = (qx in m:(m + 1) && !(qy in 2:(n ÷ 6))) || (qx in (2m):(2m + 1) && !(qy in (n - n ÷ 6):(n - 1)))
W, points = double_integrator_lattice(n, vmax, walls)
Wout = sparse(transpose(W))
L = LinearIndices(points)
stop = vmax + 1
goal = [L[qx, qy, stop, stop] for qx in (n - 3):(n - 1) for qy in 2:4]
println("double integrator lattice: $(n)² positions × $(2vmax + 1)² velocities = $(size(W, 1)) states, $(nnz(W)) edges")
V, t = timed(() -> cost_to_go(W, goal))
@printf("full value function (Dijkstra): %.3f s, %d states can reach the goal\n", t, count(isfinite, V))

rng = Random.MersenneTwister(1)
candidates = [L[qx, qy, stop, stop] for qx in 1:(m - 1) for qy in 1:n if isfinite(V[L[qx, qy, stop, stop]])]
sources = candidates[randperm(rng, length(candidates))[1:10]]
println("10 queries from the left third, stopped; exact costs $(Int.(V[sources]))")

# ---- 1. abstractions ------------------------------------------------------------
function coarsen(ψ, sizes, factors)
    return compose(ψ, aggregation_homomorphism(sizes, factors)), coarse_points(sizes, factors)
end
id = GraphHomomorphism(collect(1:prod(points)), prod(points))
abstractions = Any[]
ψ, sizes = id, points
for ℓ in 1:5
    global ψ, sizes = coarsen(ψ, sizes, (2, 2, 1, 1))
    push!(abstractions, ("space", "$(2^ℓ)×$(2^ℓ)", "1", ψ, 1))
end
ψ, sizes = id, points
for k in (2, 4, 8)
    global ψ, sizes = coarsen(ψ, sizes, (2, 2, 2, 2))
    global ψ, sizes = coarsen(ψ, sizes, (2, 2, 1, 1))
    push!(abstractions, ("space–time", "$(k^2)×$(k^2)", "$(k)", ψ, k))
end

results = []
base = [shortest_path(W, s, goal; out_edges = Wout) for s in sources]
tbase = sum(timed(() -> shortest_path(W, s, goal; out_edges = Wout))[2] for s in sources) / length(sources)
push!(results, (kind = "Dijkstra", pos = "–", vel = "–", steps = "–", blocks = size(W, 1), build = 0.0, goal_setup = 0.0,
    bound = 1.0, expanded = sum(r.expanded for r in base) / length(sources), query = tbase))
for (kind, pos, k, ψ, steps) in abstractions
    A, tbuild = timed(() -> Abstraction(W, ψ; steps))
    h, tgoal = timed(() -> lower_bound(A, goal))
    expanded, tq = 0, 0.0
    for s in sources
        r, t = timed(() -> shortest_path(W, s, goal; heuristic = h, out_edges = Wout))
        r.cost == V[s] || error("A* is not exact")
        expanded += r.expanded
        tq += t
    end
    q = length(sources)
    vel = steps == 1 ? "1" : string(steps)
    push!(results, (kind = kind, pos = pos, vel = vel, steps = k, blocks = ψ.n_target, build = tbuild, goal_setup = tgoal,
        bound = sum(h[sources] ./ V[sources]) / q, expanded = expanded / q, query = tq / q))
end

println("\n1. exact queries (averages over the 10 queries; build once per abstraction, bound once per goal)")
println("abstraction  position blocks  velocity blocks  steps  coarse states  build      bound/goal  bound/cost  expansions  query")
for r in results
    @printf("%-12s %-16s %-16s %-6s %13d %8.3f s %9.4f s %10.2f %11.0f %7.4f s\n", r.kind, r.pos, r.vel, r.steps, r.blocks,
        r.build, r.goal_setup, r.bound, r.expanded, r.query)
end
open(joinpath(@__DIR__, "..", "figures", "tropical_abstraction.csv"), "w") do io
    println(io, "abstraction,position_blocks,velocity_blocks,steps,coarse_states,build_seconds,bound_seconds,bound_over_cost,expansions,query_seconds")
    for r in results
        println(io, join((r.kind, r.pos, r.vel, r.steps, r.blocks, r.build, r.goal_setup, r.bound, r.expanded, r.query), ","))
    end
end

# ---- 2. coarse closures for many goals -------------------------------------------
println("\n2. the coarse closure for 20 goals (stopped at random free cells)")
goals = [[L[qx, qy, stop, stop]] for (qx, qy) in zip(rand(rng, 2:(n - 1), 40), rand(rng, 2:(n - 1), 40)) if !walls(qx, qy)][1:20]
println("abstraction  coarse states  dijkstra/goal  pathpack factor  pathpack/goal  same bound")
for (kind, pos, k, ψ, steps) in abstractions[[2, 3, 6, 7]]
    Ad = Abstraction(W, ψ; steps)
    Ap, tf = timed(() -> Abstraction(W, ψ; steps, method = :pathpack))
    td = sum(timed(() -> lower_bound(Ad, gl))[2] for gl in goals) / length(goals)
    tp = sum(timed(() -> lower_bound(Ap, gl))[2] for gl in goals) / length(goals)
    same = all(lower_bound(Ad, gl) == lower_bound(Ap, gl) for gl in goals)
    @printf("%-12s %13d %12.4f s %14.3f s %12.4f s  %s\n", "$(kind) $(pos)", ψ.n_target, td, tf, tp, same)
end

# ---- 3. reach-avoid ------------------------------------------------------------
println("\n3. reach-avoid: the route to the goal that keeps the most clearance from the walls")
wallcells = [(qx, qy) for qx in 1:n, qy in 1:n if walls(qx, qy)]
clear = [walls(qx, qy) ? -1.0 : minimum(hypot(qx - a, qy - b) for (a, b) in wallcells) - 0.5 for qx in 1:n, qy in 1:n]
l = [clear[I[1], I[2]] for I in CartesianIndices(points)][:]
g = fill(-Inf, length(l)); g[goal] .= Inf
R, t = timed(() -> reach_avoid(W, l, g))
@printf("reach-avoid closure: %.3f s; reachable set agrees with min-plus: %s\n", t, (R .> -Inf) == isfinite.(V))
println("source   min time  clearance of the fastest route  best clearance  time on the safest route (clearance ≥ best)")
rows, cols = findnz(W)[1:2]
for s in sources[1:5]
    fast = shortest_path(W, s, goal; out_edges = Wout)
    cfast = minimum(l[fast.path])
    c = R[s]
    keep = [l[x] >= c for x in rows]
    Wc = sparse(rows[keep], cols[keep], ones(count(keep)), size(W)...)
    safe = shortest_path(Wc, s, [t for t in goal if l[t] >= c])
    @printf("%-8d %8.0f %31.2f %15.2f %15.0f\n", s, fast.cost, cfast, c, safe.cost)
end
