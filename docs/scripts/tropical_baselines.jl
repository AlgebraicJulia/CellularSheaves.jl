# Baselines on external benchmark maps: a lattice double integrator on maps of the
# MovingAI MAPF benchmark suite (Stern et al. 2019; maps from Sturtevant 2012).
#
# Planners, all exact (A* with an admissible heuristic, reopening when needed):
#   - Dijkstra (no heuristic);
#   - kinematic: the standard double-integrator heuristic, the larger of the two
#     one-axis minimum times on an obstacle-free line (itself a pushforward along
#     the projection onto one axis, of a relaxed graph);
#   - velocity-blind abstraction: positions only (the analogue of HPA* / a
#     pattern database for this system);
#   - velocity-mode abstractions (heading × speed bins);
#   - kinematic ∨ modes: the larger of the two (still consistent).
# For each abstraction we also report the share of coarse arcs that lift (zero
# defect), on the viable subgraph (states from which some stop state is reachable).
#
#   julia -t 16 --project=docs docs/scripts/tropical_baselines.jl [vmax] [queries] [maps...]
# Maps are read from $MOVINGAI_DIR (default ~/data/movingai), downloaded from
# https://movingai.com/benchmarks/mapf/mapf-map.zip.

using CellularSheaves
using SparseArrays
using Printf
using Statistics
using Random

function read_map(path)
    lines = readlines(path)
    h = parse(Int, split(lines[2])[2])
    w = parse(Int, split(lines[3])[2])
    grid = lines[5:(4 + h)]
    blocked = falses(w, h)
    for y in 1:h, x in 1:w
        blocked[x, y] = !(grid[y][x] in ('.', 'G', 'S'))
    end
    return blocked
end

function double_integrator_on_map(blocked, vmax)
    w, h = size(blocked)
    nv = 2vmax + 1
    points = (w, h, nv, nv)
    L = LinearIndices(points)
    free(qx, qy) = 1 <= qx <= w && 1 <= qy <= h && !blocked[qx, qy]
    rows, cols = Int[], Int[]
    for I in CartesianIndices(points)
        qx, qy, ix, iy = Tuple(I)
        free(qx, qy) || continue
        vx, vy = ix - vmax - 1, iy - vmax - 1
        for ax in -1:1, ay in -1:1
            wx, wy = vx + ax, vy + ay
            (abs(wx) <= vmax && abs(wy) <= vmax) || continue
            px, py = qx + wx, qy + wy
            (free(px, py) && free(qx + wx ÷ 2, qy + wy ÷ 2)) || continue
            push!(rows, L[I]); push!(cols, L[px, py, wx + vmax + 1, wy + vmax + 1])
        end
    end
    N = prod(points)
    return sparse(rows, cols, ones(length(rows)), N, N), points
end

function heading_speed_modes(vmax, headings, speed_edges)
    M = zeros(Int, 2vmax + 1, 2vmax + 1)
    for ix in 1:(2vmax + 1), iy in 1:(2vmax + 1)
        vx, vy = ix - vmax - 1, iy - vmax - 1
        if vx == 0 && vy == 0
            M[ix, iy] = 1
            continue
        end
        band = findfirst(>=(max(abs(vx), abs(vy))), speed_edges)
        hd = mod(round(Int, atan(vy, vx) / (2π / headings)), headings)
        M[ix, iy] = 1 + (band - 1) * headings + hd + 1
    end
    return M
end

# One-axis minimum time T[Δ + R + 1, v + vmax + 1] to reach offset 0 at rest from
# offset Δ (target minus position) with velocity v, on an unbounded line (exact for
# |Δ| ≤ R − vmax²: no optimal path overshoots by more than a braking distance).
function one_axis_times(R, vmax)
    nv = 2vmax + 1
    T = fill(Inf, 2R + 1, nv)
    T[R + 1, vmax + 1] = 0
    frontier = [(0, 0)]
    t = 0
    while !isempty(frontier)
        t += 1
        next = Tuple{Int,Int}[]
        for (Δ, v) in frontier
            # predecessors (Δp, vp) with v ∈ vp ± 1 and Δ = Δp - v (position advances by v)
            for vp in (v - 1):(v + 1)
                abs(vp) <= vmax || continue
                Δp = Δ + v
                abs(Δp) <= R || continue
                if isinf(T[Δp + R + 1, vp + vmax + 1])
                    T[Δp + R + 1, vp + vmax + 1] = t
                    push!(next, (Δp, vp))
                end
            end
        end
        frontier = next
    end
    return T
end

function kinematic_heuristic(points, vmax, gx, gy, Tx, Ty, R)
    w, h, nv, _ = points
    hk = zeros(prod(points))
    L = LinearIndices(points)
    for I in CartesianIndices(points)
        qx, qy, ix, iy = Tuple(I)
        hk[L[I]] = max(Tx[gx - qx + R + 1, ix], Ty[gy - qy + R + 1, iy])
    end
    return hk
end

function restrict_map(μ, keep)
    ids = Dict{Int,Int}()
    return GraphHomomorphism([get!(ids, l, length(ids) + 1) for l in μ.vertex_map[keep]], length(ids))
end

vmax = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 3
nq = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 20
maps = length(ARGS) >= 3 ? ARGS[3:end] : ["random-32-32-20", "maze-32-32-4", "room-32-32-4", "room-64-64-8", "den312d", "warehouse-10-20-10-2-1"]
dir = get(ENV, "MOVINGAI_DIR", joinpath(homedir(), "data", "movingai"))
rng = MersenneTwister(2026)
results = []

for name in maps
    blocked = read_map(joinpath(dir, name * ".map"))
    Wfull, points = double_integrator_on_map(blocked, vmax)
    w, h, nv, _ = points
    L = LinearIndices(points)
    stop = vmax + 1
    rest = [L[qx, qy, stop, stop] for qx in 1:w, qy in 1:h if !blocked[qx, qy]]
    viable = findall(isfinite, cost_to_go(Wfull, rest))           # some stop state is reachable
    newindex = zeros(Int, size(Wfull, 1)); newindex[viable] = 1:length(viable)
    W = Wfull[viable, viable]
    out = sparse(transpose(W))
    R = max(w, h) + 2vmax^2
    Tx = one_axis_times(R, vmax)
    @printf("\n%s: %d × %d, %d free cells, vmax %d: %d states, %d viable\n", name, w, h, count(!, blocked), vmax,
        size(Wfull, 1), length(viable))

    abstractions = [("velocity-blind (1 mode)", ones(Int, nv, nv)),
                    ("8 headings × 2 speeds", heading_speed_modes(vmax, 8, [1, vmax])),
                    ("8 headings × $(vmax) speeds", heading_speed_modes(vmax, 8, collect(1:vmax)))]
    built = []
    for (aname, modes) in abstractions
        tb = @elapsed begin
            T = TruncatedPushforward(Wfull, points, (1, 2), modes)
            μ = restrict_map(T.μ, viable)
            A = Abstraction(W, μ)
        end
        ε = lifting_defect(W, μ; scope = :global, horizon = 4.0 * vmax)
        e = nonzeros(ε)
        push!(built, (aname, A, tb, count(iszero, e) / length(e)))
    end

    # queries: start and goal at rest on free cells, goal reachable
    restv = filter(!iszero, newindex[rest])
    queries = Tuple{Int,Int,Int,Int}[]                             # (source, goal state, gx, gy)
    while length(queries) < nq
        s, g = rand(rng, restv), rand(rng, restv)
        s == g && continue
        isfinite(cost_to_go(W, [g])[s]) || continue
        gx, gy = Tuple(CartesianIndices(points)[viable[g]])[1:2]
        push!(queries, (s, g, gx, gy))
    end

    methods = ["Dijkstra", "kinematic"]
    append!(methods, [b[1] for b in built])
    append!(methods, ["kinematic ∨ " * built[end][1]])
    stats = Dict(m => (Int[], Float64[]) for m in methods)
    costs = Dict(m => Float64[] for m in methods)
    # Query time includes building the heuristic for the query's goal.
    function run(m, s, g, gx, gy)
        heur = if m == "Dijkstra"
            nothing
        elseif m == "kinematic"
            kinematic_heuristic(points, vmax, gx, gy, Tx, Tx, R)[viable]
        elseif startswith(m, "kinematic ∨")
            max.(kinematic_heuristic(points, vmax, gx, gy, Tx, Tx, R)[viable], lower_bound(built[end][2], [g]))
        else
            lower_bound(built[findfirst(b -> b[1] == m, built)][2], [g])
        end
        return shortest_path(W, s, [g]; heuristic = heur, out_edges = out)
    end
    for m in methods                                               # warm up (compilation)
        run(m, queries[1]...)
    end
    for (s, g, gx, gy) in queries, m in methods
        t = @elapsed (c, _, ex) = run(m, s, g, gx, gy)
        push!(stats[m][1], ex); push!(stats[m][2], t); push!(costs[m], c)
    end
    @assert all(costs[m] ≈ costs["Dijkstra"] for m in methods) "every planner must return the optimum"
    println("method                               lift   build     expansions (mean)   time/query (mean)")
    for m in methods
        b = findfirst(x -> x[1] == m, built)
        lift = isnothing(b) ? "   –  " : @sprintf("%5.1f%%", 100built[b][4])
        tb = isnothing(b) ? "    –   " : @sprintf("%6.2f s", built[b][3])
        @printf("%-34s %s  %s   %12.0f        %10.4f s\n", m, lift, tb, mean(stats[m][1]), mean(stats[m][2]))
        push!(results, (name, length(viable), m, isnothing(b) ? NaN : built[b][4], isnothing(b) ? NaN : built[b][3],
            mean(stats[m][1]), mean(stats[m][2])))
    end
end

open(joinpath(@__DIR__, "..", "figures", "tropical_baselines_v$(vmax).csv"), "w") do io
    println(io, "map,viable_states,method,lift_fraction,build_seconds,mean_expansions,mean_seconds_per_query")
    for r in results
        println(io, join(r, ","))
    end
end
