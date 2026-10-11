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

# Velocity modes: rest, then `headings` direction sectors for each speed band
# (Chebyshev speed ≤ speed_edges[k]).
function heading_speed_modes(vmax, headings, speed_edges)
    M = zeros(Int, 2vmax + 1, 2vmax + 1)
    for ix in 1:(2vmax + 1), iy in 1:(2vmax + 1)
        vx, vy = ix - vmax - 1, iy - vmax - 1
        if vx == 0 && vy == 0
            M[ix, iy] = 1
            continue
        end
        band = findfirst(>=(max(abs(vx), abs(vy))), speed_edges)
        h = mod(round(Int, atan(vy, vx) / (2π / headings)), headings)
        M[ix, iy] = 1 + (band - 1) * headings + h + 1
    end
    return M
end

# Planners: the velocity-blind projection (one mode per stalk) and truncated
# pushforwards keeping heading × speed modes in the stalk over each cell.
ψ = projection_homomorphism(points, (1, 2))
planners = Any[]
P0, t0 = timed(() -> LayeredPlanner(W, ψ, (n, n)))
push!(planners, ("velocity-blind (1 mode)", P0, t0))
for (headings, speeds) in ((8, [2, 4, vmax]), (16, [2, 4, vmax]))
    modes = heading_speed_modes(vmax, headings, speeds)
    T = TruncatedPushforward(W, points, (1, 2), modes)
    Pm, tm = timed(() -> LayeredPlanner(T))
    push!(planners, ("$(headings) headings × $(length(speeds)) speeds ($(T.K) modes)", Pm, tm))
end
h = Dict{String,Vector{Float64}}()
for (name, Pl, tb) in planners
    hb, tbound = timed(() -> lower_bound(Pl.abstraction, goal))
    h[name] = hb
    troute = mean(timed(() -> plan_route(Pl, s, goal))[2] for s in sources)
    @printf("planner %-34s %8d coarse states, built in %5.2f s; bound for the goal %.3f s (bound/cost %.2f); route %.4f s\n",
        name, Pl.μ.n_target, tb, tbound, mean(hb[sources] ./ V[sources]), troute)
end

Wout = sparse(transpose(W))
function summarize(name, rs, ts)
    costs = [r.cost for r in rs]
    ratio = costs ./ V[sources]
    first = [hasproperty(r, :first_try) ? r.first_try : true for r in rs]
    states = [hasproperty(r, :corridor_states) ? r.corridor_states : count(isfinite, V) for r in rs]
    @printf("%-52s %6.0f%% %9.3f %7.3f %8.0f%% %12.0f %11.0f %9.4f s\n", name, 100 * mean(first), mean(ratio), maximum(ratio),
        100 * mean(ratio .== 1), mean(states), mean(r.expanded for r in rs), mean(ts))
    return (name = name, first_try = mean(first), mean_ratio = mean(ratio), max_ratio = maximum(ratio),
        optimal = mean(ratio .== 1), states = mean(states), expanded = mean(r.expanded for r in rs), seconds = mean(ts))
end

println("\nmethod                                               first try  cost/opt   worst  optimal  search states  expansions  time/query")
rows = []
q(f) = (rs = [f(s) for s in sources]; ts = [timed(() -> f(s))[2] for s in sources]; (rs, ts))
push!(rows, summarize("4-D Dijkstra (exact)", q(s -> shortest_path(W, s, goal; out_edges = Wout))...))
for (name, Pl, _) in planners
    push!(rows, summarize("4-D A*, bound of $(name)", q(s -> shortest_path(W, s, goal; heuristic = h[name], out_edges = Wout))...))
end
for (name, Pl, _) in planners, (width, mode_width) in ((0, 0), (1, 0), (2, 0), (2, 1), (4, 1), (8, 1), (16, 1))
    Pl.K == 1 && (width, mode_width) == (2, 1) && continue          # one mode: m is moot
    label = "$(split(name, " (")[1]): w=$(width)" * (Pl.K > 1 ? ", m=$(mode_width)" : "")
    push!(rows, summarize(label, q(s -> layered_path(Pl, s, goal; width, mode_width, heuristic = h[name]))...))
end

open(joinpath(@__DIR__, "..", "figures", "tropical_layered.csv"), "w") do io
    println(io, "method,first_try,mean_cost_over_optimal,max_cost_over_optimal,fraction_optimal,search_states,expansions,seconds_per_query")
    for r in rows
        println(io, join((r.name, r.first_try, r.mean_ratio, r.max_ratio, r.optimal, r.states, r.expanded, r.seconds), ","))
    end
end

# ---- certified layered control ----------------------------------------------------
# Each query returns the closed-loop trajectory, a Lyapunov certificate for the corridor
# controller (lax side), and a suboptimality bound from the optimistic refinement
# (oplax side), computed without the optimum; here the optimum is known, so the bound
# can be compared with the true ratio.
println("\ncertified layered control: Lyapunov certificate and guaranteed cost/optimum bound")
println("planner                       corridor   certified  first try  cost/opt (true)  bound (mean, worst)  time/query")
cert_rows = []
for (name, Pl, _) in planners, (width, mode_width) in ((4, 1), (16, 1))
    rs = [certified_path(Pl, s, goal; width, mode_width, lower = h[name]) for s in sources]
    ts = [timed(() -> certified_path(Pl, s, goal; width, mode_width, lower = h[name]))[2] for s in sources]
    ok = mean(r.certificate.positive_definite && r.certificate.strict && r.certificate.closed for r in rs)
    truth = [r.cost for r in rs] ./ V[sources]
    bound = [r.suboptimality_bound for r in rs]
    all(bound .>= truth .- 1e-12) || error("a bound below the true ratio")
    @printf("%-29s w=%2d m=%d %8.0f%% %9.0f%% %12.3f %12.3f %8.3f %10.4f s\n", split(name, " (")[1], width, mode_width,
        100 * ok, 100 * mean(r.first_try for r in rs), mean(truth), mean(bound), maximum(bound), mean(ts))
    push!(cert_rows, (split(name, " (")[1], width, mode_width, ok, mean(r.first_try for r in rs), mean(truth), mean(bound),
        maximum(bound), mean(ts)))
end
open(joinpath(@__DIR__, "..", "figures", "tropical_certified.csv"), "w") do io
    println(io, "planner,width,mode_width,certified_fraction,first_try,mean_cost_over_optimal,mean_bound,max_bound,seconds_per_query")
    for r in cert_rows
        println(io, join(r, ","))
    end
end
