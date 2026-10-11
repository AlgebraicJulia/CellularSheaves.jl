# Aggregation and elimination in min-plus on the double integrator lattice.
#
#   JULIA_NUM_THREADS=8 julia --project=docs docs/scripts/tropical_overlay_study.jl
#
# 1. All-pairs tables of the coarse closure (pattern databases), by PathPack's closure
#    and by one Dijkstra per column, and the per-goal bound from a table.
# 2. The overlay: exact elimination of the interiors of 16 × 16 position cells onto
#    their portals, with queries on the overlay (with and without an aggregation
#    heuristic), as the velocity range, and with it the portal set, grows.
# 3. Commutation: on a saturated separator S, aggregate-then-eliminate against
#    eliminate-then-aggregate, as operators and as cost-to-go, against the exact
#    aggregated cost.

using CellularSheaves
using PATHPACK
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

function problem(n, vmax)
    m = n ÷ 3
    walls(qx, qy) = (qx in m:(m + 1) && !(qy in 2:(n ÷ 6))) || (qx in (2m):(2m + 1) && !(qy in (n - n ÷ 6):(n - 1)))
    W, points = double_integrator_lattice(n, vmax, walls)
    L = LinearIndices(points)
    stop = vmax + 1
    goal = [L[qx, qy, stop, stop] for qx in (n - 3):(n - 1) for qy in 2:4]
    V = cost_to_go(W, goal)
    rng = Random.MersenneTwister(1)
    candidates = [L[qx, qy, stop, stop] for qx in 1:(m - 1) for qy in 1:n if isfinite(V[L[qx, qy, stop, stop]])]
    sources = candidates[randperm(rng, length(candidates))[1:10]]
    return (; W, points, L, goal, V, sources, walls)
end

function coarsen_grid(ψ, sizes, factors)
    return compose(ψ, aggregation_homomorphism(sizes, factors)), coarse_points(sizes, factors)
end
identity_hom(N) = GraphHomomorphism(collect(1:N), N)
function position_blocks(points, levels)
    ψ, sizes = identity_hom(prod(points)), points
    for _ in 1:levels
        ψ, sizes = coarsen_grid(ψ, sizes, (2, 2, 1, 1))
    end
    return ψ
end
function spacetime_blocks(points, ks)
    ψ, sizes = identity_hom(prod(points)), points
    for _ in ks
        ψ, sizes = coarsen_grid(ψ, sizes, (2, 2, 2, 2))
        ψ, sizes = coarsen_grid(ψ, sizes, (2, 2, 1, 1))
    end
    return ψ
end
timed(f) = (f(); t = @elapsed r = f(); (r, t))
mean_query(f, sources) = (rs = [f(s) for s in sources]; t = sum(timed(() -> f(s))[2] for s in sources) / length(sources);
    (expanded = mean(r.expanded for r in rs), seconds = t, costs = [r.cost for r in rs]))

println("threads: $(Threads.nthreads())")

# ---- 1. all-pairs tables ---------------------------------------------------------
P = problem(128, 6)
println("\n1. all-pairs coarse closure (pattern database), lattice 128² × 13²")
println("abstraction              coarse states  dijkstra table  pathpack table  bound/goal: solve   table   same")
rows1 = []
for (name, ψ, steps) in (("space–time 16×16 k=4", spacetime_blocks(P.points, (2, 4)), 4),
                         ("space 32×32", position_blocks(P.points, 5), 1),
                         ("space 16×16", position_blocks(P.points, 4), 1))
    A = Abstraction(P.W, ψ; steps)
    td = @elapsed Ad = Abstraction(P.W, ψ; steps, table = true)
    tp = @elapsed Ap = Abstraction(P.W, ψ; steps, method = :pathpack, table = true)
    tb = @elapsed Abstraction(P.W, ψ; steps)                     # the quotient alone
    hs, ts = timed(() -> lower_bound(A, P.goal))
    ht, tt = timed(() -> lower_bound(Ad, P.goal))
    same = hs == ht == lower_bound(Ap, P.goal) && Ad.table == Ap.table
    @printf("%-24s %13d %13.3f s %13.3f s %16.4f s %7.4f s  %s\n", name, ψ.n_target, td - tb, tp - tb, ts, tt, same)
    push!(rows1, (name, ψ.n_target, td - tb, tp - tb, ts, tt))
end

# ---- 2. the overlay -------------------------------------------------------------
println("\n2. overlay: cell interiors (16 × 16 positions × all velocities) eliminated onto portals")
println("lattice        states  portals   shortcuts  build     query: dijkstra   A* (2×2)   overlay   overlay+A*  (states expanded)")
rows2 = []
for (n, vmax) in ((128, 1), (64, 2), (64, 3))
    Q = problem(n, vmax)
    cell = [((I[1] - 1) ÷ 16) * cld(n, 16) + (I[2] - 1) ÷ 16 + 1 for I in CartesianIndices(Q.points)][:]
    O, tbuild = timed(() -> Overlay(Q.W, cell))
    Wout = sparse(transpose(Q.W))
    h = lower_bound(Abstraction(Q.W, position_blocks(Q.points, 1)), Q.goal)
    d = mean_query(s -> shortest_path(Q.W, s, Q.goal; out_edges = Wout), Q.sources)
    a = mean_query(s -> shortest_path(Q.W, s, Q.goal; heuristic = h, out_edges = Wout), Q.sources)
    o = mean_query(s -> overlay_path(O, s, Q.goal), Q.sources)
    oa = mean_query(s -> overlay_path(O, s, Q.goal; heuristic = h), Q.sources)
    all(c -> c == Q.V[Q.sources], (d.costs, a.costs, o.costs, oa.costs)) || error("a query is not exact")
    @printf("%3d² × %2d²  %9d %8d %11d %7.2f s %10.0f %10.0f %10.0f %10.0f\n", n, 2vmax + 1, size(Q.W, 1),
        length(O.portals), nnz(O.shortcuts_out), tbuild, d.expanded, a.expanded, o.expanded, oa.expanded)
    @printf("%54s %10.4f s %8.4f s %8.4f s %8.4f s  (seconds)\n", "", d.seconds, a.seconds, o.seconds, oa.seconds)
    push!(rows2, (n, vmax, size(Q.W, 1), length(O.portals), nnz(O.shortcuts_out), tbuild, d, a, o, oa))
end

# ---- 3. commutation ---------------------------------------------------------------
println("\n3. commutation of aggregation (2 × 2 positions) and elimination (cell interiors), lattice 64² × 5²")
Q = problem(64, 2)
cell = [((I[1] - 1) ÷ 16) * 4 + (I[2] - 1) ÷ 16 + 1 for I in CartesianIndices(Q.points)][:]
portal = falses(size(Q.W, 1))
for (x, y) in zip(findnz(Q.W)[1:2]...)
    cell[x] != cell[y] && (portal[x] = portal[y] = true)
end
ψ = position_blocks(Q.points, 1)
Sc = sort(unique(ψ.vertex_map[portal]))                         # blocks that meet a portal
inSc = falses(ψ.n_target); inSc[Sc] .= true
S = findall(x -> inSc[ψ.vertex_map[x]], 1:size(Q.W, 1))         # saturated separator
position = zeros(Int, ψ.n_target); position[Sc] .= eachindex(Sc)
ψS = GraphHomomorphism(position[ψ.vertex_map[S]], length(Sc))
EA, tEA = timed(() -> quotient_graph(schur_complement(Q.W, S), ψS))          # eliminate, then aggregate
AE, tAE = timed(() -> schur_complement(quotient_graph(Q.W, ψ), Sc))         # aggregate, then eliminate
ea_rows, ea_cols, ea_vals = findnz(EA)
lax_ok = all(AE[a, b] != 0 && AE[a, b] <= v + 1e-12 for (a, b, v) in zip(ea_rows, ea_cols, ea_vals) if a != b)
equal = count(AE[a, b] == v for (a, b, v) in zip(ea_rows, ea_cols, ea_vals) if a != b)
# Cost-to-go to the separator blocks of the goal cell: both orders against the exact
# aggregated cost min_{x ∈ a} V(x).
target_blocks = [position[b] for b in Sc if any(x -> cell[x] == cell[Q.goal[1]] && portal[x], findall(==(b), ψ.vertex_map))]
vAE = cost_to_go(AE, target_blocks)
vEA = cost_to_go(EA, target_blocks)
target_states = [x for x in S if position[ψ.vertex_map[x]] in target_blocks]
Vexact = tropical_pushforward(ψS, cost_to_go(Q.W, target_states)[S])
ok = isfinite.(Vexact) .& (Vexact .> 0)
@printf("separator: %d of %d blocks (%d states); operators: eliminate→aggregate %d entries (%.2f s), aggregate→eliminate %d entries (%.2f s)\n",
    length(Sc), ψ.n_target, length(S), nnz(EA), tEA, nnz(AE), tAE)
@printf("lax square AE ≤ EA on every entry of EA: %s; equal on %d of %d entries\n", lax_ok, equal, count(ea_rows .!= ea_cols))
@printf("cost-to-go / exact over %d separator blocks: aggregate→eliminate %.3f (mean), eliminate→aggregate %.3f (mean); EA exact on %.1f%%, AE exact on %.1f%%\n",
    count(ok), mean(vAE[ok] ./ Vexact[ok]), mean(vEA[ok] ./ Vexact[ok]),
    100 * mean(vEA[ok] .== Vexact[ok]), 100 * mean(vAE[ok] .== Vexact[ok]))
@printf("order checks: AE ≤ EA ≤ exact everywhere: %s\n", all(vAE .<= vEA .+ 1e-12) && all(vEA .<= Vexact .+ 1e-12))

open(joinpath(@__DIR__, "..", "figures", "tropical_overlay.csv"), "w") do io
    println(io, "n,vmax,states,portals,shortcuts,build_seconds,dijkstra_expanded,astar_expanded,overlay_expanded,overlay_astar_expanded,dijkstra_seconds,astar_seconds,overlay_seconds,overlay_astar_seconds")
    for (n, vmax, N, np, ns, tb, d, a, o, oa) in rows2
        println(io, join((n, vmax, N, np, ns, tb, d.expanded, a.expanded, o.expanded, oa.expanded, d.seconds, a.seconds, o.seconds, oa.seconds), ","))
    end
end
open(joinpath(@__DIR__, "..", "figures", "tropical_tables.csv"), "w") do io
    println(io, "abstraction,coarse_states,dijkstra_table_seconds,pathpack_table_seconds,bound_solve_seconds,bound_table_seconds")
    for r in rows1
        println(io, join(r, ","))
    end
end
