# Many right-hand sides of the coarse min-plus closure: m target sets at once.
#
#   JULIA_NUM_THREADS=8 julia --project=docs docs/scripts/tropical_multi_rhs.jl [n] [vmax]
#
# For the coarse graphs of the double integrator lattice, the cost-to-go to m goals
# is the m columns V = Wc* ⊗ B, B[:, j] the terminal cost of goal j. Dijkstra solves
# them one goal at a time (threaded over goals); PathPack factorizes Wc* = U* L* once
# and applies it to the whole block B (threaded inside the solve). m = all coarse
# states is the full closure Wc*, the cost between every pair of states.

using CellularSheaves
using PATHPACK
using LinearAlgebra
using SparseArrays
using Printf
using Random
using Base.Threads

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

function dijkstra_block(Wc, B)
    V = similar(B)
    @threads for j in axes(B, 2)
        V[:, j] = minplus_closure(Wc, view(B, :, j))
    end
    return V
end

pathpack_block(F, B) = lmul!(F, copy(B))

n = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 128
vmax = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 6
m3 = n ÷ 3
walls(qx, qy) = (qx in m3:(m3 + 1) && !(qy in 2:(n ÷ 6))) || (qx in (2m3):(2m3 + 1) && !(qy in (n - n ÷ 6):(n - 1)))
W, points = double_integrator_lattice(n, vmax, walls)
println("lattice $(n)² × $(2vmax + 1)², $(size(W, 1)) states; $(nthreads()) threads")

function coarsen_grid(ψ, sizes, factors)
    return compose(ψ, aggregation_homomorphism(sizes, factors)), coarse_points(sizes, factors)
end
id = GraphHomomorphism(collect(1:prod(points)), prod(points))
cases = Any[]
ψ, sizes = id, points
for ℓ in 1:5
    global ψ, sizes = coarsen_grid(ψ, sizes, (2, 2, 1, 1))
    ℓ >= 3 && push!(cases, ("space $(2^ℓ)×$(2^ℓ)", ψ, 1))
end
ψ, sizes = id, points
for k in (2, 4)
    global ψ, sizes = coarsen_grid(ψ, sizes, (2, 2, 2, 2))
    global ψ, sizes = coarsen_grid(ψ, sizes, (2, 2, 1, 1))
    k == 4 && push!(cases, ("space–time 16×16 k=4", ψ, k))
end

rng = Random.MersenneTwister(2)
results = []
println("coarse graph            states   edges  factor     m      dijkstra   pathpack   speed-up  same")
for (name, ψ, steps) in cases
    Wc = Abstraction(W, ψ; steps).quotient
    nc = size(Wc, 1)
    tf = @elapsed F = PATHPACK.CPU.mlu(PATHPACK.CPU.MinPlus(), Wc)
    for m in unique(min.(nc, [1, 16, 128, 1024, nc <= 11000 ? nc : 4096]))   # all pairs when it fits in memory
        targets = m == nc ? collect(1:nc) : randperm(rng, nc)[1:m]
        B = fill(Inf, nc, m)
        for (j, t) in enumerate(targets)
            B[t, j] = 0.0
        end
        Vd = dijkstra_block(Wc, B); Vp = pathpack_block(F, B)          # compile, check
        td = @elapsed dijkstra_block(Wc, B)
        tp = @elapsed pathpack_block(F, B)
        same = Vd == Vp
        @printf("%-22s %7d %7d %6.2f s %6d %10.4f s %9.4f s %8.1f×  %s\n", name, nc, nnz(Wc), tf, m, td, tp, td / tp, same)
        push!(results, (name = name, states = nc, edges = nnz(Wc), factor = tf, m = m, dijkstra = td, pathpack = tp, same = same))
    end
end

open(joinpath(@__DIR__, "..", "figures", "tropical_multi_rhs.csv"), "w") do io
    println(io, "coarse_graph,coarse_states,coarse_edges,factor_seconds,right_hand_sides,dijkstra_seconds,pathpack_seconds,threads")
    for r in results
        println(io, join((r.name, r.states, r.edges, r.factor, r.m, r.dijkstra, r.pathpack, nthreads()), ","))
    end
end
