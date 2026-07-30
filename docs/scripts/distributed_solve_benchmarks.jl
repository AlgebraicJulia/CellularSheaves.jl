# Benchmarks for the distributed sheaf solve feature guide.
#
# Everything here is measured on the *actual* harmonic-extension system Layer 1
# solves: build a sheaf Laplacian over a formation graph, pin a few vertices as
# targets (boundary data), and take H = the free–free (Dirichlet) block, which
# is symmetric positive definite. All four figures operate on H and its
# CliqueTrees.Multifrontal chordal factorization.
#
# Run with:  julia --project=docs docs/scripts/distributed_solve_benchmarks.jl
# Figures are written to docs/figures/distributed_solve/*.svg.

# Render GR figures straight to file — no on-screen gksqt window that has to be
# closed by hand. Must be set before Plots/GR initializes. ("100" = offscreen.)
get!(ENV, "GKSwstype", "100")

using CellularSheaves
using CellularSheaves.NetworkSheaves.DistributedSolve
using CliqueTrees.Multifrontal
using LinearAlgebra
using SparseArrays
using Graphs
using Plots
using LaTeXStrings
using Printf

const FIG = joinpath(@__DIR__, "..", "figures", "distributed_solve")
mkpath(FIG)

default(thickness_scaling = 1.4, legendfontsize = 8, guidefontsize = 10,
    titlefontsize = 9, tickfontsize = 8, framestyle = :box, grid = true,
    gridalpha = 0.22, markerstrokewidth = 0, size = (900, 380))
const PANELS = (layout = (1, 2), size = (1000, 400), left_margin = 5Plots.mm,
    bottom_margin = 5Plots.mm, top_margin = 4Plots.mm)
# house palette (ext/CellularSheavesPlots/src/Utils.jl)
const PAL = [:steelblue, :darkorange, :green, :crimson]

# ---------------------------------------------------------------------------
# problem construction
# ---------------------------------------------------------------------------

# Sheaf Laplacian of `g` with uniform stalk dimension `stalk` and identity
# restriction maps (the constant/consensus sheaf), as a sparse matrix.
function sheaf_laplacian(g; stalk = 2)
    s = sheaf_from_graph(g, stalk, d -> Matrix{Float64}(I, d, d); symmetric_edges = true)
    return sparse(sheaf_laplacian_matrix(s))
end

# The harmonic-extension (Dirichlet) system: pin `pinned` vertices as targets
# and return H, the free–free block — SPD whenever every component of the free
# subgraph touches a pinned vertex.
function dirichlet_block(g, pinned; stalk = 2)
    L = sheaf_laplacian(g; stalk)
    n = nv(g)
    off = (0:n) .* stalk
    freev = setdiff(1:n, pinned)
    idx = reduce(vcat, [collect((off[v] + 1):off[v + 1]) for v in freev])
    return L[idx, idx], freev
end

# four corners + center: a well-spread pin set that keeps H well-defined
function grid_pins(side)
    c = (side * side + 1) ÷ 2
    return unique([1, side, side * (side - 1) + 1, side * side, c])
end

chol(H) = cholesky!(ChordalCholesky(H), NoPivot())

function mintime(f; r = 40, w = 3)
    for _ in 1:w
        f()
    end
    return minimum(@elapsed(f()) for _ in 1:r)
end

# Write the SVG the docs use, plus a PNG copy for local inspection when
# BENCH_PNG names a directory (off by default, so the docs build is unaffected).
function savefig2(plt, name)
    savefig(plt, joinpath(FIG, name * ".svg"))
    pngdir = get(ENV, "BENCH_PNG", "")
    if !isempty(pngdir)
        mkpath(pngdir)
        savefig(plt, joinpath(pngdir, name * ".png"))
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Figure 1 — solve scaling and the factor-once / re-solve-every-step split
# ---------------------------------------------------------------------------

function figure_solve_scaling()
    sides = [8, 12, 16, 20, 24, 28, 32]
    dofs = Int[]
    t_ct = Float64[]   # CliqueTrees monolithic ldiv!
    t_rec = Float64[]  # recursive tree solve
    t_ws = Float64[]   # workspace tree solve
    t_fac = Float64[]  # one-time factorization
    maxres = 0.0

    for side in sides
        g = Graphs.grid([side, side])
        H, _ = dirichlet_block(g, grid_pins(side); stalk = 2)
        n = size(H, 1)
        F = chol(H)
        L = F.L
        ws = TreeWorkspace(L, Float64)
        b = rand(n)

        ref = copy(b)
        ldiv!(L, ref); ldiv!(L', ref)
        wsr = copy(b)
        tree_forward_ldiv!(L, wsr, ws); tree_backward_ldiv!(L, wsr, ws)
        maxres = max(maxres, norm(wsr - ref) / norm(ref))

        push!(dofs, n)
        push!(t_ct, mintime(() -> (c = copy(b); ldiv!(L, c); ldiv!(L', c))))
        push!(t_rec, mintime(() -> (c = copy(b); tree_forward_ldiv!(L, c); tree_backward_ldiv!(L, c))))
        push!(t_ws, mintime(() -> (c = copy(b); tree_forward_ldiv!(L, c, ws); tree_backward_ldiv!(L, c, ws))))
        push!(t_fac, mintime(() -> chol(H); r = 10))
    end

    p1 = plot(xscale = :log10, yscale = :log10, xlabel = L"degrees of freedom, $n$",
        ylabel = "per-solve time (s)", legend = :topleft, title = "Triangular solve cost")
    plot!(p1, dofs, t_rec, label = "tree solve (recursive)", color = PAL[4], marker = :diamond, lw = 2)
    plot!(p1, dofs, t_ct, label = "CliqueTrees monolithic", color = PAL[1], marker = :circle, lw = 2)
    plot!(p1, dofs, t_ws, label = "tree solve (workspace)", color = PAL[2], marker = :square, lw = 2)

    p2 = plot(xscale = :log10, yscale = :log10, xlabel = L"degrees of freedom, $n$",
        ylabel = "time (s)", legend = :topleft, title = "Factor once, re-solve per step")
    plot!(p2, dofs, t_fac, label = "factorization (once)", color = PAL[1], marker = :circle, lw = 2)
    plot!(p2, dofs, t_ws, label = "re-solve (per step)", color = PAL[2], marker = :square, lw = 2)

    plt = plot(p1, p2; PANELS...)
    savefig2(plt, "solve_scaling")
    @printf("[solve]  workspace closes gap to %.2f× (recursive was %.2f×); max relative residual %.1e\n",
        t_ws[end] / t_ct[end], t_rec[end] / t_ct[end], maxres)
    return maxres
end

# ---------------------------------------------------------------------------
# Communication model (shared by Figure 2)
#
# Half-duplex radios, one packet per time slot per node: in any slot a node is
# in exactly one message (as sender or receiver). A "makespan" is the number of
# slots to finish, i.e. the critical path of the message schedule. For a gather
# on a rooted tree, a parent receives its children's packets one per slot,
# earliest-ready first, and can only forward upward once it has them all — this
# `gather_makespan` computes exactly that, charging sibling serialization.
# ---------------------------------------------------------------------------

function gather_makespan(children_of, roots)
    memo = Dict{Int, Int}()
    function rd(j)
        haskey(memo, j) && return memo[j]
        ch = children_of(j)
        isempty(ch) && return (memo[j] = 0)
        t = 0
        for rc in sort!([rd(c) for c in ch])
            t = max(t, rc) + 1
        end
        return (memo[j] = t)
    end
    return maximum(rd(r) for r in roots)
end

# Tree method: gather up the supernode elimination tree then scatter back down
# (symmetric), so twice the gather makespan.
function tree_slots(L)
    S = L.S
    roots = [j for j in 1:nv(S.res) if iszero(S.pnt[j])]
    return 2 * gather_makespan(j -> collect(neighbors(S.chd, j)), roots)
end

# Centralized: a coordinator gathers every agent's state (one per slot, its
# single radio serializes) and scatters the answer back — the star tree.
central_slots(n_agents) = 2 * n_agents

# Sheaf diffusion: Richardson iteration on H with the optimal step, so the
# contraction factor is ρ = (κ−1)/(κ+1); reaching tolerance ε needs
# K = ⌈log ε / log ρ⌉ rounds, and each round every node exchanges with all its
# neighbors — Δ slots under a proper edge coloring.
function diffusion_slots(H, maxdeg; ε = 1e-6)
    ev = eigvals(Symmetric(Matrix(H)))
    κ = ev[end] / ev[1]
    ρ = (κ - 1) / (κ + 1)
    K = ceil(Int, log(ε) / log(ρ))
    return K * maxdeg, κ, K
end

# ---------------------------------------------------------------------------
# Figure 2 — communication makespan: grid vs chain
# ---------------------------------------------------------------------------

function comm_panel(build, ns; title, xlabel)
    tree = Int[]; cent = Int[]; diff = Int[]
    for n in ns
        g = build(n)
        pinned = title == "Chain formation" ? [1, nv(g)] : grid_pins(isqrt(nv(g)))
        H, freev = dirichlet_block(g, pinned; stalk = 2)
        L = chol(H).L
        push!(tree, tree_slots(L))
        push!(cent, central_slots(length(freev)))
        d, = diffusion_slots(H, Δ(g))
        push!(diff, d)
    end
    p = plot(xscale = :log10, yscale = :log10, xlabel = xlabel, ylabel = "message slots (makespan)",
        legend = :topleft, title = title)
    plot!(p, ns, diff, label = "sheaf diffusion", color = PAL[4], marker = :diamond, lw = 2)
    plot!(p, ns, cent, label = "centralized", color = PAL[1], marker = :circle, lw = 2)
    plot!(p, ns, tree, label = "distributed tree", color = PAL[2], marker = :square, lw = 2)
    return p, tree, cent, diff
end

function figure_communication()
    grid_ns = [8, 12, 16, 20, 24, 28, 32]
    pg, tg, cg, dg = comm_panel(s -> Graphs.grid([s, s]), grid_ns;
        title = "Grid formation", xlabel = "grid side")
    chain_ns = [16, 32, 64, 128, 256, 512]
    pc, tc, cc, dc = comm_panel(n -> Graphs.path_graph(n), chain_ns;
        title = "Chain formation", xlabel = L"agents, $n$")
    plt = plot(pg, pc; PANELS...)
    savefig2(plt, "communication")
    @printf("[comm]   grid 32²: tree %d vs central %d vs diffusion %d slots (%.0f× / %.0f×)\n",
        tg[end], cg[end], dg[end], cg[end] / tg[end], dg[end] / tg[end])
    @printf("[comm]   chain 512: tree %d vs central %d slots (%.1f×) — chain is the tree's worst case\n",
        tc[end], cc[end], cc[end] / tc[end])
end

# ---------------------------------------------------------------------------
# Figure 3 — low-rank edge update vs full refactorization
#
# Adding/removing one scalar coupling edge (i,j) changes H by the rank-1 term
# (e_i − e_j)(e_i − e_j)ᵀ — exactly a proximity edge for collision avoidance.
# Patching the existing factor (one lowrankupdate! + lowrankdowndate! = an edge
# on then off) versus discarding it and refactoring from scratch.
# ---------------------------------------------------------------------------

function figure_lowrank()
    sides = [12, 16, 20, 24, 28, 32, 40]
    dofs = Int[]; speedup = Float64[]; t_re = Float64[]; t_up = Float64[]
    for side in sides
        g = Graphs.grid([side, side])
        H, freev = dirichlet_block(g, grid_pins(side); stalk = 1)  # scalar edge = rank 1
        pos = Dict(v => k for (k, v) in enumerate(freev))
        e = first(e for e in edges(g) if haskey(pos, src(e)) && haskey(pos, dst(e)))
        v = zeros(size(H, 1)); v[pos[src(e)]] = 1.0; v[pos[dst(e)]] = -1.0

        tre = mintime(() -> ldlt!(ChordalLDLt(H), RowMaximum(); check = false); r = 8)
        F = ldlt!(ChordalLDLt(H), RowMaximum(); check = false)
        pair = mintime(() -> (lowrankupdate!(F, copy(v)); lowrankdowndate!(F, copy(v))); r = 40)
        tup = pair / 2   # one edge flip = one rank-1 modification

        push!(dofs, size(H, 1)); push!(t_re, tre); push!(t_up, tup); push!(speedup, tre / tup)
    end

    p1 = plot(xscale = :log10, yscale = :log10, xlabel = L"degrees of freedom, $n$",
        ylabel = "time (s)", legend = :topleft, title = "Patch vs. refactor")
    plot!(p1, dofs, t_re, label = "full refactorization", color = PAL[1], marker = :circle, lw = 2)
    plot!(p1, dofs, t_up, label = "rank-1 factor update", color = PAL[2], marker = :square, lw = 2)

    p2 = plot(xscale = :log10, xlabel = L"degrees of freedom, $n$", ylabel = "speedup (×)",
        legend = :topleft, title = "Speedup over refactoring")
    plot!(p2, dofs, speedup, label = "", color = PAL[3], marker = :square, lw = 2)

    plt = plot(p1, p2; PANELS...)
    savefig2(plt, "lowrank")
    @printf("[lowrank] rank-1 update is %.0f×–%.0f× faster than refactoring across sizes\n",
        minimum(speedup), maximum(speedup))
end

# ---------------------------------------------------------------------------
# Figure 4 — per-agent memory footprint
#
# The factor is partitioned into connected chunks (one per worker), each holding
# a disjoint slice of the factor data (no duplication). A centralized solve
# stores the whole factor on one node; a distributed solve stores only a slice
# per worker.
# ---------------------------------------------------------------------------

worker_mb(wd) = (length(wd.Dval) + length(wd.Lval)) * 8 / 1e6
factor_mb(L) = (length(L.Dval) + length(L.Lval)) * 8 / 1e6

function figure_memory()
    # panel a: per-worker footprint for one representative fleet
    side = 24
    g = Graphs.grid([side, side])
    H, _ = dirichlet_block(g, grid_pins(side); stalk = 2)
    L = chol(H).L
    W = 12
    p = partition_tree(L, W)
    nchunks = length(p.chunks)
    sizes = [worker_mb(worker_factorization(L, p, w)) for w in 1:nchunks]
    total = factor_mb(L)
    meansz = sum(sizes) / nchunks

    pa = bar(1:nchunks, sizes, label = "", color = PAL[1], xlabel = "worker",
        ylabel = "factor slice (MB)", title = @sprintf("Per-worker footprint, %d×%d grid", side, side),
        ylims = (0, maximum(sizes) * 1.35))
    hline!(pa, [meansz], color = PAL[4], ls = :dash, lw = 2, label = @sprintf("mean %.2f MB (Σ = %.2f MB = full factor)", meansz, sum(sizes)))

    # panel b: worst-worker vs centralized across sizes
    sides = [12, 16, 20, 24, 28, 32]
    ndofs = Int[]; cen = Float64[]; worst = Float64[]; mean_ = Float64[]
    for s in sides
        gg = Graphs.grid([s, s])
        HH, _ = dirichlet_block(gg, grid_pins(s); stalk = 2)
        LL = chol(HH).L
        pp = partition_tree(LL, 12)
        nc = length(pp.chunks)
        ws = [worker_mb(worker_factorization(LL, pp, w)) for w in 1:nc]
        push!(ndofs, size(HH, 1)); push!(cen, factor_mb(LL))
        push!(worst, maximum(ws)); push!(mean_, sum(ws) / nc)
    end
    pb = plot(xscale = :log10, yscale = :log10, xlabel = L"degrees of freedom, $n$", ylabel = "memory (MB)",
        legend = :topleft, title = "Distributed vs. centralized")
    plot!(pb, ndofs, cen, label = "centralized (one node)", color = PAL[1], marker = :circle, lw = 2)
    plot!(pb, ndofs, worst, label = "worst worker", color = PAL[2], marker = :square, lw = 2)
    plot!(pb, ndofs, mean_, label = "mean worker", color = PAL[3], marker = :diamond, lw = 2)

    plt = plot(pa, pb; PANELS...)
    savefig2(plt, "memory")
    @printf("[memory] %d×%d grid: worst worker %.2f MB vs centralized %.2f MB (%.0f%% of full)\n",
        side, side, maximum(sizes), total, 100 * maximum(sizes) / total)
end

# ---------------------------------------------------------------------------

function main()
    @info "Figure 1/3: solve scaling"
    figure_solve_scaling()
    @info "Figure 2/3: communication"
    figure_communication()
    @info "Figure 3/3: memory"
    figure_memory()
    @info "All figures written to $(FIG)"
    # `figure_lowrank` is kept above but intentionally not called: low-rank
    # factor updates concern a *changing* sheaf, which belongs to a future
    # dynamic-sheaf feature guide, not this fixed-system distributed solve.
end

main()
