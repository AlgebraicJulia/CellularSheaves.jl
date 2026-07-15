# The fair three-way communication comparison for the Distributed Sheaf Solve
# feature guide: distributed tree vs centralized coordinator vs sheaf diffusion,
# all charged under ONE half-duplex message-slot model, on five topologies.
# Also: the elimination-ordering study (default vs nested dissection), the
# physical routing overlay (hop distribution of tree messages), and a real
# multi-process wall-clock run.
#
#   julia --project=docs docs/scripts/distributed_solve_comparison.jl
#
# Writes docs/figures/distributed_solve/{comparison,ordering,hops}.svg and
# prints every table as markdown (paste into benchmarks.md / comparison.md).

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
using Random
using Distributed

const FIG = joinpath(@__DIR__, "..", "figures", "distributed_solve")
mkpath(FIG)

default(thickness_scaling = 1.4, legendfontsize = 8, guidefontsize = 10,
    titlefontsize = 10, tickfontsize = 8, framestyle = :box, grid = true,
    gridalpha = 0.22, markerstrokewidth = 0)
const PAL = Dict(
    :tree => :darkorange, :treend => :goldenrod,
    :cbs => :steelblue, :cin => :navy,
    :rich => :crimson, :cheb => :mediumvioletred, :cg => :seagreen)

const EPS = 1e-6

# ---------------------------------------------------------------------------
# problems: sheaf Laplacian Dirichlet blocks on five topologies (stalk 2)
# ---------------------------------------------------------------------------

function dirichlet_block(g, pinned; stalk = 2)
    s = sheaf_from_graph(g, stalk, d -> Matrix{Float64}(I, d, d); symmetric_edges = true)
    L = sparse(sheaf_laplacian_matrix(s))
    n = nv(g)
    off = (0:n) .* stalk
    freev = setdiff(1:n, pinned)
    idx = reduce(vcat, [collect((off[v] + 1):off[v + 1]) for v in freev])
    return L[idx, idx], freev
end

grid_pins(side) = unique([1, side, side * (side - 1) + 1, side * side, (side * side + 1) ÷ 2])

function rgg(n, seed)
    rng = Random.MersenneTwister(seed)
    # radius ~ sqrt(2 log n / n): connected w.h.p., mean degree ~ 2π log n / π…
    r = sqrt(2 * log(n) / n)
    pts = [(rand(rng), rand(rng)) for _ in 1:n]
    g = SimpleGraph(n)
    for i in 1:n, j in (i + 1):n
        (pts[i][1] - pts[j][1])^2 + (pts[i][2] - pts[j][2])^2 <= r^2 && add_edge!(g, i, j)
    end
    # RGG can still be disconnected at small n: bridge components along x order
    if !is_connected(g)
        ord = sortperm(pts; by = first)
        for k in 1:(n - 1)
            add_edge!(g, ord[k], ord[k + 1])
        end
    end
    return g
end

# pins that keep the free (agent) subgraph connected — the radio network must
# stay in one piece once the targets are removed as relays. Greedy: accept a
# candidate pin only if the remaining free subgraph is still connected.
function connected_pins(g, k)
    pins = Int[]
    for v in sortperm(degree(g))          # low-degree vertices first: safest to remove
        length(pins) == k && break
        cand = [pins; v]
        sg, _ = induced_subgraph(g, setdiff(1:nv(g), cand))
        is_connected(sg) && push!(pins, v)
    end
    @assert length(pins) == k "could not find $k connectivity-preserving pins"
    return pins
end

# five topologies, one accessor: name -> (graph builder, pin chooser, sizes).
# Rings are pinned at ONE vertex: pinning two opposite vertices would split the
# free radio graph into two arcs that physically cannot reach each other.
const TOPOLOGIES = [
    (:grid,  s -> Graphs.grid([s, s]),  g -> grid_pins(isqrt(nv(g))), [8, 12, 16, 20, 24]),
    (:chain, n -> path_graph(n),        g -> [1, nv(g)],              [16, 32, 64, 128, 256, 512]),
    (:ring,  n -> cycle_graph(n),       g -> [1],                     [16, 32, 64, 128, 256, 512]),
    (:star,  n -> star_graph(n),        g -> [2],                     [16, 64, 256]),
    (:rgg,   n -> rgg(n, 7),            g -> connected_pins(g, 2),    [32, 64, 128, 256]),
]

# ---------------------------------------------------------------------------
# nested-dissection elimination orders (pure Julia, no external dissector)
# ---------------------------------------------------------------------------

# 1-D: recursively eliminate the middle of each interval LAST
function nd_order_path(n)
    out = Int[]
    rec(lo, hi) = begin
        lo > hi && return
        mid = (lo + hi) ÷ 2
        rec(lo, mid - 1); rec(mid + 1, hi)
        push!(out, mid)
    end
    rec(1, n)
    return out
end

# 2-D grid (column-major vertex ids): recursively cut the longer dimension
function nd_order_grid(side)
    id(i, j) = (j - 1) * side + i
    out = Int[]
    rec(i1, i2, j1, j2) = begin
        (i1 > i2 || j1 > j2) && return
        if i2 - i1 >= j2 - j1
            m = (i1 + i2) ÷ 2
            rec(i1, m - 1, j1, j2); rec(m + 1, i2, j1, j2)
            append!(out, [id(m, j) for j in j1:j2])
        else
            m = (j1 + j2) ÷ 2
            rec(i1, i2, j1, m - 1); rec(i1, i2, m + 1, j2)
            append!(out, [id(i, m) for i in i1:i2])
        end
    end
    rec(1, side, 1, side)
    return out
end

# expand a vertex order to a dof order (stalk consecutive dofs per vertex),
# restricted to free vertices renumbered 1..|free|
function dof_order(vorder, freev, stalk)
    pos = Dict(v => k for (k, v) in enumerate(freev))
    out = Int[]
    for v in vorder
        haskey(pos, v) || continue
        k = pos[v]
        append!(out, ((k - 1) * stalk + 1):(k * stalk))
    end
    return out
end

# ---------------------------------------------------------------------------
# slot model (identical for every method)
# ---------------------------------------------------------------------------

# critical path of a combining reduce on a rooted forest, charging sibling
# serialization at each parent — used for the tree gather AND for all-reduces
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
    return maximum(rd(r) for r in roots; init = 0)
end

function tree_stats(L)
    S = L.S
    ns = Graphs.nv(S.res)
    roots = [j for j in 1:ns if iszero(S.pnt[j])]
    depth = zeros(Int, ns)
    for _ in 1:ns, j in 1:ns
        p = S.pnt[j]
        p != 0 && (depth[j] = depth[p] + 1)
    end
    mk = 2 * gather_makespan(j -> collect(Graphs.neighbors(S.chd, j)), roots)
    return (slots = mk, depth = maximum(depth; init = 0) + 1, supernodes = ns,
        fill = length(L.Lval))
end

# centralized, ideal base station: every agent one hop from the coordinator,
# whose single radio serializes the gather and the scatter
central_bs_slots(nagents) = 2 * nagents

# centralized, in-network coordinator at the graph center: n distinct packets
# converge-cast over BFS routing paths, store-and-forward, half-duplex
# (greedy schedule, deepest-first). The personalized scatter back mirrors it.
function central_innet_slots(g)
    n = nv(g)
    ecc = eccentricity(g)
    c = argmin(ecc)
    parent = bfs_parents(g, c)
    depthv = gdistances(g, c)
    queue = [Int[] for _ in 1:n]           # packets held at each vertex
    for v in 1:n
        v == c && continue
        push!(queue[v], v)
    end
    remaining = n - 1
    slots = 0
    order = sortperm(depthv; rev = true)   # deepest senders get priority
    while remaining > 0
        slots += 1
        busy = falses(n)
        for v in order
            v == c && continue
            (isempty(queue[v]) || busy[v] || busy[parent[v]]) && continue
            pkt = popfirst!(queue[v])
            busy[v] = true
            busy[parent[v]] = true
            parent[v] == c ? (remaining -= 1) : push!(queue[parent[v]], pkt)
        end
    end
    return 2 * slots
end

# diffusion / iterative methods: iterations to reach EPS, times the per-round
# neighbourhood-exchange cost Δ (edge-colored, generous: no +1)
function iterative_slots(H, g)
    ev = eigvals(Symmetric(Matrix(H)))
    κ = ev[end] / ev[1]
    ρ = (κ - 1) / (κ + 1)
    θ = (sqrt(κ) - 1) / (sqrt(κ) + 1)
    Kr = ceil(Int, log(EPS) / log(ρ))
    Kc = ceil(Int, log(2 / EPS) / log(1 / θ))
    dmax = Δ(g)
    # CG: Chebyshev-bound iterations; per iteration one neighbourhood exchange
    # plus two all-reduces = combining reduce + broadcast on the BFS tree of
    # the center (2 × its sibling-serialized critical path each)
    ecc = eccentricity(g)
    c = argmin(ecc)
    parent = bfs_parents(g, c)
    kids = [Int[] for _ in 1:nv(g)]
    for v in 1:nv(g)
        v != c && push!(kids[parent[v]], v)
    end
    ar = 2 * gather_makespan(v -> kids[v], [c])
    cg = Kc * (dmax + 2ar)
    return (rich = Kr * dmax, cheb = Kc * dmax, cg = cg, κ = κ, Kr = Kr, Kc = Kc)
end

# ---------------------------------------------------------------------------
# main sweep
# ---------------------------------------------------------------------------

chol(H; kw...) = cholesky!(ChordalCholesky(H; kw...), NoPivot())

results = Dict{Symbol, Vector{NamedTuple}}()

for (name, build, pins, sizes) in TOPOLOGIES
    rows = NamedTuple[]
    for sz in sizes
        g0 = build(sz)
        pinned = pins(g0)
        H, freev = dirichlet_block(g0, pinned)
        gfree, _ = induced_subgraph(g0, freev)
        @assert is_connected(gfree) "$name($sz): free radio graph is disconnected — bad pin choice"
        nag = length(freev)

        F = chol(H)
        ts = tree_stats(F.L)

        # nested-dissection ordering where we have a construction for it
        tsnd = if name == :chain || name == :ring
            vord = nd_order_path(nv(g0))
            Fnd = chol(H; alg = dof_order(vord, freev, 2))
            tree_stats(Fnd.L)
        elseif name == :grid
            vord = nd_order_grid(isqrt(nv(g0)))
            Fnd = chol(H; alg = dof_order(vord, freev, 2))
            tree_stats(Fnd.L)
        else
            ts
        end

        it = iterative_slots(H, gfree)
        push!(rows, (n = nag, dofs = size(H, 1),
            tree = ts.slots, treend = tsnd.slots,
            depth = ts.depth, depthnd = tsnd.depth,
            fill = ts.fill, fillnd = tsnd.fill,
            cbs = central_bs_slots(nag), cin = central_innet_slots(gfree),
            rich = it.rich, cheb = it.cheb, cg = it.cg, κ = it.κ,
            Kr = it.Kr, Kc = it.Kc))
    end
    results[name] = rows
end

# ---------------------------------------------------------------------------
# tables: emitted as `@raw html` <table> blocks (plain, so they inherit
# Documenter's ordinary table styling) with a green `td.win` on the winning
# cell(s). Styled by docs/src/assets/benchtables.css. Copy the printed blocks
# into the .md pages.
# ---------------------------------------------------------------------------

const TITLE = Dict(:grid => "grid", :chain => "chain", :ring => "ring",
    :star => "star", :rgg => "random geometric")

# a data cell, classed `win` if it ties the row's best
cell(v, best) = v == best ? "<td class=\"win\">$v</td>" : "<td>$v</td>"

function html_topology_table(name, rows)
    println("\n```@raw html")
    println("<table class=\"bench\">")
    println("<thead><tr><th>agents</th><th>tree</th><th>tree (ND)</th>",
        "<th>central (base)</th><th>central (in-net)</th><th>Richardson</th>",
        "<th>Chebyshev</th><th>CG (+all-reduce)</th><th>&kappa;(H)</th></tr></thead>")
    println("<tbody>")
    for r in rows
        best = min(r.tree, r.treend, r.cbs, r.cin, r.rich, r.cheb, r.cg)
        println("<tr><td>$(r.n)</td>",
            cell(r.tree, best), cell(r.treend, best), cell(r.cbs, best),
            cell(r.cin, best), cell(r.rich, best), cell(r.cheb, best),
            cell(r.cg, best), "<td>$(round(r.κ, digits = 1))</td></tr>")
    end
    println("</tbody></table>")
    println("```")
end

println("\n===================== MAKESPAN TABLES (raw html) =====================")
for (name, _, _, _) in TOPOLOGIES
    html_topology_table(name, results[name])
end

# winner of a default/ND metric pair (smaller is better)
pcell(a, b) = a <= b ? ("<td class=\"win\">$a</td>", "<td>$b</td>") :
                       ("<td>$a</td>", "<td class=\"win\">$b</td>")

println("\n===================== ORDERING STUDY (raw html) =====================\n")
println("```@raw html")
println("<table class=\"bench\">")
println("<thead><tr><th>topology</th><th>agents</th><th>depth</th><th>depth (ND)</th>",
    "<th>slots</th><th>slots (ND)</th><th>fill</th><th>fill (ND)</th></tr></thead>")
println("<tbody>")
for name in (:chain, :ring, :grid)
    for r in results[name]
        da, db = pcell(r.depth, r.depthnd)
        sa, sb = pcell(r.tree, r.treend)
        fa, fb = pcell(r.fill, r.fillnd)
        println("<tr><td>$(TITLE[name])</td><td>$(r.n)</td>", da, db, sa, sb, fa, fb, "</tr>")
    end
end
println("</tbody></table>")
println("```")

# ---------------------------------------------------------------------------
# figure: comparison curves, four panels
# ---------------------------------------------------------------------------

function comparison_panel(name, rows; xlabel = "agents")
    xs = [r.n for r in rows]
    p = plot(xscale = :log10, yscale = :log10, xlabel = xlabel,
        ylabel = "slots (makespan)", legend = name == :grid ? :topleft : false,
        title = string(name))
    plot!(p, xs, [r.rich for r in rows], label = "diffusion (Richardson)", color = PAL[:rich], marker = :diamond, lw = 2)
    plot!(p, xs, [r.cheb for r in rows], label = "diffusion (Chebyshev)", color = PAL[:cheb], marker = :dtriangle, lw = 2)
    plot!(p, xs, [r.cg for r in rows], label = "CG + all-reduces", color = PAL[:cg], marker = :utriangle, lw = 2)
    plot!(p, xs, [r.cin for r in rows], label = "central (in-network)", color = PAL[:cin], marker = :rect, lw = 2)
    plot!(p, xs, [r.cbs for r in rows], label = "central (base station)", color = PAL[:cbs], marker = :circle, lw = 2)
    plot!(p, xs, [min(r.tree, r.treend) for r in rows], label = "distributed tree (best order)",
        color = PAL[:tree], marker = :star5, lw = 2.6)
    return p
end

plt = plot(
    comparison_panel(:grid, results[:grid]),
    comparison_panel(:rgg, results[:rgg]),
    comparison_panel(:ring, results[:ring]),
    comparison_panel(:chain, results[:chain]),
    layout = (2, 2), size = (1050, 760),
    left_margin = 5Plots.mm, bottom_margin = 4Plots.mm)
savefig(plt, joinpath(FIG, "comparison.svg"))
@info "wrote comparison.svg"

# ---------------------------------------------------------------------------
# figure: the ordering study (chain + grid)
# ---------------------------------------------------------------------------

let rows = results[:chain], rowsg = results[:grid]
    xs = [r.n for r in rows]
    p1 = plot(xscale = :log10, yscale = :log10, xlabel = "agents (chain)",
        ylabel = "elimination-tree depth", legend = :topleft, title = "Ordering sets the depth")
    plot!(p1, xs, [r.depth for r in rows], label = "default order", color = PAL[:cbs], marker = :circle, lw = 2)
    plot!(p1, xs, [r.depthnd for r in rows], label = "nested dissection", color = PAL[:tree], marker = :star5, lw = 2.4)
    plot!(p1, xs, xs, label = L"n", color = :gray, ls = :dash, lw = 1)
    plot!(p1, xs, log2.(xs) .+ 1, label = L"\log_2 n", color = :gray, ls = :dot, lw = 1)

    p2 = plot(xscale = :log10, yscale = :log10, xlabel = "agents (chain)",
        ylabel = "slots (makespan)", legend = :topleft, title = "…and the depth sets the makespan")
    plot!(p2, xs, [r.tree for r in rows], label = "tree, default order", color = PAL[:cbs], marker = :circle, lw = 2)
    plot!(p2, xs, [r.treend for r in rows], label = "tree, nested dissection", color = PAL[:tree], marker = :star5, lw = 2.4)
    plot!(p2, xs, [r.cbs for r in rows], label = "centralized (for scale)", color = :gray, ls = :dash, lw = 1.5)

    plt2 = plot(p1, p2, layout = (1, 2), size = (1000, 400),
        left_margin = 5Plots.mm, bottom_margin = 4Plots.mm)
    savefig(plt2, joinpath(FIG, "ordering.svg"))
    @info "wrote ordering.svg"
end

# ---------------------------------------------------------------------------
# physical overlay: hop distribution of tree messages over the radio graph
# ---------------------------------------------------------------------------

function hop_distribution(name, build, pins, sz)
    g0 = build(sz)
    H, freev = dirichlet_block(g0, pins(g0))
    gfree, _ = induced_subgraph(g0, freev)
    F = chol(H)
    S = F.L.S
    ns = Graphs.nv(S.res)
    pperm = invperm(F.P.perm)
    agent_of(prow) = fld1(pperm[prow], 2)
    host = [agent_of(first(Graphs.neighbors(S.res, j))) for j in 1:ns]
    hops = Int[]
    for j in 1:ns
        p = S.pnt[j]
        p == 0 && continue
        push!(hops, gdistances(gfree, host[j])[host[p]])
    end
    return hops
end

println("\n===================== PHYSICAL ROUTING OVERLAY (markdown) =====================\n")
println("| formation | agents | tree messages | 1 hop | 2 hops | ≥3 hops | mean hops |")
println("|---|---:|---:|---:|---:|---:|---:|")
hopdata = Dict{Symbol, Vector{Int}}()
for (name, build, pins, sz) in ((:grid, TOPOLOGIES[1][2], TOPOLOGIES[1][3], 20),
                                (:rgg, TOPOLOGIES[5][2], TOPOLOGIES[5][3], 256))
    hops = hop_distribution(name, build, pins, sz)
    hopdata[name] = hops
    m = length(hops)
    n1 = count(==(1), hops); n2 = count(==(2), hops); n3 = count(>=(3), hops)
    lab = name == :grid ? "20×20 grid" : "random geometric"
    @printf("| %s | %d | %d | %d (%.0f%%) | %d (%.0f%%) | %d (%.0f%%) | %.2f |\n",
        lab, name == :grid ? 20 * 20 - 5 : 254, m,
        n1, 100n1 / m, n2, 100n2 / m, n3, 100n3 / m, sum(hops) / m)
end

let
    ps = map([:grid, :rgg]) do name
        hops = hopdata[name]
        hmax = maximum(hops)
        counts = [count(==(k), hops) for k in 1:hmax]
        bar(1:hmax, 100 .* counts ./ length(hops); label = "",
            color = name == :grid ? PAL[:tree] : PAL[:cbs],
            xlabel = "physical hops per tree message", ylabel = "share of messages (%)",
            title = name == :grid ? "20×20 grid" : "random geometric, n = 256",
            xticks = 1:hmax)
    end
    plt3 = plot(ps..., layout = (1, 2), size = (1000, 380),
        left_margin = 5Plots.mm, bottom_margin = 4Plots.mm)
    savefig(plt3, joinpath(FIG, "hops.svg"))
    @info "wrote hops.svg"
end

# ---------------------------------------------------------------------------
# wall clock: a real multi-process distributed solve (overhead check)
# ---------------------------------------------------------------------------

println("\n===================== WALL CLOCK (multi-process) =====================\n")
let side = 20
    g0 = Graphs.grid([side, side])
    H, freev = dirichlet_block(g0, grid_pins(side))
    F = chol(H)
    L = F.L
    ws = TreeWorkspace(L, Float64)
    b = rand(size(H, 1))

    t1 = minimum(@elapsed((c = copy(b); tree_forward_ldiv!(L, c, ws); tree_backward_ldiv!(L, c, ws))) for _ in 1:50)

    pids = addprocs(4; exeflags = "--project=$(Base.active_project())")
    try
        @everywhere pids @eval using CellularSheaves
        distributed_tree_solve(L, b, 4; pids)   # warm up remote compilation
        t4 = minimum(@elapsed(distributed_tree_solve(L, b, 4; pids)) for _ in 1:10)
        y1 = (c = copy(b); tree_forward_ldiv!(L, c, ws); tree_backward_ldiv!(L, c, ws); c)
        y4 = distributed_tree_solve(L, b, 4; pids)
        println("| solve | wall time | vs single process | matches |")
        println("|---|---:|---:|---|")
        @printf("| single process, workspace | %.3f ms | 1× | — |\n", 1000t1)
        @printf("| 4 worker processes, RemoteChannel | %.3f ms | %.0f× slower | %.1e |\n",
            1000t4, t4 / t1, norm(y4 - y1) / norm(y1))
        println("\nLocal processes share one machine: this measures protocol overhead and")
        println("verifies exactness, not radio latency — the slot tables above are the")
        println("communication model; this is the software correctness check.")
    finally
        rmprocs(pids)
    end
end
