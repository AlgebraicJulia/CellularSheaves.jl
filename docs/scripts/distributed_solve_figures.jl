# Structural figures for the Distributed Sheaf Solve feature guide.
#
# These are not plots: they are bespoke SVG drawings of the algorithm's
# objects — the formation, its elimination tree, one supernode's frontal
# data, the worker partition, and the message schedule — generated from the
# same computed data as docs/scripts/distributed_solve_figdata.jl (the
# 12-agent two-ring formation of the worked example). Regenerate with:
#
#   julia --project=docs docs/scripts/distributed_solve_figures.jl
#
# Output: docs/figures/distributed_solve/fig_*.svg

using CellularSheaves
using CellularSheaves.NetworkSheaves.DistributedSolve
using CliqueTrees.Multifrontal
using LinearAlgebra
using SparseArrays
using Printf
import CellularSheaves.NetworkSheaves.EuclideanSheaves: _harmonic_extension_restricted_laplacian
import Graphs: neighbors, nv

const FIG = joinpath(@__DIR__, "..", "figures", "distributed_solve")
mkpath(FIG)

# ======================= computed data (see figdata) =========================

const NA, NT = 12, 2
const TV1, TV2 = NA + 1, NA + 2
const DIM = 3
const I3 = Matrix{Float64}(I, DIM, DIM)

consensus_edges = [(1,2),(2,3),(3,4),(4,5),(5,6),(6,1),
                   (7,8),(8,9),(9,10),(10,11),(11,12),(12,7),
                   (1,7)]

sheaf = EuclideanSheaf{Float64}(fill(DIM, NA + NT))
for (i, j) in consensus_edges
    add_sheaf_edge!(sheaf, i, j, I3, I3)
end
for i in 1:6
    add_sheaf_edge!(sheaf, i, TV1, I3, I3)
end
for i in 7:12
    add_sheaf_edge!(sheaf, i, TV2, I3, I3)
end

boundary0 = Dict(TV1 => zeros(DIM), TV2 => zeros(DIM))
_, _, H, B = _harmonic_extension_restricted_laplacian(sheaf, boundary0)
F = cholesky!(ChordalCholesky(sparse(H)), NoPivot())
L = F.L
S = L.S
NS = nv(S.res)

perm = invperm(F.P.perm)
agent_of(prow) = fld1(perm[prow], DIM)

res_agents = [sort(unique(agent_of.(collect(neighbors(S.res, j))))) for j in 1:NS]
sep_agents = [sort(unique(agent_of.(collect(neighbors(S.sep, j))))) for j in 1:NS]
nn = [length(collect(neighbors(S.res, j))) for j in 1:NS]
na = [length(collect(neighbors(S.sep, j))) for j in 1:NS]
pnt = [S.pnt[j] for j in 1:NS]
chd = [collect(neighbors(S.chd, j)) for j in 1:NS]

part = partition_tree(L, 3)

# forward schedule: greedy half-duplex, sibling-serialized (matches figdata)
postord = Int[]
function _po(j)
    for c in chd[j]; _po(c); end
    push!(postord, j)
end
for j in 1:NS
    iszero(pnt[j]) && _po(j)
end
ready = Dict{Int, Int}()
sched = Tuple{Int, Int, Int}[]
for j in postord
    t = 0
    for c in sort(chd[j]; by = c -> ready[c])
        slot = max(ready[c], t) + 1
        push!(sched, (slot, c, j))
        t = slot
    end
    ready[j] = t
end
const FWD_MAKESPAN = isempty(sched) ? 0 : maximum(first, sched)

# ============================ svg helpers ====================================

const FONT = "font-family=\"Helvetica, Arial, sans-serif\""

esc(s) = replace(string(s), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;")

mutable struct SVG
    io::IOBuffer
    w::Float64
    h::Float64
end
SVG(w, h) = SVG(IOBuffer(), w, h)

raw!(s::SVG, str) = (println(s.io, str); s)

function save(s::SVG, name)
    path = joinpath(FIG, name)
    open(path, "w") do f
        println(f, "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 $(s.w) $(s.h)\" ",
            FONT, " font-size=\"13\">")
        # white card so the figure reads identically on light and dark doc themes
        println(f, "<rect x=\"0\" y=\"0\" width=\"$(s.w)\" height=\"$(s.h)\" rx=\"10\" fill=\"#ffffff\"/>")
        write(f, take!(copy(s.io)))
        println(f, "</svg>")
    end
    @info "wrote $path"
end

line!(s, x1, y1, x2, y2; stroke="#333", w=1.5, dash="", opacity=1.0, cap="round") =
    raw!(s, "<line x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\" stroke=\"$stroke\" stroke-width=\"$w\"" *
        (isempty(dash) ? "" : " stroke-dasharray=\"$dash\"") *
        " stroke-opacity=\"$opacity\" stroke-linecap=\"$cap\"/>")

circle!(s, cx, cy, r; fill="#fff", stroke="#333", w=1.5) =
    raw!(s, "<circle cx=\"$cx\" cy=\"$cy\" r=\"$r\" fill=\"$fill\" stroke=\"$stroke\" stroke-width=\"$w\"/>")

rect!(s, x, y, w_, h_; fill="#fff", stroke="none", sw=1.0, rx=0, opacity=1.0) =
    raw!(s, "<rect x=\"$x\" y=\"$y\" width=\"$w_\" height=\"$h_\" rx=\"$rx\" fill=\"$fill\" " *
        "fill-opacity=\"$opacity\" stroke=\"$stroke\" stroke-width=\"$sw\"/>")

function text!(s, x, y, str; size=13, fill="#222", anchor="middle", weight="normal", style="normal")
    raw!(s, "<text x=\"$x\" y=\"$y\" font-size=\"$size\" fill=\"$fill\" text-anchor=\"$anchor\" " *
        "font-weight=\"$weight\" font-style=\"$style\">$(esc(str))</text>")
end

# text with trailing digit runs rendered as proper subscripts via tspan
# ("D_11 y_1" etc.) — reliable in every renderer, unlike Unicode subscripts
function subtext!(s, x, y, parts::Vector{<:Tuple}; size=13, fill="#222", anchor="middle", weight="normal")
    body = join(map(parts) do (t, sub)
        esc(t) * (isempty(sub) ? "" :
            "<tspan font-size=\"$(0.7size)\" baseline-shift=\"-22%\">$(esc(sub))</tspan>")
    end)
    raw!(s, "<text x=\"$x\" y=\"$y\" font-size=\"$size\" fill=\"$fill\" text-anchor=\"$anchor\" " *
        "font-weight=\"$weight\">$body</text>")
end

# curved quadratic path (for tree edges / message arrows)
function qpath!(s, x1, y1, qx, qy, x2, y2; stroke="#333", w=1.5, dash="", opacity=1.0)
    raw!(s, "<path d=\"M $x1 $y1 Q $qx $qy $x2 $y2\" fill=\"none\" stroke=\"$stroke\" " *
        "stroke-width=\"$w\" stroke-opacity=\"$opacity\"" *
        (isempty(dash) ? "" : " stroke-dasharray=\"$dash\"") * "/>")
end

# explicit arrowhead polygon at (x,y) pointing along angle th (radians) —
# renders identically in every SVG renderer, unlike <marker>
function arrowhead!(s, x, y, th; size=8, fill="#666")
    ca, sa = cos(th), sin(th)
    p1 = (x, y)
    p2 = (x - size * ca + 0.45size * sa, y - size * sa - 0.45size * ca)
    p3 = (x - size * ca - 0.45size * sa, y - size * sa + 0.45size * ca)
    raw!(s, "<polygon points=\"$(p1[1]),$(p1[2]) $(p2[1]),$(p2[2]) $(p3[1]),$(p3[2])\" fill=\"$fill\"/>")
end

function arrow!(s, x1, y1, x2, y2; stroke="#666", w=2.0)
    th = atan(y2 - y1, x2 - x1)
    line!(s, x1, y1, x2 - 6cos(th), y2 - 6sin(th); stroke=stroke, w=w)
    arrowhead!(s, x2, y2, th; fill=stroke)
end

# ============================ shared layouts =================================

# formation layout: two hexagons + centered targets
function formation_xy(cx1, cx2, cy, R)
    pos = Dict{Int, Tuple{Float64, Float64}}()
    for k in 1:6                       # ring 1: agent 1 faces the bridge (+x)
        th = (k - 1) * pi / 3
        pos[k] = (cx1 + R * cos(th), cy - R * sin(th))
    end
    for k in 1:6                       # ring 2: agent 7 faces the bridge (−x)
        th = pi + (k - 1) * pi / 3
        pos[6 + k] = (cx2 + R * cos(th), cy - R * sin(th))
    end
    pos[TV1] = (cx1, cy)
    pos[TV2] = (cx2, cy)
    return pos
end

# house palette + worker colors
const C_AGENT   = "#4682b4"    # steelblue
const C_TARGET  = "#333333"
const C_TRACK   = "#999999"
const C_WORKER  = ["#4682b4", "#e08214", "#2f9e44", "#c2255c"]
const C_SEP     = "#c2255c"    # separator highlights
const C_RES     = "#2f9e44"

function draw_formation!(s, pos; r=13, labelsize=11, agent_fill=j -> "#eaf1f7",
        agent_stroke=j -> C_AGENT, edgecolor=(i, j) -> "#555", edgew=(i, j) -> 1.6)
    for (i, j) in consensus_edges
        (x1, y1), (x2, y2) = pos[i], pos[j]
        line!(s, x1, y1, x2, y2; stroke=edgecolor(i, j), w=edgew(i, j))
    end
    for i in 1:6
        (x1, y1), (x2, y2) = pos[i], pos[TV1]
        line!(s, x1, y1, x2, y2; stroke=C_TRACK, w=1.0, dash="3 3", opacity=0.7)
    end
    for i in 7:12
        (x1, y1), (x2, y2) = pos[i], pos[TV2]
        line!(s, x1, y1, x2, y2; stroke=C_TRACK, w=1.0, dash="3 3", opacity=0.7)
    end
    for a in 1:NA
        (x, y) = pos[a]
        circle!(s, x, y, r; fill=agent_fill(a), stroke=agent_stroke(a), w=2.0)
        text!(s, x, y + 4, a; size=labelsize, weight="bold", fill="#1a1a1a")
    end
    for tv in (TV1, TV2)
        (x, y) = pos[tv]
        d = 12
        raw!(s, "<polygon points=\"$(x),$(y-d) $(x+d),$(y) $(x),$(y+d) $(x-d),$(y)\" " *
            "fill=\"#2b2b2b\" stroke=\"#000\" stroke-width=\"1\"/>")
        text!(s, x, y + 3.5, tv == TV1 ? "T1" : "T2"; size=9, fill="#ffffff", weight="bold")
    end
end

# tree layout: x by in-order leaf walk, y by depth (root on top)
function tree_xy(x0, y0, w, dy)
    depth = zeros(Int, NS)
    for _ in 1:NS, j in 1:NS
        pnt[j] != 0 && (depth[j] = depth[pnt[j]] + 1)
    end
    leafx = Ref(0.0)
    nleaves = count(j -> isempty(chd[j]), 1:NS)
    xs = zeros(NS)
    function place(j)
        if isempty(chd[j])
            xs[j] = leafx[]
            leafx[] += 1
        else
            for c in chd[j]; place(c); end
            xs[j] = sum(xs[c] for c in chd[j]) / length(chd[j])
        end
    end
    for j in 1:NS
        iszero(pnt[j]) && place(j)
    end
    return Dict(j => (x0 + w * xs[j] / max(nleaves - 1, 1), y0 + dy * depth[j]) for j in 1:NS),
        maximum(depth)
end

function draw_tree!(s, tpos; r=17, node_fill=j -> "#ffffff", node_stroke=j -> "#333",
        strokew=j -> 2.0, label=j -> string(j), sublabel=j -> "", subanchor=j -> "start",
        skip_edges=Set{Int}())
    for j in 1:NS
        p = pnt[j]
        (p == 0 || j in skip_edges) && continue
        (x1, y1), (x2, y2) = tpos[j], tpos[p]
        line!(s, x1, y1, x2, y2; stroke="#777", w=1.6)
    end
    for j in 1:NS
        (x, y) = tpos[j]
        circle!(s, x, y, r; fill=node_fill(j), stroke=node_stroke(j), w=strokew(j))
        text!(s, x, y + 4, label(j); size=12, weight="bold")
        sl = sublabel(j)
        if !isempty(sl)
            anch = subanchor(j)
            dx = anch == "start" ? r + 5 : -(r + 5)
            text!(s, x + dx, y + 4, sl; size=9.5, fill="#555", anchor=anch)
        end
    end
end

# ====================== figure 1: formation → tree ==========================

function fig_formation_tree()
    s = SVG(1000, 470)

    text!(s, 250, 30, "Formation graph  G  (the sheaf lives here)"; size=15, weight="bold")
    text!(s, 250, 48, "12 agents, 2 escort rings, bridge 1–7; targets pinned"; size=11.5, fill="#666")

    # two example supernodes, accent-matched across both panels: the root
    # supernode 9 = agents {1,2,4} and supernode 4 = agents {9,12}
    hl = Dict(1 => C_SEP, 2 => C_SEP, 4 => C_SEP, 9 => "#e08214", 12 => "#e08214")

    pos = formation_xy(140, 370, 260, 92)
    draw_formation!(s, pos;
        agent_fill=a -> haskey(hl, a) ? (hl[a] == C_SEP ? "#fbe4ec" : "#fdeedd") : "#eaf1f7",
        agent_stroke=a -> get(hl, a, C_AGENT))

    arrow!(s, 512, 260, 585, 260; stroke="#666", w=2.2)
    text!(s, 548, 243, "factor H"; size=11.5, fill="#444", style="italic")

    text!(s, 790, 30, "Elimination tree  (the solve lives here)"; size=15, weight="bold")
    text!(s, 790, 48, "9 supernodes; {…} = agents eliminated there"; size=11.5, fill="#666")

    node_hl = Dict(9 => C_SEP, 4 => "#e08214")
    tpos, _ = tree_xy(650, 85, 280, 61)
    draw_tree!(s, tpos;
        node_fill=j -> haskey(node_hl, j) ? (node_hl[j] == C_SEP ? "#fbe4ec" : "#fdeedd") : "#eaf1f7",
        node_stroke=j -> get(node_hl, j, C_AGENT),
        strokew=j -> haskey(node_hl, j) ? 2.6 : 2.0,
        sublabel=j -> "{" * join(res_agents[j], ",") * "}",
        subanchor=j -> tpos[j][1] > 840 ? "end" : "start")

    text!(s, 790, 452, "a supernode is a group of agents eliminated together — same color, both panels";
        size=10.5, fill="#666")
    save(s, "fig_formation_tree.svg")
end

fig_formation_tree()

# ==================== figure 2: anatomy of one supernode =====================
#
# Supernode 4: res = agents {9,12} (nn = 6 dofs), sep = agent {7} (na = 3),
# children 1,2,3, parent 5. Frontal blocks drawn to scale.

function fig_supernode_anatomy()
    s = SVG(1000, 470)
    j = 4

    text!(s, 240, 30, "The bag of supernode 4 on the formation"; size=15, weight="bold")
    text!(s, 240, 48, "res = eliminated here;  sep = shared boundary with ancestors"; size=11.5, fill="#666")

    inres = a -> a in res_agents[j]
    insep = a -> a in sep_agents[j]
    pos = formation_xy(130, 350, 265, 88)
    draw_formation!(s, pos;
        agent_fill=a -> inres(a) ? "#e4f3e8" : insep(a) ? "#fbe4ec" : "#f4f4f4",
        agent_stroke=a -> inres(a) ? C_RES : insep(a) ? C_SEP : "#bbbbbb",
        edgecolor=(u, v) -> "#cccccc", edgew=(u, v) -> 1.2)

    # legend
    circle!(s, 60, 415, 8; fill="#e4f3e8", stroke=C_RES, w=2)
    text!(s, 76, 419, "res(4) = agents {9, 12} — 6 dofs"; anchor="start", size=11)
    circle!(s, 60, 440, 8; fill="#fbe4ec", stroke=C_SEP, w=2)
    text!(s, 76, 444, "sep(4) = agent {7} — 3 dofs"; anchor="start", size=11)

    # ---- frontal blocks to scale (cell = 17 px) ----
    cx, cy, c = 570, 95, 17
    text!(s, cx + 4.5c, 66, "Stored factor blocks (to scale)"; size=14, weight="bold")
    rect!(s, cx, cy, 6c, 6c; fill="#e4f3e8", stroke=C_RES, sw=2)
    rect!(s, cx, cy + 6c, 6c, 3c; fill="#fbe4ec", stroke=C_SEP, sw=2)
    for k in 1:8   # grid lines
        line!(s, cx + k * c, cy, cx + k * c, cy + 9c; stroke="#ffffff", w=1)
        line!(s, cx, cy + k * c, cx + 6c, cy + k * c; stroke="#ffffff", w=1)
    end
    # not-stored upper triangle marker
    rect!(s, cx + 6c, cy, 3c, 6c; fill="#f4f4f4", stroke="#cccccc", sw=1)
    text!(s, cx + 7.5c, cy + 3.2c, "symmetric"; size=9, fill="#999")
    rect!(s, cx + 6c, cy + 6c, 3c, 3c; fill="#fff8e6", stroke="#e08214", sw=1.6)
    text!(s, cx + 7.5c, cy + 7.4c, "update"; size=9, fill="#b06000")
    text!(s, cx + 7.5c, cy + 8.2c, "→ parent"; size=9, fill="#b06000")

    subtext!(s, cx + 3c, cy + 3.2c, [("D", "11"), ("  (6×6)", "")]; size=12, weight="bold", fill="#1c6b32")
    subtext!(s, cx + 3c, cy + 7.6c, [("L", "21"), (" (3×6)", "")]; size=11.5, weight="bold", fill="#8f1d4e")
    text!(s, cx - 10, cy + 3c, "9"; size=10, fill="#1c6b32", anchor="end")
    text!(s, cx - 10, cy + 4.6c, "12"; size=10, fill="#1c6b32", anchor="end")
    text!(s, cx - 10, cy + 7.6c, "7"; size=10, fill="#8f1d4e", anchor="end")

    # ---- dataflow ----
    fx = 830
    text!(s, fx + 40, 66, "Per-solve dataflow"; size=14, weight="bold")
    for (k, cnode) in enumerate(chd[j])
        bx = fx - 40 + (k - 1) * 58
        rect!(s, bx, 385, 44, 30; fill="#eaf1f7", stroke=C_AGENT, sw=1.6, rx=5)
        text!(s, bx + 22, 404, "child $cnode"; size=9.5)
        arrow!(s, bx + 22, 385, fx + 38 + (k - 2) * 8, 320; stroke="#888", w=1.4)
    end
    rect!(s, fx - 20, 275, 120, 44; fill="#e4f3e8", stroke=C_RES, sw=2, rx=6)
    text!(s, fx + 40, 293, "solve"; size=10.5)
    subtext!(s, fx + 40, 308, [("D", "11"), (" y", "1"), (" = c", "1")]; size=11.5, weight="bold")
    arrow!(s, fx + 40, 275, fx + 40, 195; stroke=C_SEP, w=2.2)
    subtext!(s, fx + 52, 240, [("m", "2"), (" = f", "2"), (" − L", "21"), (" y", "1")]; size=11, fill=C_SEP, anchor="start")
    text!(s, fx + 52, 256, "3 numbers, to parent"; size=9.5, fill="#8f1d4e", anchor="start")
    rect!(s, fx - 6, 150, 92, 44; fill="#ffffff", stroke="#777", sw=1.6, rx=6)
    text!(s, fx + 40, 168, "parent 5"; size=10.5)
    text!(s, fx + 40, 184, "sep = agent 7"; size=9.5, fill="#666")

    save(s, "fig_supernode_anatomy.svg")
end

fig_supernode_anatomy()

# ==================== figure 3: partition and cut edges ======================

function fig_partition()
    s = SVG(1000, 470)
    text!(s, 350, 30, "The clique tree, cut into per-worker chunks"; size=15, weight="bold")
    text!(s, 350, 48, "connected rooted subtrees; only cut edges cross the network"; size=11.5, fill="#666")

    wfill = ["#eaf1f7", "#fdeedd", "#e4f3e8"]

    tpos, _ = tree_xy(180, 90, 340, 62)

    # cut edges drawn as bold dashed crimson beneath everything
    for (child, parent) in part.cut_edges
        (x1, y1), (x2, y2) = tpos[child], tpos[parent]
        line!(s, x1, y1, x2, y2; stroke=C_SEP, w=4.0, dash="7 5", opacity=0.85)
        mx, my = (x1 + x2) / 2, (y1 + y2) / 2
        text!(s, mx + 14, my + 4, "RemoteChannel"; size=9.5, fill=C_SEP, anchor="start", style="italic")
    end

    draw_tree!(s, tpos;
        node_fill=j -> wfill[part.owner[j]],
        node_stroke=j -> C_WORKER[part.owner[j]],
        strokew=j -> 2.4,
        sublabel=j -> "{" * join(res_agents[j], ",") * "}",
        subanchor=j -> tpos[j][1] > 400 ? "start" : "end",
        skip_edges=Set(first.(part.cut_edges)))

    # per-worker slice cards
    wds = [worker_factorization(L, part, w) for w in 1:length(part.chunks)]
    for (w, wd) in enumerate(wds)
        bx, by = 660, 90 + (w - 1) * 118
        kb = (length(wd.Dval) + length(wd.Lval)) * 8 / 1e3
        rect!(s, bx, by, 300, 96; fill=wfill[w], stroke=C_WORKER[w], sw=2, rx=8)
        text!(s, bx + 16, by + 26, "worker $w"; size=13, weight="bold", anchor="start")
        text!(s, bx + 16, by + 48, "supernodes " * join(sort(part.chunks[w]), ", "); size=11, anchor="start")
        text!(s, bx + 16, by + 68, @sprintf("factor slice: %.2f KB (copied, not aliased)", kb); size=11, anchor="start")
        bup = wd.boundary_up
        bdn = wd.boundary_down
        io = String[]
        isempty(bup) || push!(io, "sends up @ " * join(bup, ","))
        isempty(bdn) || push!(io, "receives @ " * join(bdn, ","))
        text!(s, bx + 16, by + 86, isempty(io) ? "no cross-worker traffic" : join(io, " · ");
            size=10, fill="#555", anchor="start")
    end

    total_kb = (length(L.Dval) + length(L.Lval)) * 8 / 1e3
    text!(s, 350, 452, @sprintf("slices are disjoint and sum to the full factor: %.2f KB — nothing stored twice", total_kb);
        size=10.5, fill="#666")
    save(s, "fig_partition.svg")
end

fig_partition()

# ================ figure 4a: message schedule, space-time ====================

function fig_schedule()
    s = SVG(1000, 480)
    text!(s, 500, 30, "One solve = one gather up the tree, one scatter back down"; size=15, weight="bold")
    text!(s, 500, 48, "half-duplex, one packet per slot per node; siblings serialize at their parent"; size=11.5, fill="#666")

    # left: the tree, edges labeled with their forward slot
    tpos, _ = tree_xy(90, 95, 240, 56)
    slot_of = Dict((c, p) => sl for (sl, c, p) in sched)
    draw_tree!(s, tpos; r=14,
        node_fill=j -> "#eaf1f7", node_stroke=j -> C_AGENT,
        label=j -> string(j))
    for (sl, c, p) in sched
        (x1, y1), (x2, y2) = tpos[c], tpos[p]
        mx, my = (x1 + x2) / 2 + 11, (y1 + y2) / 2 + 4
        circle!(s, mx, my, 8.5; fill="#ffffff", stroke="#888", w=1.2)
        text!(s, mx, my + 3.5, sl; size=9, fill="#333", weight="bold")
    end
    text!(s, 210, 445, "edge badge = forward slot"; size=10.5, fill="#666")

    # right: space-time grid; rows = tree edges ordered by forward slot
    edges_sorted = sort(sched)
    gx, gy, cw, chh = 430, 105, 36, 34
    nslot = 2 * FWD_MAKESPAN
    for t in 1:nslot
        text!(s, gx + (t - 0.5) * cw, gy - 10, t; size=9.5, fill="#555")
    end
    text!(s, gx + FWD_MAKESPAN * cw / 2, gy - 30, "forward gather"; size=11, fill=C_AGENT, weight="bold")
    text!(s, gx + 1.5 * FWD_MAKESPAN * cw, gy - 30, "backward scatter"; size=11, fill="#b06000", weight="bold")
    for (k, (sl, c, p)) in enumerate(edges_sorted)
        y = gy + (k - 1) * chh
        text!(s, gx - 8, y + chh / 2 + 4, "$c → $p"; size=10.5, anchor="end", fill="#333")
        rect!(s, gx, y + 4, nslot * cw, chh - 8; fill="#f7f7f7")
        rect!(s, gx + (sl - 1) * cw, y + 4, cw, chh - 8; fill=C_AGENT, rx=4)
        bsl = 2 * FWD_MAKESPAN + 1 - sl        # mirrored slot, p → c
        rect!(s, gx + (bsl - 1) * cw, y + 4, cw, chh - 8; fill="#e08214", rx=4)
        text!(s, gx + (sl - 0.5) * cw, y + chh / 2 + 3.5, "↑"; size=11, fill="#fff", weight="bold")
        text!(s, gx + (bsl - 0.5) * cw, y + chh / 2 + 3.5, "↓"; size=11, fill="#fff", weight="bold")
    end
    line!(s, gx + FWD_MAKESPAN * cw, gy - 18, gx + FWD_MAKESPAN * cw, gy + length(edges_sorted) * chh + 6;
        stroke="#999", w=1.2, dash="4 4")
    text!(s, gx + nslot * cw / 2, gy + length(edges_sorted) * chh + 28,
        "makespan = $(2 * FWD_MAKESPAN) slots — set by the tree's critical path, not by the number of agents";
        size=11, fill="#333")
    save(s, "fig_schedule.svg")
end

fig_schedule()

# ================ figure 4b: the same schedule, animated =====================

function fig_schedule_anim()
    s = SVG(640, 460)
    text!(s, 320, 28, "Messages in motion: forward gather, backward scatter"; size=14.5, weight="bold")

    tpos, _ = tree_xy(120, 80, 400, 60)
    draw_tree!(s, tpos; r=15,
        node_fill=j -> "#eaf1f7", node_stroke=j -> C_AGENT)

    nslot = 2 * FWD_MAKESPAN
    slotdur = 0.55
    pause = 1.1
    cycle = nslot * slotdur + pause
    frac(t) = round(t * slotdur / cycle; digits=4)

    for (sl, c, p) in sched
        (x1, y1), (x2, y2) = tpos[c], tpos[p]
        for (t0f, xa, ya, xb, yb, col) in (
                (sl - 1, x1, y1, x2, y2, C_AGENT),                      # forward
                (2 * FWD_MAKESPAN - sl, x2, y2, x1, y1, "#e08214"))     # backward
            ta, tb = frac(t0f), frac(t0f + 1)
            raw!(s, "<circle r=\"6.5\" fill=\"$col\" opacity=\"0\">" *
                "<animateMotion dur=\"$(cycle)s\" repeatCount=\"indefinite\" " *
                "path=\"M $xa $ya L $xb $yb\" calcMode=\"linear\" " *
                "keyPoints=\"0;0;1;1\" keyTimes=\"0;$ta;$tb;1\"/>" *
                "<animate attributeName=\"opacity\" dur=\"$(cycle)s\" repeatCount=\"indefinite\" " *
                "calcMode=\"discrete\" values=\"0;1;0\" keyTimes=\"0;$ta;$tb\"/>" *
                "</circle>")
        end
    end

    # progress bar over the slot timeline
    bw, bx, by = 480, 80, 428
    rect!(s, bx, by, bw, 8; fill="#eeeeee", rx=4)
    active = round(nslot * slotdur / cycle; digits=4)
    raw!(s, "<rect x=\"$bx\" y=\"$by\" width=\"0\" height=\"8\" rx=\"4\" fill=\"#bbbbbb\">" *
        "<animate attributeName=\"width\" dur=\"$(cycle)s\" repeatCount=\"indefinite\" " *
        "calcMode=\"linear\" values=\"0;$bw;$bw\" keyTimes=\"0;$active;1\"/></rect>")
    text!(s, bx + bw / 2, by - 8, "slots 1 – $nslot"; size=9.5, fill="#666")
    save(s, "fig_schedule_anim.svg")
end

fig_schedule_anim()

# =============== figure 5: three methods, same formation =====================

adjset = Set{Tuple{Int, Int}}()
for (i, j) in consensus_edges
    push!(adjset, (i, j)); push!(adjset, (j, i))
end
host = [first(res_agents[j]) for j in 1:NS]

function fig_methods()
    s = SVG(1080, 430)
    panels = [
        ("Distributed tree (this feature)", "14 slots · exact",
         "messages follow the elimination tree"),
        ("Centralized coordinator", "24 slots · exact",
         "one radio serializes the whole fleet"),
        ("Sheaf diffusion", "≈120 slots (Richardson) · ε = 1e-6",
         "local rounds pay for conditioning: K ∝ κ(H)"),
    ]
    for (kp, (title, cost, sub)) in enumerate(panels)
        px = 30 + (kp - 1) * 350
        rect!(s, px, 55, 330, 330; fill="#fbfbfb", stroke="#dddddd", sw=1.2, rx=10)
        text!(s, px + 165, 82, title; size=13, weight="bold")
        text!(s, px + 165, 100, sub; size=10, fill="#666")
        pos = formation_xy(px + 95, px + 235, 250, 56)
        draw_formation!(s, pos; r=9, labelsize=8,
            edgecolor=(u, v) -> "#cccccc", edgew=(u, v) -> 1.0)
        if kp == 1
            for (sl, c, p) in sched
                a, b = host[c], host[p]
                (xa, ya), (xb, yb) = pos[a], pos[b]
                routed = !((a, b) in adjset)
                th = atan(yb - ya, xb - xa)
                x2, y2 = xb - 11cos(th), yb - 11sin(th)
                if routed
                    qpath!(s, xa, ya, (xa + xb) / 2, (ya + yb) / 2 - 22, x2, y2;
                        stroke=C_SEP, w=1.8, dash="5 4")
                    arrowhead!(s, x2, y2, th; fill=C_SEP, size=6)
                else
                    arrow!(s, xa, ya, x2, y2; stroke=C_AGENT, w=2.0)
                end
            end
            text!(s, px + 165, 350, "solid = radio neighbors;  dashed = routed (fill-in)";
                size=9.5, fill="#666")
        elseif kp == 2
            bx, by = px + 165, 130
            rect!(s, bx - 32, by - 14, 64, 26; fill="#2b2b2b", rx=5)
            text!(s, bx, by + 4, "base"; size=10, fill="#ffffff", weight="bold")
            for a in 1:NA
                (xa, ya) = pos[a]
                arrow!(s, xa, ya, bx + (xa - bx) * 0.12, by + 16; stroke="#888", w=1.1)
            end
            text!(s, px + 165, 350, "n up + n down every solve — 2n slots, always";
                size=9.5, fill="#666")
        else
            for (u, v) in consensus_edges
                (xa, ya), (xb, yb) = pos[u], pos[v]
                mx, my = (xa + xb) / 2, (ya + yb) / 2
                th = atan(yb - ya, xb - xa)
                arrowhead!(s, mx + 7cos(th), my + 7sin(th), th; fill=C_RES, size=6)
                arrowhead!(s, mx - 7cos(th), my - 7sin(th), th + pi; fill=C_RES, size=6)
            end
            text!(s, px + 165, 350, "every edge, every round; ×40 rounds here (×17 Chebyshev)";
                size=9.5, fill="#666")
        end
        text!(s, px + 165, 372, cost; size=12, weight="bold", fill="#333")
    end
    text!(s, 540, 412, "same 12-agent formation, same ask (one exact solve, or ε = 1e-6), same slot rules — computed, not sketched";
        size=10.5, fill="#666")
    save(s, "fig_methods.svg")
end

fig_methods()
