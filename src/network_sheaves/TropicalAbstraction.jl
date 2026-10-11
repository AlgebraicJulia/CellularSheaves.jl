"""
    TropicalAbstraction

Coarse grids for deterministic optimal control in idempotent semirings, where
they give exact solutions through bounds rather than residual corrections.

**Deterministic optimal control as a min-plus fixed point.** On a finite state
graph whose edge ``x → y`` costs ``W_{xy} ≥ 0`` (a control taking state ``x`` to
``y``), the cost-to-go to a target set ``Γ`` solves the Bellman equation

```math
V(x) = \\min_y \\big( W_{xy} + V(y) \\big), \\qquad V|_Γ = 0,
```

which is linear in the min-plus semiring ``(ℝ ∪ \\{∞\\}, \\min, +)``:
``V = W ⊗ V ⊕ b`` with ``b = 0`` on ``Γ`` and ``∞`` elsewhere, so ``V = W^* ⊗ b``
([`minplus_closure`](@ref), [`cost_to_go`](@ref)).

**Pushforward and pullback without subtraction.** For a map ``ψ`` sending each
state to a block (an aggregation such as
[`aggregation_homomorphism`](@ref CellularSheaves.NetworkSheaves.GridMultigrid.aggregation_homomorphism)),
the pullback ``E = ψ^*`` copies a block's value onto its states
([`tropical_pullback`](@ref)) and its transpose ``E^T`` takes the semiring sum
over each fibre, the minimum ([`tropical_pushforward`](@ref)). Since ``\\min`` is
idempotent, ``E^T ⊗ E = I`` with no normalization: the block average of the real
case becomes the block minimum. The Galerkin operator ``E^T ⊗ W ⊗ E`` is the
*quotient graph* ([`quotient_graph`](@ref)): one vertex per block, the cheapest
edge between two blocks, edges inside a block contracted.

**Bounds instead of corrections.** Min-plus has no subtraction and so no
residual to correct. The coarse problem is instead a *relaxation*: every path
maps to a walk of the quotient that costs no more, so the pulled-back coarse
cost-to-go is a lower bound on the fine one, and a consistent A* heuristic
([`Abstraction`](@ref), [`lower_bound`](@ref), [`shortest_path`](@ref)). The
search with it is exact and expands fewer states than Dijkstra's algorithm: the
coarse solution prunes the fine search. This is Hierarchical A* (Holte, Perez,
Zimmer & MacDonald 1996) and the abstraction heuristics of classical planning
(pattern databases, Culberson & Schaeffer 1998; merge-and-shrink, Helmert,
Haslum & Hoffmann 2007).

**Coarsening time with space.** Aggregating positions into blocks of side ``s``
while keeping one time step lets the coarse problem cross a block per step, so
its bound degrades by about ``s``. The cure is to coarsen time too: the pushforward
of the ``k``-step graph ``W^{⊗k}`` (`steps = k`), whose edges are the ``k``-step
moves. A fine path from ``x`` is cut into ``k``-step chunks from its start plus a
remainder of fewer than ``k`` steps into the target; each chunk maps to an edge
(or a contracted loop) of the quotient of ``W^{⊗k}``, and the remainder is
bounded by the coarse terminal cost ``b_c = E^T ⊗ \\bigoplus_{r<k} W^{⊗r} ⊗ b``.
The coarse cost-to-go is again a lower bound, now admissible but not always
consistent, so [`shortest_path`](@ref) reopens states. For the double integrator
the natural scaling is ``k`` steps with velocity blocks of ``k`` and position
blocks of ``k^2``.

**Other idempotent semirings.** In the bottleneck semiring ``(ℝ ∪ \\{±∞\\}, \\max,
\\min)`` the reach-avoid value is a closure ([`reach_avoid`](@ref),
[`maxmin_closure`](@ref)), and the pushforward with ``\\max`` over fibres gives an
*upper* bound (an over-approximation of the reach-avoid set; see
[`quotient_graph`](@ref) with `combine = max`).

The coarse closure can be computed once per abstraction for all targets by
PathPack's chordal semiring LU (`method = :pathpack` in [`Abstraction`](@ref),
with PATHPACK loaded). The same pushforward is the Galerkin coarse operator of
the multigrid solvers in the real semiring (see
[`GridMultigrid`](@ref CellularSheaves.NetworkSheaves.GridMultigrid)), where the
coarse solution is a correction rather than a bound.
"""
module TropicalAbstraction

export tropical_pullback, tropical_pushforward, quotient_graph, minplus_closure, cost_to_go,
    Abstraction, lower_bound, shortest_path, schur_complement, Overlay, overlay_path, maxmin_closure, reach_avoid,
    projection_homomorphism, LayeredPlanner, plan_route, corridor, layered_path,
    intertwines, TruncatedPushforward, coarsen, refine, coarsening_matrix, refinement_matrix,
    robust_quotient, upper_bound, is_control_lyapunov, greedy_policy, is_lyapunov, corridor_certificate, certified_path

using ArgCheck: @argcheck
using SparseArrays
using ..GraphHomomorphisms: GraphHomomorphism

"""
    tropical_pullback(ψ, vc) -> Vector

The pullback ``E ⊗ v_c`` along `ψ`: each state takes the value of its block.
"""
function tropical_pullback(ψ::GraphHomomorphism, vc::AbstractVector)
    @argcheck length(vc) == ψ.n_target "one value per block"
    return vc[ψ.vertex_map]
end

"""
    tropical_pushforward(ψ, v; combine = min) -> Vector

The pushforward ``E^T ⊗ v`` along `ψ` in the semiring whose addition is
`combine`: each block takes the `combine` of `v` over its states (the semiring
zero for an empty block). With `min` (min-plus) or `max` (max-min) it is a left
inverse of [`tropical_pullback`](@ref) for a surjective ``ψ``.
"""
function tropical_pushforward(ψ::GraphHomomorphism, v::AbstractVector{T}; combine = min) where {T<:Real}
    @argcheck length(v) == length(ψ.vertex_map) "one value per state"
    @argcheck combine === min || combine === max "combine must be min or max"
    S = float(T)
    vc = fill(combine === min ? typemax(S) : typemin(S), ψ.n_target)
    for (x, b) in enumerate(ψ.vertex_map)
        vc[b] = combine(vc[b], v[x])
    end
    return vc
end

function _check_costs(W)
    @argcheck size(W, 1) == size(W, 2) "W must be square"
    @argcheck all(>=(0), nonzeros(W)) "edge costs must be nonnegative"
end

# C = A ⊗ B in min-plus (Gustavson's sparse product, column by column).
function _minplus_mul(A::SparseMatrixCSC{T}, B::SparseMatrixCSC{T}) where {T}
    m, n = size(A, 1), size(B, 2)
    acc = fill(typemax(T), m)
    touched = Int[]
    colptr = [1]
    rows, vals = Int[], T[]
    Ar, Av, Br, Bv = rowvals(A), nonzeros(A), rowvals(B), nonzeros(B)
    for j in 1:n
        for q in nzrange(B, j)
            k, b = Br[q], Bv[q]
            for p in nzrange(A, k)
                i = Ar[p]
                isinf(acc[i]) && push!(touched, i)
                acc[i] = min(acc[i], Av[p] + b)
            end
        end
        sort!(touched)
        for i in touched
            push!(rows, i); push!(vals, acc[i])
            acc[i] = typemax(T)
        end
        empty!(touched)
        push!(colptr, length(rows) + 1)
    end
    return SparseMatrixCSC(m, n, colptr, rows, vals)
end

"""
    quotient_graph(W, ψ; steps = 1, combine = min) -> SparseMatrixCSC

The Galerkin operator ``E^T ⊗ W^{⊗k} ⊗ E`` of the edge weights `W` (`W[x, y]` the
weight of the edge ``x → y``; absent entries are the semiring zero) along `ψ`,
for `k = steps`: an edge between two blocks carries the best `k`-step move
between their states. Edges inside a block are dropped.

With `combine = min` (min-plus, costs ``≥ 0``) the quotient's cost-to-go is a
lower bound on the fine one; with `combine = max` (max-min, `steps = 1` only) its
bottleneck values are an upper bound.
"""
function quotient_graph(W::SparseMatrixCSC{T}, ψ::GraphHomomorphism; steps::Integer = 1, combine = min) where {T<:Real}
    @argcheck size(W, 1) == size(W, 2) "W must be square"
    @argcheck size(W, 1) == length(ψ.vertex_map) "ψ must map the states of W"
    @argcheck steps >= 1 "steps must be positive"
    @argcheck combine === min || (combine === max && steps == 1) "combine must be min, or max with one step"
    combine === min && _check_costs(W)
    nc = ψ.n_target
    if steps == 1
        rows, cols, vals = findnz(W)
    else
        S = float(T)
        r, c, v = findnz(W)
        M = sparse(r, ψ.vertex_map[c], S.(v), size(W, 1), nc, min)       # W ⊗ E
        Wf = SparseMatrixCSC{S,Int}(W)
        for _ in 2:steps
            M = _minplus_mul(Wf, M)                                         # W ⊗ (W^{k-1} ⊗ E)
        end
        rows, cols, vals = findnz(M)                                        # fine state, coarse block
        keep = [ψ.vertex_map[x] != b for (x, b) in zip(rows, cols)]
        return sparse(ψ.vertex_map[rows[keep]], cols[keep], vals[keep], nc, nc, min)
    end
    keep = [ψ.vertex_map[x] != ψ.vertex_map[y] for (x, y) in zip(rows, cols)]
    return sparse(ψ.vertex_map[rows[keep]], ψ.vertex_map[cols[keep]], vals[keep], nc, nc, combine)
end

# Binary min-heap of (key, state) pairs with lazy deletion.
function _push!(heap, item)
    push!(heap, item)
    i = length(heap)
    while i > 1
        p = i >> 1
        heap[p] <= heap[i] && break
        heap[p], heap[i] = heap[i], heap[p]
        i = p
    end
end

function _pop!(heap)
    top = heap[1]
    last = pop!(heap)
    if !isempty(heap)
        heap[1] = last
        i, n = 1, length(heap)
        while true
            l = 2i
            l > n && break
            c = (l < n && heap[l + 1] < heap[l]) ? l + 1 : l
            heap[i] <= heap[c] && break
            heap[i], heap[c] = heap[c], heap[i]
            i = c
        end
    end
    return top
end

"""
    minplus_closure(W, b) -> Vector

The least solution ``V = W^* ⊗ b`` of ``V(x) = \\min\\big(b(x), \\min_y W_{xy} + V(y)\\big)``
for costs ``W ≥ 0``: the cost-to-go with terminal cost `b` (``∞`` where there is
none), by Dijkstra's algorithm run backwards over the in-edges (columns) of `W`.
"""
function minplus_closure(W::SparseMatrixCSC{T}, b::AbstractVector{<:Real}) where {T<:Real}
    _check_costs(W)
    @argcheck length(b) == size(W, 1) "one terminal cost per state"
    S = float(promote_type(T, eltype(b)))
    g = S.(b)
    heap = Tuple{S,Int}[]
    for (x, v) in enumerate(g)
        isinf(v) || _push!(heap, (v, x))
    end
    settled = falses(length(g))
    rows, vals = rowvals(W), nonzeros(W)
    while !isempty(heap)
        _, y = _pop!(heap)
        settled[y] && continue
        settled[y] = true
        for k in nzrange(W, y)
            x = rows[k]
            d = g[y] + vals[k]
            if d < g[x]
                g[x] = d
                _push!(heap, (d, x))
            end
        end
    end
    return g
end

"""
    cost_to_go(W, targets) -> Vector

The cost-to-go to the states `targets`: [`minplus_closure`](@ref) with terminal
cost 0 on the targets.
"""
function cost_to_go(W::SparseMatrixCSC, targets)
    N = size(W, 1)
    @argcheck !isempty(targets) && all(t -> 1 <= t <= N, targets) "targets must be states of W"
    b = fill(Inf, N)
    b[targets] .= 0
    return minplus_closure(W, b)
end

# ``⊕_{r<k} W^{⊗r} ⊗ b``: the best cost to the targets in fewer than k steps,
# propagated from the frontier of exactly-r-step costs along in-edges (columns).
function _short_horizon(W, b, k)
    v = float.(b)
    frontier = Dict(x => v[x] for x in eachindex(v) if !isinf(v[x]))
    rows, vals = rowvals(W), nonzeros(W)
    for _ in 1:(k - 1)
        next = Dict{Int,Float64}()
        for (y, c) in frontier, p in nzrange(W, y)
            x, d = rows[p], c + vals[p]
            d < get(next, x, Inf) && (next[x] = d)
        end
        for (x, d) in next
            v[x] = min(v[x], d)
        end
        frontier = next
    end
    return v
end

"""
    Abstraction(W, ψ; steps = 1, method = :dijkstra, table = false)

The coarse problem of the edge costs `W` along the aggregation `ψ` with time
coarsened by `steps`: the quotient graph of the `steps`-step graph
([`quotient_graph`](@ref)), built once and reused for every target set by
[`lower_bound`](@ref).

`method = :pathpack` factorizes the coarse min-plus closure ``W_c^* = U^* L^*``
once with PathPack's chordal semiring LU (requires PATHPACK), so that each target
set costs two triangular solves. `table = true` stores the whole closure ``W_c^*``,
the cost between every pair of coarse states (all-pairs shortest paths on the
quotient, by PathPack's closure or by one Dijkstra per column; ``n_c^2`` entries):
each target set's coarse cost-to-go is then the min-plus product of the table with
the coarse terminal cost, a lookup of the columns of the target blocks. This is a
pattern database (Culberson & Schaeffer 1998) and, in control terms, the min-plus
fundamental solution of the coarse Bellman equation: the value for any terminal
cost ``φ`` is the Hopf–Lax formula ``V(B) = \\min_C W_c^*(B, C) + φ(C)``.
"""
struct Abstraction{G,C}
    W::G
    ψ::GraphHomomorphism
    steps::Int
    quotient::SparseMatrixCSC{Float64,Int}
    closure::C
    table::Union{Nothing,Matrix{Float64}}
end

struct _DijkstraClosure end

_closure(::Val{:dijkstra}, Wc) = _DijkstraClosure()
_closure(::Val{M}, Wc) where {M} = throw(ArgumentError("coarse closure method $(M) is not available" *
    (M === :pathpack ? " (load PATHPACK)" : "")))
_apply_closure(::_DijkstraClosure, Wc, bc) = minplus_closure(Wc, bc)

# The whole closure Wc*: column C is the cost from every block to C.
function _full_closure(::_DijkstraClosure, Wc)
    nc = size(Wc, 1)
    table = Matrix{Float64}(undef, nc, nc)
    Threads.@threads for c in 1:nc
        e = fill(Inf, nc)
        e[c] = 0.0
        table[:, c] = minplus_closure(Wc, e)
    end
    return table
end

function Abstraction(W::SparseMatrixCSC, ψ::GraphHomomorphism; steps::Integer = 1, method::Symbol = :dijkstra,
        table::Bool = false)
    Wc = SparseMatrixCSC{Float64,Int}(quotient_graph(W, ψ; steps))
    closure = _closure(Val(method), Wc)
    return Abstraction(W, ψ, Int(steps), Wc, closure, table ? _full_closure(closure, Wc) : nothing)
end

"""
    lower_bound(A::Abstraction, targets) -> Vector

A lower bound on the cost-to-go of every state to `targets`: the coarse
cost-to-go of its block, with terminal cost the best fewer-than-`steps` move into
the targets from each block (from the stored all-pairs table if there is one).
For `steps = 1` it is a consistent heuristic; for more steps an admissible one.
"""
function lower_bound(A::Abstraction, targets)
    N = size(A.W, 1)
    @argcheck !isempty(targets) && all(t -> 1 <= t <= N, targets) "targets must be states of W"
    b = fill(Inf, N)
    b[targets] .= 0
    bc = tropical_pushforward(A.ψ, _short_horizon(A.W, b, A.steps))
    if isnothing(A.table)
        vc = _apply_closure(A.closure, A.quotient, bc)
    else
        vc = fill(Inf, length(bc))
        for c in eachindex(bc)
            isinf(bc[c]) || (vc .= min.(vc, view(A.table, :, c) .+ bc[c]))
        end
    end
    return tropical_pullback(A.ψ, vc)
end

# A* from `source` with the out-edges given by `each_out(f, x)` (calling f(y, w) for
# every edge x → y of cost w), reopening states when a cheaper path is found.
function _astar(each_out, N, source, targets, heuristic, ::Type{S}) where {S}
    h(x) = isnothing(heuristic) ? zero(S) : S(heuristic[x])
    is_target = falses(N)
    is_target[collect(targets)] .= true
    g = fill(typemax(S), N)
    parent = zeros(Int, N)
    closed = falses(N)
    g[source] = 0
    heap = Tuple{S,Int}[]
    isinf(h(source)) || _push!(heap, (h(source), source))
    expanded = 0
    while !isempty(heap)
        key, x = _pop!(heap)
        (closed[x] || key > g[x] + h(x)) && continue
        closed[x] = true
        expanded += 1
        if is_target[x]
            path = [x]
            while path[end] != source
                push!(path, parent[path[end]])
            end
            return (cost = g[x], path = reverse!(path), expanded = expanded)
        end
        each_out(x) do y, w
            d = g[x] + w
            if d < g[y]
                g[y] = d
                parent[y] = x
                closed[y] = false                                   # reopen
                isinf(h(y)) || _push!(heap, (d + h(y), y))
            end
        end
    end
    return (cost = typemax(S), path = Int[], expanded = expanded)
end

function _check_query(N, source, targets, heuristic)
    @argcheck 1 <= source <= N "source must be a state of W"
    @argcheck !isempty(targets) && all(t -> 1 <= t <= N, targets) "targets must be states of W"
    @argcheck isnothing(heuristic) || length(heuristic) == N "one heuristic value per state"
end

"""
    shortest_path(W, source, targets; heuristic = nothing, out_edges = sparse(transpose(W)))
        -> (cost, path, expanded)

The cheapest path from `source` to `targets` by A* (Dijkstra's algorithm without
a `heuristic`), searching forward along the out-edges (`out_edges`, the transpose
of `W`, can be built once for many queries). The `heuristic` must be an
admissible lower bound on the cost-to-go to the targets, such as
[`lower_bound`](@ref); states are reopened when a cheaper path is found, so an
inconsistent heuristic still gives the exact `cost`. `expanded` counts state
expansions, reopenings included.
"""
function shortest_path(W::SparseMatrixCSC{T}, source::Integer, targets; heuristic::Union{Nothing,AbstractVector} = nothing,
        out_edges::SparseMatrixCSC = sparse(transpose(W))) where {T<:Real}
    _check_costs(W)
    N = size(W, 1)
    _check_query(N, source, targets, heuristic)
    @argcheck size(out_edges) == size(W) "out_edges must be the transpose of W"
    rows, vals = rowvals(out_edges), nonzeros(out_edges)
    each_out(f, x) = for k in nzrange(out_edges, x)
        f(rows[k], vals[k])
    end
    return _astar(each_out, N, source, targets, heuristic, float(T))
end

"""
    schur_complement(W, S) -> SparseMatrixCSC

Exact elimination of the states outside `S` in min-plus: the matrix on `S` (rows
and columns in the order of `S`)

```math
W_S = W_{SS} ⊕ W_{SI} ⊗ W_{II}^* ⊗ W_{IS},
```

whose entry ``(s, t)`` is the cheapest path from ``s`` to ``t`` through states of
``I``, the complement of `S` (or the direct edge). In any closed semiring the
closure of the Schur complement is the restriction of the closure,
``(W_S)^* = (W^*)_{SS}`` (Lehmann 1977): eliminating the interior loses nothing on
`S`. One backward search through the interior per state of `S`, threaded.
"""
function schur_complement(W::SparseMatrixCSC{T}, S::AbstractVector{<:Integer}) where {T<:Real}
    _check_costs(W)
    N = size(W, 1)
    @argcheck allunique(S) && all(s -> 1 <= s <= N, S) "S must be distinct states of W"
    local_index = zeros(Int, N)
    local_index[S] .= eachindex(S)
    rows, vals = rowvals(W), nonzeros(W)
    F = float(T)
    columns = Vector{Vector{Tuple{Int,F}}}(undef, length(S))
    Threads.@threads for b in eachindex(S)
        g = Dict{Int,F}()                  # interior labels
        best = Dict{Int,F}()               # entries (a, b) of the column
        heap = Tuple{F,Int}[]
        relax(x, d) = if local_index[x] > 0
            d < get(best, x, typemax(F)) && (best[x] = d)
        elseif d < get(g, x, typemax(F))
            g[x] = d
            _push!(heap, (d, x))
        end
        for k in nzrange(W, S[b])
            relax(rows[k], vals[k])
        end
        settled = Set{Int}()
        while !isempty(heap)
            d, y = _pop!(heap)
            (y in settled || d > g[y]) && continue
            push!(settled, y)
            for k in nzrange(W, y)
                relax(rows[k], d + vals[k])
            end
        end
        columns[b] = [(local_index[x], d) for (x, d) in best]
    end
    I_, J_, V_ = Int[], Int[], F[]
    for b in eachindex(columns), (a, d) in columns[b]
        push!(I_, a); push!(J_, b); push!(V_, d)
    end
    return sparse(I_, J_, V_, length(S), length(S))
end

"""
    Overlay(W, cell)

Exact elimination of the cell interiors (customizable route planning, Delling,
Goldberg, Pajor & Werneck 2011): `cell[x]` assigns each state to a cell, the
*portals* are the endpoints of the edges between cells, and the overlay is the
[`schur_complement`](@ref) of `W` on the portals, a shortcut for the cheapest
path between two portals through a cell's interior (plus the edges between
cells). A query searches the fine graph only inside the cells of the source and
the targets and the overlay everywhere else ([`overlay_path`](@ref)); the cells
are the cover, and the Schur complement is the min-plus harmonic extension of
the boundary data through each cell.
"""
struct Overlay{G}
    W::G
    W_out::SparseMatrixCSC{Float64,Int}
    cell::Vector{Int}
    portals::Vector{Int}
    shortcuts_out::SparseMatrixCSC{Float64,Int}    # transpose of the overlay, on all states
end

function Overlay(W::SparseMatrixCSC, cell::AbstractVector{<:Integer})
    N = size(W, 1)
    @argcheck length(cell) == N "one cell per state"
    rows, cols, _ = findnz(W)
    portal = falses(N)
    for (x, y) in zip(rows, cols)
        if cell[x] != cell[y]
            portal[x] = true
            portal[y] = true
        end
    end
    portals = findall(portal)
    WS = schur_complement(W, portals)
    a, b, v = findnz(WS)
    keep = a .!= b
    shortcuts = sparse(portals[a[keep]], portals[b[keep]], v[keep], N, N)
    return Overlay(W, SparseMatrixCSC{Float64,Int}(sparse(transpose(W))), collect(Int, cell), portals,
        sparse(transpose(shortcuts)))
end

"""
    overlay_path(O::Overlay, source, targets; heuristic = nothing) -> (cost, path, expanded)

The exact cheapest cost from `source` to `targets` on the overlay `O`: A* over the
fine edges of the cells containing the source or a target and the overlay
shortcuts between portals. The `path` lists the states visited, with a shortcut
standing for a path through a cell interior. Any admissible `heuristic` for the
fine graph, such as [`lower_bound`](@ref), is admissible here.
"""
function overlay_path(O::Overlay, source::Integer, targets; heuristic::Union{Nothing,AbstractVector} = nothing)
    N = size(O.W, 1)
    _check_query(N, source, targets, heuristic)
    active = Set(O.cell[[source; collect(targets)]])
    fr, fv = rowvals(O.W_out), nonzeros(O.W_out)
    sr, sv = rowvals(O.shortcuts_out), nonzeros(O.shortcuts_out)
    function each_out(f, x)
        for k in nzrange(O.shortcuts_out, x)
            f(sr[k], sv[k])
        end
        if O.cell[x] in active
            for k in nzrange(O.W_out, x)
                f(fr[k], fv[k])
            end
        end
    end
    return _astar(each_out, N, source, targets, heuristic, Float64)
end

"""
    projection_homomorphism(points, dims) -> GraphHomomorphism

The projection of a grid of `points` (column-major) onto the dimensions `dims`,
e.g. ``(q, v) ↦ q`` for `dims = (1, 2)` of a ``(q_x, q_y, v_x, v_y)`` grid: each
state goes to the cell of its coordinates in `dims`. Its fibres are the full
ranges of the other dimensions.
"""
function projection_homomorphism(points::NTuple{D,Integer}, dims::NTuple{K,Integer}) where {D,K}
    @argcheck allunique(dims) && all(d -> 1 <= d <= D, dims) "dims must be distinct dimensions of the grid"
    cells = ntuple(k -> Int(points[dims[k]]), K)
    Lc = LinearIndices(cells)
    vmap = [Lc[CartesianIndex(ntuple(k -> I[dims[k]], K))] for I in CartesianIndices(Tuple(points))]
    return GraphHomomorphism(vec(vmap), prod(cells))
end

# ---- min-plus sheaves over a base graph, and morphisms between them ----------

# Entrywise A ≤ B for min-plus sparse matrices (absent entries are ∞): every entry
# stored in B is stored in A with a value no larger. Rows are sorted in each column.
function _minplus_le(A::SparseMatrixCSC, B::SparseMatrixCSC; tol = 1e-12)
    size(A) == size(B) || return false
    Ar, Av, Br, Bv = rowvals(A), nonzeros(A), rowvals(B), nonzeros(B)
    for j in 1:size(B, 2)
        p, pend = first(nzrange(A, j)), last(nzrange(A, j))
        for q in nzrange(B, j)
            while p <= pend && Ar[p] < Br[q]
                p += 1
            end
            (p <= pend && Ar[p] == Br[q] && Av[p] <= Bv[q] + tol) || return false
        end
    end
    return true
end

# The sparsity pattern sorted (sparse() sorts; _minplus_mul sorts) and, for equality,
# the same stored entries with the same values.
_minplus_eq(A, B; tol = 1e-12) = _minplus_le(A, B; tol) && _minplus_le(B, A; tol)

"""
    intertwines(W_source, W_target, Φ) -> (lax, oplax, strict)

How the cochain map `Φ` (a min-plus matrix from the states of the source system to
those of the target, absent entries ``∞``) commutes with the two systems' Bellman
operators ``V ↦ W ⊗ V``: `lax` when ``W_t ⊗ Φ ≤ Φ ⊗ W_s`` entrywise (the target's
moves are no more expensive: a relaxation), `oplax` when ``Φ ⊗ W_s ≤ W_t ⊗ Φ``,
and `strict` when both hold, ``W_t ⊗ Φ = Φ ⊗ W_s``: a morphism of min-plus sheaves
in the ordinary sense. Over a base graph, with ``Φ`` acting stalk by stalk, the
intertwining relation is the commuting square of the restriction maps of every
edge; lax and oplax morphisms commute up to the order of the min-plus semiring.
"""
function intertwines(Ws::SparseMatrixCSC, Wt::SparseMatrixCSC, Φ::SparseMatrixCSC)
    @argcheck size(Φ) == (size(Wt, 1), size(Ws, 1)) "Φ maps the source states to the target states"
    left = _minplus_mul(SparseMatrixCSC{Float64,Int}(Wt), SparseMatrixCSC{Float64,Int}(Φ))     # W_t ⊗ Φ
    right = _minplus_mul(SparseMatrixCSC{Float64,Int}(Φ), SparseMatrixCSC{Float64,Int}(Ws))    # Φ ⊗ W_s
    lax = _minplus_le(left, right)
    oplax = _minplus_le(right, left)
    return (lax = lax, oplax = oplax, strict = lax && oplax)
end

"""
    TruncatedPushforward(W, points, cell_dims, fibre_modes)

The min-plus pushforward of a state graph over a grid of `points` onto the cells of
the dimensions `cell_dims`, keeping a *partition of each stalk into modes*.

**The fine sheaf.** Pushing the state graph `W` forward along the projection
``ψ`` onto the cells ([`projection_homomorphism`](@ref)) gives a min-plus sheaf
``F`` on the cell graph: its stalk over a cell ``q`` holds the values on the fibre
``ψ^{-1}(q)`` (for ``(q, v) ↦ q``, a function of the velocity), and the edge
``q → q'`` carries the block of `W` from the fibre over ``q`` to the fibre over
``q'``, moves inside a cell being loops. The value function is the greatest *lax
section* of ``F`` below the terminal cost: ``V(q) ≤ W_{q q'} ⊗ V(q')`` on every
edge. Keeping one value per stalk (its minimum) is [`quotient_graph`](@ref) along
``ψ``, which forgets the velocity.

**The coarse sheaf.** `fibre_modes` (an array over the other dimensions, the same
for every cell) partitions each stalk into ``K`` modes, e.g. heading × speed bins
of the velocity. The coarse sheaf ``G`` has stalk ``ℝ_{\\min}^K`` over every cell,
and the two maps between them are stalkwise (over the identity of the cell graph):

- *refinement* ``E : G → F`` ([`refine`](@ref)), copying each mode's value to its
  states (the pullback along the mode map ``μ``);
- *coarsening* ``P : F → G`` ([`coarsen`](@ref)), the minimum over each mode.

They form a Galois connection, ``E g ≤ f ⟺ g ≤ P f`` (``P`` is the residual of
``E``), with ``P E = I`` and ``E P ≤ I``. The coarse edge maps are the residuated
Galerkin products ``W_c = P ⊗ W ⊗ E``, so ``P`` is a *lax* morphism of min-plus
sheaves (``W_c ⊗ P ≤ P ⊗ W``) and ``E`` an *oplax* one (``E ⊗ W_c ≤ W ⊗ E``),
both checked by [`intertwines`](@ref). Neither square commutes strictly unless the
moves are constant on modes. The oplax square is what makes the coarse solution a
bound: ``E`` maps lax sections of ``G`` to lax sections of ``F``, so ``E V_c ≤ V``.
With one mode per stalk this is the velocity-blind quotient; with one mode per
state it is an isomorphism.
"""
struct TruncatedPushforward
    W::SparseMatrixCSC{Float64,Int}
    ψ::GraphHomomorphism                  # state → cell
    μ::GraphHomomorphism                  # state → coarse state (cell, mode)
    π::GraphHomomorphism                  # coarse state → cell
    K::Int                                # modes per stalk
    quotient::SparseMatrixCSC{Float64,Int}    # W_c = P ⊗ W ⊗ E, loops included
    cell_points::Tuple
    mode_neighbors::Vector{Vector{Int}}   # modes containing adjacent fibre points
end

function TruncatedPushforward(W::SparseMatrixCSC, points::NTuple{D,Integer}, cell_dims::NTuple{C,Integer},
        fibre_modes::AbstractArray{<:Integer}) where {D,C}
    _check_costs(W)
    @argcheck size(W, 1) == prod(points) "W must act on the grid of points"
    fdims = Tuple(d for d in 1:D if !(d in cell_dims))
    @argcheck size(fibre_modes) == ntuple(k -> Int(points[fdims[k]]), length(fdims)) "fibre_modes must cover the other dimensions"
    K = maximum(fibre_modes)
    @argcheck sort(unique(vec(fibre_modes))) == collect(1:K) "fibre_modes must use the modes 1:K"
    ψ = projection_homomorphism(points, cell_dims)
    LF = LinearIndices(size(fibre_modes))
    mode = [fibre_modes[LF[CartesianIndex(ntuple(k -> I[fdims[k]], length(fdims)))]] for I in CartesianIndices(Tuple(points))]
    μ = GraphHomomorphism(vec((ψ.vertex_map .- 1) .* K .+ vec(mode)), ψ.n_target * K)
    π = GraphHomomorphism([(c - 1) ÷ K + 1 for c in 1:(ψ.n_target * K)], ψ.n_target)
    rows, cols, vals = findnz(W)
    Wc = sparse(μ.vertex_map[rows], μ.vertex_map[cols], Float64.(vals), μ.n_target, μ.n_target, min)
    # Modes are adjacent when they hold Chebyshev-adjacent points of the fibre grid.
    nbrs = [Set{Int}() for _ in 1:K]
    CF = CartesianIndices(size(fibre_modes))
    for I in CF, J in CF
        if maximum(abs.(Tuple(I) .- Tuple(J))) == 1 && fibre_modes[I] != fibre_modes[J]
            push!(nbrs[fibre_modes[I]], fibre_modes[J])
        end
    end
    cells = ntuple(k -> Int(points[cell_dims[k]]), C)
    return TruncatedPushforward(SparseMatrixCSC{Float64,Int}(W), ψ, μ, π, K, Wc, cells, [sort!(collect(s)) for s in nbrs])
end

"""
    coarsen(T::TruncatedPushforward, v) -> Vector

The coarsening ``P``: each mode of each stalk takes the minimum of `v` over its
states.
"""
coarsen(T::TruncatedPushforward, v::AbstractVector) = tropical_pushforward(T.μ, v)

"""
    refine(T::TruncatedPushforward, vc) -> Vector

The refinement ``E``: each state takes the value of its mode.
"""
refine(T::TruncatedPushforward, vc::AbstractVector) = tropical_pullback(T.μ, vc)

"""
    coarsening_matrix(T) -> SparseMatrixCSC
    refinement_matrix(T) -> SparseMatrixCSC

``P`` and ``E`` as min-plus matrices (stored entries 0, the semiring one; absent
entries ``∞``), for [`intertwines`](@ref).
"""
function refinement_matrix(T::TruncatedPushforward)
    N = length(T.μ.vertex_map)
    return sparse(1:N, T.μ.vertex_map, zeros(N), N, T.μ.n_target)
end
coarsening_matrix(T::TruncatedPushforward) = sparse(transpose(refinement_matrix(T)))

"""
    LayeredPlanner(T::TruncatedPushforward)
    LayeredPlanner(W, ψ, cell_points)

A two-layer controller: a *planner* on the coarse sheaf of `T` (routes are
sequences of cells, each with a mode of its stalk) and a *low-level controller*
on the full state graph. The second form keeps one mode per stalk (the
velocity-blind projection along `ψ` onto the grid `cell_points`).

The low-level controller searches exactly over the pullback ``μ^{-1}(U)`` of a
*corridor* ``U`` of coarse states around the route: cells within Chebyshev
distance `width` of the route (with the segments between consecutive route cells
filled in), with modes within `mode_width` steps of the route's modes. ``U`` is a
union of fibres of ``μ``, so restriction to it commutes with the pushforward
(base change) and the coarse bound ([`lower_bound`](@ref) of `abstraction`) stays
admissible inside it. The corridor is a commitment: the best path inside it can
cost more than the optimum, or not exist; [`layered_path`](@ref) reports both and
can widen it.
"""
struct LayeredPlanner{G}
    abstraction::Abstraction
    μ::GraphHomomorphism
    π::GraphHomomorphism
    K::Int
    cell_points::Tuple
    mode_neighbors::Vector{Vector{Int}}
    W_out::SparseMatrixCSC{Float64,Int}
    quotient_out::SparseMatrixCSC{Float64,Int}
    W::G
end

function LayeredPlanner(T::TruncatedPushforward)
    A = Abstraction(T.W, T.μ)
    return LayeredPlanner(A, T.μ, T.π, T.K, T.cell_points, T.mode_neighbors, SparseMatrixCSC{Float64,Int}(sparse(transpose(T.W))),
        sparse(transpose(A.quotient)), T.W)
end

function LayeredPlanner(W::SparseMatrixCSC, ψ::GraphHomomorphism, cell_points::Tuple)
    @argcheck prod(cell_points) == ψ.n_target "cell_points must have one cell per block of ψ"
    A = Abstraction(W, ψ)
    return LayeredPlanner(A, ψ, GraphHomomorphism(collect(1:ψ.n_target), ψ.n_target), 1, Tuple(Int.(cell_points)),
        [Int[]], SparseMatrixCSC{Float64,Int}(sparse(transpose(W))), sparse(transpose(A.quotient)), W)
end

"""
    plan_route(P::LayeredPlanner, source, targets) -> Vector{Int}

The planner's route: the cheapest sequence of coarse states (cells with a mode)
from that of `source` to those of `targets` on the coarse graph (empty if there is
none).
"""
function plan_route(P::LayeredPlanner, source::Integer, targets)
    μ = P.μ
    return shortest_path(P.abstraction.quotient, μ.vertex_map[source], unique(μ.vertex_map[collect(targets)]);
        out_edges = P.quotient_out).path
end

# The cells of the segment from cell a to cell b, rasterized.
function _segment(C, a, b)
    ca, cb = Tuple(C[a]), Tuple(C[b])
    steps = maximum(abs.(cb .- ca))
    return [ca; [round.(Int, ca .+ (cb .- ca) .* (s / steps)) for s in 1:steps]]
end

"""
    corridor(cell_points, route, width) -> BitVector
    corridor(P::LayeredPlanner, route, width, mode_width) -> BitVector

The cells within Chebyshev distance `width` of the polyline through the cells of
`route` on the grid `cell_points`; with a planner, the coarse states over those
cells whose mode is within `mode_width` steps (in the mode adjacency) of a mode
of the route at either end of the nearby segment.
"""
function corridor(cell_points::Tuple, route::AbstractVector{<:Integer}, width::Integer)
    @argcheck width >= 0 "width must be nonnegative"
    C = CartesianIndices(cell_points)
    Lc = LinearIndices(cell_points)
    inside = falses(prod(cell_points))
    isempty(route) && return inside
    pts = length(route) == 1 ? [Tuple(C[route[1]])] :
        reduce(vcat, [_segment(C, a, b) for (a, b) in zip(route[1:(end - 1)], route[2:end])])
    for c in pts, J in CartesianIndices(ntuple(k -> max(1, c[k] - width):min(cell_points[k], c[k] + width), length(cell_points)))
        inside[Lc[J]] = true
    end
    return inside
end

function corridor(P::LayeredPlanner, route::AbstractVector{<:Integer}, width::Integer, mode_width::Integer)
    @argcheck width >= 0 && mode_width >= 0 "widths must be nonnegative"
    K = P.K
    cellof(c) = P.π.vertex_map[c]
    modeof(c) = (c - 1) % K + 1
    function dilate(modes)
        S = Set(modes)
        for _ in 1:mode_width
            S = union(S, (n for m in S for n in P.mode_neighbors[m])...)
        end
        return S
    end
    inside = falses(P.μ.n_target)
    isempty(route) && return inside
    pieces = length(route) == 1 ? [(route[1], route[1])] : collect(zip(route[1:(end - 1)], route[2:end]))
    for (a, b) in pieces
        cells = findall(corridor(P.cell_points, [cellof(a), cellof(b)], width))
        for m in dilate([modeof(a), modeof(b)]), c in cells
            inside[(c - 1) * K + m] = true
        end
    end
    return inside
end

"""
    layered_path(P::LayeredPlanner, source, targets; width = 2, mode_width = 0,
        heuristic = nothing, widen = true)
        -> (cost, path, expanded, route, width, mode_width, corridor_states, first_try)

Plan a route on the coarse sheaf, then search exactly (A* with `heuristic`, e.g.
[`lower_bound`](@ref) of `P.abstraction`) over the states above the corridor
around it. If the corridor holds no path and `widen` is set, the width doubles and
the mode width grows by one until it does (at worst the corridor is everything
and the search is exact). `first_try` reports whether the requested corridor
sufficed; `cost` is that of the path found, at least the optimum.
"""
function layered_path(P::LayeredPlanner, source::Integer, targets; width::Integer = 2, mode_width::Integer = 0,
        heuristic::Union{Nothing,AbstractVector} = nothing, widen::Bool = true)
    N = size(P.W, 1)
    _check_query(N, source, targets, heuristic)
    route = plan_route(P, source, targets)
    # No route on the relaxation: no path at all.
    isempty(route) && return (cost = Inf, path = Int[], expanded = 0, route = route, width = Int(width),
        mode_width = Int(mode_width), corridor_states = 0, first_try = false)
    rows, vals = rowvals(P.W_out), nonzeros(P.W_out)
    w, mw = Int(width), Int(mode_width)
    first_try = true
    while true
        U = corridor(P, route, w, mw)
        allowed(x) = U[P.μ.vertex_map[x]]
        each_out(f, x) = for k in nzrange(P.W_out, x)
            y = rows[k]
            allowed(y) && f(y, vals[k])
        end
        reachable_targets = [t for t in targets if allowed(t)]
        r = isempty(reachable_targets) || !allowed(source) ? (cost = Inf, path = Int[], expanded = 0) :
            _astar(each_out, N, source, reachable_targets, heuristic, Float64)
        if isfinite(r.cost) || !widen || all(U)
            return (cost = r.cost, path = r.path, expanded = r.expanded, route = route, width = w, mode_width = mw,
                corridor_states = count(x -> U[x], P.μ.vertex_map), first_try = first_try && isfinite(r.cost))
        end
        first_try = false
        w = max(1, 2w)
        mw = min(mw + 1, P.K)
    end
end

# ---- Lyapunov certificates ----------------------------------------------------

"""
    robust_quotient(W, μ) -> SparseMatrixCSC

The *pessimistic* (∀∃, "must") pushforward of the costs `W` along `μ`: an edge
``B → C`` exists when *every* state of ``B`` has a move into ``C``, and costs the
worst of their cheapest such moves,

```math
W_c^{∀∃}(B, C) = \\max_{x ∈ B} \\min_{y ∈ C} W_{xy}.
```

Against the optimistic [`quotient_graph`](@ref) (``\\min`` over both ends), which
makes the refinement ``E`` an *oplax* morphism and the coarse value a lower bound,
the pessimistic pushforward makes ``E`` a *lax* morphism of min-plus sheaves,
``W ⊗ E ≤ E ⊗ W_c^{∀∃}`` ([`intertwines`](@ref)): every coarse move can be
realized from every state of its block. Lax morphisms transport control
Lyapunov functions from the coarse system to the fine one
([`upper_bound`](@ref)). This is the "must" abstraction of modal transition
systems (Larsen & Thomsen 1988) and of abstraction-based controller synthesis
(Tabuada 2009).
"""
function robust_quotient(W::SparseMatrixCSC{T}, μ::GraphHomomorphism) where {T<:Real}
    _check_costs(W)
    @argcheck size(W, 1) == length(μ.vertex_map) "μ must map the states of W"
    N, nc = size(W, 1), μ.n_target
    r, c, v = findnz(W)
    M = sparse(r, μ.vertex_map[c], Float64.(v), N, nc, min)          # W ⊗ E: the cheapest move from x into C
    blocksize = zeros(Int, nc)
    for b in μ.vertex_map
        blocksize[b] += 1
    end
    count = Dict{Tuple{Int,Int},Int}()
    worst = Dict{Tuple{Int,Int},Float64}()
    for (x, C, w) in zip(findnz(M)...)
        B = μ.vertex_map[x]
        B == C && continue
        count[(B, C)] = get(count, (B, C), 0) + 1
        worst[(B, C)] = max(get(worst, (B, C), -Inf), w)
    end
    keys_ = [k for (k, n) in count if n == blocksize[k[1]]]
    return sparse(first.(keys_), last.(keys_), [worst[k] for k in keys_], nc, nc)
end

"""
    upper_bound(T::TruncatedPushforward, targets) -> Vector

The refinement ``E V_c`` of the coarse cost-to-go on the pessimistic quotient
([`robust_quotient`](@ref) along the mode map), with coarse terminal cost the
*largest* fine terminal cost in each mode (zero only on modes inside the
targets). Because ``E`` is a lax morphism, ``E V_c`` is a control Lyapunov
function of the fine system wherever it is finite ([`is_control_lyapunov`](@ref)):
from every such state some move lowers it by at least the move's cost. So it is
an upper bound on the cost-to-go, and the greedy feedback on it
([`greedy_policy`](@ref)) reaches the targets with total cost at most ``E V_c``
([`is_lyapunov`](@ref) certifies the closed loop).
"""
function upper_bound(T::TruncatedPushforward, targets)
    N = size(T.W, 1)
    @argcheck !isempty(targets) && all(t -> 1 <= t <= N, targets) "targets must be states of W"
    b = fill(Inf, N)
    b[targets] .= 0
    bc = tropical_pushforward(T.μ, b; combine = max)
    return refine(T, minplus_closure(robust_quotient(T.W, T.μ), bc))
end

"""
    is_control_lyapunov(W, V, b) -> Bool

Whether `V` is a control Lyapunov function of the state graph `W` for the
terminal cost `b`: wherever ``V`` is finite, ``\\min\\big(b(x), \\min_y W_{xy} + V(y)\\big) ≤ V(x)``,
i.e. ``W ⊗ V ⊕ b ≤ V``: `V` is a lax fixed point (a supersolution) of the min-plus
Bellman operator. Some move from every state off the targets lowers ``V`` by at
least its cost.
"""
function is_control_lyapunov(W::SparseMatrixCSC, V::AbstractVector, b::AbstractVector; tol = 1e-9)
    @argcheck size(W, 1) == length(V) == length(b) "one value per state"
    best = float.(b)
    for (x, y, w) in zip(findnz(W)...)
        best[x] = min(best[x], w + V[y])
    end
    return all(!isfinite(V[x]) || best[x] <= V[x] + tol for x in eachindex(V))
end

"""
    greedy_policy(W, V; out_edges = sparse(transpose(W))) -> Vector{Int}

The feedback that moves each state to the successor minimizing ``W_{xy} + V(y)``
(0 where no successor has finite value): the closed loop as a deterministic
transition system, a coalgebra ``X → \\mathcal P(X)`` with singleton successors.
"""
function greedy_policy(W::SparseMatrixCSC, V::AbstractVector; out_edges::SparseMatrixCSC = sparse(transpose(W)))
    @argcheck size(W, 1) == length(V) "one value per state"
    rows, vals = rowvals(out_edges), nonzeros(out_edges)
    next = zeros(Int, length(V))
    for x in eachindex(V)
        best = Inf
        for k in nzrange(out_edges, x)
            d = vals[k] + V[rows[k]]
            if d < best
                best, next[x] = d, rows[k]
            end
        end
    end
    return next
end

"""
    is_lyapunov(next, V, equilibria; decrease = 1) -> (positive_definite, decrescent, strict, closed, domain)

Check that `V` is a Lyapunov morphism for the closed loop `next` (a deterministic
transition system, as from [`greedy_policy`](@ref)) relative to the set
`equilibria`, on the domain where ``V`` is finite, in the discrete setting of
categorical Lyapunov theory (Ames, Moeller & Tabuada 2025; Ames, Mattenet & Moeller
2025): the stability object is ``ℝ_{≥0}`` with the stationary system
``r ↦ \\{r\\}`` and the order ``U ≤ V`` iff every element of ``U`` is at most every
element of ``V``.

- `positive_definite`: ``V = 0`` on the equilibria and ``V > 0`` elsewhere;
- `decrescent`: the square ``\\mathcal P V ∘ \\text{next} ≤ σ ∘ V`` lax commutes,
  ``V(\\text{next}(x)) ≤ V(x)``;
- `strict`: ``V(\\text{next}(x)) ≤ V(x) - `` `decrease` off the equilibria, so every
  trajectory reaches an equilibrium within ``V(x_0) /`` `decrease` steps
  (finite-time stability);
- `closed`: the domain is invariant.
"""
function is_lyapunov(next::AbstractVector{<:Integer}, V::AbstractVector, equilibria; decrease::Real = 1)
    @argcheck length(next) == length(V) "one successor per state"
    eq = falses(length(V))
    eq[collect(equilibria)] .= true
    dom = findall(isfinite, V)
    off = [x for x in dom if !eq[x]]
    positive = all(V[x] == 0 for x in dom if eq[x]) && all(V[x] > 0 for x in off)
    closed = all(next[x] != 0 && isfinite(V[next[x]]) for x in off)
    decrescent = closed && all(V[next[x]] <= V[x] for x in off)
    strict = decrescent && all(V[next[x]] <= V[x] - decrease for x in off)
    return (positive_definite = positive, decrescent = decrescent, strict = strict, closed = closed, domain = length(dom))
end

# The min-plus closure of W restricted to the `allowed` states (backward Dijkstra
# over in-edges, never entering a state outside), and the greedy feedback on it.
function _restricted_closure(W::SparseMatrixCSC, allowed::AbstractVector{Bool}, b::AbstractVector)
    g = fill(Inf, size(W, 1))
    heap = Tuple{Float64,Int}[]
    for x in eachindex(b)
        if allowed[x] && isfinite(b[x])
            g[x] = b[x]
            _push!(heap, (g[x], x))
        end
    end
    rows, vals = rowvals(W), nonzeros(W)
    settled = falses(length(g))
    while !isempty(heap)
        _, y = _pop!(heap)
        settled[y] && continue
        settled[y] = true
        for k in nzrange(W, y)
            x = rows[k]
            allowed[x] || continue
            d = g[y] + vals[k]
            if d < g[x]
                g[x] = d
                _push!(heap, (d, x))
            end
        end
    end
    return g
end

function _restricted_greedy(W_out::SparseMatrixCSC, allowed::AbstractVector{Bool}, V::AbstractVector)
    rows, vals = rowvals(W_out), nonzeros(W_out)
    next = zeros(Int, length(V))
    for x in eachindex(V)
        (allowed[x] && isfinite(V[x])) || continue
        best = Inf
        for k in nzrange(W_out, x)
            y = rows[k]
            allowed[y] || continue
            d = vals[k] + V[y]
            if d < best
                best, next[x] = d, y
            end
        end
    end
    return next
end

_above(P::LayeredPlanner, U) = BitVector(U[P.μ.vertex_map[x]] for x in 1:size(P.W, 1))

"""
    corridor_certificate(P::LayeredPlanner, source, targets; width, mode_width) -> (V, next)

The Lyapunov certificate of the layered controller's low level: the cost-to-go
``V_U`` to the targets *within* the corridor (the min-plus closure of `W`
restricted to the states above it), and the greedy feedback on it, whose closed
loop [`is_lyapunov`](@ref) certifies on ``\\{V_U < ∞\\}``. The search of
[`layered_path`](@ref) with the same corridor finds a path of cost ``V_U(\\text{source})``.
"""
function corridor_certificate(P::LayeredPlanner, source::Integer, targets; width::Integer, mode_width::Integer = 0)
    inside = _above(P, corridor(P, plan_route(P, source, targets), width, mode_width))
    b = fill(Inf, size(P.W, 1))
    b[[t for t in targets if inside[t]]] .= 0
    VU = _restricted_closure(P.W, inside, b)
    return (V = VU, next = _restricted_greedy(P.W_out, inside, VU))
end

"""
    certified_path(P::LayeredPlanner, source, targets; width = 2, mode_width = 0,
        lower = nothing, widen = true)
        -> (cost, path, route, width, mode_width, first_try, certificate, lower, suboptimality_bound)

The layered controller with its certificates, sandwiched between the two
pushforwards. The low level is the greedy feedback on the cost-to-go ``V_U``
within the corridor around the planned route (widened as in [`layered_path`](@ref)
until the source reaches the targets in it), and `path` is its closed-loop
trajectory from `source`, of cost ``V_U(\\text{source})``.

- *Stability* (the lax side): `certificate` is [`is_lyapunov`](@ref) for that
  closed loop with ``V_U``, a Lyapunov morphism into ``(ℝ_{≥0}, r ↦ \\{r\\})``
  with strict decrease, so every state above the corridor that can reach the
  targets does, within ``V_U`` steps.
- *Performance* (the oplax side): with `lower` a lower bound on the optimal
  cost-to-go, such as [`lower_bound`](@ref) of an optimistic abstraction (the
  refinement along an oplax morphism), `suboptimality_bound` ``= V_U(x) / \\text{lower}(x)``
  is a guaranteed bound on the ratio of the controller's cost to the optimum,
  computed without the optimum.
"""
function certified_path(P::LayeredPlanner, source::Integer, targets; width::Integer = 2, mode_width::Integer = 0,
        lower::Union{Nothing,AbstractVector} = nothing, widen::Bool = true)
    N = size(P.W, 1)
    _check_query(N, source, targets, lower)
    route = plan_route(P, source, targets)
    w, mw = Int(width), Int(mode_width)
    first_try = true
    while true
        U = isempty(route) ? trues(P.μ.n_target) : corridor(P, route, w, mw)
        inside = _above(P, U)
        goal_inside = [t for t in targets if inside[t]]
        b = fill(Inf, N)
        b[goal_inside] .= 0
        VU = _restricted_closure(P.W, inside, b)
        if isfinite(VU[source]) || !widen || all(U)
            next = _restricted_greedy(P.W_out, inside, VU)
            path = [source]
            while isfinite(VU[path[end]]) && VU[path[end]] > 0
                push!(path, next[path[end]])
            end
            certificate = is_lyapunov(next, VU, goal_inside)
            lb = isnothing(lower) ? NaN : Float64(lower[source])
            bound = isnothing(lower) ? NaN : VU[source] == 0 ? 1.0 : VU[source] / lb
            return (cost = VU[source], path = isfinite(VU[source]) ? path : Int[], route = route, width = w, mode_width = mw,
                first_try = first_try && isfinite(VU[source]), certificate = certificate, lower = lb,
                suboptimality_bound = bound)
        end
        first_try = false
        w = max(1, 2w)
        mw = min(mw + 1, P.K)
    end
end

"""
    maxmin_closure(A, b) -> Vector

The least solution, in the bottleneck semiring ``(\\max, \\min)``, of
``V(x) = \\max\\big(b(x), \\max_y \\min(A_{xy}, V(y))\\big)`` (absent entries of
`A` are ``-∞``): the best bottleneck value of a path from `x` to a state with
terminal value `b`. A Dijkstra-type label-setting search in decreasing order of
value, over the in-edges (columns) of `A`.
"""
function maxmin_closure(A::SparseMatrixCSC{T}, b::AbstractVector{<:Real}) where {T<:Real}
    @argcheck size(A, 1) == size(A, 2) == length(b) "A must be square with one terminal value per state"
    S = float(promote_type(T, eltype(b)))
    g = S.(b)
    heap = Tuple{S,Int}[]                                           # keys negated: a max-heap
    for (x, v) in enumerate(g)
        v == typemin(S) || _push!(heap, (-v, x))
    end
    settled = falses(length(g))
    rows, vals = rowvals(A), nonzeros(A)
    while !isempty(heap)
        _, y = _pop!(heap)
        settled[y] && continue
        settled[y] = true
        for k in nzrange(A, y)
            x = rows[k]
            d = min(vals[k], g[y])
            if d > g[x]
                g[x] = d
                _push!(heap, (-d, x))
            end
        end
    end
    return g
end

"""
    reach_avoid(W, l, g) -> Vector

The reach-avoid value on the state graph of `W` (its stored entries are the
allowed moves): with a constraint margin `l` (positive where the state is
allowed) and a target margin `g` (positive in the target),

```math
V(x) = \\max_{\\text{paths } x = x_0 → ⋯ → x_T} \\min\\Big(g(x_T), \\min_{t ≤ T} l(x_t)\\Big),
```

positive exactly on the states that can reach the target while staying allowed.
Its Bellman equation ``V = \\min(l, \\max(g, \\max_y V_y))`` distributes to
``V = \\min(l, g) ⊕ A ⊗ V`` in the bottleneck semiring with ``A_{xy} = l(x)`` on
every move, so ``V = A^* ⊗ \\min(l, g)`` ([`maxmin_closure`](@ref)).
"""
function reach_avoid(W::SparseMatrixCSC, l::AbstractVector{<:Real}, g::AbstractVector{<:Real})
    N = size(W, 1)
    @argcheck size(W, 2) == N && length(l) == N && length(g) == N "one margin per state"
    rows, cols, _ = findnz(W)
    A = sparse(rows, cols, float.(l[rows]), N, N)
    return maxmin_closure(A, min.(float.(l), float.(g)))
end

end # module
