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
    Abstraction, lower_bound, shortest_path, schur_complement, Overlay, overlay_path, maxmin_closure, reach_avoid

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
