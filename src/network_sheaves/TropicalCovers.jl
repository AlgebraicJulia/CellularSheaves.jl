"""
    TropicalCovers

Min-plus sheaves of joint state graphs over a simplicial complex: the setting of
cooperative planning for several agents that share a space.

**The complex.** The vertices are agents and a simplex is a set of agents that may
interact: an edge is a pair, a triangle a triple. Built from an *interaction graph*
(which pairs can meet), the natural complex is its clique complex
([`clique_complex`](@ref)): a set of agents is a simplex when every pair in it
interacts.

**The stalk of a simplex.** Each agent ``a`` has its own state graph ``W_a`` (arc
``s → s'`` costs ``W_a[s, s'] ≥ 0``). The agents of a simplex ``σ`` move at once,
and a *collision* rule removes the joint moves they cannot make together. The
joint state graph of ``σ`` is far too large to build, so the stalk is its
pushforward along a summary map ``μ_σ``, defined *symmetrically* in the agents:

- each agent has its own summary ``π_a`` (a [`GraphHomomorphism`](@ref) of its
  states, for instance its progress and its region of the map);
- a *contact* relation on summaries says when two agents could be close enough to
  collide;
- ``μ_σ`` keeps an agent's exact state when its summary is in contact with that of
  another agent of ``σ``, and replaces it by its summary otherwise.

No agent is privileged: permuting identical agents permutes the coarse states.
Whether an agent is kept exact depends only on summaries, so each coarse state is a
union of fibres, and each fibre is a product of single states and summary fibres.
When the contact relation is conservative (agents out of contact cannot collide in
one move), the joint moves out of a fibre are independent for the agents out of
contact, and the pushforward can be enumerated move by move without the joint graph
([`stalk_successors`](@ref)). Where every agent is out of contact this is the base
change of [`tensor_product`](@ref) and [`quotient_graph`](@ref) to the product of
the agents' quotients.

**Restriction maps.** For a face ``τ ⊂ σ``, forgetting the agents outside ``τ`` maps
the joint states of ``σ`` to those of ``τ``, and on coarse states it is again
determined: an agent kept exact in ``σ`` only because of a contact outside ``τ``
falls back to its summary ([`restrict`](@ref)). The square
``\\mathrm{restrict} ∘ μ_σ = μ_τ ∘ \\mathrm{forget}`` commutes, and the restriction
sends every arc of the stalk of ``σ`` to an arc of the stalk of ``τ`` that costs no
more (the agents of ``τ`` pay part of the cost). So the cost-to-go of a face,
pulled back along the restriction, is a subsolution on ``σ``: a lower bound. The
cost-to-go functions with these pullbacks form a cellular sheaf of min-plus
modules on the complex; its sections are the consistent bounds.

**Gluing over the cover.** A joint plan of all agents restricts to a plan of every
simplex, and its per-agent costs ``x_a`` then satisfy ``∑_{a ∈ σ} x_a ≥ Δ_σ`` for the
cost-to-go ``Δ_σ`` of every stalk. The least such integer 0-cochain is a lower bound
on the joint cost ([`simplicial_cover`](@ref), [`cover_bound`](@ref)). With edges
alone it is the weighted vertex cover of conflict-based search (Li, Felner,
Boyarski, Ma & Koenig 2019); the triangles see three-agent interactions that no
pair sees. Its linear relaxation is dual to splitting each agent's cost among the
simplices at it, a gauge transformation on the cover.

Stalks of low dimension are built explicitly and solved once
([`stalk_graph`](@ref)); higher ones are searched on demand from the coarse state
of a query, with the pulled-back face values as an admissible heuristic
([`stalk_value`](@ref)).
"""
module TropicalCovers

export SimplicialComplex, clique_complex, simplices, faces, ProductStalk, coarse_state, restrict,
    stalk_successors, StalkGraph, stalk_graph, stalk_value, CoverSheaf, section_value, cover_bound,
    simplicial_cover

using ArgCheck: @argcheck
using Graphs
using SparseArrays
using ..GraphHomomorphisms: GraphHomomorphism
using ..TropicalAbstraction: minplus_closure
import ..TropicalAbstraction: _push!, _pop!

# ---------------------------------------------------------------------------
# Simplicial complexes
# ---------------------------------------------------------------------------

"""
    SimplicialComplex(n, simplices)

An abstract simplicial complex on the vertices `1:n`: the given simplices (vectors
of vertices) and all their faces. `simplices(K, k)` lists the `k`-simplices, each a
sorted vector of `k + 1` vertices.
"""
struct SimplicialComplex
    nvertices::Int
    cells::Vector{Vector{Vector{Int}}}          # cells[k + 1]: the k-simplices, sorted
end

function SimplicialComplex(n::Integer, generators)
    @argcheck n >= 0
    sets = [Set{Vector{Int}}() for _ in 1:1]
    push!(sets[1], Int[])                        # placeholder, removed below
    function add!(σ)
        k = length(σ) - 1
        while length(sets) < k + 2
            push!(sets, Set{Vector{Int}}())
        end
        σ in sets[k + 2] && return
        push!(sets[k + 2], σ)
        k > 0 && foreach(add!, faces(σ))
    end
    for v in 1:n
        add!([v])
    end
    for σ in generators
        s = sort!(unique(collect(Int, σ)))
        @argcheck !isempty(s) && all(v -> 1 <= v <= n, s) "simplices must be nonempty sets of vertices in 1:$n"
        add!(s)
    end
    return SimplicialComplex(n, [sort!(collect(c)) for c in sets[2:end]])
end

"""
    clique_complex(g::AbstractGraph; dim = 2)

The clique complex of `g` up to dimension `dim`: a set of vertices is a simplex
when every two of them are adjacent in `g`.
"""
function clique_complex(g::AbstractGraph; dim::Integer = 2)
    @argcheck dim >= 0
    n = nv(g)
    level = [[v] for v in 1:n]
    all_cells = [level]
    for _ in 1:dim
        next = Vector{Int}[]
        for σ in level, w in (last(σ) + 1):n
            all(v -> has_edge(g, v, w), σ) && push!(next, [σ; w])
        end
        isempty(next) && break
        push!(all_cells, next)
        level = next
    end
    return SimplicialComplex(n, all_cells)
end

"""
    simplices(K, k)

The `k`-simplices of `K`, each a sorted vector of `k + 1` vertices (empty beyond its
dimension).
"""
simplices(K::SimplicialComplex, k::Integer) = k + 1 <= length(K.cells) ? K.cells[k + 1] : Vector{Int}[]

"""
    faces(σ)

The codimension-one faces of the simplex `σ` (a vector of vertices).
"""
faces(σ::AbstractVector{<:Integer}) = [deleteat!(collect(Int, σ), k) for k in eachindex(σ)]

Graphs.nv(K::SimplicialComplex) = K.nvertices
dimension(K::SimplicialComplex) = length(K.cells) - 1

function Base.show(io::IO, K::SimplicialComplex)
    print(io, "SimplicialComplex(", join((length(c) for c in K.cells), ", "), " simplices by dimension)")
end

# ---------------------------------------------------------------------------
# Product stalks
# ---------------------------------------------------------------------------

"""
    ProductStalk(members, graphs, summaries, goals; contact, collides)

The stalk of the simplex with agents `members` (vertex labels, passed to `contact`
and `collides`): the pushforward of their joint state graph along the symmetric
summary map described in [`TropicalCovers`](@ref).

- `graphs[a]`: agent `a`'s state graph, `graphs[a][s, s2]` the cost of the move
  `s → s2` (stored entries, zeros included, are the moves; costs ``≥ 0``).
- `summaries[a]`: a [`GraphHomomorphism`](@ref) from agent `a`'s states to its
  summaries.
- `goals[a]`: agent `a`'s goal state. The stalk's target is the coarse state of the
  joint goal.
- `contact(u, k, v, l)`: whether agents `u` and `v` could collide in one move from
  states with summaries `k` and `l`. It must be symmetric and *conservative*: out of
  contact, no move of the two can collide.
- `collides(u, s, s2, v, t, t2)`: whether the moves `s → s2` of `u` and `t → t2` of
  `v`, made at once, collide. `collides(u, s, s, v, t, t)` must say whether the
  states `s` and `t` conflict outright.

A coarse state is a tuple of *modes*, one per agent: a state `s` of agent `a`
(kept exact) or `n_a + k` for its summary `k`, where `n_a` is its number of states.
"""
struct ProductStalk{K,C,F}
    members::NTuple{K,Int}
    graphs::Vector{SparseMatrixCSC{Float64,Int}}
    out::Vector{SparseMatrixCSC{Float64,Int}}           # transposes: out-moves by column
    summaries::Vector{GraphHomomorphism}
    nstates::NTuple{K,Int}
    fibre_moves::Vector{Vector{Vector{Tuple{Int,Float64}}}}   # [a][k]: (s2, least cost from the fibre of k)
    fibres::Vector{Vector{Vector{Int}}}                   # [a][k]: the states with summary k
    goal::NTuple{K,Int}
    contact::C
    collides::F
end

function ProductStalk(members, graphs, summaries, goals; contact, collides)
    K = length(members)
    @argcheck K >= 1 && length(graphs) == K && length(summaries) == K && length(goals) == K
    Ws = [SparseMatrixCSC{Float64,Int}(W) for W in graphs]
    for (W, π, g) in zip(Ws, summaries, goals)
        @argcheck size(W, 1) == size(W, 2) == length(π.vertex_map) "each graph must be square and match its summary"
        @argcheck all(>=(0), nonzeros(W)) "move costs must be nonnegative"
        @argcheck 1 <= g <= size(W, 1)
    end
    fibre_moves = map(Ws, summaries) do W, π
        best = [Dict{Int,Float64}() for _ in 1:π.n_target]
        for (s, s2, c) in zip(findnz(W)...)
            d = best[π.vertex_map[s]]
            d[s2] = min(get(d, s2, Inf), c)
        end
        [sort!(collect(pairs(d))) for d in best]
    end
    fibres = map(summaries) do π
        f = [Int[] for _ in 1:π.n_target]
        for (s, k) in enumerate(π.vertex_map)
            push!(f[k], s)
        end
        f
    end
    m = NTuple{K,Int}(members)
    P = ProductStalk{K,typeof(contact),typeof(collides)}(m, Ws, [sparse(transpose(W)) for W in Ws],
        collect(GraphHomomorphism, summaries), ntuple(a -> size(Ws[a], 1), K),
        [[[(Int(s), Float64(c)) for (s, c) in v] for v in fm] for fm in fibre_moves], fibres,
        NTuple{K,Int}(goals), contact, collides)
    return P
end

Base.length(::ProductStalk{K}) where {K} = K
Base.show(io::IO, P::ProductStalk{K}) where {K} = print(io, "ProductStalk(agents ", P.members, ")")

_summary(P::ProductStalk, a, m) = m <= P.nstates[a] ? P.summaries[a].vertex_map[m] : m - P.nstates[a]
_exact(P::ProductStalk, a, m) = m <= P.nstates[a]

# Whether agent a is in contact with another agent, given every agent's summary.
_in_contact(P::ProductStalk{K}, ks, a) where {K} =
    any(b != a && P.contact(P.members[a], ks[a], P.members[b], ks[b]) for b in 1:K)

"""
    coarse_state(P::ProductStalk, states)

The coarse state ``μ(states)`` of the joint state `states` (one state per agent,
in the order of `P.members`).
"""
function coarse_state(P::ProductStalk{K}, states) where {K}
    ks = ntuple(a -> P.summaries[a].vertex_map[states[a]], K)
    return ntuple(a -> _in_contact(P, ks, a) ? Int(states[a]) : P.nstates[a] + ks[a], K)
end

"""
    restrict(Pσ::ProductStalk, Pτ::ProductStalk, x)

The restriction of the coarse state `x` of `Pσ` to the face `Pτ` (whose agents are
among those of `Pσ`): forget the other agents, and replace by its summary every
agent that is no longer in contact.
"""
function restrict(Pσ::ProductStalk{K}, Pτ::ProductStalk{L}, x) where {K,L}
    idx = ntuple(b -> findfirst(==(Pτ.members[b]), Pσ.members), L)
    @argcheck all(!isnothing, idx) "the agents of the face must be agents of the simplex"
    ks = ntuple(b -> _summary(Pσ, idx[b], x[idx[b]]), L)
    return ntuple(L) do b
        m = x[idx[b]]
        _exact(Pσ, idx[b], m) && _in_contact(Pτ, ks, b) ? m : Pτ.nstates[b] + ks[b]
    end
end

_moves(P::ProductStalk, a, m) =
    _exact(P, a, m) ? [(Int(rowvals(P.out[a])[k]), nonzeros(P.out[a])[k]) for k in nzrange(P.out[a], m)] :
    P.fibre_moves[a][m - P.nstates[a]]

"""
    stalk_successors(f, P::ProductStalk, x)

Call `f(y, cost)` for every arc `x → y` of the stalk: every joint move out of the
fibre of `x` that does not collide, with its cost, pushed forward. An arc may be
reported more than once; its cost is the least reported.
"""
function stalk_successors(f, P::ProductStalk{K}, x::NTuple{K,Int}) where {K}
    opts = ntuple(a -> _moves(P, a, x[a]), K)
    exact = ntuple(a -> _exact(P, a, x[a]), K)
    for combo in Iterators.product(opts...)
        ok = true
        for a in 1:K, b in (a + 1):K
            if exact[a] && exact[b] &&
               P.collides(P.members[a], x[a], combo[a][1], P.members[b], x[b], combo[b][1])
                ok = false
                break
            end
        end
        ok || continue
        f(coarse_state(P, ntuple(a -> combo[a][1], K)), sum(c[2] for c in combo))
    end
    return nothing
end

"""
    StalkGraph

An explicitly built stalk: `states[i]` is the `i`-th coarse state, `index` its
inverse, `W` the stalk's arcs and `Δ` its cost-to-go to the goal.
"""
struct StalkGraph{K}
    states::Vector{NTuple{K,Int}}
    index::Dict{NTuple{K,Int},Int}
    W::SparseMatrixCSC{Float64,Int}
    Δ::Vector{Float64}
end

# Every coarse state of the stalk: each tuple of summaries, with the agents in
# contact expanded to the states of their fibres (no two of them conflicting).
function _all_states(P::ProductStalk{K}) where {K}
    out = NTuple{K,Int}[]
    for ks in Iterators.product((1:P.summaries[a].n_target for a in 1:K)...)
        touching = ntuple(a -> _in_contact(P, ks, a), K)
        opts = ntuple(a -> touching[a] ? P.fibres[a][ks[a]] : [P.nstates[a] + ks[a]], K)
        any(isempty, opts) && continue
        for x in Iterators.product(opts...)
            ok = true
            for a in 1:K, b in (a + 1):K
                if touching[a] && touching[b] &&
                   P.collides(P.members[a], x[a], x[a], P.members[b], x[b], x[b])
                    ok = false
                    break
                end
            end
            ok && push!(out, x)
        end
    end
    return out
end

"""
    stalk_graph(P::ProductStalk) -> StalkGraph

Build the stalk explicitly: every coarse state, every arc, and the cost-to-go to
the goal. The number of coarse states is about the product of the agents' summary
counts, so this is meant for edges, and for higher simplices only on small maps.
"""
function stalk_graph(P::ProductStalk{K}) where {K}
    states = _all_states(P)
    index = Dict(x => i for (i, x) in enumerate(states))
    rows, cols, vals = Int[], Int[], Float64[]
    for (i, x) in enumerate(states)
        stalk_successors(P, x) do y, c
            j = index[y]
            push!(rows, i); push!(cols, j); push!(vals, c)
        end
    end
    n = length(states)
    W = sparse(rows, cols, vals, n, n, min)
    b = fill(Inf, n)
    b[index[coarse_state(P, P.goal)]] = 0.0
    return StalkGraph{K}(states, index, W, minplus_closure(W, b))
end

"""
    stalk_value(P::ProductStalk, x; heuristic = y -> 0.0, max_expanded = 10^6) -> (value, expanded, exact)

The cost-to-go of the stalk from the coarse state `x` to the goal, by A* on the
stalk built on the fly. `heuristic` must be a consistent lower bound (for instance
the pulled-back values of the faces, see [`section_value`](@ref)). If the search
stops at `max_expanded` expansions, `value` is the least key left in the queue, still
a lower bound, and `exact` is `false`.
"""
function stalk_value(P::ProductStalk{K}, x::NTuple{K,Int}; heuristic = y -> 0.0,
    max_expanded::Integer = 10^6) where {K}
    goal = coarse_state(P, P.goal)
    ids = Dict(x => 1)
    keys_ = [x]
    g = [0.0]
    h = [Float64(heuristic(x))]
    closed = falses(1)
    heap = Tuple{Float64,Int}[]
    isfinite(h[1]) || return Inf, 0, true
    _push!(heap, (h[1], 1))
    expanded = 0
    while !isempty(heap)
        f, i = _pop!(heap)
        closed[i] && continue
        f > g[i] + h[i] + 1e-9 && continue
        keys_[i] == goal && return g[i], expanded, true
        if expanded >= max_expanded
            return f, expanded, false
        end
        closed[i] = true
        expanded += 1
        stalk_successors(P, keys_[i]) do y, c
            j = get(ids, y, 0)
            if j == 0
                push!(keys_, y); push!(g, Inf); push!(h, Float64(heuristic(y))); push!(closed, false)
                j = ids[y] = length(keys_)
            end
            v = g[i] + c
            if v < g[j] && isfinite(h[j])
                g[j] = v
                closed[j] = false
                _push!(heap, (v + h[j], j))
            end
        end
    end
    return Inf, expanded, true
end

# ---------------------------------------------------------------------------
# The sheaf over the complex, and gluing
# ---------------------------------------------------------------------------

"""
    CoverSheaf(complex, graphs, summaries, goals; contact, collides, explicit_dim = 1)

The stalks of every simplex of `complex` (agents are its vertices: `graphs[v]`,
`summaries[v]`, `goals[v]` describe agent `v`; see [`ProductStalk`](@ref) for
`contact` and `collides`). Stalks of dimension up to `explicit_dim` are built
explicitly ([`stalk_graph`](@ref)); higher ones are searched on demand
([`section_value`](@ref)).
"""
struct CoverSheaf
    complex::SimplicialComplex
    stalks::Dict{Vector{Int},ProductStalk}
    explicit::Dict{Vector{Int},StalkGraph}
    explicit_dim::Int
end

function CoverSheaf(K::SimplicialComplex, graphs, summaries, goals; contact, collides, explicit_dim::Integer = 1)
    @argcheck length(graphs) == length(summaries) == length(goals) == nv(K)
    stalks = Dict{Vector{Int},ProductStalk}()
    explicit = Dict{Vector{Int},StalkGraph}()
    for k in 0:dimension(K), σ in simplices(K, k)
        P = ProductStalk(Tuple(σ), graphs[σ], summaries[σ], goals[σ]; contact, collides)
        stalks[σ] = P
        k <= explicit_dim && (explicit[σ] = stalk_graph(P))
    end
    return CoverSheaf(K, stalks, explicit, explicit_dim)
end

# The value of the face τ of σ at the restriction of σ's coarse state x.
function _face_value(C::CoverSheaf, Pσ, τ, x)
    Pτ = C.stalks[τ]
    y = restrict(Pσ, Pτ, x)
    G = get(C.explicit, τ, nothing)
    return isnothing(G) ? first(stalk_value(Pτ, y; heuristic = z -> _face_heuristic(C, Pτ, τ, z))) : G.Δ[G.index[y]]
end

# Pulled-back face values bound the cost of σ: each face pays for its own agents, and
# each agent of a k-simplex lies in k of its k + 1 faces.
function _face_heuristic(C::CoverSheaf, Pσ, σ, x)
    length(σ) <= 1 && return 0.0
    vals = [_face_value(C, Pσ, τ, x) for τ in faces(σ)]
    return max(maximum(vals), sum(vals) / (length(σ) - 1))
end

"""
    section_value(C::CoverSheaf, σ, states; max_expanded = 10^6) -> (value, expanded, exact)

The cost-to-go ``Δ_σ`` of the stalk of the simplex `σ` at the joint state `states`
(one state per vertex of the complex; only those of `σ` are read). Explicit stalks
are looked up; others are searched with the pulled-back face values as heuristic.
"""
function section_value(C::CoverSheaf, σ::AbstractVector{<:Integer}, states; max_expanded::Integer = 10^6)
    σ = sort!(collect(Int, σ))
    P = C.stalks[σ]
    x = coarse_state(P, Tuple(states[v] for v in σ))
    G = get(C.explicit, σ, nothing)
    isnothing(G) || return G.Δ[G.index[x]], 0, true
    return stalk_value(P, x; heuristic = y -> _face_heuristic(C, P, σ, y), max_expanded)
end

"""
    simplicial_cover(n, constraints) -> Float64

The least ``∑_v x_v`` over integers ``x_v ≥ 0`` (``v ∈ 1:n``) with
``∑_{v ∈ σ} x_v ≥ Δ_σ`` for every `(σ, Δ_σ)` in `constraints`. Exhaustive search
with pruning, meant for a few dozen constraints on up to about a dozen vertices.
"""
function simplicial_cover(n::Integer, constraints)
    @argcheck n <= 16 "exhaustive search: at most 16 vertices"
    cons = [(sort!(collect(Int, σ)), ceil(Int, Δ - 1e-9)) for (σ, Δ) in constraints if isfinite(Δ) && Δ > 0]
    @argcheck all(c -> all(v -> 1 <= v <= n, c[1]), cons)
    x = zeros(Int, n)
    best = Ref(sum(c[2] for c in cons; init = 0))
    closing = [[c for c in cons if last(c[1]) == v] for v in 1:n]
    function go(v, total)
        total >= best[] && return
        if v > n
            best[] = total
            return
        end
        for t in 0:maximum((c[2] for c in cons if v in c[1]); init = 0)
            x[v] = t
            all(sum(x[c[1]]) >= c[2] for c in closing[v]) && go(v + 1, total + t)
        end
        x[v] = 0
    end
    go(1, 0)
    return float(best[])
end

"""
    cover_bound(C::CoverSheaf, states; max_expanded = 10^6) -> (bound, values)

The [`simplicial_cover`](@ref) of the stalk values of every simplex of positive
dimension at the joint state `states`: a lower bound on the total cost of the
agents from `states`. `values` maps each simplex to `(value, expanded, exact)`.
"""
function cover_bound(C::CoverSheaf, states; max_expanded::Integer = 10^6)
    values = Dict{Vector{Int},Tuple{Float64,Int,Bool}}()
    for k in 1:dimension(C.complex), σ in simplices(C.complex, k)
        values[σ] = section_value(C, σ, states; max_expanded)
    end
    return simplicial_cover(nv(C.complex), [(σ, v[1]) for (σ, v) in values]), values
end

end
