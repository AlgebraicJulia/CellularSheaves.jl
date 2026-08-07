"""
    CoordinationBenchmarks

Head-to-head benchmarking of the two ways to drive a fleet onto the harmonic
extension of a pinned cellular sheaf.

Both laws share a fixed point, the harmonic extension ``q^\\star``, but they
are different dynamical systems, and this module measures what that costs.

**Diffusion**, the decentralized law of *Heterogeneous Multi-Agent Multi-Target
Tracking using Cellular Sheaves* ([arXiv:2512.24886](https://arxiv.org/abs/2512.24886),
eq. 4–5), feeds the sheaf disagreement back directly:

```math
u = -k\\, g^{+} \\eta, \\qquad \\eta = \\mathcal{H} q - B p .
```

**Direct** solves the harmonic extension first and then tracks the delivered
reference locally:

```math
u = -k\\, g^{+} (q - q^\\star), \\qquad q^\\star = \\mathcal{H}^{-1} B p .
```

Substituting each into the tracking error gives
``\\dot e_{\\mathrm{Diffusion}} = -k\\mathcal{H}e`` against ``\\dot e_{\\mathrm{Direct}} = -k e``:
Diffusion inherits the sheaf spectrum, Direct does not.

Every numerical kernel here is an existing `CellularSheaves` entry point , 
[`restricted_laplacian_blocks`](@ref) for the blocks, [`harmonic_extension`](@ref)
as the correctness oracle, and the `DistributedSolve` cached tree solve as the
hot path.

See [`benchmark_coordination`](@ref) for the top-level entry point.
"""
module CoordinationBenchmarks

using ..ControlSheaves: Tikhonov
using CellularSheaves.NetworkSheaves.DistributedSolve
using CellularSheaves.NetworkSheaves.EuclideanSheaves:
    EuclideanSheaf, add_sheaf_edge!, restricted_laplacian_blocks, harmonic_extension
using CliqueTrees.Multifrontal
using Graphs
using LinearAlgebra
using Random
using SparseArrays

export CoordinationScenario, coordination_scenario, scenario_family,
    DirectPlan, direct_plan, harmonic_reference, harmonic_reference!,
    diffusion_residual!, diffusion_control!, direct_control!,
    spectral_summary, stable_gain_ceiling,
    SettlingResult, settle_to_formation,
    StepProfile, step_profile, DIFFUSION_STEPS, DIRECT_SETUP_STEPS, DIRECT_STEP_STEPS,
    CoordinationRollout, rollout_coordination, target_trajectory, agent_groups,
    tracking_bandwidth, scenario_family_small,
    CoordinationBenchmark, benchmark_coordination,
    diffusion_round_slots, tree_makespan, tree_statistics,
    plot_coordination, animate_coordination

"""
    plot_coordination(result; metrics = (:settling, :communication), kwargs...)

Multi-panel summary of a [`CoordinationBenchmark`](@ref).

Requires `Plots` to be loaded (the recipe lives in the `CellularSheavesPlots`
package extension). Individual panels are also available directly:

```julia
plot(result, :settling)       # ticks to reach the formation, against kappa
plot(result, :compute)        # total compute for both laws
plot(result, :communication)  # total half-duplex slots for both laws
plot(result, :speedup)        # Direct advantage in both currencies
plot(result, :spectrum)       # lambda_min / lambda_max decomposition
```

A [`SettlingResult`](@ref) plots directly as its own error transient.
"""
function plot_coordination end

"""
    animate_coordination(result; kwargs...)

Reserved for animated coordination rollouts; implemented in the
`CellularSheavesPlots` package extension.
"""
function animate_coordination end

# ===========================================================================
# scenarios
# ===========================================================================

"""
    CoordinationScenario

A pinned tracking problem: an agent graph, a set of sensed targets, and the
sheaf Laplacian blocks both control laws are built from.

`H` and `Bmat` follow the tracking paper's convention, so that `η = H*q - Bmat*p`
and `q* = H \\ (Bmat*p)`.

Construct with [`coordination_scenario`](@ref).

# Fields
- `name`, `label`: identifier and display name.
- `dim`: stalk dimension (agents live in `R^dim`).
- `nagents`, `ntargets`: fleet and target counts.
- `agent_graph`: inter-agent communication graph.
- `sensing`: for each target, the agents that sense it.
- `sheaf`: the underlying [`EuclideanSheaf`](@ref).
- `H`, `Bmat`: the free–free and free–pinned Laplacian blocks.
"""
struct CoordinationScenario
    name::Symbol
    label::String
    dim::Int
    nagents::Int
    ntargets::Int
    agent_graph::SimpleGraph{Int}
    sensing::Vector{Vector{Int}}
    agent_nbrs::Vector{Vector{Int}}
    target_nbrs::Vector{Vector{Int}}
    sheaf::EuclideanSheaf{Float64}
    agent_vertices::Vector{Int}
    target_vertices::Vector{Int}
    H::SparseMatrixCSC{Float64, Int}
    Bmat::SparseMatrixCSC{Float64, Int}
end

Base.show(io::IO, s::CoordinationScenario) = print(io,
    "CoordinationScenario($(s.label): $(s.nagents) agents, $(s.ntargets) targets, stalk $(s.dim))")

function _build(name, label, agent_graph, sensing; dim = 2)
    na, nt = nv(agent_graph), length(sensing)
    idm = Matrix{Float64}(I, dim, dim)
    sheaf = EuclideanSheaf{Float64}(fill(dim, na + nt))
    for e in edges(agent_graph)
        add_sheaf_edge!(sheaf, src(e), dst(e), idm, idm)
    end
    for (k, agents) in enumerate(sensing), i in agents
        add_sheaf_edge!(sheaf, i, na + k, idm, idm)
    end

    agent_vertices = collect(1:na)
    target_vertices = collect((na + 1):(na + nt))
    L_II, L_IB = restricted_laplacian_blocks(sheaf, agent_vertices, target_vertices)

    agent_nbrs = [sort(collect(neighbors(agent_graph, i))) for i in 1:na]
    target_nbrs = [Int[] for _ in 1:na]
    for (k, agents) in enumerate(sensing), i in agents
        push!(target_nbrs[i], k)
    end

    # `restricted_laplacian_blocks` returns (L_F)_FP; the paper's B is its
    # negation, giving η = Hq - Bp and q* = H^{-1}Bp.
    return CoordinationScenario(name, label, dim, na, nt, agent_graph, sensing,
        agent_nbrs, target_nbrs, sheaf, agent_vertices, target_vertices,
        sparse(L_II), sparse(-L_IB))
end

function _random_geometric(n; seed = 7)
    rng = MersenneTwister(seed)
    pts = [(rand(rng), rand(rng)) for _ in 1:n]
    g = SimpleGraph(n)
    r2 = 2log(n) / n
    for i in 1:n, j in (i + 1):n
        sum(abs2, pts[i] .- pts[j]) <= r2 && add_edge!(g, i, j)
    end
    if !is_connected(g)
        order = sortperm(pts; by = first)
        for i in 1:(n - 1)
            add_edge!(g, order[i], order[i + 1])
        end
    end
    return g
end

function _two_ring()
    g = SimpleGraph(12)
    for (i, j) in [(1,2),(2,3),(3,4),(4,5),(5,6),(6,1),
                   (7,8),(8,9),(9,10),(10,11),(11,12),(12,7),(1,7)]
        add_edge!(g, i, j)
    end
    return g
end

"""A path spine with one pendant leaf per spine vertex: treewidth 1, but a very
different degree profile from a chain or a star."""
function _caterpillar(spine)
    g = SimpleGraph(2spine)
    for i in 1:(spine - 1); add_edge!(g, i, i + 1); end
    for i in 1:spine; add_edge!(g, i, spine + i); end
    return g
end

function _two_clique(m)
    ov = max(1, m ÷ 4)
    n = 2m - ov
    g = SimpleGraph(n)
    blockA = collect(1:m)
    blockB = vcat(collect((m - ov + 1):m), collect((m + 1):n))
    for block in (blockA, blockB), i in block, j in block
        i < j && add_edge!(g, i, j)
    end
    return g
end

"""
    coordination_scenario(name; size_parameter = 0, dim = 2, coverage = :sparse)

Build a [`CoordinationScenario`](@ref) from a named topology family.

`name` is one of `:tworing`, `:chain`, `:ring`, `:grid`, `:rgg`, `:star`,
`:twoclique`, or `:expander`. `size_parameter` sets the family's size (side
length for `:grid`, clique size for `:twoclique`, vertex count otherwise).

`coverage` controls how many agents sense a target, which is what sets
``\\lambda_{\\min}`` and therefore how hard the problem is for Diffusion:
- `:sparse` (default) pins two agents at the extremes of the formation.
- `:full` gives every agent its own target edge, the best case for Diffusion.

# Examples
```julia
prob = coordination_scenario(:chain; size_parameter = 32)
prob = coordination_scenario(:expander; size_parameter = 128, coverage = :full)
```
"""
function coordination_scenario(name::Symbol; size_parameter::Int = 0, dim = 2,
        coverage::Symbol = :sparse, seed = 11)
    g, label = if name === :tworing
        _two_ring(), "two-ring (12)"
    elseif name === :chain
        path_graph(size_parameter), "chain"
    elseif name === :ring
        cycle_graph(size_parameter), "ring"
    elseif name === :grid
        Graphs.grid([size_parameter, size_parameter]), "grid"
    elseif name === :rgg
        _random_geometric(size_parameter), "random geometric"
    elseif name === :star
        star_graph(size_parameter), "star"
    elseif name === :grid3d
        Graphs.grid([size_parameter, size_parameter, size_parameter]), "3-D grid"
    elseif name === :torus
        Graphs.grid([size_parameter, size_parameter]; periodic = true), "torus"
    elseif name === :ladder
        ladder_graph(size_parameter), "ladder"
    elseif name === :prism
        circular_ladder_graph(size_parameter), "prism"
    elseif name === :tree
        binary_tree(size_parameter), "binary tree"
    elseif name === :caterpillar
        _caterpillar(size_parameter), "caterpillar"
    elseif name === :wheel
        wheel_graph(size_parameter), "wheel"
    elseif name === :complete
        complete_graph(size_parameter), "complete"
    elseif name === :bipartite
        complete_bipartite_graph(size_parameter, size_parameter), "complete bipartite"
    elseif name === :barbell
        barbell_graph(size_parameter, size_parameter), "barbell"
    elseif name === :lollipop
        lollipop_graph(size_parameter, size_parameter), "lollipop"
    elseif name === :smallworld
        watts_strogatz(size_parameter, 4, 0.15; rng = MersenneTwister(seed)), "small world"
    elseif name === :scalefree
        barabasi_albert(size_parameter, 2; rng = MersenneTwister(seed)), "scale free"
    elseif name === :twoclique
        _two_clique(size_parameter), "two cliques"
    elseif name === :expander
        random_regular_graph(size_parameter, 3; rng = MersenneTwister(seed)), "expander"
    else
        throw(ArgumentError("unknown scenario family: $name"))
    end

    n = nv(g)
    sensing = if coverage === :full
        [collect(1:n)]
    elseif name === :tworing
        [collect(1:6), collect(7:12)]
    elseif name === :star
        [[2]]
    elseif name === :ring
        [[1], [n ÷ 2 + 1]]
    else
        [[1], [n]]
    end
    return _build(name, label, g, sensing; dim)
end

"""
    scenario_family()

The topology families and sizes swept by [`benchmark_coordination`](@ref).
Returns a vector of `(name, sizes)` pairs.
"""
scenario_family() = [
    ## trees, treewidth 1
    (:chain,       [16, 32, 64, 128]),
    (:star,        [16, 32, 64]),
    (:tree,        [4, 5, 6]),
    (:caterpillar, [8, 16, 32]),
    ## sparse, bounded treewidth
    (:ring,        [16, 32, 64, 128]),
    (:ladder,      [8, 16, 32]),
    (:prism,       [8, 16, 32]),
    (:grid,        [4, 6, 8, 10]),
    ## lattice and geometric
    (:torus,       [4, 6, 8]),
    (:grid3d,      [3, 4, 5]),
    (:rgg,         [16, 32, 64, 128]),
    ## small-world and scale-free
    (:smallworld,  [16, 32, 64]),
    (:scalefree,   [16, 32, 64]),
    ## hubs and dense blocks
    (:wheel,       [16, 32, 64]),
    (:bipartite,   [8, 16, 32]),
    (:barbell,     [8, 12, 16]),
    (:lollipop,    [8, 12, 16]),
    (:complete,    [12, 24, 36]),
    (:twoclique,   [8, 16, 32]),
    ## expanders, bounded degree, no small separators
    (:expander,    [32, 64, 128]),
]

# ===========================================================================
# spectra and gain limits
# ===========================================================================

"""
    spectral_summary(scenario) -> (; minimum, maximum, condition)

Extreme eigenvalues and condition number of ``\\mathcal{H}``.

``\\lambda_{\\min}`` is set by how much of the fleet senses a target (the
Dirichlet boundary), ``\\lambda_{\\max}`` by graph degree. Their ratio is what
governs Diffusion's settling time.
"""
function spectral_summary(s::CoordinationScenario)
    ev = eigvals(Symmetric(Matrix(s.H)))
    return (minimum = first(ev), maximum = last(ev), condition = last(ev) / first(ev))
end

"""
    stable_gain_ceiling(method, scenario, dt) -> Float64

Largest gain a zero-order-hold implementation can deploy before the discrete
closed loop diverges.

Diffusion's fastest mode decays at ``k\\lambda_{\\max}(\\mathcal H)``, limiting it to
``2/(\\Delta t\\,\\lambda_{\\max})``; Direct's every mode decays at ``k``, giving
``2/\\Delta t``. This gap is a real cost of Diffusion's `H`-weighted feedback, and it
bounds every gain used in this module.
"""
function stable_gain_ceiling(method::Symbol, s::CoordinationScenario, dt)
    lam = method === :diffusion ? spectral_summary(s).maximum : 1.0
    return 2 / (dt * lam)
end

# ===========================================================================
# the two control laws
# ===========================================================================

"""
    diffusion_residual!(eta, scenario, q, p) -> eta

Diffusion's disagreement ``\\eta_i`` (the paper's eq. 3), accumulated the way an
agent actually would: a sum over its own graph neighbours and sensed targets,
touching no global state.
"""
function diffusion_residual!(eta, s::CoordinationScenario, q, p)
    d = s.dim
    fill!(eta, 0.0)
    @inbounds for i in 1:s.nagents
        base = (i - 1) * d
        for j in s.agent_nbrs[i]
            jbase = (j - 1) * d
            for c in 1:d
                eta[base + c] += q[base + c] - q[jbase + c]
            end
        end
        for k in s.target_nbrs[i]
            for c in 1:d
                eta[base + c] += q[base + c] - p[c, k]
            end
        end
    end
    return eta
end

"""    diffusion_control!(u, eta, k) -> u

Diffusion's input ``u = -k g^{+}\\eta`` (identity ``g^{+}`` for the single-integrator
plants used here).
"""
function diffusion_control!(u, eta, k)
    @inbounds for i in eachindex(u)
        u[i] = -k * eta[i]
    end
    return u
end

"""    direct_control!(u, q, qstar, k) -> u

Direct's input ``u = -k g^{+}(q - q^\\star)``, tracking the delivered reference.
"""
function direct_control!(u, q, qstar, k)
    @inbounds for i in eachindex(u)
        u[i] = -k * (q[i] - qstar[i])
    end
    return u
end

_integrate!(q, u, dt) = (@inbounds for i in eachindex(q); q[i] += dt * u[i]; end; q)

# ===========================================================================
# Direct's harmonic solve
# ===========================================================================

"""
    DirectPlan

Everything Direct computes once per sheaf: the symbolic analysis, the numeric
factorization of ``\\mathcal{H}``, and a preallocated solve workspace.

Target motion changes only the right-hand side, so this is built once and
reused for every subsequent solve. Construct with [`direct_plan`](@ref).
"""
struct DirectPlan{F}
    scenario::CoordinationScenario
    factor::F
    workspace::TreeWorkspace
    rhs::Vector{Float64}
    scratch::Vector{Float64}
    perm::Vector{Int}
    setup_seconds::Float64
end

"""    direct_plan(scenario) -> DirectPlan

Run Direct's one-time setup and record what it cost."""
function direct_plan(s::CoordinationScenario)
    t0 = time_ns()
    factor = cholesky!(ChordalCholesky(copy(s.H)), NoPivot())
    ws = TreeWorkspace(factor.L, Float64)
    setup = (time_ns() - t0) * 1e-9
    n = size(s.H, 1)
    return DirectPlan(s, factor, ws, zeros(n), zeros(n), collect(factor.P.perm), setup)
end

"""
    harmonic_reference!(qstar, plan, p) -> qstar

Solve ``\\mathcal{H} q^\\star = B p`` for the current target positions, reusing
the cached factorization: form the right-hand side, permute, forward- and
backward-substitute over the clique tree, un-permute.
"""
function harmonic_reference!(qstar, plan::DirectPlan, p::AbstractMatrix)
    mul!(plan.rhs, plan.scenario.Bmat, vec(p))
    _permute!(plan.scratch, plan.rhs, plan.perm)         # == P' \ rhs
    tree_forward_ldiv!(plan.factor.L, plan.scratch, plan.workspace)
    tree_backward_ldiv!(plan.factor.L, plan.scratch, plan.workspace)
    _invpermute!(qstar, plan.scratch, plan.perm)         # == P \ scratch
    return qstar
end

## `P' \ v` and `P \ v` each allocate a fresh vector, putting ~2 kB per solve
## into a path whose library kernel is itself allocation-free. These write
## through preallocated buffers; both directions are pinned by tests, because
## reversing one would return a permuted-but-plausible vector that still passes
## a residual check.
function _permute!(dst, src, perm)
    @inbounds @simd for i in eachindex(dst)
        dst[i] = src[perm[i]]
    end
    return dst
end

function _invpermute!(dst, src, perm)
    @inbounds @simd for i in eachindex(src)
        dst[perm[i]] = src[i]
    end
    return dst
end

"""Formation RMS error, computed without allocating a temporary."""
function _formation_error(q, qstar, nagents)
    direct = zero(eltype(q))
    @inbounds @simd for i in eachindex(q)
        d = q[i] - qstar[i]
        direct += d * d
    end
    return sqrt(direct / nagents)
end

"""    harmonic_reference(plan, p) -> Vector

Allocating form of [`harmonic_reference!`](@ref)."""
harmonic_reference(plan::DirectPlan, p::AbstractMatrix) =
    harmonic_reference!(zeros(size(plan.scenario.H, 1)), plan, p)

"""
    harmonic_reference_oracle(scenario, p) -> Vector

The same reference computed through the library's canonical
[`harmonic_extension`](@ref), used to verify the cached path.
"""
function harmonic_reference_oracle(s::CoordinationScenario, p::AbstractMatrix)
    boundary = Dict(s.target_vertices[k] => Vector{Float64}(p[:, k]) for k in 1:s.ntargets)
    z, _ = harmonic_extension(s.sheaf, boundary)
    offsets = [0; cumsum(s.sheaf.vertex_stalks)]
    idx = reduce(vcat, (collect((offsets[v] + 1):offsets[v + 1]) for v in s.agent_vertices))
    return Vector(z[idx])
end

# ===========================================================================
# communication model
# ===========================================================================

"""Critical path of a combining reduce on a rooted forest, sibling-serialized."""
function _gather_makespan(children_of, roots)
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

"""
    tree_makespan(plan) -> Int

Direct's communication cost for one solve, in half-duplex slots: one gather up
the elimination tree and one scatter back down.

Charged under the same slot model as the
Distributed Sheaf Solve guide, one packet per node per slot, cost is
the makespan of a greedy earliest-ready schedule.
"""
function tree_makespan(plan::DirectPlan)
    S = plan.factor.L.S
    roots = [j for j in 1:nv(S.res) if iszero(S.pnt[j])]
    logical = 2 * _gather_makespan(j -> collect(neighbors(S.chd, j)), roots)

    ## `2*cp(T)` counts only messages crossing between supernodes, which
    ## presumes each supernode sits on its own worker. When the tree collapses
    ## to one or two supernodes a single worker holds essentially the whole
    ## factor: that is a centralized solve, not a distributed one, and cp(T) can
    ## even be zero. Charge the centralized pigeonhole bound instead, the
    ## coordinator's radio must receive one packet per agent and send one back.
    return nv(S.res) <= 2 ? max(logical, 2 * plan.scenario.nagents) : logical
end

"""
    tree_statistics(plan) -> (; bag_width, treewidth, fill, offdiagonal_fill, depth, supernodes)

Combinatorial structure of the elimination tree, which is what sets Direct's cost.

`bag_width` is the widest frontal matrix in degrees of freedom; `treewidth` is
that converted back to vertices (`bag_width ÷ stalk - 1`), an upper bound on the
graph's true treewidth under this ordering. `fill` is the **total** number of
stored factor entries, diagonal blocks included; `offdiagonal_fill` counts only
the off-diagonal blocks.

This is a *combinatorial* quantity, independent of the *spectral* quantity `κ`
that governs Diffusion, the two are set by different properties of the
formation and neither predicts the other.
"""
function tree_statistics(plan::DirectPlan)
    S = plan.factor.L.S
    ns = nv(S.res)
    bag_width = maximum(length(collect(neighbors(S.res, j))) +
                        length(collect(neighbors(S.sep, j))) for j in 1:ns)
    depth = zeros(Int, ns)
    for _ in 1:ns, j in 1:ns
        pj = S.pnt[j]
        pj != 0 && (depth[j] = depth[pj] + 1)
    end
    ## `Lval` holds only the off-diagonal blocks. When the elimination tree is a
    ## single supernode, a clique, say, every entry lives in the diagonal
    ## block `Dval` and `Lval` is empty, so counting `Lval` alone would report a
    ## zero-size factor for the densest graph in the study. Count both.
    return (bag_width = bag_width,
        treewidth = max(bag_width ÷ plan.scenario.dim - 1, 0),
        fill = length(plan.factor.L.Lval) + length(plan.factor.L.Dval),
        offdiagonal_fill = length(plan.factor.L.Lval),
        depth = maximum(depth; init = 0) + 1,
        supernodes = ns)
end

"""
    diffusion_round_slots(scenario) -> Int

Diffusion's communication cost for one control tick: a single edge-coloured
neighbourhood exchange, ``\\Delta`` slots. Target sensing is onboard and
charged to neither method.
"""
diffusion_round_slots(s::CoordinationScenario) = Graphs.Δ(s.agent_graph)

# ===========================================================================
# settling
# ===========================================================================

"""
    SettlingResult

What one control law spent to bring the fleet onto the harmonic formation.

# Fields
- `method`: `:diffusion` or `:direct`.
- `gain`: the gain actually deployed.
- `ticks`: control steps taken to reach the tolerance.
- `converged`: whether the tolerance was met.
- `seconds`: total compute, including any one-time setup.
- `slots`: total half-duplex communication slots.
- `peak_command`: largest speed commanded to any agent at any tick.
- `path_length`: integral of that commanded speed.
- `error_history`: formation RMS error per tick.

Ticks and seconds measure computation, not whether a vehicle could fly the
result. Each law runs at its own discrete stability ceiling, and Direct's ceiling
is `kappa` times higher, so it is free to command a more aggressive transient.
`peak_command` is what a saturating actuator would clip and is the honest
counterweight to the tick counts.
"""
struct SettlingResult
    method::Symbol
    gain::Float64
    ticks::Int
    converged::Bool
    seconds::Float64
    slots::Int
    peak_command::Float64
    path_length::Float64
    error_history::Vector{Float64}
end

function _initial_state(s::CoordinationScenario; seed = 20260801, spread = 0.55)
    rng = MersenneTwister(seed)
    d, n = s.dim, s.nagents
    q0 = zeros(d * n)
    for i in 1:n
        theta = 2pi * i / n
        base = (i - 1) * d
        q0[base + 1] = spread * cos(theta) + 0.08randn(rng)
        d >= 2 && (q0[base + 2] = spread * sin(theta) + 0.08randn(rng))
        d >= 3 && (q0[base + 3] = 0.08randn(rng))
    end
    return q0
end

"""Static target positions, spread so the formation has somewhere to go."""
function _target_state(s::CoordinationScenario)
    p = Matrix{Float64}(undef, s.dim, s.ntargets)
    for k in 1:s.ntargets
        phase = 2pi * (k - 1) / max(s.ntargets, 1)
        p[1, k] = 2.2cos(phase)
        s.dim >= 2 && (p[2, k] = 2.2sin(phase))
        s.dim >= 3 && (p[3, k] = 0.0)
    end
    return p
end

"""
    settle_to_formation(method, scenario; tolerance, dt, safety) -> SettlingResult

Cost for one law to actually *reach* the harmonic formation, from a displaced
initial state with the targets held fixed.

This is the metric the two laws should be compared on. Per-tick cost alone
misleads: Diffusion's tick is cheap precisely because it does not solve anything, so
it needs ``O(\\kappa)`` of them, while Direct is at ``q^\\star`` after one solve.
Only the total answers "what did coordination cost".

Each law runs at `safety` times its own [`stable_gain_ceiling`](@ref), which is
the fastest a deployment could actually drive it.
"""
function settle_to_formation(method::Symbol, s::CoordinationScenario;
        tolerance = 1e-3, dt = 0.005, safety = 0.4, max_ticks = 2_000_000,
        repeats = 7, seed = 20260801)
    method in (:diffusion, :direct) || throw(ArgumentError("unknown method: $method"))
    ## Wall-clock is noisy: Direct's cost is dominated by the factorization, which
    ## allocates heavily, and a GC pause landing inside it has been measured to
    ## inflate a single run by over 100x. The trajectory is deterministic, so
    ## only the timing varies, repeat and keep the fastest. Tick counts are
    ## exact and are the metric to trust; timings are indicative.
    best = _settle_once(method, s; tolerance, dt, safety, max_ticks, seed)
    for _ in 2:repeats
        candidate = _settle_once(method, s; tolerance, dt, safety, max_ticks, seed)
        candidate.seconds < best.seconds && (best = candidate)
    end
    return best
end

function _settle_once(method::Symbol, s::CoordinationScenario;
        tolerance, dt, safety, max_ticks, seed)
    k = safety * stable_gain_ceiling(method, s, dt)
    n = size(s.H, 1)

    plan = method === :direct ? direct_plan(s) : nothing
    setup = method === :direct ? plan.setup_seconds : 0.0

    p = _target_state(s)
    oracle_plan = plan === nothing ? direct_plan(s) : plan
    qstar_exact = harmonic_reference(oracle_plan, p)

    q = _initial_state(s; seed)
    u, eta, qstar = zeros(n), zeros(n), zeros(n)
    scale = _formation_error(q, qstar_exact, s.nagents)
    history = sizehint!(Float64[], 4096)

    per_round = diffusion_round_slots(s)
    direct_slots = method === :direct ? tree_makespan(plan) : 0

    ## Pass 1, instrumented, UNTIMED. Establishes the tick count, the error
    ## history and the commanded effort. The convergence check is O(n) per tick
    ## and Diffusion executes orders of magnitude more ticks than Direct, so leaving it
    ## inside the timed region charges Diffusion for the measurement apparatus and
    ## inflates the very ratio being reported (measured at 6-12% of Diffusion's total
    ## even after the check was made allocation-free).
    ticks, converged = 0, false
    peak_command, path_length = 0.0, 0.0
    while ticks < max_ticks
        err = _formation_error(q, qstar_exact, s.nagents)
        push!(history, err)
        if err <= tolerance * scale
            converged = true
            break
        end
        if method === :diffusion
            diffusion_residual!(eta, s, q, p)
            diffusion_control!(u, eta, k)
        else
            ## Fixed targets: the reference is solved once and reused. That
            ## separation of coordination from execution is the whole point.
            ticks == 0 && harmonic_reference!(qstar, plan, p)
            direct_control!(u, q, qstar, k)
        end
        commanded = maximum(abs, u)
        peak_command = max(peak_command, commanded)
        path_length += commanded * dt
        _integrate!(q, u, dt)
        ticks += 1
    end

    ## Pass 2, identical work, exactly `ticks` iterations, nothing in the loop
    ## but the algorithm. This is the number that gets reported.
    copyto!(q, _initial_state(s; seed))
    fill!(u, 0.0)
    fill!(eta, 0.0)
    t_run = time_ns()
    for tick in 0:(ticks - 1)
        if method === :diffusion
            diffusion_residual!(eta, s, q, p)
            diffusion_control!(u, eta, k)
        else
            tick == 0 && harmonic_reference!(qstar, plan, p)
            direct_control!(u, q, qstar, k)
        end
        _integrate!(q, u, dt)
    end
    run_seconds = (time_ns() - t_run) * 1e-9

    slots = method === :diffusion ? ticks * per_round : direct_slots
    return SettlingResult(method, k, ticks, converged, setup + run_seconds,
        slots, peak_command, path_length, history)
end

# ===========================================================================
# per-stage cost accounting
# ===========================================================================

"""Named stages of the Diffusion law, in execution order."""
const DIFFUSION_STEPS = [
    (:exchange,  "exchange neighbour states",    "per tick"),
    (:residual,  "accumulate eta_i",             "per tick"),
    (:control,   "u = -k g+ eta",                "per tick"),
    (:integrate, "integrate plant",              "per tick"),
]

"""One-time setup stages of the Direct law."""
const DIRECT_SETUP_STEPS = [
    (:assemble,  "assemble H, B blocks",           "once per sheaf"),
    (:symbolic,  "symbolic analysis / clique tree", "once per sheaf"),
    (:factor,    "numeric factorization",           "once per sheaf"),
    (:workspace, "allocate solve workspace",        "once per sheaf"),
]

"""Repeated stages of the Direct law."""
const DIRECT_STEP_STEPS = [
    (:rhs,       "form rhs b = Bp",                 "per solve"),
    (:permute,   "apply fill-reducing permutation", "per solve"),
    (:forward,   "forward substitution (tree)",     "per solve"),
    (:backward,  "backward substitution (tree)",    "per solve"),
    (:unpermute, "undo permutation -> q*",          "per solve"),
    (:control,   "u = -k g+ (q - q*)",              "per tick"),
    (:integrate, "integrate plant",                 "per tick"),
]

"""Per-stage cost of one control law over a fixed number of ticks."""
struct StepProfile
    method::Symbol
    scenario::CoordinationScenario
    ticks::Int
    seconds::Dict{Symbol, Float64}
    calls::Dict{Symbol, Int}
    total::Float64
    setup::Float64
end

_bench(op, reps) = (op(); minimum(begin t0 = time_ns(); op(); (time_ns() - t0) * 1e-9 end
                                  for _ in 1:reps))

"""
    step_profile(method, scenario; ticks, dt, safety, reps, solves) -> StepProfile

Time each named stage of one control law in isolation, minimum-of-`reps`, then
scale by how often the mission calls it. This attributes cost to stages; for the
end-to-end number that decides the comparison use [`settle_to_formation`](@ref).
"""
function step_profile(method::Symbol, s::CoordinationScenario;
        ticks = 1000, dt = 0.005, safety = 0.4, reps = 50, solves = ticks)
    method in (:diffusion, :direct) || throw(ArgumentError("unknown method: $method"))
    n = size(s.H, 1)
    k = safety * stable_gain_ceiling(method, s, dt)
    p = _target_state(s)
    q = _initial_state(s)
    qtmp = copy(q)
    u, eta, qstar = zeros(n), zeros(n), zeros(n)
    secs, calls = Dict{Symbol, Float64}(), Dict{Symbol, Int}()

    if method === :diffusion
        secs[:exchange] = 0.0; calls[:exchange] = ticks       # radio, charged in slots
        secs[:residual] = ticks * _bench(() -> diffusion_residual!(eta, s, q, p), reps)
        calls[:residual] = ticks
        secs[:control] = ticks * _bench(() -> diffusion_control!(u, eta, k), reps)
        calls[:control] = ticks
        secs[:integrate] = ticks * _bench(() -> _integrate!(qtmp, u, dt), reps)
        calls[:integrate] = ticks
        setup = 0.0
    else
        Hc = copy(s.H)
        secs[:assemble] = _bench(() -> copy(s.H), reps); calls[:assemble] = 1
        secs[:symbolic] = _bench(() -> ChordalCholesky(copy(Hc)), max(reps ÷ 10, 3))
        calls[:symbolic] = 1
        ## symbolic+numeric measured together, then the symbolic part removed;
        ## two independent minima can cross, so clamp at zero
        secs[:factor] = max(_bench(() -> cholesky!(ChordalCholesky(copy(Hc)), NoPivot()),
                                   max(reps ÷ 10, 3)) - secs[:symbolic], 0.0)
        calls[:factor] = 1
        plan = direct_plan(s)
        secs[:workspace] = _bench(() -> TreeWorkspace(plan.factor.L, Float64),
                                  max(reps ÷ 10, 3))
        calls[:workspace] = 1
        rhs, scratch, perm = plan.rhs, plan.scratch, plan.perm
        secs[:rhs] = solves * _bench(() -> mul!(rhs, s.Bmat, vec(p)), reps)
        calls[:rhs] = solves
        secs[:permute] = solves * _bench(() -> _permute!(scratch, rhs, perm), reps)
        calls[:permute] = solves
        secs[:forward] = solves *
            _bench(() -> tree_forward_ldiv!(plan.factor.L, scratch, plan.workspace), reps)
        calls[:forward] = solves
        secs[:backward] = solves *
            _bench(() -> tree_backward_ldiv!(plan.factor.L, scratch, plan.workspace), reps)
        calls[:backward] = solves
        secs[:unpermute] = solves * _bench(() -> _invpermute!(qstar, scratch, perm), reps)
        calls[:unpermute] = solves
        secs[:control] = ticks * _bench(() -> direct_control!(u, q, qstar, k), reps)
        calls[:control] = ticks
        secs[:integrate] = ticks * _bench(() -> _integrate!(qtmp, u, dt), reps)
        calls[:integrate] = ticks
        setup = sum(secs[sym] for (sym, _, _) in DIRECT_SETUP_STEPS)
    end
    return StepProfile(method, s, ticks, secs, calls, sum(values(secs)), setup)
end

# ===========================================================================
# moving-target rollouts
# ===========================================================================

const TARGET_PERIOD = 12.0

"""
    target_trajectory(scenario, t; period, mode, horizon)

Target positions under one of three task modes: `:static` (formation assembly),
`:orbit` (steady tracking), `:maneuver` (targets jump to new stations half-way
through, a step in the reference). Deterministic in `t`, so both laws see
identical boundary data at identical times.
"""
function target_trajectory(s::CoordinationScenario, t::Real;
        period = TARGET_PERIOD, mode::Symbol = :orbit, horizon = 2period)
    omega = 2pi / period
    phase_t = mode === :static ? 0.0 : t
    jump = (mode === :maneuver && t > horizon / 2) ? pi : 0.0
    p = Matrix{Float64}(undef, s.dim, s.ntargets)
    for k in 1:s.ntargets
        ## equally spaced on a common orbit, so they never collide
        phase = 2pi * (k - 1) / s.ntargets + jump
        p[1, k] = 2.0cos(omega * phase_t + phase)
        s.dim >= 2 && (p[2, k] = 2.0sin(omega * phase_t + phase))
        s.dim >= 3 && (p[3, k] = 0.0)
    end
    return p
end

"""    agent_groups(scenario) -> Vector{Int}

Assign each agent to the target group it escorts, for colouring figures."""
function agent_groups(s::CoordinationScenario)
    groups = fill(1, s.nagents)
    best = fill(typemax(Int), s.nagents)
    for (k, seeds) in enumerate(s.sensing)
        dist = gdistances(s.agent_graph, seeds)
        for i in 1:s.nagents
            if dist[i] < best[i]
                best[i] = dist[i]; groups[i] = k
            end
        end
    end
    return groups
end

"""A recorded closed-loop rollout against moving targets, for visualisation."""
struct CoordinationRollout
    scenario::CoordinationScenario
    method::Symbol
    mode::Symbol
    gain::Float64
    dt::Float64
    epsilon::Float64
    times::Vector{Float64}
    positions::Matrix{Float64}
    references::Matrix{Float64}
    targets::Array{Float64, 3}
    errors::Vector{Float64}
end

"""
    rollout_coordination(method, scenario; dt, period, horizon, mode, epsilon, ...)

Fly the fleet against moving targets under one control law, recording the
trajectory.

`epsilon` applies a first-order command lag `eps*vdot = -v + u` **identically to
both laws**; with `epsilon = 0` neither is filtered. It must never be applied to
one law and not the other, which is why it is a rollout parameter rather than a
property of either controller.
"""
function rollout_coordination(method::Symbol, s::CoordinationScenario;
        dt = 0.02, period = TARGET_PERIOD, horizon = 2period, safety = 0.4,
        solve_every = 1, epsilon = 0.0, mode::Symbol = :orbit, spread = 0.55,
        seed = 20260801)
    method in (:diffusion, :direct) || throw(ArgumentError("unknown method: $method"))
    epsilon >= 0 || throw(ArgumentError("epsilon must be non-negative"))
    k = safety * stable_gain_ceiling(method, s, dt)
    n = size(s.H, 1)
    times = collect(0.0:dt:horizon)
    nsteps = length(times)

    plan = direct_plan(s)
    q = _initial_state(s; seed, spread)
    u, eta, qstar, refstar, v = zeros(n), zeros(n), zeros(n), zeros(n), zeros(n)
    filtered = epsilon > 0

    positions = Matrix{Float64}(undef, n, nsteps)
    references = Matrix{Float64}(undef, n, nsteps)
    targets = Array{Float64, 3}(undef, s.dim, nsteps, s.ntargets)
    errors = Vector{Float64}(undef, nsteps)

    for (step, t) in enumerate(times)
        p = target_trajectory(s, t; period, mode, horizon)
        harmonic_reference!(refstar, plan, p)
        positions[:, step] = q
        references[:, step] = refstar
        targets[:, step, :] = p
        errors[step] = _formation_error(q, refstar, s.nagents)

        if method === :diffusion
            diffusion_residual!(eta, s, q, p); diffusion_control!(u, eta, k)
        else
            (step - 1) % solve_every == 0 && harmonic_reference!(qstar, plan, p)
            direct_control!(u, q, qstar, k)
        end
        if filtered
            @inbounds for i in eachindex(v); v[i] += dt * (u[i] - v[i]) / epsilon; end
            _integrate!(q, v, dt)
        else
            _integrate!(q, u, dt)
        end
    end
    return CoordinationRollout(s, method, mode, k, dt, epsilon, times, positions,
        references, targets, errors)
end

"""
    tracking_bandwidth(method, scenario; dt, safety) -> Float64

Closed-loop bandwidth against a moving reference. Diffusion's slowest mode responds at
`k*lambda_min`; Direct's at `k`. Targets slower than this are the quasi-static
regime, where both laws track tightly.
"""
function tracking_bandwidth(method::Symbol, s::CoordinationScenario;
        dt = 0.02, safety = 0.4)
    k = safety * stable_gain_ceiling(method, s, dt)
    return method === :diffusion ? k * spectral_summary(s).minimum : k
end

"""A short family list for quick runs, one representative size each."""
scenario_family_small() = [
    (:chain, [32]), (:star, [32]), (:ring, [32]), (:grid, [6]),
    (:rgg, [32]), (:complete, [24]), (:twoclique, [16]), (:expander, [64]),
]

# ===========================================================================
# the sweep
# ===========================================================================

"""
    CoordinationBenchmark

Result of [`benchmark_coordination`](@ref): one row per scenario, each holding
the scenario, its spectrum, and both laws' [`SettlingResult`](@ref).

Plot it with `plot(result, :settling)`, see [`plot_coordination`](@ref) for
the available metrics.
"""
struct CoordinationBenchmark
    rows::Vector{<:NamedTuple}
    tolerance::Float64
    dt::Float64
    safety::Float64
end

Base.length(b::CoordinationBenchmark) = length(b.rows)
Base.iterate(b::CoordinationBenchmark, i = 1) =
    i > length(b.rows) ? nothing : (b.rows[i], i + 1)

"""
    benchmark_coordination(; family, tolerance, dt, safety, verify) -> CoordinationBenchmark

Run both control laws to the harmonic formation across a family of topologies
and collect the cost of each.

`family` defaults to [`scenario_family`](@ref). When `verify` is true (the
default) each scenario's cached solve is checked against the library's
canonical [`harmonic_extension`](@ref) before timing.

# Example
```julia
result = benchmark_coordination()
plot(result, :settling)
```
"""
function benchmark_coordination(; family = scenario_family(), tolerance = 1e-3,
        dt = 0.005, safety = 0.4, verify = true, coverages = (:sparse, :full))
    rows = NamedTuple[]
    for (name, sizes) in family, sp in sizes, coverage in coverages
        s = coordination_scenario(name; size_parameter = sp, coverage)
        spec = spectral_summary(s)

        residual = NaN
        if verify
            plan = direct_plan(s)
            p = _target_state(s)
            cached = harmonic_reference(plan, p)
            residual = norm(cached - harmonic_reference_oracle(s, p), Inf)
        end

        diffusion = settle_to_formation(:diffusion, s; tolerance, dt, safety)
        direct = settle_to_formation(:direct, s; tolerance, dt, safety)
        tstats = tree_statistics(direct_plan(s))
        ## A run that exhausted `max_ticks` never reached the formation and its
        ## timing is a floor, not a measurement. Fail loudly rather than let it
        ## enter a table as though it were a real data point.
        (diffusion.converged && direct.converged) || error(
            "$(s.label) n=$(s.nagents) ($coverage): a law did not reach the formation")

        push!(rows, (
            name, label = s.label, size_parameter = sp, coverage,
            agents = s.nagents, dofs = size(s.H, 1),
            lambda_min = spec.minimum, lambda_max = spec.maximum,
            condition = spec.condition,
            diffusion, direct,
            bag_width = tstats.bag_width, treewidth = tstats.treewidth,
            fill = tstats.fill, tree_depth = tstats.depth,
            supernodes = tstats.supernodes,
            peak_command_ratio = direct.peak_command / diffusion.peak_command,
            speedup = diffusion.seconds / direct.seconds,
            slot_ratio = diffusion.slots / max(direct.slots, 1),
            oracle_residual = residual,
        ))
        GC.gc()
    end
    return CoordinationBenchmark(rows, tolerance, dt, safety)
end

end # module CoordinationBenchmarks
