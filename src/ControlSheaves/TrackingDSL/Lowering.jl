"""
    TrackingDSLLowering

Lowers a `ResolvedProgram` to the concrete `TrackingProblem` used by
`MultiAgentTracking.build_time_expanded_tracking_sheaf`.

The main entry points are `lower_tracking_program` (returns a `TrackingProblem`
plus a boundary dictionary) and the convenience wrapper that runs the full
parse → validate → resolve → lower pipeline.
"""
module TrackingDSLLowering

using LinearAlgebra
using ..TrackingDSLTerm
using ..TrackingDSLResolver
using ..TrackingDSLValidator: validate_tracking_program

# Import the types we're targeting (... = grandparent = ControlSheaves)
import ...MultiAgentTracking: TrackingEdge, TrackingProblem,
    agent_vertex, target_vertex

export lower_tracking_program, LoweredTrackingProblem

# ─────────────────────────────────────────────────────────────────────────────
# Output type
# ─────────────────────────────────────────────────────────────────────────────

"""
    LoweredTrackingProblem

The result of lowering a `ResolvedProgram`.

Fields:
- `problem`: a `TrackingProblem` compatible with
  `build_time_expanded_tracking_sheaf`.
- `boundary`: `Dict{Int,Vector{Float64}}` keyed by vertex IDs (as returned
  by `agent_vertex` / `target_vertex`).
"""
struct LoweredTrackingProblem
    problem::TrackingProblem
    boundary::Dict{Int,Vector{Float64}}
end

# ─────────────────────────────────────────────────────────────────────────────
# Public API
# ─────────────────────────────────────────────────────────────────────────────

"""
    lower_tracking_program(resolved::ResolvedProgram;
        consensus_weight = 1.0,
        tracking_weight  = 1.0,
        include_target_dynamics = false) -> LoweredTrackingProblem

Lower a `ResolvedProgram` to a `TrackingProblem` and a boundary dictionary.

The result `.problem` is directly usable with
`build_time_expanded_tracking_sheaf(result.problem)`.

Restrictions maps for multiple consensus constraints are combined into a
single shared `consensus_restriction` matrix (first declared wins); use the
per-constraint map from `ResolvedConsensus` if per-edge heterogeneity is
required.

For uniform consensus with a single map pair, the map is used directly.
When multiple distinct maps are present, the first is stored in
`TrackingProblem.consensus_restriction` and applied to all consensus edges.

# Example

```julia
prog = parse_tracking_program(quote
    space X = R^2
    map A : X -> X
    map B : X -> X
    map R_y : X -> X
    agent a1 dynamics (A, B) period dt
    agent a2 dynamics (A, B) period dt
    target t1
    horizon K
    time initial = 0
    time final = K
    times Tall = initial:final
    consensus c1 between (a1, a2) using (R_y, R_y) at Tall
    track tr1 agent a1 target t1 using (A, A) at final
    boundary agent a1 at initial = x0_a1
end)
ctx = Dict{Symbol,Any}(
    :K     => 5,
    :A     => [1.0 0.0; 0.0 1.0],
    :B     => [0.0; 1.0;;],
    :R_y   => [1.0 0.0; 0.0 0.0],
    :dt    => 0.05,
    :x0_a1 => [0.0, 0.0, 0.0],
)
result = lower_tracking_program(prog, ctx)
```

The special time aliases `initial = 0` and `final = K` are resolved before
lowering.
"""
function lower_tracking_program(
    resolved::ResolvedProgram;
    consensus_weight::Float64   = 1.0,
    tracking_weight::Float64    = 1.0,
    include_target_dynamics::Bool = false,
)
    n_agents  = length(resolved.agents)
    n_targets = length(resolved.targets)
    k         = resolved.k

    # All agents must share the same dynamics for TrackingProblem
    isempty(resolved.agents) &&
        throw(TrackingDimensionMismatchError("At least one agent must be declared."))

    Ad = resolved.agents[1].Ad
    Bd = resolved.agents[1].Bd
    nx = resolved.agents[1].nx
    nu = resolved.agents[1].nu

    # Validate all agents have the same dynamics
    for ra in resolved.agents[2:end]
        size(ra.Ad) == size(Ad) && ra.Ad ≈ Ad ||
            throw(TrackingDimensionMismatchError(
                "All agents must share the same Ad matrix for TrackingProblem. " *
                "Agent '$(ra.name)' differs."))
        size(ra.Bd) == size(Bd) && ra.Bd ≈ Bd ||
            throw(TrackingDimensionMismatchError(
                "All agents must share the same Bd matrix for TrackingProblem. " *
                "Agent '$(ra.name)' differs."))
    end

    # Build agent_edges from consensus constraints (unique pairs)
    agent_edges      = Tuple{Int,Int}[]
    consensus_restriction = nothing
    consensus_ts_set = Set{Int}()

    for rc in resolved.consensus
        pair = (rc.agent_indices[1], rc.agent_indices[2])
        push!(agent_edges, pair)
        union!(consensus_ts_set, rc.timesteps)
        if consensus_restriction === nothing
            # Use first consensus map as the shared restriction
            consensus_restriction = rc.restriction_maps[1]
        end
    end
    # Default consensus restriction if none was declared
    if consensus_restriction === nothing
        consensus_restriction = Matrix{Float64}(I, nx + nu, nx + nu)
    end

    consensus_timesteps = sort(collect(consensus_ts_set))

    # Build tracking edges and collect tracking timesteps
    tracking_edges = TrackingEdge[]
    tracking_ts_set = Set{Int}()
    for rt in resolved.tracking
        te = TrackingEdge(
            rt.agent_index, rt.target_index,
            rt.restriction_maps[1], rt.restriction_maps[2],
        )
        push!(tracking_edges, te)
        union!(tracking_ts_set, rt.timesteps)
    end
    tracking_timesteps = sort(collect(tracking_ts_set))

    prob = TrackingProblem(
        n_agents, n_targets, k,
        Ad, Bd,
        agent_edges, tracking_edges,
        consensus_restriction,
        consensus_timesteps, tracking_timesteps,
        include_target_dynamics,
        consensus_weight, tracking_weight,
    )

    # Build boundary dictionary
    boundary = Dict{Int,Vector{Float64}}()
    for rb in resolved.boundaries
        if rb.entity_kind == :agent
            vid = agent_vertex(prob, rb.entity_index, rb.time)
        elseif rb.entity_kind == :target
            vid = target_vertex(prob, rb.entity_index, rb.time)
        else
            throw(TrackingDimensionMismatchError("Unknown entity kind: $(rb.entity_kind)"))
        end
        boundary[vid] = rb.value
    end

    return LoweredTrackingProblem(prob, boundary)
end

"""
    lower_tracking_program(prog::TrackingProgram, ctx::AbstractDict;
        consensus_weight = 1.0,
        tracking_weight  = 1.0,
        include_target_dynamics = false) -> LoweredTrackingProblem

Convenience entry point that runs the full validate → resolve → lower pipeline
in one call.

```julia
prog = parse_tracking_program(quote
    space X = R^2
    map_decl(A, X, X)
    map_decl(B, X, X)
    agent a1 dynamics (A, B) period dt
    agent a2 dynamics (A, B) period dt
    horizon K
    time initial = 0
    time final = K
    times Tall = initial:final
    consensus c1 between (a1, a2) using (A, A) at Tall
end)
ctx = Dict{Symbol,Any}(:K => 5, :A => I(2), :B => reshape([0.0,1.0],2,1), :dt => 0.05)
result = lower_tracking_program(prog, ctx)
```
"""
function lower_tracking_program(
    prog::TrackingProgram,
    ctx::AbstractDict;
    consensus_weight::Float64   = 1.0,
    tracking_weight::Float64    = 1.0,
    include_target_dynamics::Bool = false,
)
    resolved = resolve_tracking_program(validate_tracking_program(prog), ctx)
    return lower_tracking_program(resolved;
        consensus_weight, tracking_weight, include_target_dynamics)
end

end # module TrackingDSLLowering
