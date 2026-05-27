# ---------------------------------------------------------------------------
# Result type
# ---------------------------------------------------------------------------

"""
    ScenarioResult

Bundle of trajectory data and statistics for one solved scenario.

Fields:
- `label`: scenario name.
- `times`: discrete time vector (length `k+1`).
- `agent_trajs`: `Vector` of `(k+1) × nx` state matrices, one per agent.
- `target_trajs`: `Vector` of stalk-vector sequences, one per target.
- `null_dim`: nullspace dimension of the restricted Laplacian.
- `residual`: Laplacian energy `sqrt(z' * L * z)`.
- `y_col`, `z_col`: which columns of the state matrix correspond to the y and z
  coordinates used for plotting.
"""
struct ScenarioResult
    label::String
    times::Vector{Float64}
    agent_trajs::Vector{Matrix{Float64}}
    target_trajs::Vector{Vector{Vector{Float64}}}
    null_dim::Int
    residual::Float64
    y_col::Int
    z_col::Int
end

function Base.show(io::IO, sr::ScenarioResult)
    print(io, "ScenarioResult(\"$(sr.label)\": null dim = $(sr.null_dim), ||dz|| = $(round(sr.residual; sigdigits=4)))")
end

function Base.show(io::IO, ::MIME"text/plain", sr::ScenarioResult)
    println(io, "ScenarioResult:")
    println(io, "  label    : $(sr.label)")
    println(io, "  null dim : $(sr.null_dim)")
    println(io, "  ||dz||   : $(round(sr.residual; sigdigits=4))")
    println(io, "  agents   : $(length(sr.agent_trajs))")
    print(  io, "  targets  : $(length(sr.target_trajs))")
end
