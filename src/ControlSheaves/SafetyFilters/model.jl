# The agent model the filter reads, and the error it raises when a barrier cannot reach
# the input at all.

"""
    RelativeDegreeError <: Exception

Thrown when the configuration barrier has relative degree greater than one with respect to
the filtered input, i.e. when ``P g(x) = 0``.

The barrier is a function of configuration only, so
``\\dot{h}_{ij} = 2 (P x_i - P x_j)^\\top P (f_i + g_i u_i)``. If ``P g_i`` vanishes then no
choice of ``u_i`` influences ``\\dot{h}_{ij}`` instantaneously and the row is vacuous. This is
a property of the pair (barrier, dynamics), not of where the filter sits in the pipeline.
Filtering the rotor commands of a `QuadrotorDynamics` agent hits it: the
configuration barrier there has relative degree two vertically and four horizontally.

The braking barrier of [`CBFCLFParams`](@ref) lowers the relative degree by one, which
suffices for any agent whose velocity is directly actuated — a double integrator in
particular — but not for the horizontal axes of a quadrotor, where this error is still
raised. The remaining remedies are a reference governor, filtering the commanded
acceleration of a differentially flat parameterization, or a higher-order barrier.
"""
struct RelativeDegreeError <: Exception
    norm_Pg::Float64
end

function Base.showerror(io::IO, err::RelativeDegreeError)
    print(io, "RelativeDegreeError: the configuration barrier has relative degree > 1 ",
              "with respect to the filtered input (‖P g(x)‖ = ", err.norm_Pg, "). ",
              "Filter a relative-degree-one virtual reference, filter at the acceleration ",
              "level, or use a higher-order barrier.")
end

"""
    laplacian_block_metric(H::AbstractMatrix, agent::Integer, stalk_dim::Integer)

Extract agent `agent`'s diagonal block of the agent interaction map ``H = L_Q + D`` and
return it as the CLF metric ``M``, so that the local Lyapunov function is the sheaf-weighted
norm ``V_i = \\tfrac{1}{2}\\|x_i - x_i^*\\|_H^2`` rather than the Euclidean one.

The block is symmetrized and checked for positive definiteness, since ``V_i`` is a valid
Lyapunov function only when ``M \\succ 0``.
"""
function laplacian_block_metric(H::AbstractMatrix, agent::Integer, stalk_dim::Integer)
    @argcheck agent >= 1
    @argcheck stalk_dim >= 1
    @argcheck size(H, 1) == size(H, 2)
    @argcheck agent * stalk_dim <= size(H, 1)
    idx = ((agent - 1) * stalk_dim + 1):(agent * stalk_dim)
    M = Matrix{Float64}(H[idx, idx])
    M = (M + M') / 2
    @argcheck isposdef(M)
    return M
end

"""
    ControlAffineModel(drift, input, position)

Local model ``\\dot{x} = f(x) + g(x) w`` together with the row selector ``P`` picking the
configuration coordinates on which collision barriers are defined.

`drift` and `input` are callables of the state returning ``f(x)`` and ``g(x)``.

    ControlAffineModel(n)

Single integrator on ``\\mathbb{R}^n``: ``f = 0``, ``g = I``, ``P = I``. This is the model
of a reference governor, where the filtered quantity is a reference velocity rather than a
physical actuator command.

    ControlAffineModel(dyn::AbstractAgentDynamics)

Linear agent model taken directly from [`continuous_matrices`](@ref), with ``P`` given by
[`position_indices`](@ref). No dynamics are re-derived here.
"""
struct ControlAffineModel{F, G, P <: AbstractMatrix{Float64}, V}
    drift::F
    input::G
    position::P
    velocity::V
end

ControlAffineModel(drift, input, position::AbstractMatrix) =
    ControlAffineModel(drift, input, position, nothing)

ControlAffineModel(n::Integer) = ControlAffineModel(SingleIntegrator(n))

function ControlAffineModel(dyn::AbstractAgentDynamics)
    Ac, Bc = continuous_matrices(dyn)
    nx = size(Ac, 1)
    eye = Matrix{Float64}(I, nx, nx)
    P = eye[position_indices(dyn), :]
    # An agent with no rate state gets `nothing` rather than a 0-row selector, so that the
    # barriers needing a velocity report that they cannot be used instead of silently
    # building rows out of empty vectors.
    v_idxs = velocity_indices(dyn)
    Pv = isempty(v_idxs) ? nothing : eye[v_idxs, :]
    linear_drift = let Ac = Matrix{Float64}(Ac)
        x -> Ac * x
    end
    constant_input = let Bc = Matrix{Float64}(Bc)
        x -> Bc
    end
    return ControlAffineModel(linear_drift, constant_input, P, Pv)
end

"""
    double_integrator_model(n)

Model of a double integrator on ``\\mathbb{R}^n``: state ``(p, v)``, input the acceleration,
``\\dot{p} = v``, ``\\dot{v} = u``. Position is not directly actuated, so a configuration
barrier has relative degree two and requires the braking form; see
[`CBFCLFParams`](@ref) and [`RelativeDegreeError`](@ref).

Shorthand for `ControlAffineModel(DoubleIntegrator(n))`. Use [`DoubleIntegrator`](@ref)
directly to add a damping block, a non-identity input map, or to drive an
[`AgentState`](@ref), since that is an `AbstractAgentDynamics` and this is only the
model the filter reads.
"""
double_integrator_model(n::Integer) = ControlAffineModel(DoubleIntegrator(n))

state_dimension(model::ControlAffineModel) = size(model.position, 2)

"""
    lyapunov_metric(dyn, K; Q = I)
    lyapunov_metric(A, B, K; Q = I)

Metric ``P`` for use as [`StabilityTerm`](@ref)'s `metric`, solving

```math
(A - BK)^\\top P + P (A - BK) = -Q .
```

``V = \\tfrac{1}{2} e^\\top P e`` is then a genuine control Lyapunov function for the
closed loop, because ``u = -K e`` is always a feasible witness for the dissipation
inequality: it achieves ``\\dot{V} = -\\tfrac{1}{2} e^\\top Q e``. A rate that witness
certifies at every error is ``c_3 = \\lambda_{\\min}(Q) / \\lambda_{\\max}(P)``, which is what
[`lyapunov_rate`](@ref) returns.

This is required whenever the input map is not of full row rank. The default metric ``M = I``
gives ``L_g V = g^\\top e``, which for an underactuated plant vanishes on a whole subspace —
for a double integrator it is the velocity error alone, and so is identically zero at rest —
leaving the stability row vacuous while it still demands ``-c_3 V``. Passing the stabilizing
gain the agent already uses, typically its LQR gain, removes that degeneracy and makes the
nominal command itself a witness, so the filter deviates from it only when safety requires.
"""
function lyapunov_metric(A::AbstractMatrix, B::AbstractMatrix, K::AbstractMatrix;
                         Q = I)
    @argcheck size(A, 1) == size(A, 2)
    @argcheck size(B, 1) == size(A, 1)
    @argcheck size(K) == (size(B, 2), size(A, 1))
    n = size(A, 1)
    weight = Q isa UniformScaling ? Matrix{Float64}(Q, n, n) : Matrix{Float64}(Q)
    @argcheck isposdef(weight)
    closed = Matrix{Float64}(A) - Matrix{Float64}(B) * Matrix{Float64}(K)
    maximum(real, eigvals(closed)) < 0 ||
        throw(ArgumentError("A - BK must be Hurwitz for the Lyapunov equation to have a " *
                            "positive definite solution; supply a stabilizing gain"))
    P = lyap(transpose(closed), weight)
    P = (P + transpose(P)) / 2
    @argcheck isposdef(P)
    return P
end

function lyapunov_metric(dyn::AbstractAgentDynamics, K::AbstractMatrix; Q = I)
    A, B = continuous_matrices(dyn)
    return lyapunov_metric(A, B, K; Q = Q)
end

"""
    lyapunov_rate(P, Q)

A decay rate ``c_3 = \\lambda_{\\min}(Q)/\\lambda_{\\max}(P)`` certified by the witness
``u = -Ke`` for the metric returned by [`lyapunov_metric`](@ref). It is a bound rather than
the exact best rate, which is the smallest generalized eigenvalue of ``(Q, P)``. Asking
[`StabilityTerm`](@ref) for a faster rate than this is legitimate, but the witness no longer
guarantees feasibility and the relaxation may become active.
"""
function lyapunov_rate(P::AbstractMatrix, Q = I)
    n = size(P, 1)
    weight = Q isa UniformScaling ? Matrix{Float64}(Q, n, n) : Matrix{Float64}(Q)
    return minimum(eigvals(Symmetric(weight))) / maximum(eigvals(Symmetric(P)))
end
