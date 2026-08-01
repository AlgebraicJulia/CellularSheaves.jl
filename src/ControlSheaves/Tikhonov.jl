using LinearAlgebra: Factorization

"""
    TikhonovFilter(x0; epsilon) -> TikhonovFilter

Stateful realization of the normalized Tikhonov planner

    epsilon * xdot = -x + qstar(t),

where `qstar(t)` is an already-solved harmonic reference. The normalization is
the form used by the singular-perturbation model: the boundary-layer rate is
`1 / epsilon`, independently of the spectrum of the harmonic system.
"""
mutable struct TikhonovFilter{T, V <: AbstractVector{T}}
    x::V
    epsilon::T
end

function TikhonovFilter(x0::AbstractVector{T}; epsilon::Real) where {T <: Real}
    epsilon > 0 || throw(ArgumentError("epsilon must be positive"))
    S = promote_type(float(T), typeof(float(epsilon)))
    return TikhonovFilter(Vector{S}(x0), S(epsilon))
end

"""Return the instantaneous harmonic reference solving `H * qstar = rhs`.

`H` may be an `AbstractMatrix` or a pre-computed `Factorization` (e.g.
`factorize(H)`). Passing a cached factorization avoids repeated factorization
when solving for many right-hand sides in a loop.
"""
function tikhonov_equilibrium(H::Union{AbstractMatrix, Factorization}, rhs::AbstractVector)
    size(H, 1) == size(H, 2) || throw(ArgumentError("H must be square"))
    size(H, 1) == length(rhs) || throw(DimensionMismatch("rhs must have length size(H, 1)"))
    return H \ rhs
end

"""
    tikhonov_reference_rate(H, rhs_rate)

Return `qstar_dot` by solving `H * qstar_dot = rhs_rate`. This identity applies
when the sheaf topology is fixed and target motion changes only the right-hand
side of the harmonic system.

`H` may be an `AbstractMatrix` or a pre-computed `Factorization` (e.g.
`factorize(H)`). Passing a cached factorization avoids repeated factorization
when computing reference rates in a loop.
"""
function tikhonov_reference_rate(H::Union{AbstractMatrix, Factorization}, rhs_rate::AbstractVector)
    return tikhonov_equilibrium(H, rhs_rate)
end

"""
    tikhonov_feedforward_reference(reference, reference_rate, epsilon)

Return the lag-canceling input `reference + epsilon * reference_rate` for

    epsilon * xdot = -x + input(t).

For an exact reference derivative, the error `x - reference` obeys the
homogeneous equation `epsilon * edot = -e`.
"""
function tikhonov_feedforward_reference(
        reference::AbstractVector,
        reference_rate::AbstractVector,
        epsilon::Real,
    )
    epsilon > 0 || throw(ArgumentError("epsilon must be positive"))
    length(reference) == length(reference_rate) ||
        throw(DimensionMismatch("reference and reference_rate must have equal length"))
    return reference .+ epsilon .* reference_rate
end

"""
    tikhonov_dissipation(error, epsilon)

Evaluate `Vdot = -norm(error)^2 / epsilon` for
`V(error) = 1/2 * norm(error)^2` and a stationary harmonic reference.
"""
function tikhonov_dissipation(error::AbstractVector, epsilon::Real)
    epsilon > 0 || throw(ArgumentError("epsilon must be positive"))
    return -sum(abs2, error) / epsilon
end

# Coefficient of q1-q0 in the exact first-order-hold update. The series avoids
# cancellation when the control interval is tiny relative to epsilon.
function _foh_ramp_coefficient(z)
    abs(z) < 1e-4 && return z / 2 - z^2 / 6 + z^3 / 24 - z^4 / 120
    return (z + expm1(-z)) / z
end

"""
    tikhonov_step!(filter, q0, q1, dt)
    tikhonov_step!(filter, reference_at, t, dt)

Advance the normalized Tikhonov planner exactly over one interval under a
first-order hold between harmonic references `q0` and `q1`. The callable form
evaluates `reference_at(t)` and `reference_at(t + dt)`.

Unlike an explicit ODE integrator, this update remains stable for every
`dt > 0` and `epsilon > 0`. As `epsilon` approaches zero, the returned state
approaches `q1` without a numerical singularity.
"""
function tikhonov_step!(
        filter::TikhonovFilter,
        q0::AbstractVector,
        q1::AbstractVector,
        dt::Real,
    )
    dt > 0 || throw(ArgumentError("dt must be positive"))
    length(q0) == length(filter.x) ||
        throw(DimensionMismatch("q0 must have the same length as the filter state"))
    length(q1) == length(filter.x) ||
        throw(DimensionMismatch("q1 must have the same length as the filter state"))

    z = dt / filter.epsilon
    rho = exp(-z)
    alpha = -expm1(-z)
    beta = _foh_ramp_coefficient(z)
    filter.x .= rho .* filter.x .+ alpha .* q0 .+ beta .* (q1 .- q0)
    return filter.x
end

function tikhonov_step!(filter::TikhonovFilter, reference_at, t::Real, dt::Real)
    return tikhonov_step!(filter, reference_at(t), reference_at(t + dt), dt)
end
