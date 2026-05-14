# src/simulation/dynamics.jl

function run_dynamics(::Type{D}, input) where {D <: AbstractDynamics}
    output = similar(input, Float64)
    return run_dynamics!(D, output, input)
end

function run_dynamics!(::Type{D}, output::AbstractVector{Float64}, input::AbstractVector{Float64}) where {D <: AbstractDynamics}
    error("run_dynamics! not implemented for dynamics type $(D)")
end

# Current agent/target dynamics in the Python project are zero base dynamics;
# motion comes from the control law and, for targets, the figure-8 desired velocity.
function run_dynamics!(::Type{NoDynamics}, output::AbstractVector{Float64}, input::AbstractVector{Float64})
    output .= 0.0
    return output
end

function run_dynamics!(::Type{CurrentAgentDynamics}, output::AbstractVector{Float64}, input::AbstractVector{Float64})
    output .= 0.0
    return output
end

function run_dynamics!(::Type{CurrentTargetDynamics}, output::AbstractVector{Float64}, input::AbstractVector{Float64})
    output .= 0.0
    return output
end

function get_initial_conditions!(::Type{NoDynamics}, output::AbstractVector{Float64})
    output .= 0.0
    return output
end

function get_initial_conditions!(::Type{CurrentAgentDynamics}, output::AbstractVector{Float64})
    output .= 0.0
    return output
end

function get_initial_conditions!(::Type{CurrentTargetDynamics}, output::AbstractVector{Float64})
    output .= 0.0
    return output
end

function run_dynamics!(::Type{TrophicDynamics}, output::AbstractVector{Float64}, input::AbstractVector{Float64})
    @assert length(input) == 3 "TrophicDynamics expects state length 3"
    h_pop, p_pop, t_pop = input

    r_h = 0.6
    k_cap = 100.0
    a_hp = 0.02
    a_pt = 0.01
    d_p = 0.3
    d_t = 0.1

    output[1] = r_h * h_pop * (1.0 - h_pop / k_cap) - a_hp * h_pop * p_pop
    output[2] = -d_p * p_pop + a_hp * h_pop * p_pop - a_pt * p_pop * t_pop
    output[3] = -d_t * t_pop + a_pt * p_pop * t_pop
    return output
end

function get_initial_conditions!(::Type{TrophicDynamics}, output::AbstractVector{Float64})
    @assert length(output) == 3 "TrophicDynamics initial conditions need length 3"
    output[1] = 40.0
    output[2] = 9.0
    output[3] = 2.0
    return output
end

function get_initial_conditions(::Type{D}, n::Int) where {D <: AbstractDynamics}
    out = zeros(Float64, n)
    return get_initial_conditions!(D, out)
end

function f8_dynamics(time::Real)
    A = 20.0
    B = 10.0
    a = 1.0
    b = 2.0
    delta = pi / 2.0

    C = 15.0
    c = 1.0
    delta_z = 0.0

    xdot = A * a * cos(a * time + delta)
    ydot = B * b * cos(b * time)
    zdot = C * c * cos(c * time + delta_z)

    return [xdot, ydot, zdot]
end
