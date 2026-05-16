# using DifferentialEquations
using OrdinaryDiffEq: ODEProblem, solve, Tsit5

function integrate_step(state::AbstractArray{Float64}, step::Int, dt::Float64, derivative::Function)
    shape = size(state)
    y0 = reshape(state, :)
    t0 = step * dt
    tspan = (t0, t0 + dt)

    function f!(du, y, p, t)
        du .= reshape(derivative(t, reshape(y, shape)), :)
    end

    prob = ODEProblem(f!, y0, tspan)
    # sol = solve(prob; reltol=1e-9, abstol=1e-12)
    sol = solve(prob, Tsit5(); reltol=1e-9, abstol=1e-12)
    return reshape(sol.u[end], shape)
end

function integrate_step(derivative::Function, state::AbstractArray{Float64}, step::Int, dt::Float64)
    return integrate_step(state, step, dt, derivative)
end
