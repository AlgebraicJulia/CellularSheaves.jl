using OrdinaryDiffEqTsit5: Tsit5
using SciMLBase: ODEProblem, init, reinit!, solve, solve!

const _integrator_cache = IdDict{Any,Any}()

function build_ode_problem(state::AbstractArray{Float64}, step::Int, dt::Float64, derivative::Function)
    shape = size(state)
    y0 = reshape(state, :)
    t0 = step * dt
    tspan = (t0, t0 + dt)

    function f!(du, y, p, t)
        du .= reshape(derivative(t, reshape(y, shape)), :)
    end

    return ODEProblem(f!, y0, tspan)
end

function integrate_step!(prob::ODEProblem, state::AbstractArray{Float64}, step::Int, dt::Float64, derivative::Function)
    shape = size(state)
    y0 = reshape(state, :)
    t0 = step * dt
    tf = t0 + dt

    integrator = get!(_integrator_cache, prob) do
        init(prob, Tsit5(); reltol=1e-9, abstol=1e-12)
    end

    reinit!(integrator, y0; t0=t0, tf=tf)
    solve!(integrator)
    return reshape(integrator.u, shape)
end

function integrate_step(state::AbstractArray{Float64}, step::Int, dt::Float64, derivative::Function)
    prob = build_ode_problem(state, step, dt, derivative)
    sol = solve(prob, Tsit5(); reltol=1e-9, abstol=1e-12)
    return reshape(sol.u[end], size(state))
end

function integrate_step(derivative::Function, state::AbstractArray{Float64}, step::Int, dt::Float64)
    return integrate_step(state, step, dt, derivative)
end
