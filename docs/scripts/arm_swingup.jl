# Swing-up of a torque-limited two-link arm (MechanicalHJB): solve the HJB on
# an angle² × momentum² grid by GridPolicyIteration, then fly the feedback from
# hanging down and report whether, and when, the arm comes to rest upright.
#
#   julia -t 16 --project=<env> docs/scripts/arm_swingup.jl [angle_points] [momentum_points] [τ̄₁] [τ̄₂] [backend] [forcing]
#
# backend: "cpu" (default; Julia threads) or "cuda" (one GPU, CUDA.jl); forcing: inexact
# policy evaluation (GridPolicyIteration), default 0.1.
using CellularSheaves
using CellularSheaves.ControlSheaves.MechanicalHJB
using CellularSheaves.ControlSheaves.MechanicalHJB: solve
using KernelAbstractions
using Printf

na = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 32
np = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 33
τ̄ = (length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 6.0, length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : 3.0)
backend = if length(ARGS) >= 5 && ARGS[5] == "cuda"
    @eval using CUDA
    @eval (println("device: ", CUDA.name(CUDA.device())); flush(stdout))
    @eval CUDA.CUDABackend()                   # (CUDA was loaded in this statement: a newer world)
else
    KernelAbstractions.CPU()
end

arm = TwoLinkArm()
prob = MechanicalHJBProblem(arm; angle_points = na, momentum_points = np, torque_bound = τ̄)
b1 = (arm.masses[1] * arm.centers[1] + arm.masses[2] * arm.lengths[1]) * arm.gravity
@printf("arm: gravity torque at the shoulder up to %.1f N m; torque bounds (%.1f, %.1f) N m\n", b1 + arm.masses[2] * arm.centers[2] * arm.gravity, τ̄...)
@printf("grid: %d² angles × %d² momenta = %d nodes, momentum bounds ±(%.1f, %.1f)\n", na, np, length(prob.grid),
    prob.grid.upper[3], prob.grid.upper[4])

forcing = length(ARGS) >= 6 ? parse(Float64, ARGS[6]) : 0.1
alg = GridPolicyIteration(backend = backend, maxiter = 200, forcing = forcing)
solve(MechanicalHJBProblem(arm; angle_points = 8, momentum_points = 5, torque_bound = τ̄), GridPolicyIteration(backend = backend))  # compile
t = @elapsed sol = solve(prob, alg)
@printf("solve: %d policy iterations, %d inner, %.2f s (linear %.2f s), converged %s\n", sol.iterations,
    sum(sol.linear_iterations), t, sol.seconds.linear, sol.converged)

wrap(θ) = mod(θ + π, 2π) - π
for start in ([π, 0.0, 0.0, 0.0], [π - 0.3, 0.2, 0.0, 0.0], [π / 2, 0.0, 0.0, 0.0])
    run = closed_loop(sol, start; duration = 15.0, dt = 0.002)
    θ = wrap.(run.states[1:2, :])
    near = [abs(θ[1, i]) < 0.15 && abs(θ[2, i]) < 0.15 && abs(run.states[3, i]) < 1.0 && abs(run.states[4, i]) < 0.5
            for i in axes(θ, 2)]
    # The first time after which the arm stays near upright.
    last_out = findlast(!, near)
    settled = last_out === nothing ? 0.0 : last_out < length(near) ? run.times[last_out + 1] : NaN
    @printf("from θ = (%.2f, %.2f): V = %.3f, simulated cost %.3f, %s, final θ = (%.3f, %.3f), peak |τ| = (%.2f, %.2f)\n",
        start[1], start[2], value_at(sol, start), run.cost,
        isnan(settled) ? "not upright at t = 15 s" : @sprintf("upright to stay from t = %.2f s", settled),
        θ[1, end], θ[2, end], maximum(abs, run.torques[1, :]), maximum(abs, run.torques[2, :]))
end
