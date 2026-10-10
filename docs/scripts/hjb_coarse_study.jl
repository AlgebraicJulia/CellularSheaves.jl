# Does an aggregation coarse space (GridPolicyIteration(coarse_blocks = ...),
# dense LU on the coarse grid) cut the Krylov iterations of the policy
# evaluations, and do they still grow with the resolution? One-level red–black
# SGS against two levels, for the 4-D double integrator (disc thrust) and the
# two-link arm, every evaluation solved to linear_tol (forcing = 0) so the
# counts measure the preconditioner alone.
#
#   julia -t 16 --project=<env> docs/scripts/hjb_coarse_study.jl [cpu|cuda] [small|large]
#
# Writes docs/figures/hjb/coarse_study_<backend>_<size>.csv.
using CellularSheaves
using CellularSheaves.ControlSheaves.DoubleIntegratorHJB
using CellularSheaves.ControlSheaves.MechanicalHJB
using CellularSheaves.ControlSheaves.DoubleIntegratorHJB: solve
using KernelAbstractions
using Printf
using LinearAlgebra
# One BLAS thread: the dense coarse LU is small, and BLAS threads competing with
# spinning Julia threads (JULIA_THREAD_SLEEP_THRESHOLD=infinite) slow it down by orders of magnitude.
BLAS.set_num_threads(1)

backend_name = length(ARGS) >= 1 ? ARGS[1] : "cpu"
size_name = length(ARGS) >= 2 ? ARGS[2] : "small"
backend = if backend_name == "cuda"
    @eval using CUDA
    @eval (println("device: ", CUDA.name(CUDA.device())); flush(stdout))
    @eval CUDA.CUDABackend()                   # (CUDA was loaded in this statement: a newer world)
else
    KernelAbstractions.CPU()
end

integrator(n) = HJBProblem(StateGrid(fill(-2.0, 4), fill(2.0, 4), fill(n, 4)); control_bound = 1.0, constraint = :disc)
arm(na, np) = MechanicalHJBProblem(TwoLinkArm(); angle_points = na, momentum_points = np, torque_bound = (6.0, 3.0))

# (label, problem, list of coarse_blocks to try; nothing = one level)
cases = if size_name == "small"
    [("double integrator 17^4", integrator(17), [nothing, [3, 3, 3, 3], [5, 5, 5, 5]]),
     ("double integrator 25^4", integrator(25), [nothing, [3, 3, 3, 3], [5, 5, 5, 5]]),
     ("arm 16^2x13^2", arm(16, 13), [nothing, [4, 4, 3, 3], [8, 8, 4, 4]]),
     ("arm 24^2x21^2", arm(24, 21), [nothing, [4, 4, 3, 3], [8, 8, 4, 4]])]
else
    [("double integrator 33^4", integrator(33), [nothing, [3, 3, 3, 3], [5, 5, 5, 5], [7, 7, 7, 7]]),
     ("double integrator 49^4", integrator(49), [nothing, [5, 5, 5, 5], [7, 7, 7, 7]]),
     ("arm 32^2x25^2", arm(32, 25), [nothing, [4, 4, 3, 3], [8, 8, 4, 4], [8, 8, 6, 6]]),
     ("arm 48^2x41^2", arm(48, 41), [nothing, [8, 8, 4, 4], [8, 8, 6, 6]])]
end

out = joinpath(@__DIR__, "..", "figures", "hjb", "coarse_study_$(backend_name)_$(size_name).csv")
mkpath(dirname(out))
open(out, "w") do io
    println(io, "problem,nodes,coarse_blocks,aggregates,policy_iterations,inner_total,inner_per_evaluation,seconds,setup_seconds,linear_seconds,max_value_difference")
    for (label, prob, choices) in cases
        warm = label |> contains("arm") ? arm(8, 5) : integrator(5)
        solve(warm, GridPolicyIteration(backend = backend, coarse_blocks = [2, 2, 2, 2]))   # compile
        reference = nothing
        for blocks in choices
            alg = GridPolicyIteration(backend = backend, maxiter = 100, coarse_blocks = blocks)
            t = @elapsed sol = solve(prob, alg)
            reference === nothing && (reference = sol.values)
            agg = blocks === nothing ? 0 : prod(blocks)
            inner = sum(sol.linear_iterations)
            row = (label, length(sol.values), blocks === nothing ? "none" : join(blocks, "x"), agg, sol.iterations, inner,
                round(inner / sol.iterations; digits = 1), round(t; digits = 2), round(sol.seconds.setup; digits = 2),
                round(sol.seconds.linear; digits = 2), maximum(abs, sol.values .- reference))
            println(io, join(row, ","))
            flush(io)
            @printf("%-24s coarse %-9s %5d aggregates: %2d policy its, %5d inner (%6.1f per evaluation), %7.2f s (coarse setup %.2f s)%s\n",
                label, row[3], agg, sol.iterations, inner, inner / sol.iterations, t, sol.seconds.setup,
                sol.converged ? "" : "  NOT CONVERGED")
            flush(stdout)
        end
    end
end
println("wrote ", out)
