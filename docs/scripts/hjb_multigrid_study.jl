# Multigrid for the HJB policy evaluations: Krylov iterations per evaluation
# against the resolution, for one-level red–black SGS and for V-cycles whose
# coarse operators are rediscretized or Galerkin (T A E), the coarsest grid
# solved by PathPack. Also how far the two coarse operators differ: the
# fraction of coarse points whose coefficients differ by more than 1%, under
# the converged policy. Every evaluation is solved to linear_tol (forcing = 0),
# so the counts measure the preconditioner alone.
#
#   julia -t 16 --project=<env with PATHPACK> docs/scripts/hjb_multigrid_study.jl [cpu|cuda] [small|large]
#
# Writes docs/figures/hjb/multigrid_study_<backend>_<size>.csv.
using CellularSheaves
using CellularSheaves.ControlSheaves.DoubleIntegratorHJB
using CellularSheaves.ControlSheaves.MechanicalHJB
using CellularSheaves.ControlSheaves.DoubleIntegratorHJB: solve
using KernelAbstractions
using LinearAlgebra
using PATHPACK
using Printf
BLAS.set_num_threads(1)                  # small dense work next to spinning Julia threads

const DIH = CellularSheaves.ControlSheaves.DoubleIntegratorHJB
const GS = CellularSheaves.NetworkSheaves.GridSchwarz

backend_name = length(ARGS) >= 1 ? ARGS[1] : "cpu"
size_name = length(ARGS) >= 2 ? ARGS[2] : "small"
backend = if backend_name == "cuda"
    @eval using CUDA
    @eval (println("device: ", CUDA.name(CUDA.device())); flush(stdout))
    @eval CUDA.CUDABackend()
else
    KernelAbstractions.CPU()
end

integrator(n) = HJBProblem(StateGrid(fill(-2.0, 4), fill(2.0, 4), fill(n, 4)); control_bound = 1.0, constraint = :disc)
arm(na, np) = MechanicalHJBProblem(TwoLinkArm(); angle_points = na, momentum_points = np, torque_bound = (6.0, 3.0))

# The automatic hierarchy (Multigrid(levels = :auto)): its number of levels and coarsest grid.
function hierarchy(prob, mg)
    points, periodic = Tuple(DIH._state_grid(prob).points), DIH._periodic(prob)
    levels = DIH._plan_levels(mg, points, periodic)
    n = points
    for _ in 2:levels
        n = coarse_points(n, coarsening_factors(n, periodic))
    end
    return levels, n
end

cases = size_name == "small" ?
    [("double integrator 12^4", integrator(12)), ("double integrator 25^4", integrator(25)),
     ("arm 16^2x12^2", arm(16, 12)), ("arm 24^2x21^2", arm(24, 21))] :
    [("double integrator 33^4", integrator(33)), ("double integrator 49^4", integrator(49)),
     ("double integrator 65^4", integrator(65)),
     ("arm 32^2x25^2", arm(32, 25)), ("arm 48^2x41^2", arm(48, 41)), ("arm 64^2x49^2", arm(64, 49))]

# Fraction of coarse points (level 2) where the rediscretized and Galerkin
# coefficients differ by more than 1% of the diagonal, under policy U.
function coarse_operator_difference(prob, U)
    g = DIH._state_grid(prob)
    n = Tuple(g.points)
    D = length(n)
    layout = GS.BoxLayout(n, ntuple(_ -> 1, D), 0; periodic = DIH._periodic(prob))
    op = GS.box_operator(layout, DIH._UpwindStencil(U, DIH._KernelProblem(prob, ntuple(_ -> 0, D))))
    r = coarsening_factors(n, DIH._periodic(prob))
    nc = coarse_points(n, r)
    galerkin = zeros(nc..., 2D + 1)
    galerkin_coefficients!(galerkin, op)
    Uc = zeros(nc..., size(U, D + 1))
    DIH._restrict_controls_kernel!(KernelAbstractions.CPU())(Uc, U, r, n, Val(D), Val(size(U, D + 1)); ndrange = nc)
    coarse = DIH._coarsen(prob)
    opc = GS.box_operator(coarse_layout(layout), DIH._UpwindStencil(Uc, DIH._KernelProblem(coarse, ntuple(_ -> 0, D))))
    redisc = zeros(nc..., 2D + 1)
    stencil_coefficients!(redisc, opc)
    differs = 0
    for I in CartesianIndices(nc)
        scale = abs(redisc[I, 1])
        differs += any(abs(redisc[I, k] - galerkin[I, k]) > 0.01 * scale for k in 1:(2D + 1))
    end
    return differs / prod(nc)
end

out = joinpath(@__DIR__, "..", "figures", "hjb", "multigrid_study_$(backend_name)_$(size_name).csv")
mkpath(dirname(out))
open(out, "w") do io
    println(io, "problem,nodes,method,levels,coarsest,policy_iterations,inner_total,inner_per_evaluation,seconds,setup_seconds,linear_seconds,max_value_difference,coarse_operator_difference")
    for (label, prob) in cases
        warm = occursin("arm", label) ? arm(8, 6) : integrator(8)
        solve(warm, GridPolicyIteration(backend = backend, multigrid = Multigrid(levels = 2)))       # compile
        levels, coarsest = hierarchy(prob, Multigrid())
        reference, difference = nothing, NaN
        for (method, mg) in (("one level", nothing),
                             ("rediscretize", Multigrid(coarse_operator = :rediscretize)),
                             ("galerkin", Multigrid(coarse_operator = :galerkin)))
            t = @elapsed sol = solve(prob, GridPolicyIteration(backend = backend, maxiter = 100, multigrid = mg))
            if reference === nothing
                reference = sol.values
                U = reshape(permutedims(sol.controls), DIH._state_grid(prob).points..., size(sol.controls, 1))
                difference = coarse_operator_difference(prob, U)
            end
            inner = sum(sol.linear_iterations)
            row = (label, length(sol.values), method, mg === nothing ? 1 : levels,
                mg === nothing ? "" : join(coarsest, "x"), sol.iterations, inner, round(inner / sol.iterations; digits = 1),
                round(t; digits = 2), round(sol.seconds.setup; digits = 2), round(sol.seconds.linear; digits = 2),
                maximum(abs, sol.values .- reference), round(difference; digits = 4))
            println(io, join(row, ","))
            flush(io)
            @printf("%-24s %-13s %d levels (coarsest %-9s): %2d policy its, %5d inner (%6.1f per evaluation), %7.2f s (setup %.2f s)%s\n",
                label, method, row[4], row[5], sol.iterations, inner, inner / sol.iterations, t, sol.seconds.setup,
                sol.converged ? "" : "  NOT CONVERGED")
            flush(stdout)
        end
        @printf("%-24s rediscretized vs Galerkin coarse operator: %.1f%% of level-2 points differ by > 1%%\n",
            label, 100difference)
    end
end
println("wrote ", out)
