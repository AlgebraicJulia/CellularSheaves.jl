# GPU benchmark for the implicit (matrix-free) HJB policy iteration: the same
# KernelAbstractions kernels as on CPU threads, on one NVIDIA GPU through
# CUDA.jl's CUDABackend.
#
# First checks that the GPU solve equals the CPU solve (values and iteration
# counts), then times n⁴ grids (HJB_SIZES, default 33,41,49,57,65), twice each
# (the first run includes kernel compilation for that size's types).
#
# Prepare the environment on a CPU node (so no GPU sits idle during
# installation and precompilation), e.g. with
#   julia --project=<env> -e 'using Pkg; Pkg.develop(path="<repo>"); Pkg.add(["CUDA", "KernelAbstractions"])'
#   julia --project=<env> -e 'using CUDA; CUDA.set_runtime_version!(v"12.9"; local_toolkit = true)'
# then run on a GPU node with the matching CUDA module loaded:
#   module load cuda/12.9.1 julia/1.12.6
#   julia --project=<env> docs/scripts/hjb_gpu_benchmarks.jl
using CUDA, KernelAbstractions, Printf
using CellularSheaves
using CellularSheaves.ControlSheaves.DoubleIntegratorHJB

CUDA.allowscalar(false)
const BACKEND = CUDABackend()
println(CUDA.name(CUDA.device()), ", ", Threads.nthreads(), " CPU threads for the CPU reference")
flush(stdout)

disc(n) = HJBProblem(StateGrid(fill(-2.0, 4), fill(2.0, 4), fill(n, 4)); control_bound = 1.0, constraint = :disc)

for prob in (HJBProblem(StateGrid([-2.0, -2.0], [2.0, 2.0], [41, 41]); control_bound = 1.0), disc(13))
    cpu = solve(prob, GridPolicyIteration())
    gpu = solve(prob, GridPolicyIteration(backend = BACKEND))
    @printf("check %s: max |V_gpu - V_cpu| / max|V| = %.1e, policy iterations %d vs %d, inner %d vs %d\n",
        Tuple(prob.grid.points), maximum(abs, gpu.values - cpu.values) / maximum(abs, cpu.values),
        gpu.iterations, cpu.iterations, sum(gpu.linear_iterations), sum(cpu.linear_iterations))
    flush(stdout)
end

for n in parse.(Int, split(get(ENV, "HJB_SIZES", "33,41,49,57,65"), ","))
    prob = disc(n)
    for rep in 1:2
        GC.gc()
        CUDA.reclaim()
        t = CUDA.@elapsed sol = solve(prob, GridPolicyIteration(backend = BACKEND))
        s = sol.seconds
        @printf("n=%2d (%9d unknowns) GPU run %d  %7.2f s  [assembly %.2f, linear %.2f, improvement %.2f]  PI %2d  inner %5d  conv %s\n",
            n, n^4, rep, t, s.assembly, s.linear, s.improvement, sol.iterations, sum(sol.linear_iterations), sol.converged)
        flush(stdout)
    end
end
