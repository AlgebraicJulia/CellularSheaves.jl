# Warm-up run for sysimage/build.jl: exercises the dependencies' common call
# paths so their compiled code goes into the image. CellularSheaves itself is
# not in the image; the code compiled here for its own types is discarded, but
# the dependencies' methods for standard types (sparse matrices, Float64
# vectors, KernelAbstractions launches, Krylov solves, graphs) are kept.
# PackageCompiler runs this file with the image project only; CellularSheaves
# (not an image package) comes from the full environment, stacked after it.
push!(LOAD_PATH, joinpath(get(ENV, "CS_SYSIMAGE_DIR", joinpath(homedir(), "sysimage")), "env"))
using CellularSheaves
using CellularSheaves.ControlSheaves.DoubleIntegratorHJB
using LinearAlgebra, SparseArrays, Statistics, Random, Test
using KernelAbstractions
using MPI          # loaded, not initialized: no launcher at build time
using CUDA         # loaded; no GPU on the build node

prob = HJBProblem(StateGrid(fill(-2.0, 4), fill(2.0, 4), fill(7, 4)); control_bound = 1.0, constraint = :disc)
solve(prob, GridPolicyIteration())
solve(prob, GridPolicyIteration(preconditioner = :ras))
solve(prob, PolicyIteration(evaluation = KrylovPolicyEvaluation(method = :bicgstab)))
solve(prob, PolicyIteration(evaluation = KrylovPolicyEvaluation()))
solve(prob, PolicyIteration(evaluation = SchwarzPolicyEvaluation([2, 2, 2, 2];
    local_solver = SymmetricGaussSeidelLocalSolve())))
two = HJBProblem(StateGrid([-2.0, -2.0], [2.0, 2.0], [21, 21]); control_bound = 1.0)
solve(two, PolicyIteration())

# Sparse linear algebra the Schwarz methods use.
n = 20
T = spdiagm(-1 => fill(-1.0, n - 1), 0 => fill(2.0, n), 1 => fill(-1.0, n - 1))
A = kron(sparse(1.0I, n, n), T) + kron(T, sparse(1.0I, n, n))
b = ones(n * n)
lu(A) \ b
@test A * (cholesky(Symmetric(A)) \ b) ≈ b
