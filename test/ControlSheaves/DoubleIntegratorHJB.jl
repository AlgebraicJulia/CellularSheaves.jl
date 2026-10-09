using Test
using CellularSheaves
using CellularSheaves.ControlSheaves.DoubleIntegratorHJB
using CellularSheaves.ControlSheaves.DoubleIntegratorHJB: solve
using CellularSheaves.NetworkSheaves.SchwarzMethods: SchwarzIteration, SchwarzGMRES, MulticolorSweep, ParallelSweep
const SchwarzMethods = CellularSheaves.NetworkSheaves.SchwarzMethods
using LinearAlgebra
using SparseArrays

const DIH = CellularSheaves.ControlSheaves.DoubleIntegratorHJB

# Independent reference for the bounded one-axis problem: the discounted cost
# minimized over piecewise-constant controls |u| ≤ ū on a long horizon, by
# accelerated projected gradient (FISTA) with exact adjoint gradients.
function trajectory_reference(prob, x0; dt = 0.02, horizon = 20.0, iters = 5000)
    N = round(Int, horizon / dt)
    A = [1.0 dt; 0.0 1.0]
    B = [dt^2 / 2, dt]
    w, c, r, ρ, ū = prob.position_weight, prob.velocity_weight, prob.control_weight, prob.discount, prob.control_bound
    disc = [exp(-ρ * k * dt) * dt for k in 0:(N - 1)]
    function cost_and_grad(u, x0)
        xs = Vector{Vector{Float64}}(undef, N)
        x = copy(x0)
        J = 0.0
        for k in 1:N
            xs[k] = x
            J += disc[k] * (w * x[1]^2 + c * x[2]^2 + r * u[k]^2) / 2
            x = A * x + B * u[k]
        end
        g = similar(u)
        λ = zeros(2)
        for k in N:-1:1
            g[k] = disc[k] * r * u[k] + dot(B, λ)
            λ = disc[k] .* [w * xs[k][1], c * xs[k][2]] + A' * λ
        end
        return J, g
    end
    v = randn(N)
    L = 1.0
    for _ in 1:60
        Hv = cost_and_grad(v, zeros(2))[2]
        L = norm(Hv) / norm(v)
        v = Hv / norm(Hv)
    end
    u = zeros(N); y = copy(u); t = 1.0
    for _ in 1:iters
        unew = clamp.(y - cost_and_grad(y, x0)[2] / (1.05L), -ū, ū)
        tnew = (1 + sqrt(1 + 4t^2)) / 2
        y = unew + ((t - 1) / tnew) * (unew - u)
        u, t = unew, tnew
    end
    return cost_and_grad(u, x0)[1]
end

@testset "DoubleIntegratorHJB" begin
    @testset "Riccati reference" begin
        P1 = riccati_value_matrix(1, 1.0, 0.1, 0.1, 0.5)
        A = [0.0 1.0; 0.0 0.0] - 0.25I
        B = [0.0, 1.0]
        residual = Diagonal([1.0, 0.1]) + A' * P1 + P1 * A - P1 * B * B' * P1 / 0.1
        @test norm(residual) < 1e-10
        @test isposdef(Symmetric(P1))
        # The planar problem is two copies of the one-axis problem.
        @test riccati_value_matrix(2, 1.0, 0.1, 0.1, 0.5) ≈ kron(P1, I(2))
    end

    @testset "grid and partition" begin
        g = StateGrid([-1.0, -2.0], [1.0, 2.0], [11, 21])
        @test length(g) == 231 && ndims(g) == 2
        parts = grid_partition(g, [2, 3])
        @test sort(unique(parts)) == 1:6
        doms = grid_subdomains(g, [2, 3]; overlap = 1)
        @test all(k -> k in doms[parts[k]], 1:length(g))
        @test sum(length, doms) > length(g)
        @test_throws ArgumentError HJBProblem(StateGrid([-1.0], [1.0], [5]))
        @test_throws ArgumentError HJBProblem(g; constraint = :disc)
    end

    grid1(n) = StateGrid([-2.0, -2.0], [2.0, 2.0], [n, n])
    inner(g) = [k for (k, I) in enumerate(CartesianIndices(Tuple(g.points)))
                if all(abs(g.lower[j] + (I[j] - 1) * g.spacing[j]) <= 1.0 for j in 1:2)]
    function riccati_error(n)
        prob = HJBProblem(grid1(n))
        sol = solve(prob, PolicyIteration())
        g = prob.grid
        xs = [[g.lower[j] + (I[j] - 1) * g.spacing[j] for j in 1:2] for I in CartesianIndices(Tuple(g.points))]
        exact = vec([riccati_value(prob, x) for x in xs])
        k = inner(g)
        return sol, maximum(abs, sol.values[k] - exact[k]) / maximum(exact[k])
    end

    @testset "unbounded problem converges to the Riccati value" begin
        sol41, err41 = riccati_error(41)
        sol81, err81 = riccati_error(81)
        @test sol81.converged
        @test err81 < 0.05
        @test err81 < 0.75 * err41                      # first-order convergence
        @test sol81.iterations <= 15
    end

    bounded = HJBProblem(grid1(61); control_bound = 1.0)
    sol = solve(bounded, PolicyIteration())

    @testset "bounded one-axis problem" begin
        @test sol.converged
        @test maximum(abs, sol.controls) <= 1.0 + 1e-12
        @test count(u -> abs(u) > 1 - 1e-9, sol.controls) > 0      # the bound is active
        # Policy iteration only improves on the clipped LQR policy.
        A0, b0 = DIH._assemble(bounded, DIH._lqr_controls(bounded))
        @test all(sol.values .<= A0 \ b0 .+ 1e-8)
        # A smaller control set costs more: bounded ≥ unbounded on the same grid.
        free = solve(HJBProblem(grid1(61)), PolicyIteration())
        @test all(sol.values .>= free.values .- 1e-8)
        # Agreement with direct trajectory optimization: first-order upwinding
        # overestimates, by a few percent of the value scale on this grid, and
        # the error shrinks under refinement.
        fine = solve(HJBProblem(grid1(121); control_bound = 1.0), PolicyIteration())
        points = ([1.0, 0.5], [-0.5, 1.0], [0.8, -0.8])
        refs = [trajectory_reference(bounded, x0) for x0 in points]
        scale = maximum(value_at(fine, [s1, s2]) for s1 in (-1.0, 1.0), s2 in (-1.0, 1.0))
        coarse_err = maximum(abs(value_at(sol, x0) - ref) for (x0, ref) in zip(points, refs))
        fine_err = maximum(abs(value_at(fine, x0) - ref) for (x0, ref) in zip(points, refs))
        # Regression: UMFPACK's default threshold pivoting once returned a
        # policy evaluation with residual 8.5 on this grid (values jumped to
        # 1e22 and back); partial pivoting keeps every evaluation accurate.
        @test fine.converged && all(isfinite, fine.values)
        @test maximum(fine.value_changes) < 10
        @test fine_err < 0.7 * coarse_err
        @test fine_err < 0.03 * scale
        @test value_at(sol, [0.0, 0.0]) ≈ sol.values[(length(sol.values) + 1) ÷ 2] atol = 1e-12
        @test norm(control_at(sol, [10.0, 10.0])) <= 1.0 + 1e-12          # outside: clipped LQR
    end

    @testset "Schwarz policy evaluation matches the direct solve" begin
        small = HJBProblem(grid1(41); control_bound = 1.0)
        direct = solve(small, PolicyIteration())
        stationary = solve(small, PolicyIteration(evaluation = SchwarzPolicyEvaluation([4, 4])))
        @test stationary.converged
        @test stationary.values ≈ direct.values rtol = 1e-6
        @test all(>(0), stationary.linear_iterations)
        krylov = solve(small, PolicyIteration(evaluation = SchwarzPolicyEvaluation([4, 4];
            algorithm = SchwarzGMRES(sweep = ParallelSweep(), tol = 1e-10, maxiter = 500))))
        @test krylov.values ≈ direct.values rtol = 1e-6
        # Inexact local solves: the global smoother's symmetric Gauss–Seidel
        # pass inside each subdomain, no factorization.
        sgs = SchwarzMethods.SymmetricGaussSeidelLocalSolve()
        for algorithm in (SchwarzIteration(sweep = MulticolorSweep(), tol = 1e-10, maxiter = 20_000),
                          SchwarzGMRES(sweep = ParallelSweep(), tol = 1e-10, maxiter = 500))
            inexact = solve(small, PolicyIteration(evaluation = SchwarzPolicyEvaluation([4, 4];
                algorithm, local_solver = sgs)))
            @test inexact.converged
            @test inexact.values ≈ direct.values rtol = 1e-6
        end
        # Block symmetric Gauss–Seidel: one block is the global smoother.
        for blocks in ([1, 1], [4, 4])
            blocked = solve(small, PolicyIteration(evaluation = KrylovPolicyEvaluation(
                preconditioner = :block_symmetric_gauss_seidel, blocks = blocks)))
            @test blocked.values ≈ direct.values rtol = 1e-6
        end
        @test_throws ArgumentError solve(small, PolicyIteration(evaluation = KrylovPolicyEvaluation(
            preconditioner = :block_symmetric_gauss_seidel)))
        # The parallel assembly into the fixed stencil matches a direct construction.
        U = DIH._lqr_controls(small)
        A, At, b = DIH._assemble(small, U, DIH._stencil(small.grid))
        @test At == sparse(transpose(A))
        @test all(sum(A; dims = 2) .>= small.discount - 1e-12)          # row diagonal dominance by ρ
        @test CellularSheaves.NetworkSheaves.SchwarzMethods._within_pattern(A, DIH._stencil_pattern(DIH._stencil(small.grid)))
        # The serial Newton baseline: global GMRES with a point smoother.
        for preconditioner in (:gauss_seidel, :symmetric_gauss_seidel, :jacobi)
            serial = solve(small, PolicyIteration(evaluation = KrylovPolicyEvaluation(; preconditioner)))
            @test serial.converged
            @test serial.values ≈ direct.values rtol = 1e-6
        end
        @test_throws ArgumentError solve(small, PolicyIteration(evaluation = KrylovPolicyEvaluation(preconditioner = :ilu)))
        @test all(>=(0), values(stationary.seconds)) && stationary.seconds.setup > 0
    end

    @testset "planar problem" begin
        g4 = StateGrid(fill(-2.0, 4), fill(2.0, 4), fill(9, 4))
        box = solve(HJBProblem(g4; control_bound = 1.0), PolicyIteration())
        disc = solve(HJBProblem(g4; control_bound = 1.0, constraint = :disc), PolicyIteration())
        @test box.converged && disc.converged
        @test maximum(norm, eachcol(disc.controls)) <= 1.0 + 1e-12
        # The disc is inside the box: the disc problem costs more.
        @test all(disc.values .>= box.values .- 1e-8)
        # Point symmetry x ↦ −x of the problem and the grid.
        @test disc.values ≈ reverse(disc.values) rtol = 1e-8
        # With a box bound the planar value is the sum of two one-axis values
        # (up to the boundary data, so compare near the origin).
        one = solve(HJBProblem(StateGrid([-2.0, -2.0], [2.0, 2.0], [9, 9]); control_bound = 1.0), PolicyIteration())
        for x in ([0.5, 0.0, -0.5, 0.0], [0.0, 0.5, 0.0, -0.5])
            @test value_at(box, x) ≈ value_at(one, x[[1, 3]]) + value_at(one, x[[2, 4]]) rtol = 0.1
        end
    end
end
