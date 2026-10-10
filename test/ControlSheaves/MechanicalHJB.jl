using Test
using CellularSheaves
using CellularSheaves.ControlSheaves.MechanicalHJB
using CellularSheaves.ControlSheaves.MechanicalHJB: solve
using ForwardDiff
using LinearAlgebra
using Random

const MH = CellularSheaves.ControlSheaves.MechanicalHJB
const DIH = CellularSheaves.ControlSheaves.DoubleIntegratorHJB

@testset "MechanicalHJB" begin
    @testset "two-link arm: Hamilton's equations" begin
        arm = TwoLinkArm(masses = (1.2, 0.8), lengths = (1.0, 0.7), damping = (0.0, 0.0))
        rng = Random.MersenneTwister(7)
        for _ in 1:20
            x = [2π * rand(rng) - π, 2π * rand(rng) - π, 6 * randn(rng), 3 * randn(rng)]
            M = mass_matrix(arm, x[1:2])
            @test M ≈ M' && isposdef(M)
            ∇H = ForwardDiff.gradient(y -> hamiltonian(arm, y), x)
            f = hamiltonian_vector_field(arm, x, [0.0, 0.0])
            @test f[1:2] ≈ ∇H[3:4]                         # q̇ = ∂H/∂p = M⁻¹p
            @test f[3:4] ≈ -∇H[1:2]                        # ṗ = -∂H/∂q
            τ = randn(rng, 2)
            @test hamiltonian_vector_field(arm, x, τ) ≈ f + [0, 0, τ[1], τ[2]]
            # Energy balance: dH/dt = ∇H ⋅ f = τ ⋅ q̇ (the motors' power).
            @test dot(∇H, hamiltonian_vector_field(arm, x, τ)) ≈ dot(τ, f[1:2]) atol = 1e-9
        end
        # Friction removes energy: dH/dt = -q̇ᵀβq̇.
        rough = TwoLinkArm(damping = (0.3, 0.2))
        x = [0.4, -1.1, 2.0, -0.5]
        f = hamiltonian_vector_field(rough, x, [0.0, 0.0])
        @test dot(ForwardDiff.gradient(y -> hamiltonian(rough, y), x), f) ≈ -(0.3 * f[1]^2 + 0.2 * f[2]^2)
        # Upright and hanging down are equilibria; upright has the most potential energy.
        for q in ([0.0, 0.0], [π, 0.0], [π, π], [0.0, π])
            @test norm(hamiltonian_vector_field(arm, [q; 0.0; 0.0], [0.0, 0.0])) < 1e-12
        end
        @test potential_energy(arm, [0.0, 0.0]) > potential_energy(arm, [0.3, 0.2]) > potential_energy(arm, [π, 0.0])
        @test_throws ArgumentError TwoLinkArm(masses = (0.0, 1.0))
    end

    @testset "upwind Hamiltonian minimizer with a drift" begin
        rng = Random.MersenneTwister(3)
        for _ in 1:500
            r, ū, g, Dp, Dm = 0.05 + rand(rng), 2rand(rng), 3randn(rng), randn(rng), randn(rng)
            φ(u) = r * u^2 / 2 + max(g + u, 0) * Dp + min(g + u, 0) * Dm
            u = DIH._box_argmin(r, ū, g, Dp, Dm)
            @test -ū <= u <= ū
            @test φ(u) <= minimum(φ, range(-ū, ū; length = 2001)) + 1e-9
        end
    end

    @testset "policy iteration on the torus × momentum box" begin
        arm = TwoLinkArm()
        prob = MechanicalHJBProblem(arm; angle_points = 12, momentum_points = 9, torque_bound = (6.0, 3.0))
        @test prob.grid.points == [12, 12, 9, 9]
        @test prob.grid.lower[1] ≈ -π && prob.grid.lower[1] + 12 * prob.grid.spacing[1] ≈ π
        @test_throws ArgumentError MechanicalHJBProblem(arm; angle_points = 11)
        snapshots = []
        sol = solve(prob, GridPolicyIteration(callback = s -> push!(snapshots, s)))
        @test sol.converged
        @test length(snapshots) == sol.iterations + 1
        @test all(abs.(sol.controls[1, :]) .<= 6.0 + 1e-12) && all(abs.(sol.controls[2, :]) .<= 3.0 + 1e-12)
        @test all(sol.values .>= -1e-12)
        # Policy iteration from the unactuated arm only lowers the value (monotone
        # convergence for M-matrix schemes), and the target is the cheapest state.
        @test all(sol.values .<= snapshots[2].values .+ 1e-9)
        upright = value_at(sol, [0.0, 0.0, 0.0, 0.0])
        @test abs(upright) < 1e-10 && minimum(sol.values) > -1e-10      # balanced upright costs nothing
        @test value_at(sol, [π, 0.0, 0.0, 0.0]) > upright
        # Angles wrap: θ and θ + 2π are the same state.
        x = [0.37, -2.9, 1.3, -0.4]
        @test value_at(sol, x) ≈ value_at(sol, x .+ [2π, -2π, 0, 0])
        @test control_at(sol, x) ≈ control_at(sol, x .+ [2π, 0, 0, 0])
        # Every preconditioner reaches the same discrete solution.
        for pc in (:ras, :none)
            other = solve(prob, GridPolicyIteration(preconditioner = pc))
            @test other.converged
            @test other.values ≈ sol.values rtol = 1e-7
        end
        # So do inexact evaluations (inexact Newton), with fewer Krylov iterations.
        inexact = solve(prob, GridPolicyIteration(forcing = 0.1))
        @test inexact.converged
        @test inexact.values ≈ sol.values rtol = 1e-7
        @test sum(inexact.linear_iterations) < sum(sol.linear_iterations)
        # A second level (aggregation coarse space, dense LU) changes the
        # preconditioner, not the solution.
        for correction in (:multiplicative, :additive), pc in (:red_black_sgs, :ras)
            two = solve(prob, GridPolicyIteration(preconditioner = pc, coarse_blocks = [4, 4, 3, 3],
                coarse_correction = correction))
            @test two.converged
            @test two.values ≈ sol.values rtol = 1e-7
        end
        @test_throws ArgumentError solve(prob, GridPolicyIteration(coarse_blocks = [4, 4, 3, 3], coarse_correction = :other))
        # The discrete equation holds: the improved policy is the policy evaluated,
        # and (A_u V)(x) = ℓ(x, u) at every grid point.
        g = prob.grid
        D = 4
        layout = CellularSheaves.NetworkSheaves.GridSchwarz.BoxLayout(Tuple(g.points), (1, 1, 1, 1), 0;
            overlap = 1, periodic = (true, true, false, false))
        kp = DIH._KernelProblem(prob, (0, 0, 0, 0))
        U = reshape(permutedims(sol.controls), g.points..., 2)
        op = CellularSheaves.NetworkSheaves.GridSchwarz.box_operator(layout, DIH._UpwindStencil(U, kp))
        V = CellularSheaves.NetworkSheaves.GridSchwarz.grid_zeros(op)
        CellularSheaves.NetworkSheaves.GridSchwarz.interior(V, op) .= reshape(sol.values, g.points...)
        CellularSheaves.NetworkSheaves.GridSchwarz.exchange!(CellularSheaves.NetworkSheaves.GridSchwarz.SerialBoxes(), layout, V)
        AV = CellularSheaves.NetworkSheaves.GridSchwarz.grid_zeros(op)
        CellularSheaves.NetworkSheaves.GridSchwarz.apply!(AV, op, V)
        b = CellularSheaves.NetworkSheaves.GridSchwarz.grid_zeros(op)
        DIH._launch!(DIH._grid_rhs_kernel!, op, b, U, kp, op.origin, (0, 0, 0, 0))
        residual = CellularSheaves.NetworkSheaves.GridSchwarz.interior(AV, op) .- CellularSheaves.NetworkSheaves.GridSchwarz.interior(b, op)
        @test maximum(abs, residual) <= 1e-6 * maximum(abs, sol.values)
        Unew = similar(U)
        DIH._launch!(DIH._grid_improve_kernel!, op, Unew, V, kp, op.origin)
        @test maximum(abs, Unew .- U) <= 1e-6
    end

    @testset "closed loop" begin
        arm = TwoLinkArm()
        prob = MechanicalHJBProblem(arm; angle_points = 12, momentum_points = 9)
        sol = solve(prob, GridPolicyIteration())
        run = closed_loop(sol, [0.3, -0.2, 0.0, 0.0]; duration = 2.0, dt = 0.01)
        @test size(run.states) == (4, 201) && size(run.torques) == (2, 201)
        @test all(abs.(run.torques[1, :]) .<= prob.torque_bound[1] + 1e-12)
        @test run.cost > 0 && isfinite(run.cost)
        # Without torque, a frictionless arm keeps its energy (RK4 at dt = 0.01).
        free = MechanicalHJBProblem(arm; angle_points = 12, momentum_points = 9, torque_bound = (0.0, 0.0))
        fsol = solve(free, GridPolicyIteration())
        frun = closed_loop(fsol, [2.0, 0.5, 0.0, 0.0]; duration = 2.0, dt = 0.001)
        H = [hamiltonian(arm, frun.states[:, i]) for i in axes(frun.states, 2)]
        @test maximum(abs, H .- H[1]) <= 1e-6 * abs(H[1])
    end
end
