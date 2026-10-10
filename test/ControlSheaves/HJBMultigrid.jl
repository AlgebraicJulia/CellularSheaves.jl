using Test
using CellularSheaves
using CellularSheaves.ControlSheaves.DoubleIntegratorHJB
using CellularSheaves.ControlSheaves.MechanicalHJB
using CellularSheaves.ControlSheaves.DoubleIntegratorHJB: solve
using PATHPACK

const MGH = CellularSheaves.ControlSheaves.DoubleIntegratorHJB

@testset "HJB multigrid" begin
    problems = [
        ("double integrator 64²", HJBProblem(StateGrid([-2.0, -2.0], [2.0, 2.0], [64, 64]); control_bound = 1.0), 4),
        ("double integrator 12⁴ disc", HJBProblem(StateGrid(fill(-2.0, 4), fill(2.0, 4), fill(12, 4)); control_bound = 1.0,
            constraint = :disc), 3),
        ("arm 16²×12²", MechanicalHJBProblem(TwoLinkArm(); angle_points = 16, momentum_points = 12,
            torque_bound = (6.0, 3.0)), 3)]
    for (label, prob, levels) in problems
        @testset "$label" begin
            one_level = solve(prob, GridPolicyIteration())
            @test one_level.converged
            for coarse_operator in (:rediscretize, :galerkin), coarsest in (:pathpack, :dense)
                mg = Multigrid(; levels, coarse_operator, coarsest)
                sol = solve(prob, GridPolicyIteration(multigrid = mg))
                @test sol.converged
                @test sol.values ≈ one_level.values rtol = 1e-7
                @test sum(sol.linear_iterations) < sum(one_level.linear_iterations)
            end
        end
    end

    # The coarse problems: cell-centred grids, the same physics.
    prob = HJBProblem(StateGrid([-2.0, -2.0], [2.0, 2.0], [8, 8]); control_bound = 1.0)
    coarse = MGH._coarsen(prob)
    g, gc = prob.grid, coarse.grid
    @test gc.points == [4, 4] && gc.spacing ≈ 2 .* g.spacing && gc.lower ≈ g.lower .+ g.spacing ./ 2
    @test coarse.riccati == prob.riccati && coarse.control_bound == prob.control_bound
    arm = MechanicalHJBProblem(TwoLinkArm(); angle_points = 8, momentum_points = 6)
    armc = MGH._coarsen(arm)
    @test armc.grid.points == [4, 4, 3, 3] && armc.grid.lower[1] + 4 * armc.grid.spacing[1] ≈ π + arm.grid.spacing[1] / 2

    # Invalid hierarchies.
    odd = HJBProblem(StateGrid([-2.0, -2.0], [2.0, 2.0], [9, 9]); control_bound = 1.0)
    @test_throws ArgumentError solve(odd, GridPolicyIteration(multigrid = Multigrid(levels = 2)))
    @test_throws ArgumentError solve(prob, GridPolicyIteration(multigrid = Multigrid(levels = 2, coarse_operator = :other)))
    @test_throws ArgumentError solve(prob, GridPolicyIteration(multigrid = Multigrid(levels = 2), coarse_blocks = [2, 2]))
end
