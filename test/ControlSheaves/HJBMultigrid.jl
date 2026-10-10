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

    # Odd counts: a vertex-centred coarse grid (same ends, spacing 2h).
    odd = HJBProblem(StateGrid([-2.0, -2.0], [2.0, 2.0], [9, 9]); control_bound = 1.0)
    oddc = MGH._coarsen(odd)
    @test oddc.grid.points == [5, 5] && oddc.grid.lower ≈ odd.grid.lower && oddc.grid.upper ≈ odd.grid.upper &&
          oddc.grid.spacing ≈ 2 .* odd.grid.spacing
    # A periodic dimension not divisible by 4 is left alone.
    arm6 = MechanicalHJBProblem(TwoLinkArm(); angle_points = 6, momentum_points = 7)
    @test MGH._coarsen(arm6).grid.points == [6, 6, 4, 4]

    # Automatic depth: O(log n) levels down to coarsest_size points, on grids of
    # any counts; and solves on uneven hierarchies reach the one-level solution.
    @test MGH._plan_levels(Multigrid(), (64, 64, 64, 64), ntuple(_ -> false, 4)) == 4      # 64⁴ → 32⁴ → 16⁴ → 8⁴ = 4096
    @test MGH._plan_levels(Multigrid(coarsest_size = 1), (63, 63), (false, false)) == 6    # 63 → 32 → … → 2
    @test MGH._plan_levels(Multigrid(), (48, 48, 40, 40), (true, true, false, false)) == 4  # → (24, 24, 20, 20) → (12, 12, 10, 10) → (6, 6, 5, 5)
    uneven = [
        ("double integrator 63²", HJBProblem(StateGrid([-2.0, -2.0], [2.0, 2.0], [63, 63]); control_bound = 1.0), 64),
        ("double integrator 11⁴ disc", HJBProblem(StateGrid(fill(-2.0, 4), fill(2.0, 4), fill(11, 4)); control_bound = 1.0,
            constraint = :disc), 100),
        ("arm 12²×13²", MechanicalHJBProblem(TwoLinkArm(); angle_points = 12, momentum_points = 13,
            torque_bound = (6.0, 3.0)), 200)]
    for (label, problem, size) in uneven
        @testset "$label (automatic levels)" begin
            reference = solve(problem, GridPolicyIteration())
            for coarse_operator in (:rediscretize, :galerkin)
                mg = Multigrid(; coarsest_size = size, coarse_operator)
                sol = solve(problem, GridPolicyIteration(multigrid = mg))
                @test sol.converged
                @test sol.values ≈ reference.values rtol = 1e-7
                @test sum(sol.linear_iterations) < sum(reference.linear_iterations)
            end
        end
    end

    # Invalid hierarchies.
    @test_throws ArgumentError solve(prob, GridPolicyIteration(multigrid = Multigrid(levels = 9)))
    @test_throws ArgumentError solve(prob, GridPolicyIteration(multigrid = Multigrid(levels = 1)))
    @test_throws ArgumentError solve(prob, GridPolicyIteration(multigrid = Multigrid(coarsest_size = 100)))   # 8² is already small
    @test_throws ArgumentError solve(prob, GridPolicyIteration(multigrid = Multigrid(levels = 2, coarse_operator = :other)))
    @test_throws ArgumentError solve(prob, GridPolicyIteration(multigrid = Multigrid(levels = 2), coarse_blocks = [2, 2]))
end
