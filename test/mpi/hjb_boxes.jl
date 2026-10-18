# Run under MPI by test/mpi/runtests.jl, or by hand:
#   mpiexec -n <ranks> julia --project=<test env> test/mpi/hjb_boxes.jl
using MPI
using Test
using CellularSheaves
using CellularSheaves.ControlSheaves.DoubleIntegratorHJB
using CellularSheaves.ControlSheaves.MechanicalHJB

MPI.Init()
comm = MPI.COMM_WORLD
boxes = mpi_boxes(comm)
nranks = box_count(boxes)

@testset "rank $(box_rank(boxes)) of $nranks" begin
    # Distributed policy iteration: the global red–black preconditioner is the
    # serial operator; RAS converges to the same discrete solution.
    problems = (HJBProblem(StateGrid([-2.0, -2.0], [2.0, 2.0], [41, 41]); control_bound = 1.0),
                HJBProblem(StateGrid(fill(-2.0, 4), fill(2.0, 4), fill(9, 4)); control_bound = 1.0, constraint = DiscConstraint()))
    for prob in problems
        serial = solve(prob, GridPolicyIteration())
        global_rb = solve(prob, GridPolicyIteration(communicator = boxes, smoother = SymmetricGaussSeidelSmoother()))
        @test global_rb.converged
        @test global_rb.values ≈ serial.values rtol = 1e-9
        @test global_rb.linear_iterations == serial.linear_iterations || nranks == 1 ||
              maximum(abs, global_rb.linear_iterations .- serial.linear_iterations) <= 1
        ras = solve(prob, GridPolicyIteration(communicator = boxes, smoother = RestrictedAdditiveSchwarzSmoother()))
        @test ras.converged
        @test ras.values ≈ serial.values rtol = 1e-7
        @test count(abs.(ras.controls - serial.controls) .> 1e-6) <= length(serial.values) ÷ 100
    end

    # The two-link arm: periodic joint angles, reflecting momentum bounds.
    arm = MechanicalHJBProblem(TwoLinkArm(); angle_points = 8, momentum_points = 5)
    serial = solve(arm, GridPolicyIteration())
    distributed = solve(arm, GridPolicyIteration(communicator = boxes))
    @test distributed.converged
    @test distributed.values ≈ serial.values rtol = 1e-9
    ras = solve(arm, GridPolicyIteration(communicator = boxes, smoother = RestrictedAdditiveSchwarzSmoother()))
    @test ras.converged
    @test ras.values ≈ serial.values rtol = 1e-7
    # Two levels: aggregates that do not line up with the boxes.
    for smoother in (SymmetricGaussSeidelSmoother(), RestrictedAdditiveSchwarzSmoother())
        two = solve(arm, GridPolicyIteration(; communicator = boxes, smoother, coarse = AggregateCoarseCorrection([3, 3, 2, 2])))
        @test two.converged
        @test two.values ≈ serial.values rtol = 1e-7
    end
end
