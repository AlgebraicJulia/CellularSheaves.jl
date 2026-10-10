# Run under MPI (see test/network_sheaves/GridSchwarz.jl):
#   mpiexec -n <ranks> julia --project=<test env> test/mpi/grid_boxes.jl
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
    # The halo exchange delivers the neighbours' owned values, edges and
    # corners included, wrapping around periodic dimensions: fill the owned
    # box with the global linear index.
    for (points, periodic) in (((13, 11), (false, false)), ((7, 6, 5, 4), (false, false, false, false)),
                               ((12, 11), (true, false)), ((6, 6, 5, 4), (true, true, false, false)))
        D = length(points)
        ranks = balanced_ranks(nranks, points)
        layout = BoxLayout(points, ranks, box_rank(boxes); overlap = 1, periodic)
        op = box_operator(layout, CoefficientStencil(zeros(length.(layout.owned)..., 2D + 1)))
        x = grid_zeros(op)
        L = LinearIndices(points)
        interior(x, op) .= [L[CartesianIndex(Tuple(I) .+ first.(layout.owned) .- 1)] for I in CartesianIndices(length.(layout.owned))]
        exchange!(boxes, layout, x)
        w = layout.width
        for J in CartesianIndices(x)
            G = Tuple(J) .- w .+ first.(layout.owned) .- 1          # global index of padded cell J
            G = ntuple(k -> periodic[k] ? mod1(G[k], points[k]) : G[k], D)
            inside_owned = all(first.(layout.owned) .<= G .<= last.(layout.owned))
            in_grid = all(1 .<= G .<= points)
            if inside_owned || in_grid
                @test x[J] == L[CartesianIndex(G)]
            end
        end
        @test gather_boxes(boxes, layout, x, op) == Float64.(L)
    end

    # Distributed policy iteration: the global red–black preconditioner is the
    # serial operator; RAS converges to the same discrete solution.
    problems = (HJBProblem(StateGrid([-2.0, -2.0], [2.0, 2.0], [41, 41]); control_bound = 1.0),
                HJBProblem(StateGrid(fill(-2.0, 4), fill(2.0, 4), fill(9, 4)); control_bound = 1.0, constraint = :disc))
    for prob in problems
        serial = solve(prob, GridPolicyIteration())
        global_rb = solve(prob, GridPolicyIteration(communicator = boxes, preconditioner = :red_black_sgs))
        @test global_rb.converged
        @test global_rb.values ≈ serial.values rtol = 1e-9
        @test global_rb.linear_iterations == serial.linear_iterations || nranks == 1 ||
              maximum(abs, global_rb.linear_iterations .- serial.linear_iterations) <= 1
        ras = solve(prob, GridPolicyIteration(communicator = boxes, preconditioner = :ras))
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
    ras = solve(arm, GridPolicyIteration(communicator = boxes, preconditioner = :ras))
    @test ras.converged
    @test ras.values ≈ serial.values rtol = 1e-7
    # Two levels: aggregates that do not line up with the boxes.
    for pc in (:red_black_sgs, :ras)
        two = solve(arm, GridPolicyIteration(communicator = boxes, preconditioner = pc, coarse_blocks = [3, 3, 2, 2]))
        @test two.converged
        @test two.values ≈ serial.values rtol = 1e-7
    end
end
