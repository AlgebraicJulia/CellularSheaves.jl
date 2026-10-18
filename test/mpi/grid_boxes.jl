# Run under MPI by test/mpi/runtests.jl, or by hand:
#   mpiexec -n <ranks> julia --project=<test env> test/mpi/grid_boxes.jl
using MPI
using Test
using CellularSheaves

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
end
