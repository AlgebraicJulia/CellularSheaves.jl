using Test
using CellularSheaves
using LinearAlgebra
using SparseArrays
using Random

# The dense matrix of the operator of `op` on the single box `layout` (periodic
# dimensions wrap through the halo exchange), column by column.
function dense_operator(op, layout)
    n = size(op)
    N = prod(n)
    A = zeros(N, N)
    e, y = grid_zeros(op), grid_zeros(op)
    for k in 1:N
        fill!(e, 0)
        interior(e, op)[k] = 1
        exchange!(SerialBoxes(), layout, e)
        apply!(y, op, e)
        A[:, k] = vec(interior(y, op))
    end
    return A
end

# The constant sheaf ℝ on the grid graph of `points` (wrapping periodic dimensions).
function constant_sheaf(points, periodic)
    D = length(points)
    L = LinearIndices(points)
    s = EuclideanSheaf{Float64}(ones(Int, prod(points)))
    for I in CartesianIndices(points), j in 1:D
        k = I[j] + 1
        if k > points[j]
            (periodic[j] && points[j] > 2) || continue
            k = 1
        end
        add_sheaf_edge!(s, L[I], L[Base.setindex(Tuple(I), k, j)...], ones(1, 1), ones(1, 1))
    end
    return s
end

@testset "GridMultigrid" begin
    rng = Random.MersenneTwister(5)
    # Even counts (blocks of 2^D), odd counts (a last block of one point), a
    # periodic dimension not divisible by 4 (left alone: semicoarsening), and a
    # non-periodic dimension of fewer than five points (left alone: a coarse
    # dimension keeps at least three points).
    for (n, periodic) in (((8, 6), (true, false)), ((6, 4), (false, false)), ((4, 4, 6), (true, true, false)),
                          ((7, 5), (false, false)), ((6, 9), (true, false)), ((8, 2, 5), (true, false, false)))
        D = length(n)
        r = coarsening_factors(n, periodic)
        @test r == ntuple(k -> (periodic[k] ? n[k] % 4 == 0 && n[k] >= 8 : n[k] >= 5) ? 2 : 1, D)
        nc = coarse_points(n, r)
        @test nc == ntuple(k -> r[k] == 2 ? cld(n[k], 2) : n[k], D)
        layout = BoxLayout(n, ntuple(_ -> 1, D), 0; overlap = 1, periodic)
        layoutc = coarse_layout(layout)
        @test layoutc.points == nc && layoutc.periodic == periodic
        coef = 0.1 .+ rand(rng, n..., 2D + 1)
        coef[ntuple(_ -> Colon(), D)..., 1] .+= 2D
        op = box_operator(layout, CoefficientStencil(coef))
        opc = box_operator(layoutc, CoefficientStencil(zeros(nc..., 2D + 1)))

        # The transfers as matrices: E (pullback, copy onto the block) and T
        # (pushforward transfer, block average), from the kernels.
        N, Nc = prod(n), prod(nc)
        E = zeros(N, Nc)
        T = zeros(Nc, N)
        xc, xf = grid_zeros(opc), grid_zeros(op)
        for k in 1:Nc
            fill!(xc, 0); fill!(xf, 0)
            interior(xc, opc)[k] = 1
            prolong_add!(xf, op, xc, opc)
            E[:, k] = vec(interior(xf, op))
        end
        for k in 1:N
            fill!(xc, 0); fill!(xf, 0)
            interior(xf, op)[k] = 1
            restrict_average!(xc, opc, xf, op)
            T[:, k] = vec(interior(xc, opc))
        end
        @test T * E ≈ I                                   # T is a left inverse of E
        @test T ≈ Diagonal(1 ./ vec(sum(E; dims = 1))) * E'   # the average: Eᵀ scaled by the block sizes
        ψ = aggregation_homomorphism(n, r)
        @test all(E[i, ψ.vertex_map[i]] == 1 for i in 1:N)   # E is the pullback along ψ

        # Galerkin coarse operator = T A E (and, for blocks of 2^D points, (R₀ A R₀ᵀ) / 2^D
        # of the aggregation coarse space).
        A = dense_operator(op, layout)
        coefc = zeros(nc..., 2D + 1)
        galerkin_coefficients!(coefc, op)
        Ac = dense_operator(box_operator(layoutc, CoefficientStencil(coefc)), layoutc)
        @test Ac ≈ T * A * E
        if r == ntuple(_ -> 2, D) && all(iseven, n)
            cs = AggregateCoarseSpace(layout, nc)
            @test Ac ≈ coarse_matrix(cs, op, SerialBoxes()) / 2^D
        end

        # The same pair from the sheaf machinery: the pushforward of the constant
        # sheaf along ψ has one-dimensional stalks (the constants on each block),
        # its transfer map is T, and the fibre bases are E, up to a scaling cᵥ
        # of each basis vector.
        if all(n .> 2) || !any(periodic)
            F = constant_sheaf(n, periodic)
            pushed = pushforward_sheaf(ψ, F)
            @test vertex_stalks(pushed) == ones(Int, Nc)
            bases = all_fiber_bases(ψ, F)
            Esheaf = zeros(N, Nc)
            for w in 1:Nc
                Esheaf[fiber_vertices(ψ, w), w] = bases[w]
            end
            Tsheaf = Matrix(pushforward_transfer_map(ψ, F))
            @test Tsheaf * Esheaf ≈ I
            C = Diagonal([bases[w][1] for w in 1:Nc])
            @test Esheaf ≈ E * C
            @test Tsheaf * A * Esheaf ≈ C \ (T * A * E) * C
        end
    end

    # Constant drift: the Galerkin coarse operator is the upwind operator
    # discretized again at spacing 2h (rates halved, same ρ).
    n = (8, 8, 4)
    D = 3
    ρ, a = 0.3, (1.5, -0.7, 2.0)
    coef = zeros(n..., 2D + 1)
    coef[:, :, :, 1] .= ρ + sum(abs, a)
    for j in 1:D
        coef[:, :, :, a[j] > 0 ? 1 + D + j : 1 + j] .= abs(a[j])
    end
    layout = BoxLayout(n, (1, 1, 1), 0; periodic = (true, true, true))
    coefc = zeros(coarse_points(n, (2, 2, 2))..., 2D + 1)
    galerkin_coefficients!(coefc, box_operator(layout, CoefficientStencil(coef)))
    @test all(coefc[:, :, :, 1] .≈ ρ + sum(abs, a) / 2)
    for j in 1:D
        @test all(coefc[:, :, :, (a[j] > 0 ? 1 + D + j : 1 + j)] .≈ abs(a[j]) / 2)
        @test all(iszero, coefc[:, :, :, (a[j] > 0 ? 1 + j : 1 + D + j)])
    end
    @test_throws ArgumentError coarse_points((8, 7), (2, 3))
    @test coarsening_factors((5, 4, 12, 6, 4), (false, false, true, true, true)) == (2, 1, 2, 1, 1)
end
