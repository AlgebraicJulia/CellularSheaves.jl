using Test
using CellularSheaves
using LinearAlgebra
using SparseArrays
using Random

# A random nonsymmetric axis-stencil M-matrix on a grid of `n` points (column-major),
# diagonally dominant, zero Dirichlet values outside non-periodic dimensions.
function random_stencil_matrix(rng, n, periodic)
    D = length(n)
    L = LinearIndices(n)
    rows, cols, vals = Int[], Int[], Float64[]
    for I in CartesianIndices(n)
        d = 0.1
        for j in 1:D, s in (-1, 1)
            c = rand(rng)
            d += c
            k = I[j] + s
            if !(1 <= k <= n[j])
                periodic[j] || continue
                k = mod1(k, n[j])
            end
            push!(rows, L[I]); push!(cols, L[Base.setindex(Tuple(I), k, j)...]); push!(vals, -c)
        end
        push!(rows, L[I]); push!(cols, L[I]); push!(vals, d)
    end
    return sparse(rows, cols, vals, prod(n), prod(n))
end

# The base-change identities on every level of the hierarchy, for every box.
function commutes(A, n, boxes, factors)
    ok = true
    sizes, Aℓ, cover = n, A, boxes
    for r in factors
        E, T = prolongation_matrix(sizes, r), transfer_matrix(sizes, r)
        nc = coarse_points(sizes, r)
        Ac = T * Aℓ * E
        for U in cover
            V = coarse_box(U, r)
            RU, RV = box_restriction(sizes, U), box_restriction(nc, V)
            Ei, Ti = prolongation_matrix(length.(U), r), transfer_matrix(length.(U), r)
            size(Ti, 1) == size(RV, 1) || return false
            ok &= RV * T ≈ Ti * RU && RU * E ≈ Ei * RV && T * RU' ≈ RV' * Ti && E * RV' ≈ RU' * Ei
            ok &= Ti * (RU * Aℓ * RU') * Ei ≈ RV * Ac * RV'
        end
        sizes, Aℓ, cover = nc, Ac, [coarse_box(U, r) for U in cover]
    end
    return ok
end

relative_difference(B1, B2) = norm(B1 - B2) / norm(B1)

@testset "MultilevelSchwarz" begin
    rng = Random.MersenneTwister(11)

    # The transfers are the pullback along the aggregation homomorphism and the
    # block average.
    for (n, r) in (((8, 6), (2, 2)), ((7, 5), (2, 1)), ((6, 9, 4), (2, 2, 1)))
        E, T = prolongation_matrix(n, r), transfer_matrix(n, r)
        ψ = aggregation_homomorphism(n, r)
        @test all(E[i, ψ.vertex_map[i]] == 1 for i in 1:prod(n)) && nnz(E) == prod(n)
        @test T * E ≈ I
        @test size(E, 2) == prod(coarse_points(n, r))
    end
    R = box_restriction((5, 4), (2:3, 2:4))
    @test size(R) == (6, 20) && R * collect(1.0:20) == [7, 8, 12, 13, 17, 18]
    @test coarse_box((5:12, 3:3), (2, 1)) == (3:6, 3:3)
    @test is_aligned((15, 6), (9:15, 2:5), (2, 1))       # ends at the odd last point: a block of one
    @test !is_aligned((16, 6), (9:15, 1:6), (2, 1))      # cuts the block {15, 16}
    @test !is_aligned((16, 6), (6:16, 1:6), (2, 1))      # cuts the block {5, 6}
    @test_throws ArgumentError box_restriction((5, 4), (0:3, 1:4))

    # Aligned covers on every level: an even grid with two halvings (box edges on
    # multiples of 4, overlap 4), and an odd, semicoarsened grid with a periodic
    # dimension.
    cases = [
        ((16, 12), (false, false), [(2, 2), (2, 2)],
            [(x, y) for x in (1:8, 5:12, 9:16) for y in (1:8, 5:12)]),
        ((15, 6, 4), (false, true, false), [(2, 1, 2), (2, 1, 1)],
            [(x, y, 1:4) for x in (1:8, 5:15) for y in (1:4, 3:6)]),
    ]
    for (n, periodic, factors, boxes) in cases
        A = random_stencil_matrix(rng, n, periodic)
        sizes, cover = n, boxes
        for r in factors
            @test all(U -> is_aligned(sizes, U, r), cover)
            sizes, cover = coarse_points(sizes, r), [coarse_box(U, r) for U in cover]
        end

        # Base change: restriction to the cover commutes with T and E, so the
        # Galerkin coarse operator of each local problem is the local problem of
        # the Galerkin coarse operator, on every level.
        @test commutes(A, n, boxes, factors)

        # Additive over levels and boxes: the two orders are the same operator,
        # with Jacobi smoothing and with exact local solves on every level.
        B1 = schwarz_of_multigrid(A, n, boxes, factors)
        B2 = multigrid_of_schwarz(A, n, boxes, factors)
        @test B1 ≈ B2
        exact = M -> inv(Matrix(M))
        @test schwarz_of_multigrid(A, n, boxes, factors; smoother = exact) ≈
              multigrid_of_schwarz(A, n, boxes, factors; smoother = exact)

        # V-cycles over the levels: a sum over boxes of products over levels
        # against a product over levels of sums over boxes. They differ through
        # the couplings between boxes, and agree for the cover by one box.
        V1 = schwarz_of_multigrid(A, n, boxes, factors; composition = :multiplicative)
        V2 = multigrid_of_schwarz(A, n, boxes, factors; composition = :multiplicative)
        @info "V-cycle orders on $(n): relative difference $(relative_difference(V1, V2))"
        @test !(V1 ≈ V2)
        whole = [ntuple(k -> 1:n[k], length(n))]
        @test schwarz_of_multigrid(A, n, whole, factors; composition = :multiplicative) ≈
              multigrid_of_schwarz(A, n, whole, factors; composition = :multiplicative)
    end

    # A cover that cuts blocks: the square is not cartesian, and both the operator
    # identity and the additive equality fail.
    n, factors = (16, 12), [(2, 2)]
    A = random_stencil_matrix(rng, n, (false, false))
    boxes = [(x, y) for x in (1:7, 6:16) for y in (1:8, 5:12)]
    @test !all(U -> is_aligned(n, U, factors[1]), boxes)
    @test !commutes(A, n, boxes, factors)
    B1, B2 = schwarz_of_multigrid(A, n, boxes, factors), multigrid_of_schwarz(A, n, boxes, factors)
    @info "additive orders on a misaligned cover: relative difference $(relative_difference(B1, B2))"
    @test !(B1 ≈ B2)

    @test_throws ArgumentError schwarz_of_multigrid(A, n, boxes, factors; composition = :other)
    @test_throws ArgumentError multigrid_of_schwarz(A[1:10, 1:10], n, boxes, factors)
    @test_throws ArgumentError multigrid_of_schwarz(A, n, NTuple{2,UnitRange{Int}}[], factors)
end
