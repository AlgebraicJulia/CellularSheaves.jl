using Test
using CellularSheaves
using LinearAlgebra
using SparseArrays
using Random
using MPI: mpiexec

# The sparse matrix of a stored-coefficient stencil on a box (no neighbours
# outside it), for checking the implicit kernels.
function assembled(coef)
    D = ndims(coef) - 1
    n = ntuple(k -> size(coef, k), D)
    L = LinearIndices(n)
    rows, cols, vals = Int[], Int[], Float64[]
    for I in CartesianIndices(n)
        push!(rows, L[I]); push!(cols, L[I]); push!(vals, coef[I, 1])
        for j in 1:D
            e = CartesianIndex(ntuple(k -> k == j ? 1 : 0, D))
            if I[j] > 1
                push!(rows, L[I]); push!(cols, L[I - e]); push!(vals, -coef[I, 1 + j])
            end
            if I[j] < n[j]
                push!(rows, L[I]); push!(cols, L[I + e]); push!(vals, -coef[I, 1 + D + j])
            end
        end
    end
    return sparse(rows, cols, vals, prod(n), prod(n))
end

# A random upwind M-matrix stencil: one nonzero neighbour per axis, diagonal
# dominance ρ, and no coupling across the box boundary.
function upwind_coefficients(rng, n; ρ=0.5)
    D = length(n)
    coef = zeros(n..., 2D + 1)
    for I in CartesianIndices(n)
        diagonal = ρ
        for j in 1:D
            a = rand(rng)
            diagonal += a
            if rand(rng, Bool)
                I[j] < n[j] && (coef[I, 1 + D + j] = a)
            else
                I[j] > 1 && (coef[I, 1 + j] = a)
            end
        end
        coef[I, 1] = diagonal
    end
    return coef
end

@testset "GridSchwarz" begin
    rng = Xoshiro(7)
    for n in ((17, 13), (6, 5, 7), (5, 4, 6, 3))
        coef = upwind_coefficients(rng, n)
        A = assembled(coef)
        op = GridOperator(coef; ghost=2)
        @test size(op) == n && ndims(op) == length(n)
        x = grid_zeros(op)
        interior(x, op) .= randn(rng, n...)
        y = grid_zeros(op)
        apply!(y, op, x)
        @test vec(interior(y, op)) ≈ A * vec(interior(x, op))
        # Red–black SGS equals the multicolor sparse SGS (red = even index sum first).
        r = grid_zeros(op)
        interior(r, op) .= randn(rng, n...)
        apply!(y, op, x)
        z = grid_zeros(op)
        red_black_sgs!(z, op, r)
        M = SymmetricGaussSeidel(A)
        @test length(M.colors) == 2 && 1 in M.colors[1]
        @test vec(interior(z, op)) ≈ ldiv!(zeros(prod(n)), M, vec(interior(r, op))) rtol = 1e-12
        @test all(iszero, z[1:op.origin[1], ntuple(_ -> :, length(n) - 1)...])     # zero Dirichlet ghost layer
        # Parity flips the coloring.
        z1 = grid_zeros(op)
        red_black_sgs!(z1, GridOperator(coef; ghost=2, parity=1), r)
        q = reduce(vcat, reverse(M.colors))
        Aq = A[q, q]
        Dq, Lq, Uq = Diagonal(Aq), tril(Aq, -1), triu(Aq, 1)
        @test vec(interior(z1, op))[q] ≈ Matrix(Dq + Uq) \ (Dq * (Matrix(Dq + Lq) \ vec(interior(r, op))[q])) rtol = 1e-10
        # Reductions over the interior only.
        x[1] = 1e6                                                       # a ghost value, ignored
        @test grid_dot(x, r, op, zeros(Base.tail(n)...)) ≈ dot(interior(x, op), interior(r, op))
        @test grid_reduce((a, b) -> abs(a - b), max, x, r, op, zeros(Base.tail(n)...)) ≈
            maximum(abs, interior(x, op) - interior(r, op))
        x[1] = 0
        # BiCGStab with the red–black preconditioner solves the system.
        b = grid_zeros(op)
        interior(b, op) .= randn(rng, n...)
        u = grid_zeros(op)
        ws = GridWorkspace(op)
        its, ok, residuals = grid_bicgstab!(u, op, (yy, vv) -> red_black_sgs!(yy, op, vv), b, ws; tol=1e-11, maxiter=500)
        @test ok
        @test vec(interior(u, op)) ≈ A \ vec(interior(b, op)) rtol = 1e-9
        @test last(residuals) <= 1e-11 * norm(interior(b, op))
        its0, ok0, _ = grid_bicgstab!(fill!(grid_zeros(op), 0), op, copyto!, b, ws; tol=1e-11, maxiter=2000)
        @test ok0 && its0 > its                                          # the smoother helps
    end
    @test_throws ArgumentError GridOperator(zeros(4, 4, 4))              # needs 2D + 1 coefficients
    @test_throws ArgumentError GridOperator(zeros(4, 4, 5); ghost=0)
end

@testset "GridSchwarz boxes" begin
    layout = BoxLayout((10, 7), (3, 2), 4; overlap = 1)                     # rank 4: coords (1, 1)
    @test layout.coords == (1, 1)
    @test layout.owned == (4:6, 4:7)
    @test layout.extended == (3:7, 3:7)
    @test layout.neighbors == ((3, 5), (1, -1))
    @test box_coordinates((3, 2), 5) == (2, 1)
    @test balanced_ranks(12, (40, 10)) == (6, 2)
    @test balanced_ranks(16, (33, 33, 33, 33)) == (2, 2, 2, 2)
    @test prod(balanced_ranks(6, (9, 9, 9, 9))) == 6
    @test_throws ArgumentError BoxLayout((10, 7), (3, 7), 0; overlap = 1)  # boxes of one point < width 2
    s = CoefficientStencil(zeros(3, 4, 5))
    own = box_operator(layout, s)
    @test own.origin == (2, 2) && own.padded == (7, 8) && size(own) == (3, 4)
    ext = box_operator(layout, CoefficientStencil(zeros(5, 5, 5)), :extended)
    @test ext.origin == (1, 1) && ext.padded == own.padded && size(ext) == (5, 5)
    @test own.parity == (4 + 4 - 2) & 1 && ext.parity == (3 + 3 - 2) & 1
    # One box: exchange is a no-op, gather returns the interior.
    single = BoxLayout((5, 4), (1, 1), 0)
    op = box_operator(single, CoefficientStencil(zeros(5, 4, 5)))
    x = grid_zeros(op)
    interior(x, op) .= reshape(1.0:20.0, 5, 4)
    @test exchange!(SerialBoxes(), single, x) === x
    @test gather_boxes(SerialBoxes(), single, x, op) == reshape(1.0:20.0, 5, 4)

    # The same tests under MPI on several ranks.

    script = joinpath(@__DIR__, "..", "mpi", "grid_boxes.jl")
    project = Base.active_project()
    for n in (1, 3, 4)
        ok = mpiexec() do exe                                        # the launcher environment applies inside the block
            cmd = `$exe -n $n $(Base.julia_cmd()) --project=$project --startup-file=no $script`
            success(pipeline(cmd; stdout, stderr))
        end
        @test ok
    end
end
