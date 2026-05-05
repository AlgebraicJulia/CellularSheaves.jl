using Test
using CellularSheaves
using LinearAlgebra
using BlockArrays

@testset "TrajectorySheaf" begin

    # -----------------------------------------------------------------------
    # Setup: 2-D state space, k=4 steps
    # -----------------------------------------------------------------------
    d = 2
    k = 4
    F = EuclideanSheaf{Float64}(fill(d, 1))   # base sheaf with one d-dim stalk

    # Identity dynamics
    A_id = Matrix{Float64}(I, d, d)

    # Non-trivial: 2D rotation by π/6
    θ = π / 6
    A_rot = [cos(θ) -sin(θ); sin(θ) cos(θ)]

    # -----------------------------------------------------------------------
    # 1. Constructor stores correct k and A
    # -----------------------------------------------------------------------
    @testset "constructor stores k and A" begin
        tsheaf = TrajectorySheaf(F, A_id, k)
        @test tsheaf.k == k
        @test tsheaf.A ≈ A_id
    end

    # -----------------------------------------------------------------------
    # 2. Inner sheaf has exactly k+1 vertex stalks, each of dimension d
    # -----------------------------------------------------------------------
    @testset "inner sheaf vertex stalks" begin
        tsheaf = TrajectorySheaf(F, A_id, k)
        vs = vertex_stalks(tsheaf.sheaf)
        @test length(vs) == k + 1
        @test all(v -> v == d, vs)
    end

    # -----------------------------------------------------------------------
    # 3. Inner sheaf has exactly k edges
    # -----------------------------------------------------------------------
    @testset "inner sheaf edge count" begin
        tsheaf = TrajectorySheaf(F, A_id, k)
        g = underlying_graph(tsheaf.sheaf)
        @test ne(g) == k
    end

    # -----------------------------------------------------------------------
    # 4. Restriction maps: vertex t -> A, vertex t+1 -> I
    # -----------------------------------------------------------------------
    @testset "restriction maps" begin
        tsheaf = TrajectorySheaf(F, A_rot, k)
        Id = Matrix{Float64}(I, d, d)
        for t in 1:k
            # rm from vertex t  is A (so d_e x = A*x[t] - I*x[t+1] = 0 => x[t+1] = A*x[t])
            @test get_restriction_map(tsheaf.sheaf, t, t + 1) ≈ A_rot
            # rm from vertex t+1 is I
            @test get_restriction_map(tsheaf.sheaf, t + 1, t) ≈ Id
        end
    end

    # -----------------------------------------------------------------------
    # 5. A=I: colocation returns linear interpolation between x0 and xk
    # -----------------------------------------------------------------------
    @testset "identity dynamics -> linear interpolation" begin
        tsheaf = TrajectorySheaf(F, A_id, k)
        x0 = [1.0, 2.0]
        xk = [5.0, 6.0]
        traj = colocation_trajectory(tsheaf, x0, xk)

        tol = 1e-8
        for t in 1:(k + 1)
            # Linear interpolation: x_t = x0 + (t-1)/k * (xk - x0)
            expected = x0 + ((t - 1) / k) * (xk - x0)
            @test norm(Array(traj[Block(t)]) - expected) < tol
        end
    end

    # -----------------------------------------------------------------------
    # 6. Non-trivial A: returned trajectory satisfies x[t+1] ≈ A*x[t]
    # -----------------------------------------------------------------------
    @testset "non-trivial A: dynamics consistency" begin
        tsheaf = TrajectorySheaf(F, A_rot, k)
        x0 = [1.0, 0.0]
        xk = A_rot^k * x0          # exact k-step trajectory endpoint
        traj = colocation_trajectory(tsheaf, x0, xk)

        tol = 1e-8
        for t in 1:k
            xt   = Array(traj[Block(t)])
            xt1  = Array(traj[Block(t + 1)])
            @test norm(xt1 - A_rot * xt) < tol
        end
    end

    # -----------------------------------------------------------------------
    # 7. Boundary values are exactly preserved
    # -----------------------------------------------------------------------
    @testset "boundary values preserved" begin
        tsheaf = TrajectorySheaf(F, A_rot, k)
        x0 = [1.0, 2.0]
        xk = [3.0, 4.0]
        traj = colocation_trajectory(tsheaf, x0, xk)

        @test Array(traj[Block(1)])     ≈ x0
        @test Array(traj[Block(k + 1)]) ≈ xk
    end

    # -----------------------------------------------------------------------
    # 8. Returned trajectory minimises Laplacian energy among all 0-cochains
    #    with the given boundary values (energy of perturbed cochain is higher)
    # -----------------------------------------------------------------------
    @testset "minimum Laplacian energy" begin
        tsheaf = TrajectorySheaf(F, A_rot, k)
        x0 = [1.0, 0.0]
        xk = [0.0, 1.0]
        traj = colocation_trajectory(tsheaf, x0, xk)

        L = sheaf_laplacian_matrix(tsheaf.sheaf)
        x_arr = Array(traj)
        E_opt = x_arr' * L * x_arr

        # Perturb an interior DOF and check that energy increases
        x_pert = copy(x_arr)
        x_pert[d + 1] += 0.5          # perturb first component of block 2
        E_pert = x_pert' * L * x_pert

        @test E_opt <= E_pert + 1e-10
    end

    # -----------------------------------------------------------------------
    # 9. Passing a wrong-length x0 or xk throws ArgumentError
    # -----------------------------------------------------------------------
    @testset "wrong-length x0 / xk throws ArgumentError" begin
        tsheaf = TrajectorySheaf(F, A_id, k)
        @test_throws ArgumentError colocation_trajectory(tsheaf, [1.0], [1.0, 2.0])
        @test_throws ArgumentError colocation_trajectory(tsheaf, [1.0, 2.0], [1.0])
        @test_throws ArgumentError colocation_trajectory(tsheaf, [1.0, 2.0, 3.0], [1.0, 2.0])
    end

    # -----------------------------------------------------------------------
    # 10. Passing k < 1 throws ArgumentError
    # -----------------------------------------------------------------------
    @testset "k < 1 throws ArgumentError" begin
        @test_throws ArgumentError TrajectorySheaf(F, A_id, 0)
        @test_throws ArgumentError TrajectorySheaf(F, A_id, -1)
    end

    # -----------------------------------------------------------------------
    # 11. Non-square A or A of wrong size throws ArgumentError
    # -----------------------------------------------------------------------
    @testset "non-square or wrong-size A throws ArgumentError" begin
        # Non-square
        @test_throws ArgumentError TrajectorySheaf(F, zeros(Float64, 2, 3), k)
        # Square but wrong size (d=2 but A is 3×3)
        @test_throws ArgumentError TrajectorySheaf(F, Matrix{Float64}(I, 3, 3), k)
    end

end
