using Test
using CellularSheaves
using LinearAlgebra
using BlockArrays

@testset "NullspaceTrajectoryFamily" begin
    Ac = reshape([0.0], 1, 1)
    Bc = reshape([1.0], 1, 1)
    h = 1.0
    F = EuclideanSheaf{Float64}(fill(1, 1))

    @testset "shape and block layout" begin
        k = 4
        ts = ControlledTrajectorySheaf(F, Ac, Bc, h, k)
        x1 = [0.0]
        xk1 = [1.0]
        z_p, N = feasible_control_trajectory_basis(ts, x1, xk1)
        @test size(N, 2) > 0

        family = nullspace_trajectory_family(ts, z_p, N; amplitude=0.5)
        @test length(family) == size(N, 2)

        for z in family
            @test z isa BlockVector
            @test length(blocks(z)) == (k + 1) + k
            for t in 1:k+1
                @test length(Array(z[Block(t)])) == ts.state_dim
            end
            for t in 1:k
                @test length(Array(z[Block(k + 1 + t)])) == ts.control_dim
            end
        end
    end

    @testset "basis-direction reconstruction" begin
        k = 3
        ts = ControlledTrajectorySheaf(F, Ac, Bc, h, k)
        x1 = [0.0]
        xk1 = [1.0]
        z_p, N = feasible_control_trajectory_basis(ts, x1, xk1)
        amp = 0.25

        family_plus = nullspace_trajectory_family(ts, z_p, N; amplitude=amp)
        for j in 1:size(N, 2)
            @test Array(family_plus[j]) ≈ z_p + amp * N[:, j] atol=1e-10
        end

        family_pm = nullspace_trajectory_family(ts, z_p, N; amplitude=amp, include_negative=true)
        @test length(family_pm) == 2 * size(N, 2)
        for j in 1:size(N, 2)
            @test Array(family_pm[2j - 1]) ≈ z_p + amp * N[:, j] atol=1e-10
            @test Array(family_pm[2j]) ≈ z_p - amp * N[:, j] atol=1e-10
        end
    end

    @testset "feasibility preservation" begin
        k = 3
        ts = ControlledTrajectorySheaf(F, Ac, Bc, h, k)
        x1 = [0.0]
        xk1 = [1.0]
        z_p, N = feasible_control_trajectory_basis(ts, x1, xk1)
        family = nullspace_trajectory_family(ts, z_p, N; amplitude=1.0)

        for z in family
            @test Array(z[Block(1)]) ≈ x1 atol=1e-10
            @test Array(z[Block(k + 1)]) ≈ xk1 atol=1e-10
            for t in 1:k
                xt = Array(z[Block(t)])
                xt1 = Array(z[Block(t + 1)])
                ut = Array(z[Block(k + 1 + t)])
                @test norm(ts.Ad * xt + ts.Bd * ut - xt1) < 1e-10
            end
        end
    end

    @testset "empty basis behavior" begin
        k = 1
        ts = ControlledTrajectorySheaf(F, Ac, Bc, h, k)
        z_p, N = feasible_control_trajectory_basis(ts, [0.0], [1.0])
        @test size(N, 2) == 0
        @test isempty(nullspace_trajectory_family(ts, z_p, N))
    end

    @testset "argument validation" begin
        k = 2
        ts = ControlledTrajectorySheaf(F, Ac, Bc, h, k)
        z_p, N = feasible_control_trajectory_basis(ts, [0.0], [1.0])
        p = (k + 1) * ts.state_dim + k * ts.control_dim

        @test_throws ArgumentError nullspace_trajectory_family(ts, zeros(p + 1), N)
        @test_throws ArgumentError nullspace_trajectory_family(ts, z_p, zeros(p + 1, max(1, size(N, 2))))
        @test_throws ArgumentError nullspace_trajectory_family(ts, z_p, N; amplitude=-0.1)
    end
end
