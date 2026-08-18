using Test
using CellularSheaves
using CellularSheaves.ControlSheaves.Tikhonov
using LinearAlgebra
import CellularSheaves.NetworkSheaves.EuclideanSheaves:
    _harmonic_extension_restricted_laplacian

@testset "TikhonovFilter" begin
    H = [2.0 -1.0; -1.0 2.0]
    rhs = [1.0, -0.5]
    qstar = tikhonov_equilibrium(H, rhs)

    @test H * qstar ≈ rhs
    @test tikhonov_reference_rate(H, [0.2, -0.1]) ≈ H \ [0.2, -0.1]
    @test tikhonov_dissipation([1.0, -2.0], 0.1) ≈ -50.0
    @test_throws ArgumentError TikhonovFilter(zeros(2); epsilon = 0.0)
    @test_throws ArgumentError tikhonov_step!(TikhonovFilter(zeros(2); epsilon = 0.1), zeros(2), zeros(2), 0.0)

    @testset "Factorization overloads" begin
        F = factorize(H)
        rhs_rate = [0.2, -0.1]
        @test tikhonov_equilibrium(F, rhs) ≈ tikhonov_equilibrium(H, rhs)
        @test tikhonov_reference_rate(F, rhs_rate) ≈ tikhonov_reference_rate(H, rhs_rate)
        @test tikhonov_equilibrium(F, rhs) ≈ H \ rhs
    end

    @testset "constant reference exact solution" begin
        epsilon = 0.2
        x0 = [3.0, -2.0]
        filter = TikhonovFilter(x0; epsilon)
        dt = 0.07
        for _ in 1:10
            tikhonov_step!(filter, qstar, qstar, dt)
        end
        expected = qstar + exp(-10dt / epsilon) * (x0 - qstar)
        @test filter.x ≈ expected atol = 2e-14
    end

    @testset "affine reference and feedforward" begin
        epsilon = 0.3
        a = [0.4, -0.2]
        velocity = [0.5, 0.15]
        reference(t) = a + t * velocity
        duration = 2.0
        dt = 0.05

        baseline = TikhonovFilter([1.2, -0.7]; epsilon)
        x0 = copy(baseline.x)
        for k in 0:Int(duration / dt)-1
            tikhonov_step!(baseline, reference, k * dt, dt)
        end
        expected = reference(duration) - epsilon * velocity +
                   exp(-duration / epsilon) * (x0 - a + epsilon * velocity)
        @test baseline.x ≈ expected atol = 2e-14

        corrected = TikhonovFilter(a; epsilon)
        for k in 0:Int(duration / dt)-1
            q0 = reference(k * dt)
            q1 = reference((k + 1) * dt)
            u0 = tikhonov_feedforward_reference(q0, velocity, epsilon)
            u1 = tikhonov_feedforward_reference(q1, velocity, epsilon)
            tikhonov_step!(corrected, u0, u1, dt)
        end
        @test corrected.x ≈ reference(duration) atol = 2e-14
    end

    @testset "small epsilon approaches direct delivery" begin
        filter = TikhonovFilter([1.0, -1.0]; epsilon = 1e-12)
        q0 = [0.0, 0.0]
        q1 = [2.0, 3.0]
        tikhonov_step!(filter, q0, q1, 0.05)
        @test all(isfinite, filter.x)
        @test filter.x ≈ q1 atol = 1e-9
    end

    @testset "harmonic extension supplies the reference" begin
        sheaf = EuclideanSheaf{Float64}(fill(1, 3))
        I1 = ones(1, 1)
        add_sheaf_edge!(sheaf, 1, 2, I1, I1)
        add_sheaf_edge!(sheaf, 2, 3, I1, I1)
        boundary = Dict(3 => [2.0])
        _, interior, Hsheaf, LIB = _harmonic_extension_restricted_laplacian(sheaf, boundary)
        rhs_sheaf = -Matrix(LIB) * [2.0]
        qstar_sheaf = tikhonov_equilibrium(Matrix(Hsheaf), rhs_sheaf)

        harmonic, _ = harmonic_extension(sheaf, boundary)
        offsets = [0; cumsum(sheaf.vertex_stalks)]
        interior_dofs = vcat([collect(offsets[v] + 1:offsets[v + 1]) for v in interior]...)
        @test qstar_sheaf ≈ Array(harmonic)[interior_dofs]

        filter = TikhonovFilter(zeros(length(qstar_sheaf)); epsilon = 0.1)
        for _ in 1:120
            tikhonov_step!(filter, qstar_sheaf, qstar_sheaf, 0.01)
        end
        @test norm(filter.x - qstar_sheaf) < 1e-4
    end
end
