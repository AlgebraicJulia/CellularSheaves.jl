using Test
using CellularSheaves
using LinearAlgebra
using BlockArrays

@testset "ControlledTrajectoryExamples" begin

    # -----------------------------------------------------------------------
    # 1. Double integrator
    #    Standard point-mass model: ẍ = u
    #    Ac = [0 1; 0 0], Bc = [0; 1]
    # -----------------------------------------------------------------------
    @testset "Double integrator" begin
        Ac = [0.0  1.0; 0.0  0.0]
        Bc = reshape([0.0, 1.0], 2, 1)
        h  = 0.25
        k  = 8
        F  = EuclideanSheaf{Float64}(fill(2, 1))
        ts = ControlledTrajectorySheaf(F, Ac, Bc, h, k)

        n = ts.state_dim    # 2
        m = ts.control_dim  # 1
        p = (k + 1) * n + k * m   # expected trajectory length

        @test p == (k + 1) * n + k * m

        x1  = [0.0, 0.0]
        xk1 = [1.0, 0.0]

        Q  = Matrix{Float64}(I, n, n)
        Ru = Matrix{Float64}(I, m, m)
        Qf = 10.0 * Matrix{Float64}(I, n, n)

        H, f, _ = lqr_objective(ts, Q, Ru; Qf=Qf)
        z_opt, _, _, _ = optimal_control_trajectory(ts, x1, xk1, H, f)

        # Total coordinate length
        @test length(Array(z_opt)) == p

        # Initial and terminal state blocks match requested endpoints
        @test Array(z_opt[Block(1)])     ≈ x1
        @test Array(z_opt[Block(k + 1)]) ≈ xk1
    end

    # -----------------------------------------------------------------------
    # 2. Vehicle platoon
    #    Two stacked double integrators; one per vehicle.
    # -----------------------------------------------------------------------
    @testset "Vehicle platoon" begin
        Ac = [0.0  1.0  0.0  0.0;
              0.0  0.0  0.0  0.0;
              0.0  0.0  0.0  1.0;
              0.0  0.0  0.0  0.0]
        Bc = [0.0  0.0;
              1.0  0.0;
              0.0  0.0;
              0.0  1.0]
        h  = 0.25
        k  = 8
        F  = EuclideanSheaf{Float64}([2, 2])
        ts = ControlledTrajectorySheaf(F, Ac, Bc, h, k)

        n = ts.state_dim    # 4
        m = ts.control_dim  # 2

        x1  = [0.0, 0.0, 2.0, 0.0]
        xk1 = [1.0, 0.0, 3.0, 0.0]

        Q  = Matrix{Float64}(I, n, n)
        Ru = Matrix{Float64}(I, m, m)
        Qf = 10.0 * Matrix{Float64}(I, n, n)

        H, f, _ = lqr_objective(ts, Q, Ru; Qf=Qf)
        z_opt, _, _, _ = optimal_control_trajectory(ts, x1, xk1, H, f)

        # Optimization succeeds (no exception) and endpoints are satisfied
        @test Array(z_opt[Block(1)])     ≈ x1
        @test Array(z_opt[Block(k + 1)]) ≈ xk1

        # Endpoint constraints per vehicle: vehicle 1 position and velocity
        @test Array(z_opt[Block(1)])[1:2]     ≈ x1[1:2]
        @test Array(z_opt[Block(k+1)])[1:2]   ≈ xk1[1:2]
        # Vehicle 2
        @test Array(z_opt[Block(1)])[3:4]     ≈ x1[3:4]
        @test Array(z_opt[Block(k+1)])[3:4]   ≈ xk1[3:4]
    end

    # -----------------------------------------------------------------------
    # 3. Planar quadrotor
    #    Hover-linearized model: 6-state, 2-control
    # -----------------------------------------------------------------------
    @testset "Planar quadrotor" begin
        g_quad = 9.81
        m_quad = 0.5
        I_quad = 0.01
        ℓ_quad = 0.25

        Ac = [0.0  0.0    0.0    1.0  0.0  0.0;
              0.0  0.0    0.0    0.0  1.0  0.0;
              0.0  0.0    0.0    0.0  0.0  1.0;
              0.0  0.0   -g_quad 0.0  0.0  0.0;
              0.0  0.0    0.0    0.0  0.0  0.0;
              0.0  0.0    0.0    0.0  0.0  0.0]

        Bc = [0.0             0.0;
              0.0             0.0;
              0.0             0.0;
              0.0             0.0;
              1.0/m_quad      1.0/m_quad;
              ℓ_quad/(2.0*I_quad)  -ℓ_quad/(2.0*I_quad)]

        h = 0.05
        k = 12
        F = EuclideanSheaf{Float64}(fill(6, 1))
        ts = ControlledTrajectorySheaf(F, Ac, Bc, h, k)

        n     = ts.state_dim    # 6
        m_dim = ts.control_dim  # 2

        x1  = zeros(6)
        xk1 = [0.5, 0.0, 0.0, 0.0, 0.0, 0.0]

        Q  = Diagonal([1.0, 1.0, 2.0, 0.5, 0.5, 0.5])
        Ru = Matrix{Float64}(I, m_dim, m_dim)
        Qf = Diagonal([10.0, 10.0, 20.0, 1.0, 1.0, 1.0])

        H, f, _ = lqr_objective(ts, Matrix(Q), Ru; Qf=Matrix(Qf))
        z_opt, _, _, _ = optimal_control_trajectory(ts, x1, xk1, H, f)

        # Control blocks have 2 coordinates per time step
        for t in 1:k
            ctrl_block = Array(z_opt[Block(k + 1 + t)])
            @test length(ctrl_block) == 2
        end

        # All entries are finite
        @test all(isfinite, Array(z_opt))
    end

    # -----------------------------------------------------------------------
    # 4. Mass-spring-damper chain
    #    Two masses connected by springs/dampers, wall on the left.
    # -----------------------------------------------------------------------
    @testset "Mass-spring-damper chain" begin
        k_spring = 1.0
        c_damp   = 0.2
        mass_val = 1.0

        # State ordering: [x₁, ẋ₁, x₂, ẋ₂] — matches the two-vertex EuclideanSheaf
        Ac = [0.0           1.0         0.0          0.0;
              -2k_spring   -2c_damp     k_spring     c_damp;
               0.0          0.0         0.0          1.0;
               k_spring     c_damp     -k_spring    -c_damp]

        Bc = [0.0          0.0;
              1.0/mass_val 0.0;
              0.0          0.0;
              0.0          1.0/mass_val]

        h = 0.5
        k = 10
        F = EuclideanSheaf{Float64}([2, 2])
        ts = ControlledTrajectorySheaf(F, Ac, Bc, h, k)

        n     = ts.state_dim    # 4
        m_dim = ts.control_dim  # 2

        x1  = [0.0, 0.0, 0.0, 0.0]
        xk1 = [0.5, 0.0, 1.0, 0.0]   # mass 1 at 0.5 m, mass 2 at 1.0 m, both at rest

        Q  = Matrix{Float64}(I, n, n)
        Ru = Matrix{Float64}(I, m_dim, m_dim)
        Qf = 5.0 * Matrix{Float64}(I, n, n)

        H, f, _ = lqr_objective(ts, Q, Ru; Qf=Qf)
        z_opt, _, _, _ = optimal_control_trajectory(ts, x1, xk1, H, f)

        # Optimal trajectory satisfies the requested endpoints
        @test Array(z_opt[Block(1)])     ≈ x1
        @test Array(z_opt[Block(k + 1)]) ≈ xk1

        # State dimension is larger than the single-agent double integrator (2)
        @test ts.state_dim > 2
    end

end
