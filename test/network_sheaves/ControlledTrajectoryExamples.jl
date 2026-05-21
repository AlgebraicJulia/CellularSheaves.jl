using Test
using CellularSheaves
using LinearAlgebra
using BlockArrays

function public_dynamics_constraints(ts)
    n = ts.state_dim
    m = ts.control_dim
    k = ts.k
    p = (k + 1) * n + k * m

    Aeq = zeros(k * n + 2n, p)
    row = 1

    Aeq[row:row+n-1, 1:n] .= Matrix{Float64}(I, n, n)
    row += n

    for t in 1:k
        x_t  = (t - 1) * n + 1 : t * n
        x_tp = t * n + 1 : (t + 1) * n
        u_t  = (k + 1) * n + (t - 1) * m + 1 : (k + 1) * n + t * m

        Aeq[row:row+n-1, x_t]  .= ts.Ad
        Aeq[row:row+n-1, u_t]  .= ts.Bd
        Aeq[row:row+n-1, x_tp] .-= Matrix{Float64}(I, n, n)
        row += n
    end

    Aeq[row:row+n-1, k*n+1:(k+1)*n] .= Matrix{Float64}(I, n, n)
    return Aeq
end

function full_space_kkt_residual(ts, z, H, f)
    Z = nullspace(public_dynamics_constraints(ts))
    return norm(Z' * (H * Array(z) + f))
end

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
        z_opt, _, _, N = optimal_control_trajectory(ts, x1, xk1, H, f)

        # Total coordinate length
        @test length(Array(z_opt)) == p

        @test size(N, 2) == expected_feasible_dimension(ts)

        # Initial and terminal state blocks match requested endpoints
        @test Array(z_opt[Block(1)])     ≈ x1
        @test Array(z_opt[Block(k + 1)]) ≈ xk1
        @test full_space_kkt_residual(ts, z_opt, H, f) < 1e-8
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
        z_opt, _, _, N = optimal_control_trajectory(ts, x1, xk1, H, f)

        # Optimization succeeds (no exception) and endpoints are satisfied
        @test size(N, 2) == expected_feasible_dimension(ts)
        @test Array(z_opt[Block(1)])     ≈ x1
        @test Array(z_opt[Block(k + 1)]) ≈ xk1
        @test full_space_kkt_residual(ts, z_opt, H, f) < 1e-8

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
        z_opt, _, _, null_basis = optimal_control_trajectory(ts, x1, xk1, H, f)

        # Control blocks have 2 coordinates per time step
        @test size(null_basis, 2) == expected_feasible_dimension(ts)
        for t in 1:k
            ctrl_block = Array(z_opt[Block(k + 1 + t)])
            @test length(ctrl_block) == 2
        end

        # All entries are finite
        @test all(isfinite, Array(z_opt))
        # Tolerance is 1e-7 (vs 1e-8 elsewhere): the quadrotor's LDL null-basis uses
        # pivots near the sqrt(eps) threshold, producing a slightly less orthogonal
        # basis than other examples; the optimality gap is still well within 1e-7.
        @test full_space_kkt_residual(ts, z_opt, H, f) < 1e-7
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
        z_opt, _, _, null_basis = optimal_control_trajectory(ts, x1, xk1, H, f)

        # Optimal trajectory satisfies the requested endpoints
        @test size(null_basis, 2) == expected_feasible_dimension(ts)
        @test Array(z_opt[Block(1)])     ≈ x1
        @test Array(z_opt[Block(k + 1)]) ≈ xk1

        # State dimension is larger than the single-agent double integrator (2)
        @test ts.state_dim > 2
        @test full_space_kkt_residual(ts, z_opt, H, f) < 1e-8
    end

    # -----------------------------------------------------------------------
    # 5. SVD vs LDL nullspace dimensionality agreement
    #    Verifies that harmonic_extension_svd_diagnostics and
    #    harmonic_extension_ldl_diagnostics agree on the nullity of the
    #    restricted Laplacian for the controlled-trajectory boundary problem,
    #    and that both match expected_feasible_dimension.
    # -----------------------------------------------------------------------
    @testset "SVD and LDL nullspace dimensionality agree" begin
        systems = [
            # (name, Ac, Bc, h, k, F)
            ("Double integrator",
             [0.0 1.0; 0.0 0.0],
             reshape([0.0, 1.0], 2, 1),
             0.25, 8,
             EuclideanSheaf{Float64}(fill(2, 1))),

            ("Vehicle platoon",
             [0.0 1.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 1.0; 0.0 0.0 0.0 0.0],
             [0.0 0.0; 1.0 0.0; 0.0 0.0; 0.0 1.0],
             0.25, 8,
             EuclideanSheaf{Float64}([2, 2])),

            ("Planar quadrotor",
             [0.0 0.0 0.0 1.0 0.0 0.0;
              0.0 0.0 0.0 0.0 1.0 0.0;
              0.0 0.0 0.0 0.0 0.0 1.0;
              0.0 0.0 -9.81 0.0 0.0 0.0;
              0.0 0.0 0.0 0.0 0.0 0.0;
              0.0 0.0 0.0 0.0 0.0 0.0],
             [0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0; 2.0 2.0; 12.5 -12.5],
             0.05, 12,
             EuclideanSheaf{Float64}(fill(6, 1))),

            ("Mass-spring-damper chain",
             [0.0 1.0 0.0 0.0; -2.0 -0.4 1.0 0.2; 0.0 0.0 0.0 1.0; 1.0 0.2 -1.0 -0.2],
             [0.0 0.0; 1.0 0.0; 0.0 0.0; 0.0 1.0],
             0.5, 10,
             EuclideanSheaf{Float64}([2, 2])),
        ]

        for (name, Ac, Bc, h, k, F) in systems
            @testset "$name" begin
                ts = ControlledTrajectorySheaf(F, Ac, Bc, h, k)
                n  = ts.state_dim
                x1  = zeros(n)
                xk1 = zeros(n); xk1[1] = 1.0

                boundary = Dict{Int,Vector{Float64}}(
                    1     => x1,
                    k + 2 => xk1,
                )

                ldl_diag = harmonic_extension_ldl_diagnostics(ts.sheaf, boundary)
                svd_diag = harmonic_extension_svd_diagnostics(ts.sheaf, boundary)

                # LDL and SVD must agree on nullity
                @test ldl_diag.nullity_estimate == svd_diag.nullity_estimate

                # Both must match the controllability-based expected dimension
                expected = expected_feasible_dimension(ts)
                @test ldl_diag.nullity_estimate == expected
                @test svd_diag.nullity_estimate == expected
            end
        end
    end

end
