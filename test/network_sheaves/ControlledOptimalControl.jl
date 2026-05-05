using Test
using CellularSheaves
using LinearAlgebra
using BlockArrays
using SparseArrays

@testset "ControlledOptimalControl" begin

    # -----------------------------------------------------------------------
    # Primary test system: scalar integrator  ẋ = u
    #   Ac = [0], Bc = [1], h = 1  →  Ad = [1], Bd = [1]
    # -----------------------------------------------------------------------
    Ac_int = reshape([0.0], 1, 1)
    Bc_int = reshape([1.0], 1, 1)
    h      = 1.0
    F1     = EuclideanSheaf{Float64}(fill(1, 1))

    k  = 4
    ts = ControlledTrajectorySheaf(F1, Ac_int, Bc_int, h, k)

    n = ts.state_dim    # 1
    m = ts.control_dim  # 1
    p = (k + 1) * n + k * m  # 9

    Q  = Matrix{Float64}(I, n, n) * 0.0   # no state cost
    Ru = Matrix{Float64}(I, m, m)          # unit control cost
    Qf = Matrix{Float64}(I, n, n) * 0.0   # no terminal cost

    # -----------------------------------------------------------------------
    # 1. lqr_objective returns H, f, c of the correct sizes and sparse type
    # -----------------------------------------------------------------------
    @testset "lqr_objective sizes and sparse type" begin
        H, f, c = lqr_objective(ts, Q, Ru; Qf=Qf)
        @test size(H) == (p, p)
        @test H isa SparseMatrixCSC
        @test length(f) == p
        @test c isa Float64
    end

    # -----------------------------------------------------------------------
    # 2. With zero references, f is zero for the standard LQR cost
    # -----------------------------------------------------------------------
    @testset "zero references give zero f" begin
        H, f, c = lqr_objective(ts, Q, Ru; Qf=Qf)
        @test norm(f) ≈ 0.0
        @test c ≈ 0.0
    end

    # -----------------------------------------------------------------------
    # 3. optimal_control_trajectory returns trajectory satisfying endpoint states
    # -----------------------------------------------------------------------
    @testset "endpoint states satisfied" begin
        x1  = [0.0]
        xk1 = [1.0]
        H, f, _ = lqr_objective(ts, Q, Ru; Qf=Qf)
        z_opt, α_opt, z_p, N = optimal_control_trajectory(ts, x1, xk1, H, f)

        # Extract first and last state from BlockVector
        @test Array(z_opt[Block(1)]) ≈ x1
        @test Array(z_opt[Block(k + 1)]) ≈ xk1
    end

    # -----------------------------------------------------------------------
    # 4. Returned trajectory satisfies the discrete controlled dynamics
    # -----------------------------------------------------------------------
    @testset "satisfies discrete dynamics" begin
        x1  = [0.0]
        xk1 = [1.0]
        H, f, _ = lqr_objective(ts, Q, Ru; Qf=Qf)
        z_opt, _, _, _ = optimal_control_trajectory(ts, x1, xk1, H, f)

        Ad = ts.Ad
        Bd = ts.Bd
        for t in 1:k
            xt  = Array(z_opt[Block(t)])
            xt1 = Array(z_opt[Block(t + 1)])
            ut  = Array(z_opt[Block(k + 1 + t)])
            @test norm(Ad * xt + Bd * ut - xt1) < 1e-10
        end
    end

    # -----------------------------------------------------------------------
    # 5. Scalar integrator minimum-energy: constant controls equal 1/k,
    #    states are linearly interpolated
    # -----------------------------------------------------------------------
    @testset "minimum-energy scalar integrator" begin
        x1  = [0.0]
        xk1 = [1.0]
        H, f, _ = lqr_objective(ts, Q, Ru; Qf=Qf)
        z_opt, _, _, _ = optimal_control_trajectory(ts, x1, xk1, H, f)

        # All controls equal 1/k
        for t in 1:k
            ut = Array(z_opt[Block(k + 1 + t)])
            @test ut ≈ [1.0 / k] atol=1e-10
        end

        # States are linearly interpolated: x_t = (t-1)/k
        for t in 1:k+1
            xt = Array(z_opt[Block(t)])
            @test xt ≈ [(t - 1) / k] atol=1e-10
        end
    end

    # -----------------------------------------------------------------------
    # 6. Optimizer lies in the affine space: z_opt ≈ z_p + N * α_opt
    # -----------------------------------------------------------------------
    @testset "optimizer lies in affine space" begin
        x1  = [0.0]
        xk1 = [1.0]
        H, f, _ = lqr_objective(ts, Q, Ru; Qf=Qf)
        z_opt, α_opt, z_p, N = optimal_control_trajectory(ts, x1, xk1, H, f)

        @test Array(z_opt) ≈ Array(z_p) + N * α_opt atol=1e-10
    end

    # -----------------------------------------------------------------------
    # 7. If null_basis has 0 columns, function returns z_p
    # -----------------------------------------------------------------------
    @testset "unique trajectory (k=1 has no free controls)" begin
        # For k=1, scalar integrator: z = (x1, x2, u1) with x1 fixed, x2 fixed,
        # dynamics give x2 = x1 + u1, so u1 = xk1 - x1 is determined uniquely.
        k1 = 1
        ts1 = ControlledTrajectorySheaf(F1, Ac_int, Bc_int, h, k1)
        x1  = [0.0]
        xk1 = [1.0]

        p1 = (k1 + 1) * n + k1 * m
        H1 = Matrix{Float64}(I, p1, p1)
        f1 = zeros(p1)

        z_opt1, α_opt1, z_p1, N1 = optimal_control_trajectory(ts1, x1, xk1, H1, f1)

        @test size(N1, 2) == 0
        @test length(α_opt1) == 0
        @test Array(z_opt1) ≈ Array(z_p1)
    end

    # -----------------------------------------------------------------------
    # 8. Wrong-sized Q, Ru, Qf, x_ref, u_ref, H, or f throw ArgumentError
    # -----------------------------------------------------------------------
    @testset "argument validation" begin
        H_ok, f_ok, _ = lqr_objective(ts, Q, Ru; Qf=Qf)

        # Wrong Q size
        @test_throws ArgumentError lqr_objective(ts, zeros(2, 2), Ru; Qf=Qf)
        # Non-symmetric Q (n=1 is always symmetric, so use a 2-state system)
        F2_test = EuclideanSheaf{Float64}(fill(2, 1))
        Ac2_test = zeros(2, 2); Bc2_test = reshape([1.0, 0.0], 2, 1)
        ts2_test = ControlledTrajectorySheaf(F2_test, Ac2_test, Bc2_test, h, 2)
        Q_asym = [1.0 0.5; 0.0 1.0]   # 2×2 but asymmetric
        Ru_1x1 = Matrix{Float64}(I, 1, 1)
        @test_throws ArgumentError lqr_objective(ts2_test, Q_asym, Ru_1x1)
        # Wrong Ru size
        @test_throws ArgumentError lqr_objective(ts, Q, zeros(2, 2); Qf=Qf)
        # Non-symmetric Ru (need m >= 2: MIMO system)
        Bc_mimo = reshape([1.0, 1.0], 1, 2)  # 1 state, 2 controls
        ts_mimo = ControlledTrajectorySheaf(F1, Ac_int, Bc_mimo, h, 2)
        Ru_asym = [1.0 0.5; 0.0 1.0]   # 2×2 but asymmetric
        @test_throws ArgumentError lqr_objective(ts_mimo, Q, Ru_asym)
        # Wrong Qf size
        @test_throws ArgumentError lqr_objective(ts, Q, Ru; Qf=zeros(2, 2))
        # Non-symmetric Qf (use same 2-state system)
        Qf_asym = [1.0 0.3; 0.0 1.0]
        @test_throws ArgumentError lqr_objective(ts2_test, Matrix{Float64}(I, 2, 2), Ru_1x1; Qf=Qf_asym)
        # Non-PD Ru
        Ru_sing = zeros(m, m)
        @test_throws ArgumentError lqr_objective(ts, Q, Ru_sing)
        # Wrong x_ref shape
        @test_throws ArgumentError lqr_objective(ts, Q, Ru; x_ref=zeros(n, k))
        # Wrong u_ref shape
        @test_throws ArgumentError lqr_objective(ts, Q, Ru; u_ref=zeros(m, k + 1))

        # Wrong H size for optimal_control_trajectory
        @test_throws ArgumentError optimal_control_trajectory(ts, [0.0], [1.0], zeros(p + 1, p + 1))
        # Wrong f length
        @test_throws ArgumentError optimal_control_trajectory(ts, [0.0], [1.0], H_ok, zeros(p + 1))
    end

    # -----------------------------------------------------------------------
    # 9. Deliberately unbounded reduced problem throws ArgumentError
    # -----------------------------------------------------------------------
    @testset "unbounded reduced problem throws" begin
        # Use k=3 so null_basis has 2 columns.
        k3 = 3
        ts3 = ControlledTrajectorySheaf(F1, Ac_int, Bc_int, h, k3)
        n3 = ts3.state_dim
        m3 = ts3.control_dim
        p3 = (k3 + 1) * n3 + k3 * m3

        # H = 0 (zero Hessian) with non-zero f creates an unbounded problem.
        H_zero = zeros(Float64, p3, p3)
        f_nonzero = ones(Float64, p3)

        @test_throws ArgumentError optimal_control_trajectory(ts3, [0.0], [1.0], H_zero, f_nonzero)
    end

    # -----------------------------------------------------------------------
    # 10. Reduced first-order optimality condition holds
    # -----------------------------------------------------------------------
    @testset "reduced first-order optimality" begin
        x1  = [0.0]
        xk1 = [1.0]
        H, f, _ = lqr_objective(ts, Q, Ru; Qf=Qf)
        z_opt, α_opt, z_p, N = optimal_control_trajectory(ts, x1, xk1, H, f)

        q_p = Array(z_p)
        Rred = N' * H * N
        rred = N' * (H * q_p + f)
        @test norm(Rred * α_opt + rred) < 1e-10
    end

    # -----------------------------------------------------------------------
    # Additional: lqr_objective with non-zero references shifts f and c
    # -----------------------------------------------------------------------
    @testset "non-zero references shift f and c" begin
        x_ref = 0.5 * ones(n, k + 1)
        u_ref = 0.1 * ones(m, k)
        H, f, c = lqr_objective(ts, Q, Ru; Qf=Qf, x_ref=x_ref, u_ref=u_ref)

        # f = -H * z̄
        z_ref = vcat(vec(x_ref), vec(u_ref))
        @test f ≈ -(H * z_ref)
        @test c ≈ 0.5 * dot(z_ref, H * z_ref)
    end

    # -----------------------------------------------------------------------
    # Additional: 2D double integrator smoke test
    # -----------------------------------------------------------------------
    @testset "2D double integrator smoke test" begin
        Ac2 = [0.0 1.0; 0.0 0.0]
        Bc2 = reshape([0.0, 1.0], 2, 1)
        F2  = EuclideanSheaf{Float64}(fill(2, 1))
        k2  = 4
        h2  = 0.1
        ts2 = ControlledTrajectorySheaf(F2, Ac2, Bc2, h2, k2)

        n2 = ts2.state_dim   # 2
        m2 = ts2.control_dim # 1

        Q2  = Matrix{Float64}(I, n2, n2)
        Ru2 = Matrix{Float64}(I, m2, m2)
        Qf2 = 10.0 * Matrix{Float64}(I, n2, n2)

        H2, f2, _ = lqr_objective(ts2, Q2, Ru2; Qf=Qf2)
        p2 = (k2 + 1) * n2 + k2 * m2

        @test size(H2) == (p2, p2)

        x1_2  = [0.0, 0.0]
        xk1_2 = [1.0, 0.0]
        z_opt2, _, _, _ = optimal_control_trajectory(ts2, x1_2, xk1_2, H2, f2)

        # Endpoint states
        @test Array(z_opt2[Block(1)]) ≈ x1_2
        @test Array(z_opt2[Block(k2 + 1)]) ≈ xk1_2

        # Dynamics satisfied
        Ad2, Bd2 = ts2.Ad, ts2.Bd
        for t in 1:k2
            xt  = Array(z_opt2[Block(t)])
            xt1 = Array(z_opt2[Block(t + 1)])
            ut  = Array(z_opt2[Block(k2 + 1 + t)])
            @test norm(Ad2 * xt + Bd2 * ut - xt1) < 1e-10
        end
    end

    # -----------------------------------------------------------------------
    # Additional: infeasible endpoint pair throws ArgumentError
    # -----------------------------------------------------------------------
    @testset "infeasible endpoints throw ArgumentError" begin
        # The scalar integrator can reach any x_{k+1} given x_1, so use a system
        # with a terminal velocity constraint that cannot be satisfied.
        # Double integrator with k=1 step: x_{k+1} must satisfy
        # x2 = Ad*x1 + Bd*u1. With k=1 there is 1 control dof (u1), so
        # x_pos can be anything but x_vel is constrained.
        # Force an impossible velocity at the terminal state by using k=1
        # and asking for x_vel change that is impossible in 1 step.
        Ac2 = [0.0 1.0; 0.0 0.0]; Bc2 = reshape([0.0, 1.0], 2, 1)
        F2  = EuclideanSheaf{Float64}(fill(2, 1))
        # k=1, h=0.1 gives Ad≈I+h*Ac, Bd≈h*Bc (small h)
        ts_infeas = ControlledTrajectorySheaf(F2, Ac2, Bc2, 0.1, 1)
        # Only 1 control dof, but 2 state constraints. The system IS reachable
        # for compatible endpoints. For k=1, position AND velocity at x_{k+1} are
        # both prescribed but only u1 is free → overdetermined → infeasible for
        # a generic xk1 where the two constraints conflict.
        # Ad*x1 + Bd*u1 = xk1 with n=2, m=1, k=1 is a 2-eq 1-unknown system.
        # For a generic xk1 this will be infeasible.
        x1_inf  = [0.0, 0.0]
        xk1_inf = [1.0, 1.0]   # overdetermined: conflicting target
        @test_throws ArgumentError feasible_control_trajectory_basis(ts_infeas, x1_inf, xk1_inf)
    end

end
