using Test
using CellularSheaves.AgentControllers
using LinearAlgebra

# Tests assert structural and dynamical properties -- conservation, sign,
# convergence, the certificate inequalities -- rather than pinning numbers.
# The point of the certificate is that it is checkable, so the interesting
# tests check it against a real trajectory instead of against a stored value.

x0_for_rate() = se3_pack([0.15, -0.1, 1.42], zeros(3),
                         expm_so3([0.03, -0.02, 0.0]), zeros(3))

@testset "GeometricSE3" begin

    @testset "SO(3) utilities" begin
        a = [0.3, -1.2, 0.7]
        b = [1.0, 2.0, -0.5]
        @test hat_so3(a) * b ≈ cross(a, b)
        @test vee_so3(hat_so3(a)) ≈ a
        @test hat_so3(a)' ≈ -hat_so3(a)

        for w in ([0.0, 0.0, 0.0], [1e-12, 0.0, 0.0], [0.4, -0.2, 1.1], [3.0, 0.0, 0.0])
            R = expm_so3(w)
            @test R' * R ≈ I atol = 1e-12
            @test det(R) ≈ 1.0 atol = 1e-12
        end
        # Rotation about z by pi/2 sends x to y.
        @test expm_so3([0.0, 0.0, pi / 2]) * [1.0, 0.0, 0.0] ≈ [0.0, 1.0, 0.0] atol = 1e-12

        # Projection repairs a perturbed rotation and fixes a clean one.
        R = expm_so3([0.2, 0.3, -0.1])
        @test project_SO3(R) ≈ R atol = 1e-12
        Rp = project_SO3(R .+ 1e-3 .* randn(3, 3))
        @test Rp' * Rp ≈ I atol = 1e-12
        @test det(Rp) ≈ 1.0 atol = 1e-12
    end

    @testset "deriv_unit_vector matches finite differences" begin
        # A smooth curve A(t), differentiated analytically and numerically.
        A(t) = [2.0 + sin(t), 0.5 * t^2, 3.0 - cos(2t)]
        Ad(t) = [cos(t), t, 2 * sin(2t)]
        Add(t) = [-sin(t), 1.0, 4 * cos(2t)]

        t, h = 0.7, 1e-6
        q, qd, qdd = deriv_unit_vector(A(t), Ad(t), Add(t))
        @test norm(q) ≈ 1.0
        qp, _, _ = deriv_unit_vector(A(t + h), Ad(t + h), Add(t + h))
        qm, _, _ = deriv_unit_vector(A(t - h), Ad(t - h), Add(t - h))
        @test qd ≈ (qp .- qm) ./ (2h) atol = 1e-6
        @test qdd ≈ (qp .- 2 .* q .+ qm) ./ h^2 atol = 1e-3

        @test_throws ArgumentError deriv_unit_vector(zeros(3), zeros(3), zeros(3))
    end

    @testset "nonlinear plant" begin
        dyn = SE3QuadrotorDynamics(m = 0.6, J = Matrix(Diagonal([0.012, 0.012, 0.02])))
        @test state_dim(dyn) == 18
        @test position_indices(dyn) == 1:3
        @test velocity_indices(dyn) == 4:6

        x0 = se3_pack([0.0, 0.0, 1.5], zeros(3), Matrix{Float64}(I, 3, 3), zeros(3))
        p, v, R, Om = se3_unpack(x0)
        @test p == [0.0, 0.0, 1.5] && v == zeros(3) && R == I && Om == zeros(3)

        # Hover: thrust exactly cancels weight, nothing moves.
        x1 = se3_step(dyn, x0, dyn.m * dyn.g, zeros(3), 0.05)
        @test x1 ≈ x0 atol = 1e-12

        # Free fall: z-up convention means gravity pulls -z.
        xf = se3_step(dyn, x0, 0.0, zeros(3), 0.1)
        @test xf[3] < 1.5
        @test xf[6] ≈ -dyn.g * 0.1 atol = 1e-10

        # Excess thrust climbs.
        xu = se3_step(dyn, x0, 1.5 * dyn.m * dyn.g, zeros(3), 0.1)
        @test xu[3] > 1.5

        # Attitude stays on the manifold under a sustained moment.
        x = x0
        for _ in 1:200
            x = se3_step(dyn, x, dyn.m * dyn.g, [0.002, -0.001, 0.0005], 0.01)
        end
        Rt = reshape(x[7:15], 3, 3)
        @test Rt' * Rt ≈ I atol = 1e-10
        @test det(Rt) ≈ 1.0 atol = 1e-10

        # Free rigid body conserves rotational kinetic energy and, with a
        # diagonal inertia and no moment, the body-frame angular momentum
        # magnitude.
        xr = se3_pack(zeros(3), zeros(3), Matrix{Float64}(I, 3, 3), [0.5, -0.3, 0.2])
        E0 = 0.5 * dot(xr[16:18], dyn.J * xr[16:18])
        h0 = norm(dyn.J * xr[16:18])
        x = xr
        for _ in 1:500
            x = se3_step(dyn, x, 0.0, zeros(3), 0.001)
        end
        @test 0.5 * dot(x[16:18], dyn.J * x[16:18]) ≈ E0 rtol = 1e-6
        @test norm(dyn.J * x[16:18]) ≈ h0 rtol = 1e-6
    end

    @testset "control law at hover" begin
        dyn = SE3QuadrotorDynamics(m = 0.6)
        ctrl = GeometricSE3Controller()
        x0 = se3_pack([0.2, -0.1, 1.5], zeros(3), Matrix{Float64}(I, 3, 3), zeros(3))
        ref = SE3Reference(x = [0.2, -0.1, 1.5])

        f, M, err = geometric_control(ctrl, dyn, x0, ref)
        @test f ≈ dyn.m * dyn.g atol = 1e-10
        @test norm(M) < 1e-10
        @test err.Rc ≈ I atol = 1e-10
        @test err.Psi < 1e-12
        @test se3_lyapunov(ctrl, dyn, err) < 1e-12

        # Displaced to +x: the commanded attitude must tilt so that body z
        # leans toward -x, i.e. the thrust pushes back.
        xd = se3_pack([0.5, 0.0, 1.5], zeros(3), Matrix{Float64}(I, 3, 3), zeros(3))
        _, _, errd = geometric_control(ctrl, dyn, xd, ref)
        @test errd.Rc[1, 3] < 0
    end

    @testset "closed loop converges from a displaced start" begin
        dyn = SE3QuadrotorDynamics(m = 0.6, J = Matrix(Diagonal([0.012, 0.012, 0.02])))
        ctrl = GeometricSE3Controller(kx = 12.0, kv = 8.0, kR = 10.0, kOmega = 2.0)
        ref = SE3Reference(x = [0.0, 0.0, 1.5])

        x = se3_pack([1.0, -0.8, 0.9], zeros(3), Matrix{Float64}(I, 3, 3), zeros(3))
        dt = 0.002
        e0 = norm(x[1:3] .- ref.x)
        for _ in 1:5000                      # 10 s
            f, M, _ = geometric_control(ctrl, dyn, x, ref)
            x = se3_step(dyn, x, f, M, dt)
        end
        @test norm(x[1:3] .- ref.x) < 0.02 * e0
        @test norm(x[4:6]) < 0.05
    end

    @testset "closed loop tracks a moving reference and actually tilts" begin
        dyn = SE3QuadrotorDynamics(m = 0.6, J = Matrix(Diagonal([0.012, 0.012, 0.02])))
        ctrl = GeometricSE3Controller(kx = 12.0, kv = 8.0, kR = 10.0, kOmega = 2.0)

        # A circle brisk enough to demand real bank angle: r * omega^2 gives a
        # centripetal acceleration of 1.6 m/s^2, about 9 degrees of tilt.
        r, w = 1.0, 1.26
        pos(t) = [r * cos(w * t), r * sin(w * t), 1.5]
        vel(t) = [-r * w * sin(w * t), r * w * cos(w * t), 0.0]
        acc(t) = [-r * w^2 * cos(w * t), -r * w^2 * sin(w * t), 0.0]
        jrk(t) = [r * w^3 * sin(w * t), -r * w^3 * cos(w * t), 0.0]
        snp(t) = [r * w^4 * cos(w * t), r * w^4 * sin(w * t), 0.0]

        dt = 0.002
        x = initial_state(dyn, pos(0.0), vel(0.0), acc(0.0))
        tilt = 0.0
        errs = Float64[]
        for k in 0:6000
            t = k * dt
            ref = SE3Reference(x = pos(t), v = vel(t), a = acc(t), jerk = jrk(t), snap = snp(t))
            f, M, e = geometric_control(ctrl, dyn, x, ref)
            x = se3_step(dyn, x, f, M, dt)
            t > 4 && push!(errs, norm(e.ex))
            R = reshape(x[7:15], 3, 3)
            tilt = max(tilt, acos(clamp(R[3, 3], -1.0, 1.0)))
        end
        @test maximum(errs) < 0.05
        # The vehicle is genuinely banked: a hover linearization would not be
        # exercised by a run that never leaves R[3,3] ~ 1.
        @test tilt > deg2rad(5)
    end

    @testset "certificate is feasible and is respected by a trajectory" begin
        dyn = SE3QuadrotorDynamics(m = 0.6, J = Matrix(Diagonal([0.012, 0.012, 0.02])))
        region, psi1 = 0.3, 0.05
        base = GeometricSE3Controller(kx = 12.0, kv = 8.0, kR = 10.0, kOmega = 2.0)

        tuned = se3_tune_lyapunov(base, dyn; psi1, ex_max = region, ev_max = region, n = 30)
        cert = se3_certificate(tuned, dyn; psi1, ex_max = region, ev_max = region)
        @test cert.feasible
        @test isempty(cert.violations)
        @test cert.alpha1 > 0
        @test cert.alpha2 >= cert.alpha1
        @test cert.c3 > 0
        # Tuning does not touch the control law, only the certificate.
        @test tuned.kx == base.kx && tuned.kR == base.kR && tuned.kOmega == base.kOmega

        # A soft attitude loop is the failure mode worth naming: the vehicle
        # still flies, but no region is certifiable, so the composite argument
        # has nothing to stand on.
        soft = GeometricSE3Controller(kx = 12.0, kv = 8.0, kR = 1.5, kOmega = 0.35)
        softtuned = se3_tune_lyapunov(soft, dyn; psi1, ex_max = region, ev_max = region, n = 20)
        softcert = se3_certificate(softtuned, dyn; psi1, ex_max = region, ev_max = region)
        @test !softcert.feasible
        @test !isempty(softcert.violations)
        @test softcert.c3 == 0.0

        # Lee's closed-form bounds on c1, c2 are advisory, not binding: a
        # certificate outside them is still a certificate.
        wide = se3_tune_lyapunov(base, dyn; psi1, ex_max = region, ev_max = region,
                                 n = 20, widen = 8.0)
        widecert = se3_certificate(wide, dyn; psi1, ex_max = region, ev_max = region)
        @test widecert.feasible == isempty(widecert.violations)

        # An infeasible region is reported as such rather than silently
        # returning a meaningless rate.
        toobig = se3_certificate(tuned, dyn; psi1, ex_max = 50.0, ev_max = 50.0)
        @test !toobig.feasible
        @test toobig.c3 == 0.0

        # The inequalities the certificate asserts, checked along a real
        # trajectory started inside the region.
        ref = SE3Reference(x = [0.0, 0.0, 1.5])
        x = se3_pack([0.15, -0.1, 1.42], zeros(3), expm_so3([0.03, -0.02, 0.0]), zeros(3))
        dt = 0.001
        Vs = Float64[]
        for _ in 1:4000
            f, M, e = geometric_control(tuned, dyn, x, ref)
            z = se3_error_vector(e)
            V = se3_lyapunov(tuned, dyn, e)
            if e.Psi < cert.psi1 && norm(e.ex) < cert.region_ex_max
                @test cert.alpha1 * dot(z, z) <= V + 1e-9
                @test V <= cert.alpha2 * dot(z, z) + 1e-9
            end
            push!(Vs, V)
            x = se3_step(dyn, x, f, M, dt)
        end
        # V decays at least as fast as the certified rate: the certificate is
        # a sufficient bound, so the realized decay must be no slower.
        T = (length(Vs) - 1) * dt
        @test Vs[end] <= Vs[1] * exp(-cert.c3 * T) + 1e-8
        @test Vs[end] < Vs[1]

        # The measured rate is positive and, this bound being loose, well
        # above the certified one.
        rate, Vm = se3_measured_rate(tuned, dyn, x0_for_rate(), ref, dt, 3000)
        @test length(Vm) == 3000
        @test rate > 0
        @test rate > cert.c3
    end

    @testset "SE3AgentState in the layered stack" begin
        dyn = SE3QuadrotorDynamics(m = 0.6, J = Matrix(Diagonal([0.012, 0.012, 0.02])))
        ctrl = GeometricSE3Controller(kx = 12.0, kv = 8.0, kR = 10.0, kOmega = 2.0)

        # Start on the "airstrip", well away from the commanded station, so
        # the transient into formation is what the rollout shows.
        x0 = initial_state(dyn, [1.5, 0.0, 1.5], zeros(3), zeros(3))
        agent = SE3AgentState(x0, dyn, ctrl, 0.05)

        # A 4-dimensional stalk (SE(3) homogeneous coordinates) passes through
        # unmodified; the trailing 1 must be ignored, not tracked.
        qstar = [0.0, 0.6, 1.5, 1.0]
        qdot = [0.0, 0.0, 0.0, 0.0]
        qddot = [0.0, 0.0, 0.0, 0.0]

        local xk, rk
        dt = 0.005
        for _ in 1:3000
            xk, rk = step_agent!(agent, qstar, qdot, qddot, dt)
        end
        @test length(xk) == 18
        @test length(rk) == 3
        @test rk ≈ qstar[1:3] atol = 1e-6          # filter converged to the reference
        @test norm(xk[1:3] .- qstar[1:3]) < 0.02   # vehicle converged to the filter
        @test agent.last isa SE3Errors

        # The shorter call forms work.
        a2 = SE3AgentState(x0, dyn, ctrl, 0.05)
        @test length(first(step_agent!(a2, qstar, dt))) == 18
        @test length(first(step_agent!(a2, qstar, qdot, dt))) == 18
    end

    @testset "agrees with the hover linearization near hover" begin
        # Same vehicle, both plants, a gentle reference. The nonlinear model
        # and the linearized one should not disagree materially when the
        # linearization's own hypotheses hold; if they do, the frame or sign
        # convention in the nonlinear model is wrong.
        m, Ixx = 0.6, 0.012
        nl = SE3QuadrotorDynamics(m = m, J = Matrix(Diagonal([Ixx, Ixx, 0.02])))
        ctrl = GeometricSE3Controller(kx = 12.0, kv = 8.0, kR = 10.0, kOmega = 2.0)

        lin = QuadrotorDynamics(m = m, Ixx = Ixx, Iyy = Ixx)
        Q = Matrix(Diagonal([500.0, 500.0, 500.0, 150.0, 150.0,
                             100.0, 100.0, 100.0, 5.0, 5.0]))
        Rw = Matrix(Diagonal([0.005, 0.005, 0.005]))
        lqr = LQRController(lin, 0.005, Q, Rw)

        goal = [0.15, -0.1, 1.5]
        start = [0.0, 0.0, 1.5]

        xn = initial_state(nl, start, zeros(3), zeros(3))
        xl = initial_state(lin, start, zeros(3))
        xl_ref = initial_state(lin, goal, zeros(3))
        Ad, Bd = AgentControllers.discrete_matrices(lin, 0.005)

        for _ in 1:2000
            f, M, _ = geometric_control(ctrl, nl, xn, SE3Reference(x = goal))
            xn = se3_step(nl, xn, f, M, 0.005)
            xl = Ad * xl .+ Bd * (-lqr.K * (xl .- xl_ref))
        end
        @test norm(xn[1:3] .- goal) < 0.02
        @test norm(xl[1:3] .- goal) < 0.02
        @test norm(xn[1:3] .- xl[1:3]) < 0.05
    end

    @testset "per-vehicle gains keep a heterogeneous fleet certifiable" begin
        # The six escort vehicles: mass and inertia both vary. Fixed attitude
        # gains detune the heavy end of the fleet out of the certifiable set;
        # gains specified by closed-loop attitude response do not.
        fleet = [SE3QuadrotorDynamics(m = 0.5 + 0.05i,
                                      J = Matrix(Diagonal([0.01 + 0.002i,
                                                           0.01 + 0.002i,
                                                           0.02 + 0.004i])))
                 for i in 1:6]
        region, psi1 = 0.3, 0.05

        # Same omega_n for everyone => same closed-loop attitude response.
        scaled = [GeometricSE3Controller(d) for d in fleet]
        certs = [last(se3_autotune(c, d; psi1, ex_max = region, ev_max = region, n = 24))
                 for (c, d) in zip(scaled, fleet)]
        @test all(c -> c.feasible, certs)
        @test all(c -> c.c3 > 0, certs)
        # Attitude gains track inertia, as the constructor's derivation says.
        @test scaled[end].kR > scaled[1].kR

        # Assumption 1's fleet constants are the min/max across agents. That
        # they exist and are finite is the whole content of the hypothesis for
        # a heterogeneous fleet.
        alpha1 = minimum(c.alpha1 for c in certs)
        alpha2 = maximum(c.alpha2 for c in certs)
        c3 = minimum(c.c3 for c in certs)
        @test 0 < alpha1 <= alpha2
        @test c3 > 0

        # One shared attitude gain, tuned for the lightest vehicle, does not
        # survive the heaviest.
        fixed = GeometricSE3Controller(kx = 12.0, kv = 8.0, kR = 10.0, kOmega = 2.0)
        fixed_certs = [last(se3_autotune(fixed, d; psi1, ex_max = region, ev_max = region, n = 24))
                       for d in fleet]
        @test !all(c -> c.feasible, fixed_certs)
    end

    @testset "se3_autotune reports rather than hides infeasibility" begin
        dyn = SE3QuadrotorDynamics(m = 0.6, J = Matrix(Diagonal([0.012, 0.012, 0.02])))
        soft = GeometricSE3Controller(kx = 12.0, kv = 8.0, kR = 1.0, kOmega = 0.2)
        _, cert = se3_autotune(soft, dyn; psi1 = 0.05, ex_max = 1.0, ev_max = 1.0, n = 16)
        @test !cert.feasible
        @test !isempty(cert.violations)
    end

    @testset "analytic feedforward removes the filter's lag" begin
        dyn = SE3QuadrotorDynamics(m = 0.6, J = Matrix(Diagonal([0.012, 0.012, 0.02])))
        ctrl = GeometricSE3Controller(dyn)
        eps, dt, speed = 0.2, 0.005, 0.35
        station(t) = [speed * t, 0.0, 1.5]
        statvel(t) = [speed, 0.0, 0.0]

        function tail_error(ff)
            x0 = initial_state(dyn, station(0.0), statvel(0.0), zeros(3))
            a = SE3AgentState(x0, dyn, ctrl, eps; feedforward = ff)
            errs = Float64[]
            for k in 0:2999
                t = k * dt
                xk, _ = step_agent!(a, station(t), statvel(t), zeros(3), dt)
                t > 10 && push!(errs, norm(xk[1:3] .- station(t)))
            end
            return sum(errs) / length(errs)
        end

        plain, ff = tail_error(false), tail_error(true)
        # Without feedforward the lag is set by the filter time constant and
        # the reference speed; with it, the forcing term cancels.
        @test plain > 0.5 * eps * speed
        @test ff < 0.2 * plain
    end
end
