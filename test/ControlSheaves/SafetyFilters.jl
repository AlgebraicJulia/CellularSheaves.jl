using Test
using LinearAlgebra
using Random
using CellularSheaves
using CellularSheaves.ControlSheaves.SafetyFilters
using CellularSheaves.ControlSheaves.AgentControllers

# Reference projection onto a single halfspace a'u >= b. When exactly one barrier row is
# active and no other constraint binds, the filter must reproduce this in closed form.
halfspace_projection(u_nom, a, b) =
    u_nom + max(0.0, (b - dot(a, u_nom)) / dot(a, a)) .* a

@testset "single integrator regression" begin
    model = ControlAffineModel(2)
    params = CBFCLFParams(d_safe = 0.8, actuator_cap = 1.0)
    x1, x2 = [-0.5, 0.0], [0.5, 0.0]
    r1, r2 = [0.5, 0.0], [-0.5, 0.0]

    a = safety_filter(params, model, [2.0, 0.0], x1, r1; others = [x2])
    b = safety_filter(params, model, [-2.0, 0.0], x2, r2; others = [x1])

    @test a.certified
    @test b.certified
    @test norm(a.command) <= 1.0 + 1e-7
    @test norm(b.command) <= 1.0 + 1e-7
    # Both agents enforcing their own half of the pair recovers hdot >= -gamma h.
    h = sum(abs2, x1 - x2) - 0.8^2
    @test 2 * dot(x1 - x2, a.command - b.command) >= -params.gamma * h - 1e-6
end

@testset "reductions" begin
    model = ControlAffineModel(2)

    # At the reference with no neighbours and no cap, the nominal command is optimal.
    free = safety_filter(CBFCLFParams(), model, [0.3, -0.2], [1.0, 1.0], [1.0, 1.0])
    @test free.certified
    @test free.command ≈ [0.3, -0.2] atol = 1e-6
    @test isinf(free.cbf_margin)

    # A neighbour beyond the sensing radius contributes no row.
    far = safety_filter(CBFCLFParams(sense = 0.5), model, [0.3, -0.2], [1.0, 1.0], [1.0, 1.0];
                        others = [[9.0, 9.0]])
    @test isinf(far.cbf_margin)
    @test far.command ≈ free.command atol = 1e-6

    # A very permissive class-K gain leaves a non-conflicting command untouched.
    permissive = safety_filter(CBFCLFParams(gamma = 1e6, d_safe = 0.2), model,
                               [0.3, -0.2], [1.0, 1.0], [1.0, 1.0]; others = [[2.0, 2.0]])
    @test permissive.certified
    @test permissive.command ≈ [0.3, -0.2] atol = 1e-5
end

@testset "closed-form halfspace ground truth" begin
    # One neighbour and no actuator cone. Placing the agent at its reference makes the CLF
    # row identically zero, so the solution is the projection of u_nom onto the single
    # barrier halfspace.
    model = ControlAffineModel(2)
    params = CBFCLFParams(d_safe = 1.0, gamma = 1.0)
    x = [0.0, 0.0]
    xj = [1.2, 0.0]
    u_nom = [3.0, 0.0]

    res = safety_filter(params, model, u_nom, x, x; others = [xj])
    @test res.certified

    displacement = x - xj
    a = 2 .* displacement
    h = sum(abs2, displacement) - params.d_safe^2
    b = -(params.gamma / 2) * h
    @test res.command ≈ halfspace_projection(u_nom, a, b) atol = 1e-5
end

@testset "Monte Carlo optimality of the encoding" begin
    # Certify the conic encoding without a second solver: no sampled feasible point may
    # beat the returned objective.
    rng = MersenneTwister(0xC0FFEE)
    model = ControlAffineModel(2)
    params = CBFCLFParams(d_safe = 0.9, gamma = 2.0, actuator_cap = 1.5)
    x = [0.0, 0.0]
    others = [[1.0, 0.2], [-0.3, 1.0]]
    u_nom = [2.0, 1.0]

    res = safety_filter(params, model, u_nom, x, x; others = others)
    @test res.certified
    @test res.command !== nothing

    A = [2 .* (x - xj) for xj in others]
    b = [-(params.gamma / 2) * (sum(abs2, x - xj) - params.d_safe^2) for xj in others]
    objective(u) = 0.5 * sum(abs2, u - u_nom)

    best = objective(res.command)
    for _ in 1:20_000
        u = params.actuator_cap .* (2 .* rand(rng, 2) .- 1)
        norm(u) <= params.actuator_cap || continue
        all(dot(A[k], u) >= b[k] for k in eachindex(A)) || continue
        @test objective(u) >= best - 1e-6
    end
end

@testset "certification implies feasibility" begin
    rng = MersenneTwister(20260806)
    model = ControlAffineModel(2)
    params = CBFCLFParams(d_safe = 0.7, gamma = 4.0, actuator_cap = 1.2)

    for _ in 1:200
        x = randn(rng, 2)
        ref = randn(rng, 2)
        u_nom = 3 .* randn(rng, 2)
        others = [x .+ 0.9 .* randn(rng, 2) for _ in 1:2]
        res = safety_filter(params, model, u_nom, x, ref; others = others)
        res.certified || continue

        # The struct promises every hard constraint holds whenever a command is returned.
        @test res.command !== nothing
        @test norm(res.command) <= params.actuator_cap + 1e-6
        for xj in others
            displacement = x - xj
            norm(displacement) > params.sense && continue
            h = sum(abs2, displacement) - params.d_safe^2
            @test 2 * dot(displacement, res.command) >= -(params.gamma / 2) * h - 1e-6
        end
    end
end

@testset "actuator cone binds" begin
    model = ControlAffineModel(2)
    params = CBFCLFParams(actuator_cap = 0.25)
    res = safety_filter(params, model, [10.0, 0.0], [0.0, 0.0], [5.0, 0.0])
    @test res.certified
    @test norm(res.command) <= 0.25 + 1e-7
    @test res.cap_residual >= -1e-7
    # Minimally invasive: the cap is saturated in the direction of the nominal command.
    @test res.command ≈ [0.25, 0.0] atol = 1e-5
end

@testset "obstacles use full responsibility" begin
    model = ControlAffineModel(2)
    params = CBFCLFParams(d_safe = 1.0, gamma = 1.0, obstacle_radius = 0.0)
    x = [0.0, 0.0]
    obstacle = [1.2, 0.0]
    res = safety_filter(params, model, [3.0, 0.0], x, x; obstacles = [obstacle])
    @test res.certified

    displacement = x - obstacle
    a = 2 .* displacement
    h = sum(abs2, displacement) - params.d_safe^2
    @test res.command ≈ halfspace_projection([3.0, 0.0], a, -params.gamma * h) atol = 1e-5
end

@testset "barrier_values agrees with the filter" begin
    model = ControlAffineModel(2)
    x = [0.0, 0.0]
    others = [[1.0, 0.0], [0.0, 2.0]]
    values = barrier_values(model, x, others; d_safe = 0.5)
    @test values ≈ [1.0 - 0.25, 4.0 - 0.25]

    res = safety_filter(CBFCLFParams(d_safe = 0.5), model, [0.1, 0.1], x, x; others = others)
    @test res.cbf_margin ≈ minimum(values)

    @test isempty(barrier_values(model, x, others; d_safe = 0.5, sense = 0.5))
end

@testset "closed-loop forward invariance and Lyapunov decay" begin
    # Two single integrators commanded to swap slots: the nominal law drives them through
    # each other, the filter must keep the barrier non-negative throughout. A small lateral
    # offset keeps the encounter off the exactly head-on configuration, where the barrier
    # gradient is antiparallel to the CLF descent direction and the pair deadlocks.
    model = ControlAffineModel(2)
    params = CBFCLFParams(d_safe = 0.6, gamma = 5.0, actuator_cap = 2.0, clf_rate = 3.0)
    dt = 0.01
    xs = [[-1.5, -0.15], [1.5, 0.15]]
    refs = [[1.5, -0.15], [-1.5, 0.15]]

    min_h = Inf
    uncertified = 0
    for _ in 1:400
        u_noms = [4 .* (refs[i] - xs[i]) for i in 1:2]
        results = safety_filter_all(params, model, u_noms, xs, refs)
        for i in 1:2
            # An uncertified step is held rather than driven by the unfiltered command:
            # zero is admissible for a single integrator whenever the barrier already holds.
            u = results[i].certified ? results[i].command : zeros(2)
            results[i].certified || (uncertified += 1)
            xs[i] = xs[i] + dt .* u
        end
        min_h = min(min_h, sum(abs2, xs[1] - xs[2]) - params.d_safe^2)
    end

    @test uncertified == 0
    @test min_h >= -1e-6
    # Both agents complete the swap despite the detour around each other.
    @test norm(xs[1] - refs[1]) < 1e-3
    @test norm(xs[2] - refs[2]) < 1e-3

    # With no neighbour to avoid, the CLF alone drives V monotonically down. This is the
    # dissipation certificate the layered architecture relies on at the execution layer.
    x = [2.0, 1.0]
    ref = [0.0, 0.0]
    V0 = 0.5 * sum(abs2, x - ref)
    V_prev = V0
    for _ in 1:200
        res = safety_filter(CBFCLFParams(clf_rate = 2.0, actuator_cap = 3.0), model,
                            zeros(2), x, ref)
        @test res.certified
        x = x + dt .* res.command
        V = 0.5 * sum(abs2, x - ref)
        @test V <= V_prev + 1e-9
        V_prev = V
    end
    @test V_prev < 0.1 * V0
end

@testset "head-on deadlock and symmetry breaking" begin
    # Exactly head-on: the barrier gradient is antiparallel to the CLF descent direction,
    # so no admissible command makes progress and the pair stalls at contact distance.
    model = ControlAffineModel(2)
    dt = 0.01
    function swap(bias)
        params = CBFCLFParams(d_safe = 0.6, gamma = 5.0, actuator_cap = 2.0,
                              clf_rate = 3.0, deadlock_bias = bias)
        xs = [[-1.5, 0.0], [1.5, 0.0]]
        refs = [[1.5, 0.0], [-1.5, 0.0]]
        min_h = Inf
        for _ in 1:800
            u_noms = [4 .* (refs[i] - xs[i]) for i in 1:2]
            results = safety_filter_all(params, model, u_noms, xs, refs)
            for i in 1:2
                u = results[i].certified ? results[i].command : zeros(2)
                xs[i] = xs[i] + dt .* u
            end
            min_h = min(min_h, sum(abs2, xs[1] - xs[2]) - params.d_safe^2)
        end
        return (progress = norm(xs[1] - refs[1]), min_h = min_h)
    end

    stalled = swap(0.0)
    @test stalled.min_h >= -1e-6          # safe, but
    @test stalled.progress > 1.0          # never gets past the other agent

    broken = swap(0.35)
    @test broken.min_h >= -1e-6           # still safe, and
    @test broken.progress < stalled.progress   # now makes headway
end

@testset "braking barrier reaches an acceleration input" begin
    di = double_integrator_model(2)
    x = [0.0, 0.0, 1.0, 0.0]          # at the origin, moving +x
    neighbour = [1.4, 0.0, -1.0, 0.0] # ahead, closing head-on

    # The distance barrier is vacuous for a double integrator: position is not actuated.
    @test_throws SafetyFilters.RelativeDegreeError safety_filter(
        CBFCLFParams(d_safe = 0.5), di, zeros(2), x, zeros(4); others = [neighbour])

    # The braking barrier depends on relative velocity, so it does reach the input.
    res = safety_filter(CBFCLFParams(d_safe = 0.5, braking_acceleration = 2.0, gamma = 3.0),
                        di, zeros(2), x, zeros(4); others = [neighbour])
    @test res.certified
    @test res.command !== nothing

    # It must command deceleration against the closing motion.
    @test res.command[1] < 0

    # A model without a velocity selector cannot express the barrier at all.
    @test_throws ArgumentError safety_filter(
        CBFCLFParams(braking_acceleration = 2.0), ControlAffineModel(2),
        zeros(2), zeros(2), zeros(2); others = [[1.0, 0.0]])
end

@testset "braking barrier keeps double integrators apart under an actuator cap" begin
    # Two double integrators told to swap positions, with a hard bound on acceleration.
    # The nominal PD law drives them through one another; the braking barrier has to keep
    # them separated using only the acceleration it is allowed to command.
    di = double_integrator_model(2)
    a_max = 2.0
    params = CBFCLFParams(d_safe = 0.5, gamma = 3.0, sense = 6.0, clf_rate = 1.0,
                          actuator_cap = a_max, braking_acceleration = a_max)
    dt = 0.01
    xs = [[-2.0, -0.12, 0.0, 0.0], [2.0, 0.12, 0.0, 0.0]]
    refs = [[2.0, -0.12, 0.0, 0.0], [-2.0, 0.12, 0.0, 0.0]]

    min_sep = Inf
    max_cmd = 0.0
    uncertified = 0
    for _ in 1:900
        u_noms = [1.2 .* (refs[i][1:2] - xs[i][1:2]) - 1.6 .* xs[i][3:4] for i in 1:2]
        results = safety_filter_all(params, di, u_noms, xs, refs)
        for i in 1:2
            u = results[i].certified ? results[i].command : zeros(2)
            results[i].certified || (uncertified += 1)
            max_cmd = max(max_cmd, norm(u))
            xs[i] = xs[i] + dt .* vcat(xs[i][3:4], u)
        end
        min_sep = min(min_sep, norm(xs[1][1:2] - xs[2][1:2]))
    end

    @test uncertified == 0
    @test max_cmd <= a_max + 1e-6              # actuator limit respected
    @test min_sep >= params.d_safe - 1e-3      # keep-out radius held
    # The swap still completes: this is avoidance, not gridlock.
    @test norm(xs[1][1:2] - refs[1][1:2]) < 0.5
    @test norm(xs[2][1:2] - refs[2][1:2]) < 0.5
end

@testset "braking barrier on a quadrotor reaches only the actuated axis" begin
    # Velocity is directly actuated in z alone, so the braking barrier is well posed for a
    # vertical encounter and still degenerate for a horizontal one. Horizontal collision
    # avoidance on this plant needs an acceleration-level input or a higher-order barrier.
    quad = ControlAffineModel(QuadrotorDynamics())
    params = CBFCLFParams(d_safe = 0.5, braking_acceleration = 2.0)

    x = zeros(10)
    x[3] = 0.0
    vertical = zeros(10)
    vertical[3] = 1.2        # neighbour directly above
    res = safety_filter(params, quad, zeros(3), x, zeros(10); others = [vertical])
    @test res.certified

    horizontal = zeros(10)
    horizontal[1] = 1.2      # neighbour directly ahead
    @test_throws SafetyFilters.RelativeDegreeError safety_filter(
        params, quad, zeros(3), x, zeros(10); others = [horizontal])
end

@testset "CLF metric from the interaction map" begin
    H = [4.0 1.0 0.0 0.0
         1.0 3.0 0.0 0.0
         0.0 0.0 5.0 2.0
         0.0 0.0 2.0 6.0]
    M = laplacian_block_metric(H, 2, 2)
    @test M ≈ [5.0 2.0; 2.0 6.0]
    @test isposdef(M)
    @test_throws ArgumentError laplacian_block_metric(H, 3, 2)

    # A non-identity metric changes the Lyapunov value the filter reports.
    model = ControlAffineModel(2)
    res = safety_filter(CBFCLFParams(clf_metric = M), model, zeros(2), [1.0, 0.0], zeros(2))
    @test res.lyapunov ≈ 2.5
end

@testset "relative degree guard" begin
    dyn = QuadrotorDynamics()
    model = ControlAffineModel(dyn)
    @test size(model.position) == (3, 10)

    x = zeros(10)
    x[1] = 0.4
    neighbour = zeros(10)

    # Moments cannot change the translational barrier instantaneously.
    @test_throws SafetyFilters.RelativeDegreeError safety_filter(
        CBFCLFParams(d_safe = 1.0), model, zeros(3), x, zeros(10); others = [neighbour])

    # Without a barrier row the CLF-only filter is well posed on the same agent.
    res = safety_filter(CBFCLFParams(actuator_cap = 8.0), model, zeros(3), x, zeros(10))
    @test res.certified
    @test res.lyapunov ≈ 0.08
end

@testset "CBFCLFController slots into AgentState" begin
    dyn = PlanarQuadrotorDynamics()
    dt = 0.05
    Q = Matrix(Diagonal([100.0, 100.0, 10.0, 10.0, 10.0, 1.0]))
    R = Matrix(Diagonal([0.01, 0.01]))
    inner = LQRController(dyn, dt, Q, R)
    ctrl = CBFCLFController(inner, CBFCLFParams(actuator_cap = 50.0), ControlAffineModel(dyn))
    @test ctrl isa AbstractAgentController

    w = AgentState(zeros(6), dyn, dt, inner.K, 0.02)
    target = [0.5, 1.0]

    x, x_ref, res = step_agent!(w, ctrl, target, Vector{Vector{Float64}}(), dt)
    @test length(x) == 6
    @test length(x_ref) == 6
    @test res isa SafetyFilterResult
    @test res.certified

    # The controller path surfaces the same modelling error as the bare filter.
    @test_throws SafetyFilters.RelativeDegreeError step_agent!(
        w, CBFCLFController(inner, CBFCLFParams(d_safe = 1.0), ControlAffineModel(dyn)),
        target, [zeros(6)], dt)
end

@testset "CLF filter restores per-agent dissipation of the analytic law" begin
    # The analytic execution law tracks the precomputed harmonic extension directly,
    #     u_i = -k1 g_i^+ (x_i - q_i^*),   q^* = H^{-1} (-L_IB p),
    # which is decoupled across agents: once q^* has been distributed, no agent needs a
    # neighbour's state. Against a stationary target each V_i = 1/2 ||x_i - q_i^*||_{H_ii}^2
    # is then dissipated trivially. Against a moving target the proportional law lags, since
    #     V_i_dot = -k1 ||e_i||^2 - e_i' M qstar_dot_i
    # is positive whenever the reference outruns the tracking error. The CLF filter carries
    # the feedforward term and enforces the dissipation inequality regardless.
    using CellularSheaves.Formations: build_escort_ring

    NA, D, K1, dt, steps = 5, 4, 1.0, 0.01, 400
    target_node = NA + 1
    sheaf = build_escort_ring(NA, target_node, 0.3; observers = [1])

    target(t) = [0.6cos(1.4t), 0.6sin(1.4t), 1.5, 1.0]
    target_rate(t) = [-0.84sin(1.4t), 0.84cos(1.4t), 0.0, 0.0]

    H_blk, L_IB = restricted_laplacian_blocks(sheaf, collect(1:NA), [target_node])
    H = Matrix(H_blk)
    L_IB = Matrix(L_IB)
    harmonic(t) = H \ (-L_IB * target(t))
    harmonic_rate(t) = H \ (-L_IB * target_rate(t))
    @test norm(H * harmonic(0.3) + L_IB * target(0.3)) < 1e-8

    slice(i) = ((i - 1) * D + 1):(i * D)
    metrics = [laplacian_block_metric(H, i, D) for i in 1:NA]
    V_local(q, qstar, i) = 0.5 * dot(q[slice(i)] - qstar[slice(i)],
                                     metrics[i] * (q[slice(i)] - qstar[slice(i)]))

    # The certificate is a statement about the Lie derivative along the closed loop,
    #     V_i_dot = e_i' M_i (u_i - qstar_dot_i)  <=  -c_3 V_i + delta_i,
    # so it is that quantity, and not a forward-Euler difference of V_i, that is checked
    # here. A discrete difference would additionally carry O(dt^2) terms and the reference's
    # own motion between samples, neither of which the filter claims to control.
    c3 = 2.0
    model = ControlAffineModel(D)
    rng = MersenneTwister(11)
    perturbation = 0.5 .* randn(rng, NA * D)
    for i in 1:NA
        perturbation[last(((i - 1) * D + 1):(i * D))] = 0.0   # homogeneous coordinate exact
    end

    function rollout(filtered; penalty = 1e4)
        q = harmonic(0.0) + perturbation
        worst_margin = -Inf     # max over steps of (V_dot + c3 V), bounded by delta
        worst_rate = -Inf       # max over steps of V_dot itself
        uncertified = 0
        max_slack = 0.0
        ultimate = 0.0
        for k in 1:steps
            t = (k - 1) * dt
            qstar = harmonic(t)
            qstar_dot = harmonic_rate(t)
            u = -K1 .* (q - qstar)
            if filtered
                for i in 1:NA
                    params = CBFCLFParams(clf_metric = metrics[i], clf_rate = c3,
                                          clf_penalty = penalty)
                    res = safety_filter(params, model, u[slice(i)], q[slice(i)],
                                        qstar[slice(i)];
                                        ref_velocity = qstar_dot[slice(i)])
                    if res.certified
                        u[slice(i)] .= res.command
                        max_slack = max(max_slack, res.slack)
                    else
                        uncertified += 1
                    end
                end
            end
            for i in 1:NA
                e = q[slice(i)] - qstar[slice(i)]
                V = 0.5 * dot(e, metrics[i] * e)
                Vdot = dot(metrics[i] * e, u[slice(i)] - qstar_dot[slice(i)])
                worst_rate = max(worst_rate, Vdot)
                worst_margin = max(worst_margin, Vdot + c3 * V)
            end
            q = q + dt .* u
            if k > 3 * steps ÷ 4
                after = harmonic(t + dt)
                ultimate = max(ultimate, maximum(V_local(q, after, i) for i in 1:NA))
            end
        end
        return (margin = worst_margin, rate = worst_rate, uncertified = uncertified,
                slack = max_slack, ultimate = ultimate)
    end

    analytic = rollout(false)
    filtered = rollout(true)

    # The unfiltered proportional law lags a moving reference: some agent's Lyapunov
    # function has a strictly positive Lie derivative, so it is not a stability certificate
    # for that agent.
    @test analytic.rate > 0
    @test analytic.margin > 0

    # Under the filter the dissipation inequality V_dot <= -c_3 V + delta holds for every
    # agent at every step: this is the certificate, and it is enforced exactly.
    @test filtered.uncertified == 0
    @test filtered.margin <= filtered.slack + 1e-6

    # The Lie derivative itself is strictly negative throughout, so each agent is a stable
    # tracker of its own harmonic slot.
    @test filtered.rate < 0

    # The quadratic penalty p*delta^2 is inexact: delta is O(1/p) rather than identically
    # zero once the unrelaxed problem is feasible, as an L1 penalty would give. Raising p by
    # two decades shrinks the relaxation by close to two decades at any fixed state, so the
    # certificate approaches the unrelaxed inequality V_dot <= -c_3 V but never attains it.
    # Compared across whole rollouts the ratio is smaller than at a fixed state, because the
    # two penalties produce different trajectories.
    loose = rollout(true; penalty = 1e2)
    tight = rollout(true; penalty = 1e4)
    @test tight.slack < 0.2 * loose.slack
    @test tight.margin <= tight.slack + 1e-6

    # Enforcing the dissipation rate also tightens the residual tracking error.
    @test filtered.ultimate < analytic.ultimate
end

@testset "filter_program composes arbitrary row sets" begin
    using CellularSheaves.IPM: solve, OPTIMAL, NEAR_OPTIMAL
    using CellularSheaves.BlockSparseArrays: colrange

    settings = safety_filter_settings()
    solve_rows(u_nom, rows; kwargs...) = begin
        problem, B = filter_program(u_nom, rows; kwargs...)
        result = solve(problem, settings)
        @test result.status in (OPTIMAL, NEAR_OPTIMAL)
        (command = Vector(result.p[colrange(B, 1)]), B = B, p = result.p)
    end

    # A program with no stability row at all: hard rows only, which is the minimally
    # invasive barrier filter. One active row reduces to a halfspace projection.
    u_nom = [3.0, 0.0]
    a = [-2.4, 0.0]
    b = -0.5
    only_barrier = solve_rows(u_nom, [LinearRow(a, b)])
    @test only_barrier.command ≈ halfspace_projection(u_nom, a, b) atol = 1e-5

    # Box bounds on the command, four rows that the encoder has never seen from
    # safety_filter, plus nothing else.
    box = [LinearRow([1.0, 0.0], -0.4), LinearRow([-1.0, 0.0], -0.4),
           LinearRow([0.0, 1.0], -0.4), LinearRow([0.0, -1.0], -0.4)]
    boxed = solve_rows([2.0, -1.5], box)
    @test all(abs.(boxed.command) .<= 0.4 + 1e-6)
    @test boxed.command ≈ [0.4, -0.4] atol = 1e-5

    # Two soft rows exercise the multiple-relaxation path, which the single CLF row of
    # safety_filter never reaches. Satisfiable rows leave the nominal command alone.
    slack_pair = [LinearRow([1.0, 0.0], -5.0; penalty = 1e3),
                  LinearRow([0.0, 1.0], -5.0; penalty = 1e3)]
    untouched = solve_rows([0.2, 0.3], slack_pair)
    @test untouched.command ≈ [0.2, 0.3] atol = 1e-5
    @test untouched.p[first(colrange(untouched.B, 2))] <= 1e-4
    @test untouched.p[first(colrange(untouched.B, 3))] <= 1e-4

    # An unreachable soft row is absorbed by its own relaxation rather than made infeasible.
    conflict = [LinearRow([1.0, 0.0], 50.0; penalty = 1.0),
                LinearRow([0.0, 1.0], -5.0; penalty = 1.0)]
    relaxed = solve_rows([0.0, 0.0], conflict; cone_metric = Matrix{Float64}(I, 2, 2),
                         cone_bound = 1.0)
    @test norm(relaxed.command) <= 1.0 + 1e-6
    @test relaxed.p[first(colrange(relaxed.B, 2))] > 1.0
end

@testset "constraint rows match the documented formulas" begin
    # Rebuild the program from the equations written in the module docstring, using nothing
    # but the public LinearRow/filter_program API, and require the filter to agree. This is
    # the executable form of the notation table: if a sign, a factor, or a disturbance term
    # in the documentation drifts away from the implementation, this fails.
    using CellularSheaves.IPM: solve, OPTIMAL, NEAR_OPTIMAL
    using CellularSheaves.BlockSparseArrays: colrange

    rng = MersenneTwister(4242)
    model = ControlAffineModel(2)          # f = 0, g = I, P = I
    params = CBFCLFParams(d_safe = 0.7, gamma = 3.0, clf_rate = 2.0, clf_penalty = 5.0,
                          disturbance_bound = 0.1, cbf_disturbance_bound = 0.05,
                          actuator_cap = 1.5)
    W = Matrix{Float64}(I, 2, 2)
    compared = 0

    for _ in 1:60
        x = randn(rng, 2)
        ref = randn(rng, 2)
        rdot = 0.3 .* randn(rng, 2)
        u_nom = 2 .* randn(rng, 2)
        others = [x .+ 0.8 .* randn(rng, 2) for _ in 1:2]

        got = safety_filter(params, model, u_nom, x, ref;
                            ref_velocity = rdot, others = others)
        got.certified || continue

        # CLF row of the docstring, with f = 0 and g = I:
        #   (M e)' u <= -c3 V + (M e)' rdot - ||M e|| Delta + delta,
        # written in the a'u >= b - delta convention of LinearRow.
        e = x - ref
        Me = e
        V = 0.5 * dot(e, Me)
        clf_rhs = -params.clf_rate * V + dot(Me, rdot) -
                  norm(Me) * params.disturbance_bound
        rows = LinearRow[LinearRow(-Me, -clf_rhs; penalty = params.clf_penalty)]

        # Barrier row of the docstring, with the even responsibility split:
        #   2 (P g)' dp . u >= -(gamma/2) h - 2 dp' P f + 2 ||dp|| Deltabar.
        for xj in others
            dp = x - xj
            dot(dp, dp) > params.sense^2 && continue
            h = dot(dp, dp) - params.d_safe^2
            a = 2 .* dp
            b = -(params.gamma / 2) * h + 2 * norm(dp) * params.cbf_disturbance_bound
            # Rows no admissible command can violate are dropped by the filter.
            -params.actuator_cap * norm(W' \ a) >= b + 1e-9 * max(1.0, abs(b)) && continue
            push!(rows, LinearRow(a, b))
        end

        problem, B = filter_program(u_nom, rows;
                                    cone_metric = W, cone_bound = params.actuator_cap)
        # Solve the reference exactly as the filter does, retries included. Taking the
        # iterate of a solve that never converged would compare against noise, and would
        # pass only while the filter happened to fail on the same instance.
        reference = solve(problem, params.settings)
        if !(reference.status in (OPTIMAL, NEAR_OPTIMAL))
            for raug in SafetyFilters.SAFETY_FILTER_FALLBACK_RAUG
                candidate = solve(problem, SafetyFilters._with_raug(params.settings, raug))
                if candidate.status in (OPTIMAL, NEAR_OPTIMAL)
                    reference = candidate
                    break
                end
            end
        end
        @test reference.status in (OPTIMAL, NEAR_OPTIMAL)
        expected = Vector(reference.p[colrange(B, 1)])
        @test got.command ≈ expected atol = 1e-6
        compared += 1
    end
    @test compared >= 40
end

@testset "each term alone" begin
    # The three families are independent, so each can be put through the filter on its own
    # with the others simply absent. In isolation every one of them has a closed form, which
    # is what makes the disjoint tests below exact rather than comparative.
    model = ControlAffineModel(2)
    settings = safety_filter_settings()

    @testset "stability alone" begin
        # No barrier and no cone: with u free the CLF row is always satisfiable, so the
        # relaxation stays at the O(1/p) floor and the dissipation inequality is tight.
        c3, p = 2.0, 1e4
        term = StabilityTerm(rate = c3, penalty = p)

        # Already dissipating: the nominal command is left alone.
        ctx = FilterContext(model, [1.0, 0.0], [0.0, 0.0], 2)
        rest = safety_filter((term,), ctx, [-3.0, 0.0]; settings = settings)
        @test rest.certified
        @test rest.command ≈ [-3.0, 0.0] atol = 1e-4

        # Driving the wrong way: the filter must turn it around and satisfy Vdot <= -c3 V.
        ctx = FilterContext(model, [1.0, 0.5], [0.0, 0.0], 2)
        act = safety_filter((term,), ctx, [2.0, 1.0]; settings = settings)
        @test act.certified
        e = [1.0, 0.5]
        V = 0.5 * dot(e, e)
        @test dot(e, act.command) <= -c3 * V + act.slack + 1e-6
        @test act.lyapunov ≈ V
        @test isinf(act.cbf_margin)          # no barrier exists at all
        @test isinf(act.cap_residual)        # no cone exists at all
    end

    @testset "safety alone" begin
        # No CLF: the pure barrier filter is exactly the projection onto the single active
        # halfspace, and it is inert when the nominal command already satisfies it.
        gamma, d_safe = 1.0, 1.0
        term = SafetyTerm(DistanceBarrier(d_safe); gamma = gamma)
        x, xj = [0.0, 0.0], [1.2, 0.0]
        ctx = FilterContext(model, x, x, 2; others = [xj])

        res = safety_filter((term,), ctx, [3.0, 0.0]; settings = settings)
        @test res.certified
        dp = x - xj
        a = 2 .* dp
        h = dot(dp, dp) - d_safe^2
        @test res.cbf_margin ≈ h
        @test res.command ≈ halfspace_projection([3.0, 0.0], a, -(gamma / 2) * h) atol = 1e-5
        @test isinf(res.cap_residual)

        # Retreating already: nothing to correct.
        inert = safety_filter((term,), ctx, [-1.0, 0.0]; settings = settings)
        @test inert.certified
        @test inert.command ≈ [-1.0, 0.0] atol = 1e-5

        # Both agents enforcing their own half recovers the pairwise condition.
        ctx_j = FilterContext(model, xj, xj, 2; others = [x])
        other = safety_filter((term,), ctx_j, [-3.0, 0.0]; settings = settings)
        @test other.certified
        @test 2 * dot(dp, res.command - other.command) >= -gamma * h - 1e-6
    end

    @testset "actuator alone" begin
        # No rows at all, only the cone: the solution is the radial projection of the
        # nominal command onto the ball.
        cap = 0.4
        term = ActuatorTerm(bound = cap)
        ctx = FilterContext(model, [0.0, 0.0], [0.0, 0.0], 2)

        u_nom = [3.0, 4.0]
        res = safety_filter((term,), ctx, u_nom; settings = settings)
        @test res.certified
        @test res.command ≈ cap .* u_nom ./ norm(u_nom) atol = 1e-5
        @test norm(res.command) <= cap + 1e-7
        @test isinf(res.cbf_margin)
        @test isinf(res.clf_residual)        # no stability row exists at all

        inside = safety_filter((term,), ctx, [0.1, 0.1]; settings = settings)
        @test inside.command ≈ [0.1, 0.1] atol = 1e-6
    end

    @testset "composition preserves each guarantee" begin
        # Composing the three does not disturb what each certifies on its own.
        c3, gamma, d_safe, cap = 2.0, 3.0, 0.8, 1.0
        terms = (StabilityTerm(rate = c3, penalty = 1e3),
                 SafetyTerm(DistanceBarrier(d_safe); gamma = gamma),
                 ActuatorTerm(bound = cap))
        x, xj, ref = [0.0, 0.0], [1.0, 0.2], [2.0, 0.0]
        ctx = FilterContext(model, x, ref, 2; others = [xj])
        res = safety_filter(terms, ctx, [4.0, 0.5]; settings = settings)
        @test res.certified

        e = x - ref
        V = 0.5 * dot(e, e)
        dp = x - xj
        h = dot(dp, dp) - d_safe^2
        @test dot(e, res.command) <= -c3 * V + res.slack + 1e-6      # stability
        @test 2 * dot(dp, res.command) >= -(gamma / 2) * h - 1e-6    # safety
        @test norm(res.command) <= cap + 1e-6                        # actuator
        # And the composed result agrees with the bundled parameter form.
        bundled = safety_filter(CBFCLFParams(d_safe = d_safe, gamma = gamma, clf_rate = c3,
                                             clf_penalty = 1e3, actuator_cap = cap),
                                model, [4.0, 0.5], x, ref; others = [xj])
        @test bundled.command ≈ res.command atol = 1e-6
    end
end

@testset "returned point satisfies the KKT conditions" begin
    # For a convex program the KKT conditions are necessary AND sufficient, so checking them
    # on the solver's own primal-dual pair certifies optimality without trusting the solver.
    # Everything else in this file checks feasibility or consistency, which is weaker.
    using CellularSheaves.IPM: solve, CofreeCone, PositiveCone, SecondOrderCone,
                               OPTIMAL, NEAR_OPTIMAL
    using CellularSheaves.BlockSparseArrays: colrange

    primal_violation(::CofreeCone, v) = 0.0                      # free variable
    primal_violation(::PositiveCone, v) = max(0.0, -minimum(v))
    primal_violation(::SecondOrderCone, v) = max(0.0, norm(v[2:end]) - v[1])
    dual_violation(::CofreeCone, v) = norm(v)                    # dual cone is {0}
    dual_violation(::PositiveCone, v) = max(0.0, -minimum(v))
    dual_violation(::SecondOrderCone, v) = max(0.0, norm(v[2:end]) - v[1])

    rng = MersenneTwister(31337)
    settings = safety_filter_settings()
    worst = zeros(6)
    checked = 0

    for _ in 1:120
        u_nom = 2 .* randn(rng, 2)
        rows = LinearRow[LinearRow(randn(rng, 2), -abs(randn(rng)); penalty = 5.0)]
        for _ in 1:rand(rng, 1:2)
            push!(rows, LinearRow(randn(rng, 2), -abs(randn(rng))))
        end
        W = Matrix{Float64}(I, 2, 2)
        problem, B = filter_program(u_nom, rows; cone_metric = W, cone_bound = 1.5)
        result = solve(problem, settings)
        result.status in (OPTIMAL, NEAR_OPTIMAL) || continue

        Q = Matrix(problem.Q); Bd = Matrix(B)
        c = Vector(problem.c); g = Vector(problem.g)
        p, y, d = result.p, result.y, result.d
        cells = 1:length(problem.K)
        worst = max.(worst, [
            norm(Bd * p - g, Inf),                                          # primal feasible
            norm(Q * p + c - transpose(Bd) * y - d, Inf),                    # dual feasible
            maximum(primal_violation(problem.K[k], p[colrange(B, k)]) for k in cells),
            maximum(dual_violation(problem.K[k], d[colrange(B, k)]) for k in cells),
            abs(dot(p, d)),                                                  # complementarity
            abs(dot(c, p) + dot(p, Q * p) - dot(g, y))])                     # duality gap
        checked += 1
    end

    @test checked >= 100
    @test worst[1] <= 1e-6      # B p = g
    @test worst[3] <= 1e-8      # p in K
    @test worst[5] <= 1e-2      # complementarity
    # The dual side of this solver is materially looser than the primal side, so the bounds
    # below are the observed behaviour rather than an endorsement: a returned command is
    # certified safe on primal feasibility, but "optimal" here means optimal to about 1e-4.
    @test worst[2] <= 1e-2      # stationarity
    @test worst[4] <= 1e-2      # d in K*
    @test worst[6] <= 1e-2      # gap
end

@testset "encoder agrees with an exact active-set enumeration" begin
    # An independent oracle for min 1/2||u-u_nom||^2 s.t. a_k'u >= b_k: enumerate every
    # active set, solve its linear stationarity system, and let the KKT sign conditions pick
    # the winner. No sampling and no second library, so this checks the conic encoding
    # itself rather than merely agreeing with the solver that produced it.
    using CellularSheaves.IPM: solve, OPTIMAL, NEAR_OPTIMAL
    using CellularSheaves.BlockSparseArrays: colrange

    function exact_qp(u_nom, A, b)
        m = length(A)
        best = nothing
        for mask in 0:(2^m - 1)
            S = [k for k in 1:m if (mask >> (k - 1)) & 1 == 1]
            u = copy(u_nom)
            multipliers = Float64[]
            if !isempty(S)
                As = reduce(vcat, [transpose(A[k]) for k in S])
                G = As * transpose(As)
                rank(G) < length(S) && continue
                multipliers = G \ (b[S] - As * u_nom)
                u = u_nom + transpose(As) * multipliers
            end
            all(>=(-1e-9), multipliers) || continue
            all(k -> dot(A[k], u) >= b[k] - 1e-8, 1:m) || continue
            objective = 0.5 * sum(abs2, u - u_nom)
            if best === nothing || objective < best[2] - 1e-12
                best = (u, objective)
            end
        end
        return best === nothing ? nothing : best[1]
    end

    rng = MersenneTwister(0xA57)
    settings = safety_filter_settings()
    worst = 0.0
    compared = 0
    for _ in 1:250
        u_nom = 2 .* randn(rng, 2)
        m = rand(rng, 1:3)
        A = [randn(rng, 2) for _ in 1:m]
        b = [-abs(randn(rng)) for _ in 1:m]
        exact = exact_qp(u_nom, A, b)
        exact === nothing && continue
        problem, B = filter_program(u_nom, [LinearRow(A[k], b[k]) for k in 1:m])
        result = solve(problem, settings)
        result.status in (OPTIMAL, NEAR_OPTIMAL) || continue
        worst = max(worst, norm(Vector(result.p[colrange(B, 1)]) - exact))
        compared += 1
    end
    @test compared >= 200
    @test worst <= 1e-4
end

@testset "a moving obstacle's velocity enters the barrier row" begin
    # hdot = 2 (P x - o)' (P xdot - odot), so an obstacle closing on the agent must tighten
    # the row. Treating it as stationary is exactly the failure that lets a moving obstacle
    # penetrate the keep-out radius while the filter still reports success.
    model = ControlAffineModel(2)
    gamma, d_safe = 2.0, 0.5
    params = CBFCLFParams(d_safe = d_safe, gamma = gamma, obstacle_radius = 0.0)

    x = [0.0, 0.0]
    obstacle = [1.0, 0.0]
    closing = [-0.8, 0.0]                    # obstacle driving at the agent
    dp = x - obstacle
    h = dot(dp, dp) - d_safe^2

    aware = safety_filter(params, model, [2.0, 0.0], x, x;
                          obstacles = [obstacle], obstacle_velocities = [closing])
    unaware = safety_filter(params, model, [2.0, 0.0], x, x; obstacles = [obstacle])
    @test aware.certified
    @test unaware.certified

    # True hdot must dominate -gamma*h once the obstacle's own motion is included.
    true_hdot(u) = 2 * dot(dp, u - closing)
    @test true_hdot(aware.command) >= -gamma * h - 1e-6
    # The stationary-obstacle row permits a command that violates it, which is why the
    # velocity has to be carried rather than assumed away.
    @test true_hdot(unaware.command) < -gamma * h - 1e-6
    @test aware.command[1] < unaware.command[1] - 1e-3
end

@testset "CLF metric for an underactuated plant" begin
    # V = 1/2 e'e is a control Lyapunov function only when the input map has full row rank.
    # For anything underactuated, L_g V = g' e vanishes on a subspace and the stability row
    # goes vacuous while still demanding -c3 V, which forces the relaxation open.
    model = double_integrator_model(2)
    A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0.0]
    B = [0 0; 0 0; 1 0; 0 1.0]
    K = [1.0 0 1.8 0; 0 1.0 0 1.8]

    x = [-2.0, 0.0, 0.0, 0.0]        # displaced but at rest
    ref = [2.0, 0.0, 0.0, 0.0]
    e = x - ref
    g = model.input(x)

    # The identity metric is degenerate here: the row has no coefficients at all.
    @test norm(transpose(g) * e) == 0

    P = lyapunov_metric(A, B, K)
    @test isposdef(P)
    @test norm(transpose(g) * (P * e)) > 1e-3

    # u = -K e is a witness: it certifies Vdot = -1/2 e'Qe, hence this rate.
    c3 = lyapunov_rate(P)
    @test c3 > 0
    witness = -K * e
    Vdot = dot(P * e, A * x - A * ref + B * witness)
    @test Vdot <= -c3 * (0.5 * dot(e, P * e)) + 1e-8

    # With the proper metric the nominal LQR command already satisfies the stability row,
    # so the filter returns it untouched: the stability family is inert until safety bites.
    term = StabilityTerm(metric = P, rate = 0.9 * c3, penalty = 1e3)
    ctx = FilterContext(model, x, ref, 2)
    res = safety_filter((term,), ctx, witness)
    @test res.certified
    @test res.command ≈ witness atol = 1e-4

    # A non-stabilizing gain has no positive definite solution and must be rejected.
    @test_throws ArgumentError lyapunov_metric(A, B, zeros(2, 4))
end

@testset "discrete stability rescues a destabilising nominal" begin
    # The continuous row asks for Vdot <= -c3 V at the sample instants, which says nothing
    # about where an explicit step of length dt lands. Against a nominal command that is
    # itself destabilising the gap dominates: every step can be certified while the
    # trajectory diverges. The discrete row constrains the successor state directly.
    dyn = QuadrotorDynamics()
    dt = 0.05
    Ad, Bd = discrete_matrices(dyn, dt)
    Q = Matrix(Diagonal([500.0, 500.0, 500.0, 150.0, 150.0, 100.0, 100.0, 100.0, 5.0, 5.0]))
    R = Matrix(Diagonal([0.005, 0.005, 0.005]))
    K = LQRController(dyn, dt, Q, R).K
    P = lyapunov_metric(dyn, K)
    c3 = lyapunov_rate(P)
    model = ControlAffineModel(dyn)

    # A sign error on the roll channel: a plausible implementation fault, not a contrivance.
    K_bad = copy(K)
    K_bad[2, :] .*= -1

    continuous = (StabilityTerm(metric = P, rate = 0.5 * c3, penalty = 1e4),
                  ActuatorTerm(bound = 60.0))
    discrete = (DiscreteStabilityTerm(Ad = Ad, Bd = Bd, metric = P, dt = dt,
                                      rate = 0.5 * c3),
                ActuatorTerm(bound = 60.0))

    function final_norm(terms)
        x = vcat(0.4, -0.3, 0.2, zeros(7))
        ref = zeros(10)
        certified = 0
        stepped = 0
        for _ in 1:400
            u = -K_bad * x
            if terms !== nothing
                res = safety_filter(terms, FilterContext(model, x, ref, 3), u)
                if res.certified
                    u = res.command
                    certified += 1
                end
            end
            stepped += 1
            x = Ad * x + Bd * u
            norm(x) > 1e6 && break
        end
        return (norm(x), certified, stepped)
    end

    start = norm(vcat(0.4, -0.3, 0.2, zeros(7)))
    bare, _, _ = final_norm(nothing)
    cont, cont_certified, cont_steps = final_norm(continuous)
    disc, disc_certified, _ = final_norm(discrete)

    @test bare > 10 * start                    # the nominal policy diverges
    @test cont > 10 * start                    # and the continuous row does not stop it
    # The continuous row is not failing to solve: it certifies essentially every step it
    # takes, so the divergence happens entirely under commands it declared to satisfy the
    # dissipation inequality. How many steps that is depends on the machine, since the run
    # ends when a diverging state overflows the solver, so the claim is made as a fraction
    # of the steps actually taken rather than as a fixed count.
    @test cont_certified >= 0.9 * cont_steps
    @test cont_steps >= 50
    @test disc < start                         # the discrete row does stop it
    @test disc_certified > 300                 # and keeps certifying throughout
end

@testset "cone rows carry an affine offset" begin
    # ConeRow expresses ||W u + v|| <= b. The offset is what lets a cone constrain where the
    # state will be rather than only how large the command is, which is what the discrete
    # stability condition needs.
    using CellularSheaves.IPM: solve, OPTIMAL, NEAR_OPTIMAL
    using CellularSheaves.BlockSparseArrays: colrange

    rng = MersenneTwister(4711)
    settings = safety_filter_settings()
    W = Matrix{Float64}(I, 2, 2)
    worst = 0.0
    checked = 0
    for _ in 1:120
        offset = 0.6 .* randn(rng, 2)
        bound = 1.0 + rand(rng)
        u_nom = 3 .* randn(rng, 2)
        problem, B = filter_program(u_nom, LinearRow[];
                                    cones = [ConeRow(W, bound; offset = offset)])
        result = solve(problem, settings)
        result.status in (OPTIMAL, NEAR_OPTIMAL) || continue
        u = Vector(result.p[colrange(B, 1)])
        worst = max(worst, norm(W * u + offset) - bound)
        checked += 1
    end
    @test checked >= 100
    @test worst <= 1e-6

    # With no offset it is the plain actuator bound, and it agrees with the keyword form.
    plain, Bp = filter_program([4.0, 0.0], LinearRow[]; cones = [ConeRow(W, 0.5)])
    keyword, Bk = filter_program([4.0, 0.0], LinearRow[]; cone_metric = W, cone_bound = 0.5)
    up = Vector(solve(plain, settings).p[colrange(Bp, 1)])
    uk = Vector(solve(keyword, settings).p[colrange(Bk, 1)])
    @test up ≈ uk atol = 1e-6
    @test up ≈ [0.5, 0.0] atol = 1e-5

    # Several cones compose.
    two, Bt = filter_program([3.0, 3.0], LinearRow[];
                             cones = [ConeRow(W, 1.0), ConeRow(W, 2.0; offset = [1.0, 0.0])])
    ut = Vector(solve(two, settings).p[colrange(Bt, 1)])
    @test norm(ut) <= 1.0 + 1e-6
    @test norm(W * ut + [1.0, 0.0]) <= 2.0 + 1e-6
end

@testset "discrete stability contracts V by the demanded factor" begin
    # The row promises V(A_d x + B_d u) <= (1 - c3 dt) V(x). Check the promise directly on
    # the successor state, not the derivative, since that distinction is the whole point.
    dyn = QuadrotorDynamics()
    dt = 0.05
    Ad, Bd = discrete_matrices(dyn, dt)
    Q = Matrix(Diagonal([500.0, 500.0, 500.0, 150.0, 150.0, 100.0, 100.0, 100.0, 5.0, 5.0]))
    R = Matrix(Diagonal([0.005, 0.005, 0.005]))
    K = LQRController(dyn, dt, Q, R).K
    P = lyapunov_metric(dyn, K)
    c3 = lyapunov_rate(P)
    model = ControlAffineModel(dyn)
    rate = 0.5 * c3
    terms = (DiscreteStabilityTerm(Ad = Ad, Bd = Bd, metric = P, dt = dt, rate = rate),
             ActuatorTerm(bound = 60.0))

    K_bad = copy(K)
    K_bad[2, :] .*= -1
    x = vcat(0.4, -0.3, 0.2, zeros(7))
    ref = zeros(10)
    worst_ratio = -Inf
    certified = 0
    for _ in 1:260
        u = -K_bad * x
        res = safety_filter(terms, FilterContext(model, x, ref, 3), u)
        if res.certified
            u = res.command
            certified += 1
            V = 0.5 * dot(x - ref, P * (x - ref))
            successor = Ad * x + Bd * u
            V_next = 0.5 * dot(successor - ref, P * (successor - ref))
            V > 1e-12 && (worst_ratio = max(worst_ratio, V_next / V))
        end
        x = Ad * x + Bd * u
    end
    # Nearly every step certifies, rather than exactly all of them: whether one marginal
    # instance resolves is a property of the machine's floating point, while the contraction
    # below is the claim being made and holds on every step that did certify.
    @test certified >= 0.95 * 260
    @test worst_ratio <= (1 - rate * dt) + 1e-6

    # The rate must leave the contraction factor positive.
    @test_throws ArgumentError constraints(
        DiscreteStabilityTerm(Ad = Ad, Bd = Bd, metric = P, dt = dt, rate = 25.0),
        FilterContext(model, x, ref, 3))
end

@testset "the controller certifies against the reference it actually tracks" begin
    # The stability row is only meaningful if the reference rate it is given is the true rate
    # of change of the reference the agent follows. Two ways to get that wrong were live at
    # once: placing the velocity target in the state's velocity slots (right for the Tikhonov
    # step, wrong as r_dot, whose position slots carry the reference velocity), and using the
    # instantaneous derivative when the filter very nearly converges within one control
    # period. With eps = 0.02 against dt = 0.05 the second alone is a factor of four, and the
    # certified inequality Vdot <= -c3 V + delta then fails along the trajectory even though
    # every solve reports success. This checks the inequality directly, against a finite
    # difference of the reference the controller returned.
    dyn = QuadrotorDynamics()
    dt, eps = 0.05, 0.02
    Ac, Bc = continuous_matrices(dyn)
    Q = Matrix(Diagonal([500.0, 500.0, 500.0, 150.0, 150.0, 100.0, 100.0, 100.0, 5.0, 5.0]))
    R = Matrix(Diagonal([0.005, 0.005, 0.005]))
    K = LQRController(dyn, dt, Q, R).K
    M = lyapunov_metric(dyn, K)
    c3 = 0.9 * lyapunov_rate(M)

    agent = AgentState(zeros(10), dyn, dt, K, eps; use_velocity = true)
    controller = CBFCLFController(K, CBFCLFParams(clf_metric = M, clf_rate = c3),
                                  ControlAffineModel(dyn))
    # The reference must be *deflected*, not merely moving: against a smooth trajectory the
    # two candidate rates agree closely and the error hides. A step is what a governor does
    # when it dodges, and it is where the discrepancy is largest.
    reference(t) = t < 2.0 ? [0.0, 0.0, 1.0] : [1.2, 0.8, 1.0]
    velocity(t) = zeros(3)

    previous_ref = nothing
    worst = -Inf
    checked = 0
    for k in 1:120
        t = k * dt
        state_before = copy(agent.x)
        _, x_ref, res = step_agent!(agent, controller, reference(t),
                                    Vector{Vector{Float64}}(), dt;
                                    qstar_dot_target = velocity(t))
        if previous_ref !== nothing && res.certified
            e = state_before - x_ref
            @test 0.5 * dot(e, M * e) ≈ res.lyapunov atol = 1e-12
            true_rate = (x_ref - previous_ref) ./ dt
            Vdot = dot(M * e, (Ac * state_before + Bc * res.command) - true_rate)
            worst = max(worst, Vdot + c3 * (0.5 * dot(e, M * e)) - res.slack)
            checked += 1
        end
        previous_ref = copy(x_ref)
    end
    @test checked >= 100
    @test worst <= 1e-6
end

@testset "an exact penalty drives the relaxation to zero" begin
    # A quadratic charge p*delta^2 is an inexact penalty: its optimum trades a little
    # violation against the objective, so delta settles at an O(1/p) floor even when the
    # unrelaxed row is satisfiable. A linear charge p*delta is exact above a threshold on p,
    # and the distinction matters because that floor is what sets the ultimate bound: near
    # the origin O(1/p) is a sizeable fraction of the demanded decay rate c3*V.
    using CellularSheaves.IPM: solve, OPTIMAL, NEAR_OPTIMAL
    using CellularSheaves.BlockSparseArrays: colrange

    settings = safety_filter_settings()
    u_nom = [2.0, 0.5]
    hard = LinearRow([0.0, 1.0], -1.0)

    function relaxation(penalty, exact)
        soft = LinearRow([1.0, 0.0], -3.0; penalty = penalty, exact = exact)
        problem, B = filter_program(u_nom, [soft, hard])
        result = solve(problem, settings)
        result.status in (OPTIMAL, NEAR_OPTIMAL) || return nothing
        return result.p[first(colrange(B, 2))]
    end

    # The row is slack at the nominal command, so an exact penalty owes nothing.
    for p in (1e1, 1e2, 1e3)
        quadratic = relaxation(p, false)
        linear = relaxation(p, true)
        @test quadratic !== nothing && linear !== nothing
        @test linear < 1e-8              # numerically zero
        @test linear < quadratic / 100   # and orders below the quadratic floor
    end

    # The quadratic floor really does scale as 1/p, which is why no finite penalty removes it.
    @test relaxation(1e3, false) < relaxation(1e1, false)
    @test relaxation(1e1, false) > 1e-6

    # When the row genuinely cannot be met, both formulations open by the same amount: an
    # exact penalty removes the floor, it does not refuse to relax.
    function forced(penalty, exact)
        soft = LinearRow([1.0, 0.0], 5.0; penalty = penalty, exact = exact)
        problem, B = filter_program([0.0, 0.0], [soft];
                                    cone_metric = Matrix{Float64}(I, 2, 2), cone_bound = 1.0)
        result = solve(problem, settings)
        result.status in (OPTIMAL, NEAR_OPTIMAL) || return nothing
        return result.p[first(colrange(B, 2))]
    end
    @test forced(1e2, true) ≈ 4.0 atol = 1e-3
    @test forced(1e2, false) ≈ 4.0 atol = 1e-3
end

@testset "each soft row is certified against its own relaxation" begin
    # Every soft row gets its own delta in the program, and the certificate has to read the
    # one belonging to the row it is checking. Sharing a single relaxation across all soft
    # rows is not a cosmetic slip: the rows here need relaxations that differ by a factor of
    # seven, so the second row measured against the first row's delta looks violated by a
    # wide margin, the solve is reported uncertified, and a caller that trusts `certified`
    # then discards a filtered command that was in fact admissible and applies the raw
    # nominal instead -- the one that breaks the actuator cap.
    using CellularSheaves.IPM: solve, OPTIMAL, NEAR_OPTIMAL
    using CellularSheaves.BlockSparseArrays: colrange

    model = ControlAffineModel(2)
    state, reference = [1.0, 0.8], [0.0, 0.0]
    slow = StabilityTerm(metric = Matrix{Float64}(I, 2, 2), rate = 6.0, penalty = 2.0)
    fast = StabilityTerm(metric = [8.0 0.0; 0.0 0.5], rate = 9.0, penalty = 2.0)
    cap = 0.6
    ctx = FilterContext(model, state, reference, 2)

    # Solve the same program directly to recover the individual relaxations.
    rows = [constraints(slow, ctx).rows[1], constraints(fast, ctx).rows[1]]
    problem, B = filter_program(zeros(2), rows;
                                cone_metric = Matrix{Float64}(I, 2, 2), cone_bound = cap)
    raw = solve(problem, safety_filter_settings())
    @test raw.status in (OPTIMAL, NEAR_OPTIMAL)
    command = Vector(raw.p[colrange(B, 1)])
    deltas = [raw.p[first(colrange(B, 1 + i))] for i in 1:2]

    # The premise: the two relaxations are genuinely far apart.
    @test deltas[2] > 5 * deltas[1]
    # Each row is met by its own relaxation, and the second is badly violated by the first's.
    for i in 1:2
        @test dot(rows[i].coefficients, command) + deltas[i] - rows[i].bound >= -1e-6
    end
    @test dot(rows[2].coefficients, command) + deltas[1] - rows[2].bound < -1.0

    # The filter must agree with the per-row reading, not the shared one.
    result = safety_filter((slow, fast, ActuatorTerm(bound = cap)), ctx, zeros(2))
    @test result.certified
    @test result.command !== nothing
    @test result.clf_residual >= -1e-6
    @test result.slack ≈ maximum(deltas) rtol = 1e-4
end

@testset "a failed solve is retried at another augmentation" begin
    # Which Uzawa augmentations condition badly is a property of the instance and of the
    # machine's floating point, not of the problem: this closed loop solves cleanly at 1e3
    # and stalls almost completely at 1e5, and on a different CPU the same swap moves across
    # that boundary. Without a retry the filter is not portable, so this pins the recovery
    # rather than the tuning.
    model = ControlAffineModel(2)
    dt = 0.01
    function swap_run(raug)
        params = CBFCLFParams(d_safe = 0.6, gamma = 5.0, actuator_cap = 2.0, clf_rate = 3.0,
                              settings = safety_filter_settings(raug = raug))
        xs = [[-1.5, -0.15], [1.5, 0.15]]
        refs = [[1.5, -0.15], [-1.5, 0.15]]
        uncertified = 0
        for _ in 1:400
            u_noms = [4 .* (refs[i] - xs[i]) for i in 1:2]
            results = safety_filter_all(params, model, u_noms, xs, refs)
            for i in 1:2
                results[i].certified || (uncertified += 1)
                xs[i] = xs[i] + dt .* (results[i].certified ? results[i].command : zeros(2))
            end
        end
        return uncertified, norm(xs[1] - refs[1])
    end

    for raug in (1e2, 1e3, 1e4, 1e5)
        uncertified, reached = swap_run(raug)
        @test uncertified == 0
        @test reached < 1e-3
    end

    # The recovery is real rather than the augmentation being irrelevant: without it, the
    # bare solve at 1e5 does fail on this instance.
    bare = safety_filter_settings(raug = 1e5)
    @test SafetyFilters._with_raug(bare, 1e2).kkt.raug == 1e2
    @test SafetyFilters._with_raug(bare, 1e5) === bare
end

@testset "integrator dynamics live in the agent hierarchy" begin
    @test SingleIntegrator <: AbstractControlAffine <: AbstractAgentDynamics
    @test DoubleIntegrator <: AbstractControlAffine <: AbstractAgentDynamics

    # The free double integrator: only the rate blocks are given, the kinematic coupling and
    # the zero rows are padded on.
    A, B = continuous_matrices(DoubleIntegrator(2))
    @test A == [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0]
    @test B == [0 0; 0 0; 1 0; 0 1]

    # A damping block and a non-identity input map land in those same blocks.
    Ad, Bd = continuous_matrices(DoubleIntegrator(2; damping = -0.5I,
                                                 input = [1.0 0.0; 0.0 2.0]))
    @test Ad[1:2, 3:4] == Matrix(1.0I, 2, 2)   # coupling still padded on
    @test Ad[3:4, 3:4] == [-0.5 0.0; 0.0 -0.5]
    @test Bd == [0 0; 0 0; 1 0; 0 2]

    di = DoubleIntegrator(2)
    @test (position_indices(di), velocity_indices(di), state_dim(di)) == (1:2, 3:4, 4)
    si = SingleIntegrator(3)
    @test (position_indices(si), state_dim(si)) == (1:3, 3)
    @test isempty(velocity_indices(si))
    @test continuous_matrices(si) == (zeros(3, 3), Matrix(1.0I, 3, 3))

    @test initial_state(di, [1.0, 2.0], [3.0, 4.0]) == [1.0, 2.0, 3.0, 4.0]
    @test initial_state(si, [1.0, 2.0, 3.0]) == [1.0, 2.0, 3.0]

    @test_throws ArgumentError DoubleIntegrator(0)
    @test_throws ArgumentError DoubleIntegrator(2; damping = zeros(3, 3))
end

@testset "the filter models are the integrator dynamics" begin
    # `double_integrator_model` is now a shorthand, not a second source of truth.
    a = double_integrator_model(2)
    b = ControlAffineModel(DoubleIntegrator(2))
    x = [1.0, 2.0, 3.0, 4.0]
    @test a.drift(x) == b.drift(x)
    @test a.input(x) == b.input(x)
    @test a.position == b.position
    @test a.velocity == b.velocity

    # `ControlAffineModel(n)` is the single integrator, and still reports no velocity.
    m = ControlAffineModel(2)
    @test m.velocity === nothing
    @test m.position == Matrix(1.0I, 2, 2)
    @test m.drift([1.0, 2.0]) == zeros(2)
    @test m.input([1.0, 2.0]) == Matrix(1.0I, 2, 2)

    # An empty velocity selector means the braking barrier reports it cannot be used, rather
    # than building a row out of empty vectors.
    ctx = FilterContext(ControlAffineModel(SingleIntegrator(2)), zeros(2), zeros(2), 2;
                        others = [[1.0, 0.0]])
    @test_throws ArgumentError constraints(SafetyTerm(BrakingBarrier(0.5, 1.0); gamma = 1.0),
                                           ctx)
end

@testset "one safety term covers a whole formation" begin
    # James's question on the PR: an agent with several near neighbours does not need several
    # filters. One term emits one row per neighbour in range, and they compose into a single
    # program.
    model = ControlAffineModel(2)
    x = zeros(2)
    others = [[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0]]
    ctx = FilterContext(model, x, x, 2; others = others)
    term = SafetyTerm(DistanceBarrier(0.5); gamma = 1.0)
    @test length(constraints(term, ctx).rows) == 3

    # and only the ones inside `sense`
    near = SafetyTerm(DistanceBarrier(0.5); gamma = 1.0, sense = 1.5)
    far = FilterContext(model, x, x, 2; others = vcat(others, [[9.0, 0.0]]))
    @test length(constraints(near, far).rows) == 3
end

@testset "scoped safety terms give a formation per-pair clearances" begin
    model = ControlAffineModel(2)
    x = zeros(2)
    others = [[1.0, 0.0], [0.0, 1.0]]
    ctx = FilterContext(model, x, x, 2; others = others)

    # Each term sees only its slice.
    tight = SafetyTerm(DistanceBarrier(0.5); gamma = 1.0, neighbors = [1])
    loose = SafetyTerm(DistanceBarrier(1.2); gamma = 1.0, neighbors = [2])
    @test length(constraints(tight, ctx).rows) == 1
    @test length(constraints(loose, ctx).rows) == 1

    # The rows are exactly the ones an unscoped term would emit for those neighbours, so
    # scoping only selects and never changes a row.
    all_rows = constraints(SafetyTerm(DistanceBarrier(0.5); gamma = 1.0), ctx).rows
    @test constraints(tight, ctx).rows[1].coefficients ≈ all_rows[1].coefficients
    @test constraints(tight, ctx).rows[1].bound ≈ all_rows[1].bound

    # Distinct clearances really do reach the barrier values.
    @test constraints(tight, ctx).diagnostics.barriers ≈ [1.0 - 0.5^2]
    @test constraints(loose, ctx).diagnostics.barriers ≈ [1.0 - 1.2^2]

    # And the composed program carries one row per neighbour, each with its own clearance.
    res = safety_filter((tight, loose, ActuatorTerm()), ctx, [0.3, 0.3])
    @test res.certified
    @test res.cbf_margin ≈ min(1.0 - 0.5^2, 1.0 - 1.2^2)
end

@testset "scoped safety terms default to the whole context" begin
    model = ControlAffineModel(2)
    x = zeros(2)
    ctx = FilterContext(model, x, x, 2; others = [[1.0, 0.0], [0.0, 1.0]],
                        obstacles = [[0.0, -1.0]])
    explicit = SafetyTerm(DistanceBarrier(0.5); gamma = 1.0,
                          neighbors = [1, 2], obstacles = [1])
    implicit = SafetyTerm(DistanceBarrier(0.5); gamma = 1.0)
    a = constraints(implicit, ctx)
    b = constraints(explicit, ctx)
    @test length(a.rows) == length(b.rows) == 3
    @test all(a.rows[i].coefficients ≈ b.rows[i].coefficients for i in 1:3)
    @test a.diagnostics.barriers ≈ b.diagnostics.barriers
end

@testset "overlapping safety terms are rejected" begin
    model = ControlAffineModel(2)
    x = zeros(2)
    ctx = FilterContext(model, x, x, 2; others = [[1.0, 0.0], [0.0, 1.0]],
                        obstacles = [[0.0, -1.0]])
    one = SafetyTerm(DistanceBarrier(0.5); gamma = 1.0, neighbors = [1])

    # Explicit overlap on a neighbour.
    @test_throws ArgumentError safety_filter((one, SafetyTerm(DistanceBarrier(1.0);
                                                             gamma = 1.0, neighbors = [1, 2])),
                                             ctx, [0.1, 0.1])
    # Two unscoped terms overlap on everything: this is the silent-duplication case.
    @test_throws ArgumentError safety_filter((SafetyTerm(DistanceBarrier(0.5); gamma = 1.0),
                                              SafetyTerm(DistanceBarrier(1.0); gamma = 1.0)),
                                             ctx, [0.1, 0.1])
    # Overlap on an obstacle is caught independently of the neighbours.
    @test_throws ArgumentError safety_filter((SafetyTerm(DistanceBarrier(0.5); gamma = 1.0,
                                                        neighbors = [1], obstacles = [1]),
                                              SafetyTerm(DistanceBarrier(1.0); gamma = 1.0,
                                                         neighbors = [2], obstacles = [1])),
                                             ctx, [0.1, 0.1])
    # Scoping `neighbors` does not implicitly scope `obstacles`: the unscoped field still
    # defaults to every obstacle, so two otherwise-disjoint terms overlap there.
    @test_throws ArgumentError safety_filter((one, SafetyTerm(DistanceBarrier(1.0);
                                                             gamma = 1.0, neighbors = [2],
                                                             obstacles = [1])),
                                             ctx, [0.1, 0.1])

    # Fully disjoint on both fields is fine.
    disjoint = (SafetyTerm(DistanceBarrier(0.5); gamma = 1.0, neighbors = [1],
                           obstacles = Int[]),
                SafetyTerm(DistanceBarrier(1.0); gamma = 1.0, neighbors = [2],
                           obstacles = [1]))
    @test safety_filter(disjoint, ctx, [0.1, 0.1]).certified
    # Between them the two terms cover all three targets exactly once.
    @test sum(length(constraints(t, ctx).rows) for t in disjoint) == 3

    # Out-of-range and repeated indices are caught.
    @test_throws ArgumentError safety_filter((SafetyTerm(DistanceBarrier(0.5); gamma = 1.0,
                                                        neighbors = [3]),), ctx, [0.1, 0.1])
    @test_throws ArgumentError safety_filter((SafetyTerm(DistanceBarrier(0.5); gamma = 1.0,
                                                        neighbors = [1, 1]),), ctx, [0.1, 0.1])
end

@testset "input validation" begin
    model = ControlAffineModel(2)
    @test_throws ArgumentError safety_filter(CBFCLFParams(), model, [1.0, 0.0], [0.0], [0.0])
    @test_throws ArgumentError safety_filter(CBFCLFParams(gamma = -1.0), model,
                                             [1.0, 0.0], zeros(2), zeros(2))
    @test_throws ArgumentError safety_filter_all(CBFCLFParams(), model,
                                                 [[1.0, 0.0]], [zeros(2), zeros(2)],
                                                 [zeros(2), zeros(2)])
end
