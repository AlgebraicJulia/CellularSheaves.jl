using Test
using CellularSheaves
using CellularSheaves.ControlSheaves.CoordinationBenchmarks
using LinearAlgebra

@testset "CoordinationBenchmarks" begin

    @testset "scenario construction" begin
        for (name, sp) in ((:tworing, 0), (:chain, 12), (:ring, 12), (:grid, 3),
                           (:rgg, 12), (:star, 12), (:twoclique, 8), (:expander, 12))
            s = coordination_scenario(name; size_parameter = sp)
            @test s.nagents > 0
            @test size(s.H, 1) == s.dim * s.nagents
            @test size(s.Bmat) == (s.dim * s.nagents, s.dim * s.ntargets)
            # pinning makes the free block positive definite (paper's Lemma 2)
            @test isposdef(Symmetric(Matrix(s.H)))
        end
        @test_throws ArgumentError coordination_scenario(:nonsense)
    end

    @testset "both laws share the harmonic fixed point" begin
        # This is the paper's Lemma 1: eta vanishes exactly at q*, so Diffusion is at
        # equilibrium precisely where Direct's solve lands. If this fails, the two
        # laws are not solving the same problem and no comparison is meaningful.
        for (name, sp) in ((:chain, 16), (:ring, 16), (:grid, 3), (:twoclique, 8))
            s = coordination_scenario(name; size_parameter = sp)
            plan = direct_plan(s)
            p = hcat([[2.2cos(2pi * k / s.ntargets), 2.2sin(2pi * k / s.ntargets)]
                      for k in 1:s.ntargets]...)
            qstar = harmonic_reference(plan, p)

            eta = zeros(length(qstar))
            diffusion_residual!(eta, s, qstar, p)
            @test norm(eta, Inf) < 1e-10

            # and q* really solves H q* = B p
            @test norm(s.H * qstar - s.Bmat * vec(p), Inf) < 1e-10
        end
    end

    @testset "cached solve matches harmonic_extension" begin
        # The hot path is a cached clique-tree solve; the library's canonical
        # harmonic_extension is the oracle it must reproduce.
        for (name, sp) in ((:chain, 16), (:grid, 3), (:star, 12))
            s = coordination_scenario(name; size_parameter = sp)
            plan = direct_plan(s)
            p = hcat([[1.5k, -0.7k] for k in 1:s.ntargets]...)
            cached = harmonic_reference(plan, p)
            oracle = CoordinationBenchmarks.harmonic_reference_oracle(s, p)
            @test norm(cached - oracle, Inf) < 1e-9
        end
    end

    @testset "gain ceilings follow the spectrum" begin
        s = coordination_scenario(:chain; size_parameter = 16)
        spec = spectral_summary(s)
        dt = 0.005
        @test spec.minimum > 0
        @test spec.maximum >= spec.minimum
        @test spec.condition ≈ spec.maximum / spec.minimum
        # Direct's ceiling is topology-independent; Diffusion's is throttled by lambda_max
        @test stable_gain_ceiling(:direct, s, dt) ≈ 2 / dt
        @test stable_gain_ceiling(:diffusion, s, dt) ≈ 2 / (dt * spec.maximum)
        @test stable_gain_ceiling(:diffusion, s, dt) < stable_gain_ceiling(:direct, s, dt)
    end

    @testset "target coverage sets lambda_min" begin
        # Full coverage gives H = L + I, so lambda_min is exactly 1 and kappa
        # collapses. This is what makes Diffusion competitive, and it is about the
        # Dirichlet boundary rather than graph density.
        sparse_cov = coordination_scenario(:ring; size_parameter = 32)
        full_cov = coordination_scenario(:ring; size_parameter = 32, coverage = :full)
        @test spectral_summary(full_cov).minimum ≈ 1.0 atol = 1e-8
        @test spectral_summary(full_cov).condition < spectral_summary(sparse_cov).condition
    end

    @testset "settling" begin
        s = coordination_scenario(:chain; size_parameter = 16)
        diffusion = settle_to_formation(:diffusion, s)
        direct = settle_to_formation(:direct, s)

        @test diffusion.converged
        @test direct.converged
        @test diffusion.method === :diffusion && direct.method === :direct
        @test diffusion.ticks > 0 && direct.ticks > 0
        @test length(diffusion.error_history) == diffusion.ticks + 1
        # error decreases monotonically enough to have actually converged
        @test last(diffusion.error_history) < first(diffusion.error_history)
        @test last(direct.error_history) < first(direct.error_history)
        # Direct is exact per solve, so it settles in far fewer ticks than a
        # diffusion on an ill-conditioned sheaf
        @test direct.ticks < diffusion.ticks

        @test_throws ArgumentError settle_to_formation(:nonsense, s)
    end

    @testset "communication model" begin
        s = coordination_scenario(:chain; size_parameter = 16)
        plan = direct_plan(s)
        @test diffusion_round_slots(s) == 2          # a path has max degree 2
        @test tree_makespan(plan) > 0
        @test iseven(tree_makespan(plan))      # one gather plus one scatter
    end

    @testset "benchmark sweep" begin
        # one row per (family, coverage) pair
        result = benchmark_coordination(family = [(:chain, [12]), (:ring, [12])])
        @test length(result) == 4
        @test Set(r.coverage for r in result) == Set((:sparse, :full))
        for row in result
            @test row.agents == 12
            @test row.condition > 1
            @test row.diffusion.converged && row.direct.converged
            @test row.oracle_residual < 1e-9
            @test row.speedup > 0
            @test row.slot_ratio > 0
        end

        single = benchmark_coordination(family = [(:chain, [12])], coverages = (:sparse,))
        @test length(single) == 1
    end

    @testset "every family in the atlas builds and solves" begin
        # The example plots all twenty families; a family that fails to build or
        # whose plan disagrees with the oracle would ship a wrong figure.
        families = unique(first(f) for f in scenario_family())
        @test length(families) == 20
        for name in families
            sp = last(first(filter(f -> first(f) == name, scenario_family())))
            s = coordination_scenario(name; size_parameter = first(sp))
            plan = direct_plan(s)
            p = hcat([[1.3k, -0.9k] for k in 1:s.ntargets]...)
            @test norm(harmonic_reference(plan, p) -
                       CoordinationBenchmarks.harmonic_reference_oracle(s, p), Inf) < 1e-8
        end
    end

    @testset "settling matches the analytic rate" begin
        # Guards against a settling loop that reports plausible-looking but wrong
        # tick counts. The slowest mode contracts by 2*safety/kappa per tick for
        # Diffusion and by 2*safety for Direct, so the tick count is predictable
        # from the spectrum alone. An initial condition that does not fully excite
        # the slowest mode converges sooner, hence the one-sided band.
        tol, dt, safety = 1e-3, 0.005, 0.4
        predict(rate) = ceil(Int, log(tol) / log1p(-rate))
        for (name, sp) in ((:chain, 24), (:ring, 24), (:star, 16), (:grid, 4))
            s = coordination_scenario(name; size_parameter = sp)
            kappa = spectral_summary(s).condition
            d = settle_to_formation(:diffusion, s; tolerance = tol, dt, safety)
            r = settle_to_formation(:direct, s; tolerance = tol, dt, safety)
            pd = predict(2safety / kappa)
            @test 0.2pd <= d.ticks <= 1.3pd
            @test r.ticks <= 1.3predict(2safety)
        end
    end

    @testset "settling records control effort" begin
        s = coordination_scenario(:chain; size_parameter = 24)
        d = settle_to_formation(:diffusion, s)
        r = settle_to_formation(:direct, s)
        @test d.peak_command > 0 && r.peak_command > 0
        @test d.path_length > 0 && r.path_length > 0
        # Direct runs at a kappa-times-higher gain ceiling, so it commands a more
        # aggressive transient. The example reports this as a fairness caveat; if
        # it ever stopped being true the caveat would be wrong.
        @test r.peak_command > d.peak_command
    end

    @testset "solve is allocation free on the hot path" begin
        # Measured inside a function: at global scope, boxing of non-const globals
        # is attributed to the call and reports phantom allocations.
        s = coordination_scenario(:grid; size_parameter = 4)
        plan = direct_plan(s)
        p = hcat([[1.0k, -1.0k] for k in 1:s.ntargets]...)
        qstar = harmonic_reference(plan, p)
        probe(pl, q, pp) = @allocated harmonic_reference!(q, pl, pp)
        probe(plan, qstar, p)
        @test probe(plan, qstar, p) <= 64

        err(a, b, n) = @allocated CoordinationBenchmarks._formation_error(a, b, n)
        v = copy(qstar)
        err(v, qstar, s.nagents)
        @test err(v, qstar, s.nagents) == 0
    end

    @testset "permutations match the factorization convention" begin
        # harmonic_reference! hand-rolls P'\v and P\v to avoid allocating. If the
        # two directions were swapped the solve would silently return a permuted
        # vector, so pin both against the library operators.
        s = coordination_scenario(:grid; size_parameter = 3)
        plan = direct_plan(s)
        v = collect(1.0:length(plan.rhs))
        fwd, bwd = similar(v), similar(v)
        CoordinationBenchmarks._permute!(fwd, v, plan.perm)
        CoordinationBenchmarks._invpermute!(bwd, v, plan.perm)
        @test fwd == v[plan.perm]
        @test bwd[plan.perm] == v
    end

    @testset "tree statistics" begin
        for (name, sp) in ((:chain, 24), (:grid, 4), (:complete, 10), (:twoclique, 8))
            s = coordination_scenario(name; size_parameter = sp)
            t = tree_statistics(direct_plan(s))
            @test t.treewidth >= 1
            @test t.depth >= 1
            @test t.supernodes >= 1
            # fill counts the diagonal blocks too, so it is never below the
            # off-diagonal count and is never zero for a nonempty factor
            @test t.fill >= t.offdiagonal_fill
            @test t.fill > 0
        end
        # a path is a tree
        @test tree_statistics(direct_plan(coordination_scenario(:chain; size_parameter = 24))).treewidth == 1
        # the complete graph collapses to one supernode: all fill is diagonal
        ct = tree_statistics(direct_plan(coordination_scenario(:complete; size_parameter = 10)))
        @test ct.offdiagonal_fill == 0
        @test ct.fill > 0
    end

    @testset "degenerate trees are charged the centralized bound" begin
        # A single supernode means one worker holds the whole factor, which is a
        # centralized solve. Counting only inter-supernode messages would report a
        # near-zero slot count and fabricate an enormous communication advantage.
        s = coordination_scenario(:complete; size_parameter = 10)
        plan = direct_plan(s)
        @test tree_statistics(plan).supernodes <= 2
        @test tree_makespan(plan) >= 2 * s.nagents
        @test iseven(tree_makespan(plan))
    end

    @testset "step profile" begin
        s = coordination_scenario(:grid; size_parameter = 4)
        d = step_profile(:diffusion, s; ticks = 50, reps = 5, solves = 50)
        r = step_profile(:direct, s; ticks = 50, reps = 5, solves = 50)
        @test d.total > 0 && r.total > 0
        @test d.setup == 0                       # a diffusion has no setup stage
        @test r.setup > 0
        @test r.setup < r.total
        # every advertised stage is actually timed
        for (sym, _, _) in DIFFUSION_STEPS
            @test haskey(d.seconds, sym)
        end
        for (sym, _, _) in vcat(DIRECT_SETUP_STEPS, DIRECT_STEP_STEPS)
            @test haskey(r.seconds, sym)
        end
    end

    @testset "rollouts" begin
        s = coordination_scenario(:chain; size_parameter = 12)
        for mode in (:static, :orbit, :maneuver)
            d = rollout_coordination(:diffusion, s; dt = 0.02, horizon = 6.0, mode)
            r = rollout_coordination(:direct, s; dt = 0.02, horizon = 6.0, mode)
            @test length(d.times) == length(d.errors) == length(r.times)
            @test all(isfinite, d.errors) && all(isfinite, r.errors)
            # same targets, same initial condition, only the control line differs
            @test d.times == r.times
            @test d.errors[1] ≈ r.errors[1]
        end

        # a static target means both laws converge; a moving one means neither does
        st = rollout_coordination(:direct, s; dt = 0.02, horizon = 30.0, mode = :static)
        @test st.errors[end] < 1e-3st.errors[1]

        # the maneuver mode must actually produce a step in the reference
        mv = rollout_coordination(:direct, s; dt = 0.02, horizon = 30.0, mode = :maneuver)
        @test argmax(mv.errors) > 1
        @test maximum(mv.errors) > 10mv.errors[end]

        # the command filter is applied identically to both laws
        fd = rollout_coordination(:diffusion, s; dt = 0.02, horizon = 6.0, epsilon = 0.1)
        @test all(isfinite, fd.errors)

        @test_throws ArgumentError rollout_coordination(:nonsense, s)
    end

    @testset "tracking bandwidth" begin
        s = coordination_scenario(:chain; size_parameter = 16)
        dt = 0.02
        spec = spectral_summary(s)
        @test tracking_bandwidth(:direct, s; dt) ≈
              0.4stable_gain_ceiling(:direct, s, dt)
        @test tracking_bandwidth(:diffusion, s; dt) ≈
              0.4stable_gain_ceiling(:diffusion, s, dt) * spec.minimum
        @test tracking_bandwidth(:direct, s; dt) > tracking_bandwidth(:diffusion, s; dt)
    end
end
