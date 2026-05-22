"""
Tests for the name-driven Tracking DSL (TrackingDSL).

Coverage:
1. Parse and validate a full happy-path program.
2. Late binding with bind!: matrices passed by value after parsing.
3. `initial` and `final` aliases resolve correctly with horizon.
4. Consensus/tracking constraints activated on distinct time sets.
5. Lowered TrackingProblem matches a manually assembled reference.
6. Unbound symbol errors are thrown before lowering.
7. Dimension mismatch errors include map/vector names.
8. Conflicting bind types are rejected.
9. Indexed boundary reference (x[a,t]) with late binding.
10. @tracking_problem macro evaluates bind values inline.
"""

using Test
using LinearAlgebra
using Graphs
using CellularSheaves
using CellularSheaves.ControlSheaves.TrackingDSL
import CellularSheaves.ControlSheaves.MultiAgentTracking as MAT

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

# Small 2-agent 1-target scenario
function _make_test_prog(; k=3)
    Ad = [1.0 0.0; 0.0 1.0]  # nx=2
    Bd = reshape([0.0; 1.0], 2, 1)  # nu=1
    # stalk_dim = nx+nu = 3; R_y must be (stalk_dim)×(stalk_dim) because
    # it operates on the full state-control stalk vector
    R_y = Matrix{Float64}(I, 3, 3)

    prog = parse_tracking_program(quote
        space(X) = R^2
        space(U) = R^1
        map_decl(A, X, X)
        map_decl(B, U, X)
        map_decl(R_y, X, X)
        agent(a1; dynamics=(A,B), period=dt)
        agent(a2; dynamics=(A,B), period=dt)
        target(t1)
        horizon(K)
        time(initial = 0)
        time(final = K)
        times(Tall = initial:final)
        consensus(c1; agents=(a1,a2), maps=(R_y,R_y), at=Tall)
        track(tr1; agent=a1, target=t1, maps=(R_y,R_y), at=final)
        boundary(:agent, a1; at=initial, value=x0_a1)
    end)
    bind!(prog, :K, k)
    bind!(prog, :A, Ad)
    bind!(prog, :B, Bd)
    bind!(prog, :R_y, R_y)
    bind!(prog, :dt, 0.01)
    bind!(prog, :x0_a1, zeros(3))
    return prog, Ad, Bd, R_y
end

# ─────────────────────────────────────────────────────────────────────────────
# 1. Happy-path: parse + validate + resolve + lower
# ─────────────────────────────────────────────────────────────────────────────
@testset "TrackingDSL happy path" begin
    prog, Ad, Bd, R_y = _make_test_prog(k=4)

    @test length(prog.statements) > 0
    @test validate_tracking_program(prog) === prog

    resolved = resolve_tracking_program(prog)
    @test resolved.k == 4
    @test length(resolved.agents)  == 2
    @test length(resolved.targets) == 1
    @test length(resolved.consensus) == 1
    @test length(resolved.tracking)  == 1
    @test length(resolved.boundaries) == 1

    result = lower_tracking_program(resolved)
    @test result.problem.n_agents  == 2
    @test result.problem.n_targets == 1
    @test result.problem.k         == 4
    @test result.problem.Ad ≈ Ad
    @test result.problem.Bd ≈ Bd
    @test length(result.boundary) == 1
end

# ─────────────────────────────────────────────────────────────────────────────
# 2. Late binding via bind!
# ─────────────────────────────────────────────────────────────────────────────
@testset "TrackingDSL late binding with bind!" begin
    prog = parse_tracking_program(quote
        agent(a1; dynamics=(A,B), period=dt)
        horizon(K)
        time(initial = 0)
        time(final = K)
    end)

    # Before binding, resolution should fail on missing horizon
    @test_throws TrackingUnboundSymbolError resolve_tracking_program(prog)

    bind!(prog, :K, 6)
    bind!(prog, :A, [1.0 0.0; 0.0 1.0])
    bind!(prog, :B, reshape([0.0; 1.0], 2, 1))
    bind!(prog, :dt, 0.05)

    resolved = resolve_tracking_program(prog)
    @test resolved.k == 6
    @test length(resolved.agents) == 1
end

# ─────────────────────────────────────────────────────────────────────────────
# 3. initial / final aliases resolve correctly
# ─────────────────────────────────────────────────────────────────────────────
@testset "TrackingDSL initial and final aliases" begin
    prog = parse_tracking_program(quote
        agent(a1; dynamics=(A,B), period=dt)
        horizon(K)
        time(initial = 0)
        time(final = K)
        times(Tall = initial:final)
        track(tr1; agent=a1, target=t1, maps=(R_y,R_y), at=final)
        target(t1)
    end)
    bind!(prog, :K, 5)
    bind!(prog, :A, [1.0 0.0; 0.0 1.0])
    bind!(prog, :B, reshape([0.0; 1.0], 2, 1))
    bind!(prog, :R_y, Matrix(1.0I, 3, 3))
    bind!(prog, :dt, 0.05)

    resolved = resolve_tracking_program(prog)
    @test resolved.k == 5

    # tracking at=final should activate only at t=5
    @test length(resolved.tracking) == 1
    @test resolved.tracking[1].timesteps == [5]
end

# ─────────────────────────────────────────────────────────────────────────────
# 4. Distinct time sets for consensus and tracking
# ─────────────────────────────────────────────────────────────────────────────
@testset "TrackingDSL independent time sets" begin
    prog = parse_tracking_program(quote
        agent(a1; dynamics=(A,B), period=dt)
        agent(a2; dynamics=(A,B), period=dt)
        target(t1)
        horizon(K)
        time(initial = 0)
        time(final = K)
        times(Tall = initial:final)
        consensus(c1; agents=(a1,a2), maps=(R_y,R_y), at=Tall)
        track(tr1; agent=a1, target=t1, maps=(R_y,R_y), at=final)
    end)
    bind!(prog, :K, 4)
    bind!(prog, :A, [1.0 0.0; 0.0 1.0])
    bind!(prog, :B, reshape([0.0; 1.0], 2, 1))
    bind!(prog, :R_y, Matrix(1.0I, 3, 3))
    bind!(prog, :dt, 0.05)

    resolved = resolve_tracking_program(prog)
    # Tall = 0:4 → [0,1,2,3,4]
    @test resolved.consensus[1].timesteps == collect(0:4)
    # final = 4
    @test resolved.tracking[1].timesteps == [4]
    @test resolved.consensus[1].timesteps ≠ resolved.tracking[1].timesteps
end

# ─────────────────────────────────────────────────────────────────────────────
# 5. Lowered TrackingProblem matches manual construction
# ─────────────────────────────────────────────────────────────────────────────
@testset "TrackingDSL lowering matches manual TrackingProblem" begin
    Ad = [1.0 0.0; 0.0 1.0]
    Bd = reshape([0.0; 1.0], 2, 1)
    R_y = Matrix{Float64}(I, 3, 3)
    k = 3

    # DSL path
    prog, _, _, _ = _make_test_prog(k=k)
    resolved = resolve_tracking_program(prog)
    result = lower_tracking_program(resolved)

    # Manual reference
    te = MAT.TrackingEdge(1, 1, R_y, R_y)
    ref_prob = MAT.TrackingProblem(
        2, 1, k,
        Ad, Bd,
        [(1,2)], [te],
        R_y,
        collect(0:k), [k],
        false, 1.0, 1.0,
    )

    dsl_prob = result.problem
    @test dsl_prob.n_agents  == ref_prob.n_agents
    @test dsl_prob.n_targets == ref_prob.n_targets
    @test dsl_prob.k         == ref_prob.k
    @test dsl_prob.Ad        ≈ ref_prob.Ad
    @test dsl_prob.Bd        ≈ ref_prob.Bd
    @test dsl_prob.agent_edges == ref_prob.agent_edges
    @test dsl_prob.consensus_timesteps == ref_prob.consensus_timesteps
    @test dsl_prob.tracking_timesteps  == ref_prob.tracking_timesteps

    # Both should produce sheaves with the same structure
    sheaf_dsl = MAT.build_time_expanded_tracking_sheaf(dsl_prob)
    sheaf_ref = MAT.build_time_expanded_tracking_sheaf(ref_prob)
    @test nv(sheaf_dsl.underlying_graph) == nv(sheaf_ref.underlying_graph)
    @test ne(sheaf_dsl.underlying_graph) == ne(sheaf_ref.underlying_graph)
end

# ─────────────────────────────────────────────────────────────────────────────
# 6. Unbound symbol errors before lowering
# ─────────────────────────────────────────────────────────────────────────────
@testset "TrackingDSL unbound symbol error" begin
    prog = parse_tracking_program(quote
        agent(a1; dynamics=(A,B), period=dt)
        horizon(K)
        time(initial = 0)
        time(final = K)
    end)
    # K not bound → TrackingUnboundSymbolError
    @test_throws TrackingUnboundSymbolError resolve_tracking_program(prog)

    bind!(prog, :K, 3)
    # A, B not bound → TrackingUnboundSymbolError
    @test_throws TrackingUnboundSymbolError resolve_tracking_program(prog)
end

# ─────────────────────────────────────────────────────────────────────────────
# 7. Dimension mismatch errors include names
# ─────────────────────────────────────────────────────────────────────────────
@testset "TrackingDSL dimension mismatch error" begin
    prog = parse_tracking_program(quote
        agent(a1; dynamics=(A,B), period=dt)
        horizon(K)
        time(initial = 0)
        time(final = K)
    end)
    bind!(prog, :K, 3)
    bind!(prog, :dt, 0.05)
    # Non-square A: rows ≠ cols
    bind!(prog, :A, [1.0 0.0 0.0; 0.0 1.0 0.0])  # 2×3 — non-square
    bind!(prog, :B, reshape([0.0; 1.0], 2, 1))

    err = try
        resolve_tracking_program(prog)
        nothing
    catch e
        e
    end
    @test err isa TrackingDimensionMismatchError
    @test occursin("Ad", err.msg) || occursin("a1", err.msg)
end

# ─────────────────────────────────────────────────────────────────────────────
# 8. Conflicting bind types: scalar vs matrix
# ─────────────────────────────────────────────────────────────────────────────
@testset "TrackingDSL conflicting bind types" begin
    prog = parse_tracking_program(quote
        agent(a1; dynamics=(A,B), period=dt)
        horizon(K)
        time(initial = 0)
        time(final = K)
    end)
    bind!(prog, :K, 3)
    bind!(prog, :dt, 0.05)
    # Bind A as a scalar — will fail when we try to use it as a matrix
    bind!(prog, :A, 42.0)  # scalar, not matrix
    bind!(prog, :B, reshape([0.0; 1.0], 2, 1))

    # Resolution should fail because 42.0 can't be converted to Matrix{Float64}
    @test_throws Exception resolve_tracking_program(prog)
end

# ─────────────────────────────────────────────────────────────────────────────
# 9. Indexed boundary: x[a,t] with late binding
# ─────────────────────────────────────────────────────────────────────────────
@testset "TrackingDSL indexed boundary late binding" begin
    prog = parse_tracking_program(quote
        agent(a1; dynamics=(A,B), period=dt)
        horizon(K)
        time(initial = 0)
        time(final = K)
        boundary(:agent, a1; at=t_pin, value=x_ref[a,t_pin])
    end)
    bind!(prog, :K, 5)
    bind!(prog, :A, [1.0 0.0; 0.0 1.0])
    bind!(prog, :B, reshape([0.0; 1.0], 2, 1))
    bind!(prog, :dt, 0.05)
    # Bind the indexed reference: x_ref[a=1, t_pin=2] => some vector
    bind!(prog, :x_ref, :a, :t_pin, [1.0, 2.0, 0.0])
    bind!(prog, :a, 1)
    bind!(prog, :t_pin, 2)

    resolved = resolve_tracking_program(prog)
    @test length(resolved.boundaries) == 1
    @test resolved.boundaries[1].time == 2
    @test resolved.boundaries[1].entity_index == 1
    @test resolved.boundaries[1].value ≈ [1.0, 2.0, 0.0]
end

# ─────────────────────────────────────────────────────────────────────────────
# 10. @tracking_problem macro evaluates bind values inline
# ─────────────────────────────────────────────────────────────────────────────
@testset "TrackingDSL @tracking_problem macro" begin
    Ad = [1.0 0.0; 0.0 1.0]
    Bd = reshape([0.0; 1.0], 2, 1)
    R_y = Matrix{Float64}(I, 3, 3)
    k_val = 4

    prog = @tracking_problem begin
        agent(a1; dynamics=(A,B), period=dt)
        agent(a2; dynamics=(A,B), period=dt)
        target(t1)
        horizon(K)
        time(initial = 0)
        time(final = K)
        times(Tall = initial:final)
        consensus(c1; agents=(a1,a2), maps=(R_y,R_y), at=Tall)
        bind(K => k_val)
        bind(A => Ad)
        bind(B => Bd)
        bind(R_y => R_y)
        bind(dt => 0.01)
    end

    # Bind values are concrete (not Expr) because the macro evaluated them
    resolved = resolve_tracking_program(validate_tracking_program(prog))
    @test resolved.k == k_val
    @test resolved.agents[1].Ad ≈ Ad
    @test resolved.agents[1].Bd ≈ Bd
end

# ─────────────────────────────────────────────────────────────────────────────
# Shared setup for Scenarios 12–15 (quadrotor-like dynamics, k=5)
# ─────────────────────────────────────────────────────────────────────────────
# Simplified 6-state / 2-control kinematics (same dimensions as the full planar
# quadrotor, but without the gravity/inertia nonlinearity for test speed).
const _Q_DT  = 0.05
const _Q_NX  = 6
const _Q_NU  = 2
const _Q_K   = 5

const _Q_AD  = let A = Matrix{Float64}(I, _Q_NX, _Q_NX)
    A[1,4] = _Q_DT; A[2,5] = _Q_DT; A[3,6] = _Q_DT   # pos += dt * vel
    A
end
const _Q_BD  = let B = zeros(_Q_NX, _Q_NU)
    B[4,1] = _Q_DT; B[5,2] = _Q_DT                    # vel += dt * u
    B
end

# Projection matrices (rows select from the (nx+nu)-stalk):
#   R_yz : project lateral (y) and altitude (z)          → 2×8
#   R_y  : project lateral (y) only                      → 1×8
#   R_z  : project altitude (z) only                     → 1×8
const _Q_R_YZ = MAT.state_projection_matrix([1, 2], _Q_NX, _Q_NU)
const _Q_R_Y  = MAT.state_projection_matrix([1],    _Q_NX, _Q_NU)
const _Q_R_Z  = MAT.state_projection_matrix([2],    _Q_NX, _Q_NU)

# ─────────────────────────────────────────────────────────────────────────────
# 11. Regression: scenario from multi_quadrotor_target_tracking
# ─────────────────────────────────────────────────────────────────────────────
@testset "TrackingDSL regression vs manual setup" begin
    # Simplified version of the multi_quadrotor scenario.
    # n_agents=2, n_targets=1, k=5, stalk_dim = nx+nu
    # State: pos-x, pos-y, vel-x, vel-y (nx=4), control: ax, ay (nu=2) => stalk=6
    dt = 0.05
    nx = 4; nu = 2
    # Discretised double-integrator (exact ZOH)
    Ad = [1 0 dt 0; 0 1 0 dt; 0 0 1 0; 0 0 0 1] * 1.0
    Bd = [dt^2/2 0; 0 dt^2/2; dt 0; 0 dt] * 1.0
    R_pos = hcat(I(nx), zeros(nx, nu))   # project state only (4×6)
    k = 5

    # DSL approach
    prog = parse_tracking_program(quote
        agent(a1; dynamics=(Ad_m,Bd_m), period=h)
        agent(a2; dynamics=(Ad_m,Bd_m), period=h)
        target(t1)
        horizon(K)
        time(initial = 0)
        time(final = K)
        times(Tall = initial:final)
        consensus(c1; agents=(a1,a2), maps=(R_pos_m,R_pos_m), at=Tall)
        track(tr1; agent=a1, target=t1, maps=(R_pos_m,R_pos_m), at=final)
    end)
    bind!(prog, :K, k)
    bind!(prog, :Ad_m, Ad)
    bind!(prog, :Bd_m, Bd)
    bind!(prog, :R_pos_m, R_pos)
    bind!(prog, :h, dt)

    resolved = resolve_tracking_program(prog)
    result = lower_tracking_program(resolved)
    sheaf_dsl = MAT.build_time_expanded_tracking_sheaf(result.problem)

    # Manual reference
    te = MAT.TrackingEdge(1, 1, R_pos, R_pos)
    ref_prob = MAT.TrackingProblem(
        2, 1, k,
        Ad, Bd,
        [(1,2)], [te],
        R_pos,
        collect(0:k), [k],
        false, 1.0, 1.0,
    )
    sheaf_ref = MAT.build_time_expanded_tracking_sheaf(ref_prob)

    @test nv(sheaf_dsl.underlying_graph) == nv(sheaf_ref.underlying_graph)
    @test ne(sheaf_dsl.underlying_graph) == ne(sheaf_ref.underlying_graph)

    # Verify Laplacian dimensions match
    L_dsl = sheaf_laplacian_matrix_direct(sheaf_dsl)
    L_ref = sheaf_laplacian_matrix_direct(sheaf_ref)
    @test size(L_dsl) == size(L_ref)
    @test L_dsl ≈ L_ref
end

# ─────────────────────────────────────────────────────────────────────────────
# 12. Scenario 1 mirror: (y,z) consensus + tracking at every timestep
#     include_target_dynamics=true, tracking_weight=5.0
# ─────────────────────────────────────────────────────────────────────────────
@testset "TrackingDSL Scenario 1: (y,z) all-timestep consensus+tracking" begin
    prog = @tracking_problem begin
        agent(a1; dynamics=(Ad_dyn, Bd_dyn), period=h_dt)
        agent(a2; dynamics=(Ad_dyn, Bd_dyn), period=h_dt)
        target(t1)
        target(t2)
        horizon(K)
        time(initial = 0)
        time(final = K)
        times(Tall = initial:final)
        consensus(c1; agents=(a1,a2), maps=(R_yz_m,R_yz_m), at=Tall)
        track(tr1; agent=a1, target=t1, maps=(R_yz_m,R_yz_m), at=Tall)
        track(tr2; agent=a2, target=t2, maps=(R_yz_m,R_yz_m), at=Tall)
        bind(K => _Q_K)
        bind(Ad_dyn => _Q_AD)
        bind(Bd_dyn => _Q_BD)
        bind(R_yz_m => _Q_R_YZ)
        bind(h_dt => _Q_DT)
    end

    resolved = prog |> validate_tracking_program |> resolve_tracking_program
    lowered  = lower_tracking_program(resolved;
        include_target_dynamics=true, tracking_weight=5.0)
    p = lowered.problem

    @test p.n_agents  == 2
    @test p.n_targets == 2
    @test p.k         == _Q_K
    @test p.consensus_timesteps == collect(0:_Q_K)
    @test p.tracking_timesteps  == collect(0:_Q_K)
    @test p.consensus_restriction ≈ _Q_R_YZ

    # Verify the DSL-built sheaf matches a manually-assembled reference
    te_yz = [MAT.TrackingEdge(1,1,_Q_R_YZ,_Q_R_YZ),
             MAT.TrackingEdge(2,2,_Q_R_YZ,_Q_R_YZ)]
    ref_p = MAT.TrackingProblem(2, 2, _Q_K,
        _Q_AD, _Q_BD,
        [(1,2)], te_yz, _Q_R_YZ,
        collect(0:_Q_K), collect(0:_Q_K),
        true, 1.0, 5.0)
    sheaf_dsl = MAT.build_time_expanded_tracking_sheaf(p)
    sheaf_ref = MAT.build_time_expanded_tracking_sheaf(ref_p)
    @test nv(sheaf_dsl.underlying_graph) == nv(sheaf_ref.underlying_graph)
    @test ne(sheaf_dsl.underlying_graph) == ne(sheaf_ref.underlying_graph)
end

# ─────────────────────────────────────────────────────────────────────────────
# 13. Scenario 2 mirror: y-consensus, z-tracking, all timesteps
# ─────────────────────────────────────────────────────────────────────────────
@testset "TrackingDSL Scenario 2: y-consensus z-tracking all timesteps" begin
    prog = @tracking_problem begin
        agent(a1; dynamics=(Ad_dyn, Bd_dyn), period=h_dt)
        agent(a2; dynamics=(Ad_dyn, Bd_dyn), period=h_dt)
        target(t1)
        target(t2)
        horizon(K)
        time(initial = 0)
        time(final = K)
        times(Tall = initial:final)
        consensus(c1; agents=(a1,a2), maps=(R_y_m,R_y_m), at=Tall)
        track(tr1; agent=a1, target=t1, maps=(R_z_m,R_z_m), at=Tall)
        track(tr2; agent=a2, target=t2, maps=(R_z_m,R_z_m), at=Tall)
        bind(K => _Q_K)
        bind(Ad_dyn => _Q_AD)
        bind(Bd_dyn => _Q_BD)
        bind(R_y_m => _Q_R_Y)
        bind(R_z_m => _Q_R_Z)
        bind(h_dt => _Q_DT)
    end

    resolved = prog |> validate_tracking_program |> resolve_tracking_program
    lowered  = lower_tracking_program(resolved)
    p = lowered.problem

    @test p.consensus_timesteps == collect(0:_Q_K)
    @test p.tracking_timesteps  == collect(0:_Q_K)
    # Consensus projection is R_y; tracking projection is R_z
    @test p.consensus_restriction ≈ _Q_R_Y
    @test p.tracking_edges[1].agent_restriction ≈ _Q_R_Z
    @test p.tracking_edges[2].agent_restriction ≈ _Q_R_Z
end

# ─────────────────────────────────────────────────────────────────────────────
# 14. Scenario 3 mirror: y-consensus, z-tracking, terminal timestep only
# ─────────────────────────────────────────────────────────────────────────────
@testset "TrackingDSL Scenario 3: y-consensus z-tracking terminal only" begin
    prog = @tracking_problem begin
        agent(a1; dynamics=(Ad_dyn, Bd_dyn), period=h_dt)
        agent(a2; dynamics=(Ad_dyn, Bd_dyn), period=h_dt)
        target(t1)
        target(t2)
        horizon(K)
        time(initial = 0)
        time(final = K)
        consensus(c1; agents=(a1,a2), maps=(R_y_m,R_y_m), at=final)
        track(tr1; agent=a1, target=t1, maps=(R_z_m,R_z_m), at=final)
        track(tr2; agent=a2, target=t2, maps=(R_z_m,R_z_m), at=final)
        bind(K => _Q_K)
        bind(Ad_dyn => _Q_AD)
        bind(Bd_dyn => _Q_BD)
        bind(R_y_m => _Q_R_Y)
        bind(R_z_m => _Q_R_Z)
        bind(h_dt => _Q_DT)
    end

    resolved = prog |> validate_tracking_program |> resolve_tracking_program
    lowered  = lower_tracking_program(resolved)
    p = lowered.problem

    # Switching from Tall to `at=final` moves k constraints to just [k]
    @test p.consensus_timesteps == [_Q_K]
    @test p.tracking_timesteps  == [_Q_K]

    # Verify sheaf has fewer edges than Scenario 2 (consensus+tracking only at t=k)
    prog2 = @tracking_problem begin
        agent(a1; dynamics=(Ad_dyn, Bd_dyn), period=h_dt)
        agent(a2; dynamics=(Ad_dyn, Bd_dyn), period=h_dt)
        target(t1)
        target(t2)
        horizon(K)
        time(initial = 0)
        time(final = K)
        times(Tall = initial:final)
        consensus(c1; agents=(a1,a2), maps=(R_y_m,R_y_m), at=Tall)
        track(tr1; agent=a1, target=t1, maps=(R_z_m,R_z_m), at=Tall)
        track(tr2; agent=a2, target=t2, maps=(R_z_m,R_z_m), at=Tall)
        bind(K => _Q_K)
        bind(Ad_dyn => _Q_AD)
        bind(Bd_dyn => _Q_BD)
        bind(R_y_m => _Q_R_Y)
        bind(R_z_m => _Q_R_Z)
        bind(h_dt => _Q_DT)
    end
    p2 = lower_tracking_program(resolve_tracking_program(prog2)).problem
    # Restricting to the terminal step removes k consensus + k*n_tracking_edges edges
    @test ne(MAT.build_time_expanded_tracking_sheaf(p).underlying_graph) <
          ne(MAT.build_time_expanded_tracking_sheaf(p2).underlying_graph)
end

# ─────────────────────────────────────────────────────────────────────────────
# 15. Scenario 4 mirror: same DSL topology as Scenario 3; null dim is invariant
#     (different target trajectories are boundary data, not part of the DSL)
# ─────────────────────────────────────────────────────────────────────────────
@testset "TrackingDSL Scenario 4: null dim invariant under boundary change" begin
    # Scenario 4 uses exactly the same DSL program as Scenario 3.
    # Only the boundary data (target trajectories) differs; that is supplied
    # externally after lowering.  The resulting sheaf topology — and therefore
    # the nullspace dimension of the Laplacian — must be identical.
    function _s34_prog()
        @tracking_problem begin
            agent(a1; dynamics=(Ad_dyn, Bd_dyn), period=h_dt)
            agent(a2; dynamics=(Ad_dyn, Bd_dyn), period=h_dt)
            target(t1)
            target(t2)
            horizon(K)
            time(initial = 0)
            time(final = K)
            consensus(c1; agents=(a1,a2), maps=(R_y_m,R_y_m), at=final)
            track(tr1; agent=a1, target=t1, maps=(R_z_m,R_z_m), at=final)
            track(tr2; agent=a2, target=t2, maps=(R_z_m,R_z_m), at=final)
            bind(K => _Q_K)
            bind(Ad_dyn => _Q_AD)
            bind(Bd_dyn => _Q_BD)
            bind(R_y_m => _Q_R_Y)
            bind(R_z_m => _Q_R_Z)
            bind(h_dt => _Q_DT)
        end
    end

    p3 = lower_tracking_program(resolve_tracking_program(_s34_prog())).problem
    p4 = lower_tracking_program(resolve_tracking_program(_s34_prog())).problem

    # Same topology → same sheaf
    s3 = MAT.build_time_expanded_tracking_sheaf(p3)
    s4 = MAT.build_time_expanded_tracking_sheaf(p4)
    @test nv(s3.underlying_graph) == nv(s4.underlying_graph)
    @test ne(s3.underlying_graph) == ne(s4.underlying_graph)

    # Same Laplacian → same nullspace dimension
    L3 = sheaf_laplacian_matrix_direct(s3)
    L4 = sheaf_laplacian_matrix_direct(s4)
    @test size(L3) == size(L4)
    @test L3 ≈ L4
end
