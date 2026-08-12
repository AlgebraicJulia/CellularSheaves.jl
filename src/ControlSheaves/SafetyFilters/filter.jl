# One solve attempt. A solver exception is the same outcome as a bad status, so both are
# funnelled to `nothing` and handled by the caller; an interrupt is not and is rethrown.
function _attempt(problem, settings)
    return try
        solve(problem, settings)
    catch err
        err isa InterruptException && rethrow()
        nothing
    end
end

# Composition. Terms are independent; the two operations below are not, and so live here
# rather than inside any one term:
#
#   * a barrier row no admissible command can violate is dropped, which needs the cone;
#   * a head-on deadlock is broken by perturbing the nominal command, which needs the rows.

"""
    SafetyFilterResult

Outcome of one filter solve.

`command` is `nothing` unless the returned control passed an independent primal check of
every hard constraint, so a `nothing` command means the execution layer must treat the
instant as infeasible or numerically unresolved rather than applying an uncertified control.
`cbf_margin` is the smallest barrier value seen, `Inf` when no barrier was in range. `slack`
is the largest relaxation taken by any soft row, so it reports the worst case in which
stability was yielded; each row is nonetheless certified against its own relaxation.
"""
struct SafetyFilterResult
    command::Union{Nothing, Vector{Float64}}
    slack::Float64
    status::Any
    lyapunov::Float64
    cbf_margin::Float64
    cbf_residual::Float64
    cap_residual::Float64
    clf_residual::Float64
    safety_certified::Bool
    clf_certified::Bool
    certified::Bool
end

"""
    safety_filter(terms, ctx, u_nom; settings, deadlock_bias = 0.0)

Compose the [`AbstractFilterTerm`](@ref)s over the context, solve the resulting conic
quadratic program, and certify the result against every hard constraint.

Each term is independent and may be tested on its own; this function only concatenates their
rows, drops the ones the actuator cone renders unreachable, and hands the rest to
[`filter_program`](@ref). Correctness therefore composes — but feasibility does not, since
hard barrier rows and a hard cone can have empty intersection.
"""
function safety_filter(terms, ctx::FilterContext, u_nom::AbstractVector{<:Real};
                       settings = safety_filter_settings(), deadlock_bias::Real = 0.0)
    command_nom = Vector{Float64}(u_nom)
    @argcheck length(command_nom) == ctx.nu
    @argcheck deadlock_bias >= 0

    soft = LinearRow[]
    hard = LinearRow[]
    cone = nothing
    extra_cones = ConeRow[]
    lyapunov = 0.0
    margin = Inf

    for term in terms
        contribution = constraints(term, ctx)
        for row in contribution.rows
            row.penalty === nothing ? push!(hard, row) : push!(soft, row)
        end
        if haskey(contribution, :cones)
            append!(extra_cones, contribution.cones)
        end
        if haskey(contribution, :cone) && contribution.cone !== nothing
            cone === nothing || throw(ArgumentError("at most one actuator term is supported"))
            cone = contribution.cone
        end
        diagnostics = contribution.diagnostics
        haskey(diagnostics, :lyapunov) && (lyapunov = diagnostics.lyapunov)
        if haskey(diagnostics, :barriers) && !isempty(diagnostics.barriers)
            margin = min(margin, minimum(diagnostics.barriers))
        end
    end

    W = cone === nothing ? Matrix{Float64}(I, ctx.nu, ctx.nu) : cone.metric
    cap = cone === nothing ? Inf : cone.bound
    # Rows no admissible command can violate are slack at every solution, and their large
    # slacks degrade the conditioning of the KKT system badly enough that the solve itself
    # can fail, so they are dropped. This is not merely an optimisation: retaining them
    # costs certified solves. The test is exact rather than heuristic, since over the cone
    # the least attainable value of a row is -cap ||W^-T a||.
    retained = filter(row -> !_row_is_redundant(row.coefficients, row.bound, W, cap), hard)

    if deadlock_bias > 0 && !isempty(retained)
        command_nom = _break_deadlock(command_nom,
                                      [row.coefficients for row in retained],
                                      [row.bound for row in retained], Float64(deadlock_bias))
    end

    rows = vcat(soft, retained)
    problem, B = filter_program(command_nom, rows; cones = extra_cones,
                                cone_metric = cone === nothing ? nothing : W,
                                cone_bound = cap)
    # A solver failure is a reportable outcome, not an error: it becomes an uncertified
    # result below. An interrupt is not, and must not be swallowed by that convention.
    # A failed solve is retried at the augmentations of `SAFETY_FILTER_FALLBACK_RAUG`, since
    # which augmentations condition badly is instance- and machine-dependent.
    result = _attempt(problem, settings)
    if result === nothing || !(result.status in (OPTIMAL, NEAR_OPTIMAL))
        for raug in SAFETY_FILTER_FALLBACK_RAUG
            retry = _with_raug(settings, raug)
            retry === settings && continue
            candidate = _attempt(problem, retry)
            if candidate !== nothing && candidate.status in (OPTIMAL, NEAR_OPTIMAL)
                result = candidate
                break
            end
        end
    end
    if result === nothing || !(result.status in (OPTIMAL, NEAR_OPTIMAL))
        return SafetyFilterResult(nothing, Inf,
                                  result === nothing ? :solver_error : result.status,
                                  lyapunov, margin, -Inf, -Inf, Inf, false, false, false)
    end

    command = Vector(result.p[colrange(B, 1)])
    # Soft rows lead `rows`, and `filter_program` gives the i-th of them its own relaxation in
    # column cell 1 + i. Each row must be certified against *its own* relaxation: sharing one
    # would certify a row whose relaxation is smaller than the one it was checked against.
    relaxations = [result.p[first(colrange(B, 1 + i))] for i in eachindex(soft)]
    raw_slack = isempty(relaxations) ? 0.0 : maximum(relaxations)
    slack = max(0.0, raw_slack)
    # Every row is certified in the same convention it was written: a'u >= b, relaxed by
    # delta where the row is soft.
    cbf_residual = isempty(retained) ? Inf :
        minimum(dot(row.coefficients, command) - row.bound for row in retained)
    cap_residual = isfinite(cap) ? cap - norm(W * command) : Inf
    clf_residual = isempty(soft) ? Inf :
        minimum(dot(soft[i].coefficients, command) + max(0.0, relaxations[i]) - soft[i].bound
                for i in eachindex(soft))

    tol = 10 * settings.feas_tol
    safety_certified = cbf_residual >= -tol && cap_residual >= -tol
    clf_certified = all(>=(-tol), relaxations) && clf_residual >= -tol
    certified = safety_certified && clf_certified
    return SafetyFilterResult(certified ? command : nothing, slack, result.status, lyapunov,
                              margin, cbf_residual, cap_residual, clf_residual,
                              safety_certified, clf_certified, certified)
end
