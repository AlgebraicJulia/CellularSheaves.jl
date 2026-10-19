# The conic encoder, and the two composition-level helpers that need the assembled rows.
# Nothing here knows what a row means.

"""
    LinearRow(coefficients, bound; penalty = nothing)

One scalar inequality ``a^\\top u \\ge b`` of a local filter program.

A row with `penalty === nothing` is hard and is never violated. A row carrying a penalty
``p`` is soft: it is relaxed to ``a^\\top u \\ge b - \\delta`` for a variable
``\\delta \\ge 0`` charged in the objective.

`exact` selects how that charge is levied. The default ``p \\delta^2`` is an *inexact*
penalty: its optimum trades a little violation against the objective, so ``\\delta`` settles
at an ``O(1/p)`` floor even when the row is satisfiable, and the certified inequality is only
``a^\\top u \\ge b - O(1/p)``. Setting `exact = true` charges ``p \\delta`` instead, which
is an exact penalty: above a threshold on ``p`` the relaxation is driven to zero whenever the
unrelaxed row can be met, so the certificate is the unrelaxed one. The floor is worth
removing because near the origin it is a sizeable fraction of a decay rate proportional to
the state.

The exact form needs only a moderate ``p`` — that is the point of it — and the solver is
observed to fail above roughly ``10^4``, so it should not be paired with the very large
penalties an inexact charge would require.

Stability rows and barrier rows lower to the same conic block, so the encoder does not
distinguish them; a filter is composed by building the list of rows it needs and handing it
to [`filter_program`](@ref). A CLF row is the relaxable one, a barrier row is hard.
"""
struct LinearRow
    coefficients::Vector{Float64}
    bound::Float64
    penalty::Union{Nothing, Float64}
    exact::Bool
end

LinearRow(coefficients::AbstractVector, bound::Real; penalty = nothing, exact::Bool = false) =
    LinearRow(Vector{Float64}(coefficients), Float64(bound),
              penalty === nothing ? nothing : Float64(penalty), exact)

"""
    ConeRow(metric, bound; offset = nothing)

One second-order cone constraint ``\\|W u + v\\|_2 \\le b`` of a local filter program.

The offset is what lets a cone say something about where the state will *be* rather than
only about the size of the command: the discrete stability condition of
[`DiscreteStabilityTerm`](@ref) is a cone on ``L(A_d x + B_d u)``, whose constant part is
``L A_d x``. With `offset === nothing` this reduces to the actuator bound
``\\|W u\\|_2 \\le b``.
"""
struct ConeRow
    metric::Matrix{Float64}
    bound::Float64
    offset::Vector{Float64}
end

function ConeRow(metric::AbstractMatrix, bound::Real; offset = nothing)
    W = Matrix{Float64}(metric)
    v = offset === nothing ? zeros(size(W, 1)) : Vector{Float64}(offset)
    @argcheck length(v) == size(W, 1)
    @argcheck bound > 0
    return ConeRow(W, Float64(bound), v)
end

"""
    filter_program(nominal, rows; cones = ConeRow[], cone_metric = nothing, cone_bound = Inf)

Assemble the conic quadratic program

```math
\\min_{u, \\delta} \\ \\tfrac{1}{2}\\|u - u_{\\text{nom}}\\|^2 + \\sum_k p_k \\delta_k^2
\\quad\\text{s.t.}\\quad a_k^\\top u \\ge b_k - \\delta_k, \\quad \\|W u\\|_2 \\le u_{\\max},
```

over the [`LinearRow`](@ref)s in `rows`, and return `(problem, B)` where `B` is the
block-sparse constraint matrix. Column cell 1 holds ``u`` and the cells that follow hold one
relaxation per soft row, then one slack per row, then the cone variable when `cone_bound` is
finite, so `colrange(B, 1)` recovers the command and `colrange(B, 2)` the first relaxation.

Separating this from the constraint mathematics lets the same encoder serve any set of rows:
a horizon of them for a predictive filter, a different objective, or additional constraint
families, without touching the barrier and Lyapunov definitions that produce them.
"""
function filter_program(nominal::AbstractVector{<:Real}, rows::AbstractVector{LinearRow};
                        cones::AbstractVector{ConeRow} = ConeRow[],
                        cone_metric = nothing, cone_bound::Real = Inf)
    command_nom = Vector{Float64}(nominal)
    nu = length(command_nom)
    n = length(rows)
    @argcheck all(length(r.coefficients) == nu for r in rows)

    # The two-keyword form is the common actuator case; it is one cone with no offset.
    cone_list = collect(cones)
    if cone_metric !== nothing && isfinite(cone_bound)
        cone_list = vcat(cone_list, ConeRow(cone_metric, cone_bound))
    end
    @argcheck all(size(c.metric, 2) == nu for c in cone_list)

    soft = findall(r -> r.penalty !== nothing, rows)
    relaxation_of = Dict(k => 1 + i for (i, k) in enumerate(soft))
    slack_of(k) = 1 + length(soft) + k
    cone_cell(j) = 1 + length(soft) + n + j

    row_ids = Int[]
    col_ids = Int[]
    blocks = Matrix{Float64}[]
    for (k, row) in enumerate(rows)
        push!(row_ids, k); push!(col_ids, 1)
        push!(blocks, reshape(row.coefficients, 1, nu))
        if haskey(relaxation_of, k)
            push!(row_ids, k); push!(col_ids, relaxation_of[k]); push!(blocks, ones(1, 1))
        end
        push!(row_ids, k); push!(col_ids, slack_of(k)); push!(blocks, -ones(1, 1))
    end
    for (j, cone) in enumerate(cone_list)
        # The cone variable is (t; z) with ||z|| <= t; pinning t = bound and
        # z = W u + offset expresses ||W u + offset|| <= bound.
        push!(row_ids, n + j); push!(col_ids, 1)
        push!(blocks, -vcat(zeros(1, nu), cone.metric))
        push!(row_ids, n + j); push!(col_ids, cone_cell(j))
        push!(blocks, Matrix{Float64}(I, size(cone.metric, 1) + 1, size(cone.metric, 1) + 1))
    end

    B = blocksparse(row_ids, col_ids, blocks)
    Q = allocblockdiag(B)
    fill!(Q, 0)
    block(Q, 1, 1, 1) .= Matrix{Float64}(I, nu, nu)
    c = zeros(size(B, 2))
    c[1:nu] .= -command_nom
    for (k, cell) in relaxation_of
        if rows[k].exact
            # An exact (l1) penalty: the relaxation is charged linearly, so above a
            # threshold on p it is driven to exactly zero whenever the unrelaxed row is
            # satisfiable, rather than to the O(1/p) floor a quadratic charge leaves.
            c[first(colrange(B, cell))] = rows[k].penalty
        else
            block(Q, cell, cell, cell) .= reshape([2.0 * rows[k].penalty], 1, 1)
        end
    end
    g = Float64[row.bound for row in rows]
    for cone in cone_list
        append!(g, vcat(cone.bound, cone.offset))
    end
    cone_types = AbstractCone[CofreeCone()]
    append!(cone_types, (PositiveCone() for _ in 1:length(soft)))
    append!(cone_types, (PositiveCone() for _ in 1:n))
    append!(cone_types, (SecondOrderCone() for _ in 1:length(cone_list)))
    return IPMProblem(Q, B, c, g, cone_types), B
end

# A barrier row a'w >= b is unreachable, hence redundant, when the smallest value a'w can
# take over the actuator cone ||W u|| <= cap already exceeds b. That minimum is
# -cap * ||W^-T a||. Returns false whenever the cap is absent or W is not invertible, so
# the test never drops a row it cannot prove redundant.
function _row_is_redundant(row::AbstractVector{Float64}, rhs::Float64,
                           W::AbstractMatrix{Float64}, cap::Float64)
    isfinite(cap) || return false
    size(W, 1) == size(W, 2) || return false
    reach = try
        cap * norm(W' \ row)
    catch
        return false
    end
    isfinite(reach) || return false
    return -reach >= rhs + 1e-9 * max(1.0, abs(rhs))
end

# Rotate `u` by a quarter turn within each successive pair of coordinates. The map is skew,
# so the image is orthogonal to `u` and, crucially, reverses sign with `u`. Two agents of a
# pair see opposite barrier normals, so this hands them opposite evasion directions and the
# encounter resolves; an evasion built from each agent's own command would hand both the
# same direction and preserve the head-on symmetry.
function _quarter_turn(u::AbstractVector{Float64})
    t = zeros(length(u))
    for i in 1:2:(length(u) - 1)
        t[i] = -u[i + 1]
        t[i + 1] = u[i]
    end
    return t
end

# Perturb the nominal command tangentially to the most violated barrier row. Motion along
# the row normal is exactly the motion that changes the barrier, so its orthogonal
# complement slips past an obstacle without approaching it. Applied only when the nominal
# command actually violates a row, i.e. when the filter must act.
function _break_deadlock(u_nom::Vector{Float64}, A::Vector{Vector{Float64}},
                         b::Vector{Float64}, bias::Float64)
    residuals = [dot(A[i], u_nom) - b[i] for i in eachindex(A)]
    k = argmin(residuals)
    residuals[k] >= 0 && return u_nom

    scale = norm(A[k])
    scale <= 0 && return u_nom
    tangent = _quarter_turn(A[k] ./ scale)

    tangent_norm = norm(tangent)
    tangent_norm <= 0 && return u_nom
    return u_nom + bias * norm(u_nom) .* (tangent ./ tangent_norm)
end
