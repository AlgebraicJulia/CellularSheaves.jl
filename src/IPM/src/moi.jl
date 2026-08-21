#
# MathOptInterface wrapper.
#
# The optimizer's native input is the IPM primal form
#
#   min ½ pᵀQp − fᵀp   s.t.   Bp = g,   p ∈ K,
#
# transcribed directly onto MOI:
#
#   min ½ pᵀQp − fᵀp   ↔   ScalarQuadraticFunction objective
#   Bp = g             ↔   VectorAffineFunction-in-Zeros
#   p ∈ K              ↔   VectorOfVariables-in-{cone}
#
# Everything outside this form (affine-in-cone, scalar constraints, bounds,
# …) is delivered to us by MOI's bridges.
#
import MathOptInterface as MOI

############################################################################################
# constants
############################################################################################

# cones we accept a VectorOfVariables in — one cone block each
const SupportedCone = Union{
    MOI.Nonnegatives,
    MOI.SecondOrderCone,
    MOI.ExponentialCone,
    MOI.PowerCone,
    MOI.Scaled{MOI.PositiveSemidefiniteConeTriangle},
}

const _TERM = Dict(
    OPTIMAL                 => MOI.OPTIMAL,
    NEAR_OPTIMAL            => MOI.ALMOST_OPTIMAL,
    STALLED                 => MOI.SLOW_PROGRESS,
    NUMERICAL_FAILURE       => MOI.NUMERICAL_ERROR,
    ITERATION_LIMIT         => MOI.ITERATION_LIMIT,
    PRIMAL_INFEASIBLE       => MOI.INFEASIBLE,
    DUAL_INFEASIBLE         => MOI.DUAL_INFEASIBLE,
    ILL_POSED               => MOI.INVALID_MODEL,
    NEAR_PRIMAL_INFEASIBLE  => MOI.ALMOST_INFEASIBLE,
    NEAR_DUAL_INFEASIBLE    => MOI.ALMOST_DUAL_INFEASIBLE,
    NEAR_ILL_POSED          => MOI.ALMOST_INFEASIBLE,
)

############################################################################################
# Optimizer
############################################################################################

mutable struct Optimizer{T} <: MOI.AbstractOptimizer
    settings::IPMSettings{T}
    # the assembled primal-form problem (built by copy_to, solved by optimize!)
    problem::Union{Nothing, IPMProblem{T}}
    result::Union{Nothing, IPMResult{T}}
    jmap::Dict{MOI.VariableIndex, Int}                  # MOI variable → IPM column j
    # CSC-style cone layout: vmap gives each cone's slice (colptr), vidx the flat
    # IPM columns in MOI coordinate order (row indices) read back through it
    vmap::Dict{MOI.ConstraintIndex, UnitRange{Int}}     # cone (block-column v) → its slice in vidx
    vidx::Vector{Int}
    umap::Dict{MOI.ConstraintIndex, UnitRange{Int}}     # Zeros (block-row u) → IPM (B, y) rows (irng)
    sense::MOI.OptimizationSense
    objconst::T
    solvetime::Float64

    function Optimizer{T}() where {T}
        return new{T}(
            IPMSettings{T}(), nothing, nothing,
            Dict{MOI.VariableIndex, Int}(),
            Dict{MOI.ConstraintIndex, UnitRange{Int}}(),
            Int[],
            Dict{MOI.ConstraintIndex, UnitRange{Int}}(),
            MOI.MIN_SENSE, zero(T), 0.0,
        )
    end
end

function Optimizer()
    return Optimizer{Float64}()
end

# +1 for MIN (IPM's native sense), −1 for MAX
function _sensesign(::Type{T}, sense::MOI.OptimizationSense) where {T}
    if sense == MOI.MAX_SENSE
        return -one(T)
    else
        return one(T)
    end
end

function _sensesign(opt::Optimizer{T}) where {T}
    return _sensesign(T, opt.sense)
end

############################################################################################
# ipmcone
############################################################################################

# IPM cone for a MOI set (dimension comes from the block width)
function ipmcone(::MOI.Nonnegatives)
    return PositiveCone()
end

# MOI's SecondOrderCone is (t, x) with t ≥ ‖x‖; the IPM cone's coordinate 1 is
# the same scalar bound, so the block maps across with no reordering
function ipmcone(::MOI.SecondOrderCone)
    return SecondOrderCone()
end

# MOI's ExponentialCone (x,y,z) is y·exp(x/y) ≤ z; the IPM cone (a,b,c) is
# b·log(a/b) ≥ c, i.e. (a,b,c) = (z,y,x) — the reversed triple (see coneperm)
function ipmcone(::MOI.ExponentialCone)
    return ExponentialCone()
end

# MOI's PowerCone(α) (x,y,z): x^α·y^(1-α) ≥ |z| matches the IPM cone directly
function ipmcone(set::MOI.PowerCone)
    return PowerCone(set.exponent)
end

# MOI's Scaled PSD-triangle already carries the √2 off-diagonal scaling the IPM
# svec uses; only the ordering differs (upper vs lower col-major, see coneperm)
function ipmcone(::MOI.Scaled{MOI.PositiveSemidefiniteConeTriangle})
    return SemidefiniteCone()
end

# coneperm(set) gives the order in which MOI coordinates are appended to build
# the IPM cone's coordinate layout: identity where the conventions agree, a
# reversal for the exponential cone.
function coneperm(set::MOI.Nonnegatives)
    return 1:MOI.dimension(set)
end

function coneperm(set::MOI.SecondOrderCone)
    return 1:MOI.dimension(set)
end

function coneperm(::MOI.ExponentialCone)
    return (3, 2, 1)
end

function coneperm(set::MOI.PowerCone)
    return 1:MOI.dimension(set)
end

# MOI packs the PSD triangle upper-column-major (entry (i,j), i≤j, at index
# j(j-1)/2 + i); the IPM svec packs lower-column-major (column j: (j,j) then
# (i,j) for i>j). Emit the MOI index of each IPM slot's entry.
function coneperm(set::MOI.Scaled{MOI.PositiveSemidefiniteConeTriangle})
    d = MOI.side_dimension(set)
    order = Int[]
    for j in 1:d
        for i in j:d
            push!(order, div((i - 1) * i, 2) + j)
        end
    end
    return order
end

############################################################################################
# solver attributes
############################################################################################

function MOI.get(::Optimizer, ::MOI.SolverName)
    return "CellularSheaves IPM"
end

function MOI.supports_incremental_interface(::Optimizer)
    return false
end

function MOI.is_empty(opt::Optimizer)
    return isnothing(opt.result) && isempty(opt.jmap)
end

function MOI.empty!(opt::Optimizer{T}) where {T}
    opt.problem = nothing
    opt.result = nothing
    opt.jmap = Dict{MOI.VariableIndex, Int}()
    opt.vmap = Dict{MOI.ConstraintIndex, UnitRange{Int}}()
    opt.vidx = Int[]
    opt.umap = Dict{MOI.ConstraintIndex, UnitRange{Int}}()
    opt.sense = MOI.MIN_SENSE
    opt.objconst = zero(T)
    opt.solvetime = 0.0
    return
end

function MOI.supports(::Optimizer, ::MOI.Silent)
    return true
end

# Silent is just a view onto verbose: silent ⇔ verbose == 0; un-silencing
# restores the iteration log (verbose = 1)
function MOI.get(opt::Optimizer, ::MOI.Silent)
    return opt.settings.verbose == 0
end

function MOI.set(opt::Optimizer, ::MOI.Silent, v::Bool)
    opt.settings = IPMSettings(opt.settings; verbose = v ? 0 : 1)
    return
end

function MOI.supports(::Optimizer, ::MOI.RawOptimizerAttribute)
    return true
end

function MOI.get(opt::Optimizer, a::MOI.RawOptimizerAttribute)
    return getfield(opt.settings, Symbol(a.name))
end

function MOI.set(opt::Optimizer, a::MOI.RawOptimizerAttribute, v)
    opt.settings = IPMSettings(opt.settings; Symbol(a.name) => v)
    return
end

############################################################################################
# supported objective / constraints
############################################################################################

function MOI.supports(::Optimizer, ::MOI.ObjectiveSense)
    return true
end

function MOI.supports(::Optimizer{T}, ::MOI.ObjectiveFunction{<:Union{MOI.ScalarAffineFunction{T}, MOI.ScalarQuadraticFunction{T}}}) where {T}
    return true
end

function MOI.supports_constraint(::Optimizer{T}, ::Type{MOI.VectorAffineFunction{T}}, ::Type{MOI.Zeros}) where {T}
    return true
end

function MOI.supports_constraint(::Optimizer, ::Type{MOI.VectorOfVariables}, ::Type{<:SupportedCone})
    return true
end

############################################################################################
# copy_to
############################################################################################

#
# lay out the columns: each VoV-in-cone block contiguous, then the free
# variables as a trailing CofreeCone block. Returns (jmap, K, s, vmap, vidx,
# jdup, ncol), where jmap maps each MOI variable to its (first) IPM column,
# (vmap, vidx) are the CSC-style cone layout — vmap gives each cone's slice and
# vidx the flat IPM columns in MOI coordinate order — and ncol is the total
# column count.
#
# A variable already claimed by one cone that reappears in a second cone can't
# share a column, so the second block gets a fresh synthesized column z and we
# record jdup[z] = x: _equalities pins z = x with an internal row. This is the
# slack reformulation MOI's VectorSlackBridge applies to affine-in-cone, done
# here for the variable-in-two-cones case.
#
function _layout(src::MOI.ModelLike)
    vars = MOI.get(src, MOI.ListOfVariableIndices())

    K = AbstractCone[]
    s = Int[]
    jmap = Dict{MOI.VariableIndex, Int}()
    vmap = Dict{MOI.ConstraintIndex, UnitRange{Int}}()
    vidx = Int[]
    jdup = Dict{Int, Int}()   # synthesized column z → the claimed column x it copies
    j = 0

    for (F, S) in MOI.get(src, MOI.ListOfConstraintTypesPresent())
        if F === MOI.VectorOfVariables && S <: SupportedCone
            for vx in MOI.get(src, MOI.ListOfConstraintIndices{F, S}())
                fun = MOI.get(src, MOI.ConstraintFunction(), vx)
                set = MOI.get(src, MOI.ConstraintSet(), vx)

                jstrt = j
                d = length(fun.variables)
                resize!(vidx, jstrt + d)

                # lay members out in the IPM cone's coordinate order; write each
                # MOI coordinate k's column into its vidx slot for ordered readback
                for k in coneperm(set)
                    j += 1
                    jx = fun.variables[k]

                    if haskey(jmap, jx)
                        jdup[j] = jmap[jx]
                    else
                        jmap[jx] = j
                    end

                    vidx[jstrt + k] = j
                end

                push!(K, ipmcone(set))
                push!(s, d)
                vmap[vx] = jstrt + 1:j
            end
        end
    end

    jstrt = j

    for jx in vars
        if !haskey(jmap, jx)
            j += 1
            jmap[jx] = j
        end
    end

    if jstrt < j
        push!(K, CofreeCone())
        push!(s, j - jstrt)
    end

    return jmap, K, s, vmap, vidx, jdup, j
end

#
# equalities Bp = g from VectorAffineFunction-in-Zeros, plus the internal z = x
# rows pinning each synthesized column to the variable it copies. Returns
# (B, g, umap), where umap maps each Zeros constraint (block-row u) to its rows
# in (B, y); the coupling rows carry no umap entry (they aren't MOI constraints).
#
function _equalities(::Type{T}, src::MOI.ModelLike, jmap, ncol, jdup) where {T}
    I = Int[]
    J = Int[]
    V = T[]
    g = T[]
    umap = Dict{MOI.ConstraintIndex, UnitRange{Int}}()
    i = 0

    for ux in MOI.get(src, MOI.ListOfConstraintIndices{MOI.VectorAffineFunction{T}, MOI.Zeros}())
        fun = MOI.get(src, MOI.ConstraintFunction(), ux)
        d = MOI.output_dimension(fun)

        for t in fun.terms
            push!(I, i + t.output_index)
            push!(J, jmap[t.scalar_term.variable])
            push!(V, t.scalar_term.coefficient)
        end

        for k in 1:d                 # Ax + b = 0  ⇒  Ax = -b
            push!(g, -fun.constants[k])
        end

        umap[ux] = i + 1:i + d
        i += d
    end

    for (j, k) in jdup               # column j copies column k:  j − k = 0
        i += 1
        push!(I, i); push!(J, j); push!(V,  one(T))
        push!(I, i); push!(J, k); push!(V, -one(T))
        push!(g, zero(T))
    end

    return sparse(I, J, V, i, ncol), g, umap
end

#
# objective ½pᵀQp − fᵀp. Returns (Q, f, sense, objconst).
#
function _objective(::Type{T}, src::MOI.ModelLike, jmap, ncol) where {T}
    sense = MOI.get(src, MOI.ObjectiveSense())
    F = MOI.get(src, MOI.ObjectiveFunctionType())
    obj = MOI.get(src, MOI.ObjectiveFunction{F}())

    I = Int[]
    J = Int[]
    V = T[]
    f = zeros(T, ncol)

    σ = _sensesign(T, sense)   # IPM minimizes

    if obj isa MOI.ScalarQuadraticFunction
        for term in obj.quadratic_terms
            a = jmap[term.variable_1]
            b = jmap[term.variable_2]

            push!(I, a)
            push!(J, b)
            push!(V, σ * term.coefficient)

            if a != b
                push!(I, b)
                push!(J, a)
                push!(V, σ * term.coefficient)
            end
        end

        aff = obj.affine_terms
    else
        aff = obj.terms
    end

    for term in aff
        f[jmap[term.variable]] -= σ * term.coefficient   # min σ·(aᵀp) = −fᵀp
    end

    return sparse(I, J, V, ncol, ncol), f, sense, MOI.constant(obj)
end

#
# copy_to: read the model in primal form into an IPMProblem
#
function MOI.copy_to(dest::Optimizer{T}, src::MOI.ModelLike) where {T}
    MOI.empty!(dest)

    jmap, K, s, vmap, vidx, jdup, ncol = _layout(src)
    B, g, umap = _equalities(T, src, jmap, ncol, jdup)
    Q, f, sense, objconst = _objective(T, src, jmap, ncol)

    dest.jmap = jmap
    dest.vmap = vmap
    dest.vidx = vidx
    dest.umap = umap
    dest.sense = sense
    dest.objconst = objconst
    dest.problem = IPMProblem(Q, B, f, g, K, s)
    return MOI.Utilities.identity_index_map(src)
end

############################################################################################
# optimize!
############################################################################################

function MOI.optimize!(opt::Optimizer)
    t0 = time()
    opt.result = solve(opt.problem, opt.settings)
    opt.solvetime = time() - t0
    return
end

############################################################################################
# results
############################################################################################

function MOI.get(opt::Optimizer, ::MOI.RawStatusString)
    return string(opt.result.status)
end

function MOI.get(opt::Optimizer, ::MOI.SolveTimeSec)
    return opt.solvetime
end

function MOI.get(opt::Optimizer, ::MOI.ResultCount)
    return isnothing(opt.result) ? 0 : 1
end

function MOI.get(opt::Optimizer, ::MOI.TerminationStatus)
    isnothing(opt.result) && return MOI.OPTIMIZE_NOT_CALLED
    return get(_TERM, opt.result.status, MOI.OTHER_ERROR)
end

function MOI.get(opt::Optimizer, attr::MOI.PrimalStatus)
    ok = !isnothing(opt.result) && opt.result.status in (OPTIMAL, NEAR_OPTIMAL)
    (attr.result_index == 1 && ok) || return MOI.NO_SOLUTION
    return MOI.FEASIBLE_POINT
end

function MOI.get(opt::Optimizer, attr::MOI.DualStatus)
    ok = !isnothing(opt.result) && opt.result.status in (OPTIMAL, NEAR_OPTIMAL)
    (attr.result_index == 1 && ok) || return MOI.NO_SOLUTION
    return MOI.FEASIBLE_POINT
end

function MOI.get(opt::Optimizer, attr::MOI.VariablePrimal, jx::MOI.VariableIndex)
    MOI.check_result_index_bounds(opt, attr)
    return opt.result.p[opt.jmap[jx]]
end

# ConstraintPrimal: VoV-in-cone → the variable values p; Zeros → Bp − g (≈ 0)
function MOI.get(opt::Optimizer, attr::MOI.ConstraintPrimal, vx::MOI.ConstraintIndex{MOI.VectorOfVariables, <:SupportedCone})
    MOI.check_result_index_bounds(opt, attr)
    return opt.result.p[opt.vidx[opt.vmap[vx]]]
end

function MOI.get(opt::Optimizer{T}, attr::MOI.ConstraintPrimal, ux::MOI.ConstraintIndex{MOI.VectorAffineFunction{T}, MOI.Zeros}) where {T}
    MOI.check_result_index_bounds(opt, attr)
    prob = opt.problem
    # prob.B, prob.g are stored permuted (C cols, R rows); result.p is original.
    # match coordinates: original cols → internal, apply B − g, internal → original
    r = prob.R \ (prob.B * (prob.C * opt.result.p) - prob.g)
    return r[opt.umap[ux]]
end

# ConstraintDual: VoV-in-cone → the conic dual d ∈ K*; Zeros → the equality
# multiplier y. A ConstraintDual lives in the dual cone and its sign is pinned
# by that geometry, so it does not flip with objective sense (unlike the
# recovered ObjectiveValue below).
function MOI.get(opt::Optimizer, attr::MOI.ConstraintDual, vx::MOI.ConstraintIndex{MOI.VectorOfVariables, <:SupportedCone})
    MOI.check_result_index_bounds(opt, attr)
    return opt.result.d[opt.vidx[opt.vmap[vx]]]
end

function MOI.get(opt::Optimizer{T}, attr::MOI.ConstraintDual, ux::MOI.ConstraintIndex{MOI.VectorAffineFunction{T}, MOI.Zeros}) where {T}
    MOI.check_result_index_bounds(opt, attr)
    return opt.result.y[opt.umap[ux]]
end

function MOI.get(opt::Optimizer{T}, attr::MOI.ObjectiveValue) where {T}
    MOI.check_result_index_bounds(opt, attr)
    return _sensesign(opt) * opt.result.pobj + opt.objconst
end

function MOI.get(opt::Optimizer{T}, attr::MOI.DualObjectiveValue) where {T}
    MOI.check_result_index_bounds(opt, attr)
    return _sensesign(opt) * opt.result.dobj + opt.objconst
end
