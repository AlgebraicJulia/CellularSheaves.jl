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

# cones we accept a VectorOfVariables in — one cone block each
const SupportedCone = Union{
    MOI.Nonnegatives,
    MOI.SecondOrderCone,
    MOI.PositiveSemidefiniteConeTriangle,
    MOI.ExponentialCone,
    MOI.PowerCone,
}

# IPM cone for a MOI set (dimension comes from the block width)
function ipmcone(::MOI.Nonnegatives)
    return PositiveCone()
end

function ipmcone(::MOI.SecondOrderCone)
    return SecondOrderCone()
end

function ipmcone(::MOI.PositiveSemidefiniteConeTriangle)
    return SemidefiniteCone()
end

function ipmcone(::MOI.ExponentialCone)
    return ExponentialCone()
end

function ipmcone(set::MOI.PowerCone)
    return PowerCone(set.exponent)
end

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

mutable struct Optimizer{T} <: MOI.AbstractOptimizer
    settings::IPMSettings{T}
    silent::Bool
    # the assembled primal-form problem (built by copy_to, solved by optimize!)
    problem::Union{Nothing, IPMProblem}
    result::Union{Nothing, IPMResult{T}}
    vmap::Dict{MOI.VariableIndex, Int}  # MOI variable → IPM column
    sense::MOI.OptimizationSense
    objconst::T
    solvetime::Float64

    function Optimizer{T}() where {T}
        return new{T}(IPMSettings{T}(), false, nothing, nothing, Dict{MOI.VariableIndex, Int}(), MOI.MIN_SENSE, zero(T), 0.0)
    end
end

function Optimizer()
    return Optimizer{Float64}()
end

function MOI.get(::Optimizer, ::MOI.SolverName)
    return "CellularSheaves IPM"
end

function MOI.supports_incremental_interface(::Optimizer)
    return false
end

function MOI.is_empty(opt::Optimizer)
    return isnothing(opt.result) && isempty(opt.vmap)
end

function MOI.empty!(opt::Optimizer{T}) where {T}
    opt.problem = nothing
    opt.result = nothing
    opt.vmap = Dict{MOI.VariableIndex, Int}()
    opt.sense = MOI.MIN_SENSE
    opt.objconst = zero(T)
    opt.solvetime = 0.0
    return
end

function MOI.supports(::Optimizer, ::MOI.Silent)
    return true
end

function MOI.get(opt::Optimizer, ::MOI.Silent)
    return opt.silent
end

function MOI.set(opt::Optimizer, ::MOI.Silent, v::Bool)
    opt.silent = v
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

#
# lay out the columns: each VoV-in-cone block contiguous, then the free
# variables as a trailing CofreeCone block. Returns (vmap, K, s), where
# vmap maps each MOI variable to its IPM column. (vinv, the column → variable
# ordering, is the build-time intermediate that fixes the layout.)
#
function _layout(src::MOI.ModelLike)
    vars = MOI.get(src, MOI.ListOfVariableIndices())

    K = AbstractCone[]
    s = Int[]
    vinv = MOI.VariableIndex[]
    incone = Set{MOI.VariableIndex}()

    for (F, S) in MOI.get(src, MOI.ListOfConstraintTypesPresent())
        if !(F === MOI.VectorOfVariables && S <: SupportedCone)
            continue
        end

        for ci in MOI.get(src, MOI.ListOfConstraintIndices{F, S}())
            fun = MOI.get(src, MOI.ConstraintFunction(), ci)
            set = MOI.get(src, MOI.ConstraintSet(), ci)
            for v in fun.variables
                v in incone && error("variable in more than one cone block (unsupported)")
                push!(incone, v)
            end
            push!(K, ipmcone(set))
            push!(s, length(fun.variables))
            append!(vinv, fun.variables)
        end
    end

    free = [v for v in vars if v ∉ incone]
    if !isempty(free)
        push!(K, CofreeCone())
        push!(s, length(free))
        append!(vinv, free)
    end

    vmap = Dict(v => c for (c, v) in enumerate(vinv))
    return vmap, K, s
end

#
# equalities Bp = g from VectorAffineFunction-in-Zeros. Returns (B, g).
#
function _equalities(::Type{T}, src::MOI.ModelLike, vmap) where {T}
    n = length(vmap)

    Bi = Int[]
    Bj = Int[]
    Bv = T[]
    g = T[]
    row = 0
    for ci in MOI.get(src, MOI.ListOfConstraintIndices{MOI.VectorAffineFunction{T}, MOI.Zeros}())
        fun = MOI.get(src, MOI.ConstraintFunction(), ci)
        d = MOI.output_dimension(fun)
        for t in fun.terms
            push!(Bi, row + t.output_index)
            push!(Bj, vmap[t.scalar_term.variable])
            push!(Bv, t.scalar_term.coefficient)
        end
        for k in 1:d                 # Ax + b = 0  ⇒  Ax = -b
            push!(g, -fun.constants[k])
        end
        row += d
    end

    return sparse(Bi, Bj, Bv, row, n), g
end

#
# objective ½pᵀQp − fᵀp. Returns (Q, f, sense, objconst).
#
function _objective(::Type{T}, src::MOI.ModelLike, vmap) where {T}
    n = length(vmap)

    sense = MOI.get(src, MOI.ObjectiveSense())
    F = MOI.get(src, MOI.ObjectiveFunctionType())
    obj = MOI.get(src, MOI.ObjectiveFunction{F}())

    Qi = Int[]
    Qj = Int[]
    Qv = T[]
    f = zeros(T, n)

    if sense == MOI.MAX_SENSE   # IPM minimizes
        σ = -one(T)
    else
        σ = one(T)
    end

    if obj isa MOI.ScalarQuadraticFunction
        for qt in obj.quadratic_terms
            a = vmap[qt.variable_1]
            b = vmap[qt.variable_2]
            push!(Qi, a)
            push!(Qj, b)
            push!(Qv, σ * qt.coefficient)
            if a != b
                push!(Qi, b)
                push!(Qj, a)
                push!(Qv, σ * qt.coefficient)
            end
        end
        aff = obj.affine_terms
    else
        aff = obj.terms
    end

    for at in aff
        f[vmap[at.variable]] -= σ * at.coefficient   # min σ·(aᵀp) = −fᵀp
    end

    return sparse(Qi, Qj, Qv, n, n), f, sense, MOI.constant(obj)
end

#
# copy_to: read the model in primal form into an IPMProblem
#
function MOI.copy_to(dest::Optimizer{T}, src::MOI.ModelLike) where {T}
    MOI.empty!(dest)

    vmap, K, s = _layout(src)
    B, g = _equalities(T, src, vmap)
    Q, f, sense, objconst = _objective(T, src, vmap)

    dest.vmap = vmap
    dest.sense = sense
    dest.objconst = objconst
    dest.problem = IPMProblem(Q, B, f, g, K, s)
    return MOI.Utilities.identity_index_map(src)
end

#
# optimize!: solve the assembled problem
#
function MOI.optimize!(opt::Optimizer)
    settings = opt.silent ? IPMSettings(opt.settings; verbose = 0) : opt.settings

    t0 = time()
    opt.result = solve(opt.problem, settings)
    opt.solvetime = time() - t0
    return
end

#
# results
#
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

function _solved_ok(opt::Optimizer)
    return !isnothing(opt.result) && opt.result.status in (OPTIMAL, NEAR_OPTIMAL)
end

function MOI.get(opt::Optimizer, attr::MOI.PrimalStatus)
    (attr.result_index == 1 && _solved_ok(opt)) || return MOI.NO_SOLUTION
    return MOI.FEASIBLE_POINT
end

function MOI.get(opt::Optimizer, attr::MOI.DualStatus)
    (attr.result_index == 1 && _solved_ok(opt)) || return MOI.NO_SOLUTION
    return MOI.FEASIBLE_POINT
end

function MOI.get(opt::Optimizer, attr::MOI.VariablePrimal, vi::MOI.VariableIndex)
    MOI.check_result_index_bounds(opt, attr)
    return opt.result.p[opt.vmap[vi]]
end

function MOI.get(opt::Optimizer{T}, attr::MOI.ObjectiveValue) where {T}
    MOI.check_result_index_bounds(opt, attr)
    σ = opt.sense == MOI.MAX_SENSE ? -one(T) : one(T)
    return σ * opt.result.pobj + opt.objconst
end
