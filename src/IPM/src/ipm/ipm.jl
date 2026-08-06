@enum IPMStatus CONTINUE OPTIMAL NEAR_OPTIMAL STALLED NUMERICAL_FAILURE ITERATION_LIMIT PRIMAL_INFEASIBLE DUAL_INFEASIBLE ILL_POSED NEAR_PRIMAL_INFEASIBLE NEAR_DUAL_INFEASIBLE NEAR_ILL_POSED
@enum RefStatus REACHED REFINE_STALLED REFINE_ITMAX

include("history/history.jl")
include("settings/settings.jl")
include("problem.jl")
include("result/result.jl")
include("workspace/workspace.jl")
include("solver/solver.jl")
