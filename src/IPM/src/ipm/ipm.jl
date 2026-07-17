@enum IPMStatus CONTINUE OPTIMAL NEAR_OPTIMAL STALLED NUMERICAL_FAILURE ITERATION_LIMIT PRIMAL_INFEASIBLE DUAL_INFEASIBLE ILL_POSED
@enum RefStatus REACHED_FORCE REACHED_FLOOR REFINE_STALLED

include("history/history.jl")
include("settings/settings.jl")
include("problem.jl")
include("result.jl")
include("workspace/workspace.jl")
include("solver/solver.jl")
