module ControlSheaves

export AbstractDynamics,
       NoDynamics,
       TrophicDynamics,
       CurrentAgentDynamics,
       CurrentTargetDynamics,
       AbstractEntity,
       Entity,
       Agent,
       Target,
       run_dynamics,
       run_dynamics!,
       update_dynamics!,
       compute_control_output!,
       build_ode_problem,
       build_integrator,
       integrate_step,
       integrate_step!

abstract type AbstractDynamics end
struct NoDynamics <: AbstractDynamics end
struct TrophicDynamics <: AbstractDynamics end
struct CurrentAgentDynamics <: AbstractDynamics end
struct CurrentTargetDynamics <: AbstractDynamics end

include("integrate.jl")
include("entity.jl")
include("MultiAgentTracking.jl")

import .MultiAgentTracking
export MultiAgentTracking

end # module ControlSheaves
