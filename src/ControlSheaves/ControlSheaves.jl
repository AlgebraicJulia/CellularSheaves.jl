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
include("Tikhonov.jl")
include("entity.jl")
include("MultiAgentTracking.jl")
include("TrackingDSL/TrackingDSL.jl")

import .Tikhonov
export Tikhonov

import .MultiAgentTracking
export MultiAgentTracking

import .TrackingDSL
export TrackingDSL

include("AgentControllers.jl")
import .AgentControllers
export AgentControllers

include("DistributedLayeredControl.jl")
import .DistributedLayeredControl
export DistributedLayeredControl

include("CoordinationBenchmarks.jl")
import .CoordinationBenchmarks
export CoordinationBenchmarks

include("CoordinationProfiling.jl")
import .CoordinationProfiling
export CoordinationProfiling

end # module ControlSheaves
