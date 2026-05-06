"""
    HerdingPlatoon

Module for vehicle platoon models with linear herding behavior emerging from
sheaf Laplacian structure on ladder graphs.

This module provides functions to construct ControlledTrajectorySheaf instances
that model vehicle platoons where the herding behavior (vehicles attracting
to each other) emerges from the sheaf structure rather than being explicitly
encoded in the dynamics matrices.

# Exported functions
- `controlled_vehicle_platoon_with_herding`
"""
module HerdingPlatoon

using ..EuclideanSheaves: EuclideanSheaf, add_sheaf_edge!
using ..SheafInterface: vertex_stalks
using ..TrajectorySheaves: ControlledTrajectorySheaf, continuous_to_discrete_zoh
using LinearAlgebra
using SparseArrays

export controlled_vehicle_platoon_with_herding

"""
    controlled_vehicle_platoon_with_herding(
        Ac::AbstractMatrix{T},
        Bc::AbstractMatrix{T},
        h::Real,
        k::Int,
        herding_gain::Real;
        k_herding::Real = herding_gain
    ) -> ControlledTrajectorySheaf{T}

Construct a ControlledTrajectorySheaf for a vehicle platoon where herding behavior
emerges from the sheaf Laplacian rather than being explicitly encoded in the dynamics.

The resulting sheaf is defined on a ladder graph with:
- Two parallel paths (backbones) representing each vehicle's state trajectory over time
- Each backbone has k+2 vertices following the ControlledTrajectorySheaf pattern
- Rungs connect corresponding time steps between vehicles
- Herding gain controls the strength of attraction between vehicles

# Arguments
- `Ac`: n×n continuous-time state matrix for each vehicle's individual dynamics
- `Bc`: n×m continuous-time input matrix for each vehicle
- `h`: sample period for discretization
- `k`: number of time steps in the trajectory
- `herding_gain`: positive scalar controlling herding strength (default: 1.0)
- `k_herding`: alternative parameter name for herding_gain (for backward compatibility)

# Returns
- `ControlledTrajectorySheaf{T}` representing the herding vehicle platoon system

# Notes
The herding behavior emerges from the sheaf structure:
- When herding_gain = 0, vehicles evolve independently (decoupled)
- When herding_gain > 0, vehicles are attracted to each other with force
  proportional to their separation
- The sheaf Laplacian automatically generates the appropriate coupling terms
"""
function controlled_vehicle_platoon_with_herding(
    Ac::AbstractMatrix{T},
    Bc::AbstractMatrix{T},
    h::Real,
    k::Int,
    herding_gain::Real;
    k_herding::Real = herding_gain
) where T
    # Use herding_gain (with fallback to k_herding for backward compatibility)
    γ = herding_gain
    
    # Discretize the continuous-time system
    Ad, Bd = continuous_to_discrete_zoh(Ac, Bc, h)
    n = size(Ac, 1)  # state dimension
    m = size(Bc, 2)  # control dimension
    
    # Construct ladder graph sheaf following ControlledTrajectorySheaf pattern for each vehicle
    # We'll create a sheaf with 2*(k+2) vertices:
    #   Vertices 1 to k+2: vehicle 1 trajectory (following CTL pattern)
    #   Vertices k+3 to 2(k+2): vehicle 2 trajectory (following CTL pattern)
    
    # Stalks following CTL pattern for each vehicle:
    #   First vertex: x₁ (n-dimensional, dummy initial)
    #   Middle k vertices: (xₜ, uₜ) ((n+m)-dimensional)
    #   Last vertex: x_{k+1} (n-dimensional, terminal)
    
    stalks = vcat(
        [n],                           # v1: x₁⁽¹⁾ (dummy initial)
        fill(n+m, k),                  # v2 to v_{k+1}: (xₜ⁽¹⁾, uₜ⁽¹⁾)
        [n],                           # v_{k+2}: x_{k+1}⁽¹⁾ (terminal)
        [n],                           # v_{k+3}: x₁⁽²⁾ (dummy initial)
        fill(n+m, k),                  # v_{k+4} to v_{2k+3}: (xₜ⁽²⁾, uₜ⁽²⁾)
        [n]                            # v_{2k+4}: x_{k+1}⁽²⁾ (terminal)
    )
    
    sheaf = EuclideanSheaf{T}(stalks)
    
    # Backbone edges for vehicle 1 (following CTL pattern)
    # Edge (1,2): connects x₁⁽¹⁾ dummy (state-only) to (x₁⁽¹⁾, u₁⁽¹⁾) (state-control)
    add_sheaf_edge!(sheaf, 1, 2, 
        Matrix{T}(I, n, n),                           # ρ_{1→e} = I_n
        hcat(Matrix{T}(I, n, n), zeros(T, n, m))      # ρ_{2→e} = [I_n | 0]
    )
    
    # Dynamics edges for vehicle 1: (t+1, t+2) for t=1..k-1
    for t in 1:k-1
        add_sheaf_edge!(sheaf, 1+t, 1+t+1, 
            hcat(Ad, Bd),                           # ρ_{t+1→e} = [A_d | B_d]
            hcat(Matrix{T}(I, n, n), zeros(T, n, m)) # ρ_{t+2→e} = [I | 0]
        )
    end
    
    # Final dynamics edge for vehicle 1: (k+1, k+2)
    add_sheaf_edge!(sheaf, k+1, k+2, 
        hcat(Ad, Bd),                           # ρ_{k+1→e} = [A_d | B_d]
        Matrix{T}(I, n, n)                      # ρ_{k+2→e} = I_n
    )
    
    # Backbone edges for vehicle 2 (shifted by k+2)
    # Edge (k+2, k+3): connects x₁⁽²⁾ dummy (state-only) to (x₁⁽²⁾, u₁⁽²⁾) (state-control)
    add_sheaf_edge!(sheaf, k+2, k+3, 
        Matrix{T}(I, n, n),                           # ρ_{k+2→e} = I_n
        hcat(Matrix{T}(I, n, n), zeros(T, n, m))      # ρ_{k+3→e} = [I_n | 0]
    )
    
    # Dynamics edges for vehicle 2: (t+1, t+2) for t=1..k-1 (shifted by k+2)
    for t in 1:k-1
        add_sheaf_edge!(sheaf, (k+2)+t, (k+2)+t+1, 
            hcat(Ad, Bd),                           # ρ→e = [A_d | B_d]
            hcat(Matrix{T}(I, n, n), zeros(T, n, m)) # ρ→e = [I | 0]
        )
    end
    
    # Final dynamics edge for vehicle 2: (2(k+2)-1, 2(k+2))
    add_sheaf_edge!(sheaf, 2(k+2)-1, 2(k+2), 
        hcat(Ad, Bd),                           # ρ→e = [A_d | B_d]
        Matrix{T}(I, n, n)                      # ρ→e = I_n
    )
    
    # Rung edges (herding coupling) at each time step
    # For each t=1..k+1: connect corresponding state vertices between vehicles
    for t in 1:k+1
        if t == 1
            # Initial state: x₁ for both vehicles (vertices 1 and k+3)
            v1 = 1
            v2 = k + 3
            ρ1 = Matrix{T}(I, n, n)                           # I_n
            ρ2 = -γ * Matrix{T}(I, n, n)                      # -γ*I_n
        elseif t == k+1
            # Final state: x_{k+1} for both vehicles (vertices k+2 and 2(k+2))
            v1 = k + 2
            v2 = 2*(k+2)
            ρ1 = Matrix{T}(I, n, n)                           # I_n
            ρ2 = -γ * Matrix{T}(I, n, n)                      # -γ*I_n
        else
            # Intermediate state: xₜ for both vehicles (1 < t < k+1)
            # Vehicle 1: vertex t+1 ((xₜ⁽¹⁾, uₜ⁽¹⁾))
            # Vehicle 2: vertex (k+2)+t ((xₜ⁽²⁾, uₜ⁽²⁾))
            v1 = t + 1
            v2 = (k+2) + t
            # We want to extract only the state portion (first n components) from the state-control pairs
            ρ1 = hcat(Matrix{T}(I, n, n), zeros(T, n, m))   # [I_n | 0_{n×m}] - selects state
            ρ2 = -γ * hcat(Matrix{T}(I, n, n), zeros(T, n, m)) # -γ*[I_n | 0_{n×m}] - selects state
        end
        
        # Add the rung edge for herding coupling at time step t
        add_sheaf_edge!(sheaf, v1, v2, ρ1, ρ2)
    end
    
    return ControlledTrajectorySheaf{T}(k, T(h), Ac, Bc, Ad, Bd, sheaf, n, m)
end

end # HerdingPlatoon