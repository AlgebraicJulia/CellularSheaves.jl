# # Multilayer Hierarchical Control for Team-of-Teams Escort Formation

# This example implements a three-layer control architecture using **sheaf pushforwards** to coordinate a "team of teams." 

# ## Scenario
# We have 15 quadrotor agents organized into three rings:
# - **Ring 1 (6 agents):** Escorting Target 1.
# - **Ring 2 (6 agents):** Escorting Target 2.
# - **Ring 3 (3 agents):** Providing communication support, tracking the midpoint of the two escort rings.

# ## Hierarchical Architecture
# 1. **High-Level Planner (Pushforward Sheaf $f_\star F$):** Operates on a coarse graph $H$ (3 nodes: T1, T2, Support). It computes the optimal "center" positions for the three rings.
# 2. **Mid-Level Planner (Sheaf $F$):** Operates on the full agent graph $G$. It takes the ring centers from the high-level and computes the specific agent offsets via harmonic extension.
# 3. **Low-Level Controller (LQR + Feedforward):** Executes the 3-stage (Pos, Vel, Accel) tracking for each agent.

using CellularSheaves
using CellularSheaves.Formations
using CellularSheaves.AgentControllers
using CellularSheaves.DistributedLayeredControl
using CellularSheaves.Pushforwards
using CellularSheaves.GraphHomomorphisms
using LinearAlgebra
using Statistics
using Distributed
using Plots
using Printf
using Graphs

# --- Setup ---
default(framestyle = :box, grid = true, gridalpha = 0.18, gridstyle = :dot,
    titlefontsize = 10, guidefontsize = 9, legendfontsize = 8, tickfontsize = 8,
    markerstrokewidth = 0, size = (800, 400))

const RING_COLOR = :steelblue
const TARGET_COLOR = :black

dyn = QuadrotorDynamics()
DT = 0.05
nx = 10
epsilon = 0.02

# LQR Costs
Q_diag_std = [500.0, 500.0, 500.0, 150.0, 150.0, 100.0, 100.0, 100.0, 5.0, 5.0]
Q_std = Matrix(Diagonal(Q_diag_std))
R_lqr = Matrix(Diagonal([0.005, 0.005, 0.005]))

# --- 1. Define the Mid-Level Sheaf F over Graph G ---
# Agents: 1-6 (Ring 1), 7-12 (Ring 2), 13-15 (Ring 3)
# Targets: 16 (T1), 17 (T2), 18 (S)
const NA1 = 6
const NA2 = 6
const NA3 = 3
const NT = 3
const T1_NODE = 16
const T2_NODE = 17
const TS_NODE = 18
const r_ring = 0.3

function build_multilayer_sheaf(na1, na2, na3, r1, r2, r3)
    ## Total nodes: (na1+1) + (na2+1) + (na3+1) = 15 + 3 = 18
    ## Ring 1: 1..na1 agents, node na1+1 is T1
    ## Ring 2: na1+2..na1+na2+1 agents, node na1+na2+2 is T2
    ## Ring 3: na1+na2+3..na1+na2+na3+2 agents, node na1+na2+na3+3 is S
    F = EuclideanSheaf{Float64}(fill(4, 18))
    
    function add_ring!(sheaf, agent_range, target_node, radius)
        n = length(agent_range)
        for i in 1:n
            u = agent_range[i]
            v = agent_range[i % n + 1]
            angle_u = (i - 1) * 2π / n
            angle_v = (i % n) * 2π / n
            du = [cos(angle_u), sin(angle_u), 0.0] * radius
            dv = [cos(angle_v), sin(angle_v), 0.0] * radius
            add_sheaf_edge!(sheaf, u, v, se3_translation_matrix(du), se3_translation_matrix(dv))
        end
        u_first = agent_range[1]
        d_first = [cos(0), sin(0), 0.0] * radius
        add_sheaf_edge!(sheaf, u_first, target_node, se3_translation_matrix(d_first), Matrix{Float64}(I, 4, 4))
    end

    add_ring!(F, 1:na1, na1 + 1, r1)
    add_ring!(F, na1 + 2 : na1 + na2 + 1, na1 + na2 + 2, r2)
    add_ring!(F, na1 + na2 + 3 : na1 + na2 + na3 + 2, na1 + na2 + na3 + 3, r3)
    
    # Add coordination edges between observer agents to create cross-edges in the pushforward
    # Observer 1 (Node 1) <-> Observer 3 (Node na1 + na2 + 3)
    # Observer 2 (Node na1 + 2) <-> Observer 3 (Node na1 + na2 + 3)
    obs1 = 1
    obs2 = na1 + 2
    obs3 = na1 + na2 + 3
    
    add_sheaf_edge!(F, obs1, obs3, Matrix{Float64}(I, 4, 4), Matrix{Float64}(I, 4, 4))
    add_sheaf_edge!(F, obs2, obs3, Matrix{Float64}(I, 4, 4), Matrix{Float64}(I, 4, 4))
    
    return F
end

F = build_multilayer_sheaf(NA1, NA2, NA3, r_ring, r_ring * 1.2, r_ring * 0.8)

# --- 2. Define the High-Level Graph H and Homomorphism f ---
# H has 3 nodes: T1, T2, S
# f maps agents and targets to their respective ring centers
function build_homomorphism(na1, na2, na3)
    ## vertex_map: 
    ## Ring 1 (agents + T1) -> 1
    ## Ring 2 (agents + T2) -> 2
    ## Ring 3 (agents + S)   -> 3
    v_map = vcat(fill(1, na1 + 1), fill(2, na2 + 1), fill(3, na3 + 1))
    
    # The GraphHomomorphism constructor takes the vertex map and the number of target nodes.
    # We don't need to pass the SimpleGraph object itself.
    return GraphHomomorphism(v_map, 3)
end

f = build_homomorphism(NA1, NA2, NA3)

# --- 3. Construct and Validate the Pushforward Sheaf ---
# The pushforward sheaf f_*F represents the "coarse-grained" coordination 
# of the three rings.
PfF = pushforward_sheaf(f, F)

println("Pushforward Sheaf Stalk Dimensions: ", vertex_stalks(PfF))
# We expect the stalk dimension at each node in H to be the dimension of 
# the global sections of the fiber. For an escort ring, this is 4 (SE(3) translation).

# Inspect restriction maps in PfF to verify the "midpoint" geometry
# The edges in H are (1,3) and (2,3).
for e in edges(underlying_graph(PfF))
    u, v = src(e), dst(e)
    rm_u = get_restriction_map(PfF, u, v)
    rm_v = get_restriction_map(PfF, v, u)
    println("Edge ($u, $v) restriction map size: $(size(rm_u))")
end

# Validation: Check that the global sections of the pushforward sheaf 
# represent the 3 ring centers. We expect the global section space to be 
# spanned by 3 vectors (each 4D) that are collinear in the sense that 
# they represent the same translation in the world frame.
B_pf = nullspace_ldlt(coboundary_map(PfF)' * coboundary_map(PfF))
println("Pushforward global section basis size: $(size(B_pf))")

# The columns of B_pf are the basis for the global sections of PfF.
# We verify that the dimensions match our expectation (3 rings * 4D = 12D total space, 
# but the global sections should be triples of 4D vectors that are colinear.).
@assert size(B_pf, 2) == 4 "Expected 4-dimensional global section space for SE(3) translation"

