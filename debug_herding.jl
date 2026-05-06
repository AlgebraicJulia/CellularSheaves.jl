using CellularSheaves
using LinearAlgebra

# Simple test: double integrator system
Ac = [0.0 1.0; 0.0 0.0]
Bc = [0.0 1.0; 0.0 0.0]
h = 0.25
k = 2  # Small for easier debugging
herding_gain = 0.5

println("Input: Ac size ", size(Ac), ", Bc size ", size(Bc))
println("h = ", h, ", k = ", k, ", herding_gain = ", herding_gain)

n = size(Ac, 1)
m = size(Bc, 2)
println("n = ", n, ", m = ", m)

# Let's manually trace through our function construction
stalks = vcat(
    [n],                           # v1: x₁⁽¹⁾ (dummy initial)
    fill(n+m, k),                  # v2 to v_{k+1}: (xₜ⁽¹⁾, uₜ⁽¹⁾)
    [n],                           # v_{k+2}: x_{k+1}⁽¹⁾ (terminal)
    [n],                           # v_{k+3}: x₁⁽²⁾ (dummy initial)
    fill(n+m, k),                  # v_{k+4} to v_{2k+3}: (xₜ⁽²⁾, uₜ⁽²⁾)
    [n]                            # v_{2k+4}: x_{k+1}⁽²⁾ (terminal)
)

println("Constructed stalks: ", stalks)

sheaf = EuclideanSheaf{Float64}(stalks)
println("Sheaf vertex stalks: ", vertex_stalks(sheaf))

# Helper functions
v1_state(t) = t
v1_ctrl(t) = k+2 + t
v2_state(t) = (k+2) + t
v2_ctrl(t) = 2*(k+2) + t

println()
println("Vehicle info:")
println("  k+2 = ", k+2)
println("  2*(k+2) = ", 2*(k+2))

for t in 1:k+1
    if t == 1
        v1 = 1
        v2 = k + 3
        println("t=$t (initial): v1=$v1, v2=$v2")
    elseif t == k+1
        v1 = k + 2
        v2 = 2*(k+2)
        println("t=$t (final): v1=$v1, v2=$v2")
    else
        v1 = t + 1
        v2 = (k+2) + t
        println("t=$t (intermediate): v1=$v1, v2=$v2")
    end
    
    println("  vertex_stalks[v1] = ", get_vertex_stalk(sheaf, v1))
    println("  vertex_stalks[v2] = ", get_vertex_stalk(sheaf, v2))
end