# to run at the robotarium, multi-agent multi-target tracking
# we are going to be able to run this code on the robotarium, but some simulations will have to be done first
# to test the tracking algorithm before deploying it on the robotarium, we will need to simulate the environment and the targets. 
# The robotarium is a rectangular arena where the agents can move around and track the targets.
# the targets are going to be able to move around the arena (in 2D) and the agents will be attach to a line (1D), 
# which have a beginning and an end. and will be able to move along the lines to tracki the targets as close as possible.
# We can use a simple 2D simulation where the targets move in a random walk and the agents try to track them using a simple nearest neighbor algorithm.

using LinearAlgebra
using Random
using Plots
using CellularSheaves

mutable struct Target
    position::Vector{Float64}
    velocity::Vector{Float64}
end

mutable struct MovingAgent
    position::Vector{Float64}
    velocity::Vector{Float64}
end

# Define rectangular arena dimensions
# arena_width = 10.0
# arena_height = 10.0

# Define the number of agents and targets
# num_agents = 1
# num_targets = 1

# scenario 1: one agent moving along an horizontal line, one target moving in a random walk
"""
Vertex 1: horizontal agent position in R² but only moving along a horizontal line (y = constant)
Vertex 2: target position in R²

The edge stalk is R and compares only the x-coordinates.
"""

# Function to initialize a random target
function initialize_random_target(;
    arena_width = 10.0,
    arena_height = 10.0,
    speed = 0.5
    )
    
    position = [
        rand() * arena_width,
        rand() * arena_height
    ]

    direction = 2π * rand()

    velocity = speed .* [
        cos(direction),
        sin(direction)
    ]

    return Target(position, velocity)
    # ex this has a form of
end


function update_random_target!(
    target::Target,
    dt::Float64;
    arena_width = 10.0,
    arena_height = 10.0,
    turning_noise = 0.4
    )
    
    angle_noise = turning_noise * sqrt(dt) * randn()

    rotation = [
        cos(angle_noise) -sin(angle_noise)
        sin(angle_noise)  cos(angle_noise)
    ]

    target.velocity .= rotation * target.velocity

    # Update the position
    target.position .+= dt .* target.velocity

    # rebound at the horizontal boundaries
    if target.position[1] < 0.0
        target.position[1] = -target.position[1]
        target.velocity[1] *= -1
    elseif target.position[1] > arena_width
        target.position[1] = 2 * arena_width - target.position[1]
        target.velocity[1] *= -1
    end

    # rebound at the vertical boundaries
    if target.position[2] < 0.0
        target.position[2] = -target.position[2]
        target.velocity[2] *= -1
    elseif target.position[2] > arena_height
        target.position[2] = 2 * arena_height - target.position[2]
        target.velocity[2] *= -1
    end

    return nothing
end


function build_horizontal_tracking_sheaf()
    # Two vertices, each with a 2-dimensional stalk
    sheaf = EuclideanSheaf{Float64}([2, 2])

    # R² -> R
    # Extract only the x-coordinate
    R_agent = [1.0 0.0]
    R_target = [1.0 0.0]

    add_sheaf_edge!(
        sheaf,
        1,          # agent vertex
        2,          # target vertex
        R_agent,
        R_target
    )

    return sheaf
end

# z = [
#     agent.position[1],
#     agent.position[2],
#     target.position[1],
#     target.position[2]
# ]

function horizontal_sheaf_control(
    sheaf::EuclideanSheaf,
    agent::MovingAgent,
    target::Target;
    gain::Float64 = 1.5
)
    # 0-cochain: agent stalk followed by target stalk
    z = vcat(agent.position, target.position)
    # ex: z = [agent_x, agent_y, target_x, target_y]

    # CellularSheaves.jl coboundary operator
    δ = coboundary_map(sheaf)

    # Sheaf Laplacian:
    # Lz = δ'δz
    Lz = δ' * (δ * z)

    # First two entries correspond to the agent vertex stalk
    agent_laplacian = Lz[1:2]

    # Negative gradient of the sheaf disagreement energy
    control = -gain .* agent_laplacian

    return collect(control)
end

function update_horizontal_agent!(
    agent,
    target,
    tracking_sheaf,
    dt;
    y_track = y_track,
    gain = 1.5,
    max_speed = 1.2,
    arena_width = arena_width
)
    # Control produced by the actual sheaf Laplacian
    unconstrained_control = horizontal_sheaf_control(
        sheaf,
        agent,
        target;
        gain = gain
    )

    # Allowed tangent direction for y = constant
    P_horizontal = [
        1.0 0.0
        0.0 0.0
    ]

    constrained_control =
        P_horizontal * unconstrained_control

    # Limit speed
    speed = norm(constrained_control)

    if speed > max_speed
        constrained_control .*= max_speed / speed
    end

    agent.velocity .= constrained_control
    agent.position .+= dt .* agent.velocity

    # Arena and track constraints
    agent.position[1] = clamp(
        agent.position[1],
        0.0,
        arena_width
    )

    agent.position[2] = y_track
    agent.velocity[2] = 0.0

    return nothing
end

# Simulation parameters
arena_width = 10.0
arena_height = 10.0
y_track = 3.0

dt = 0.1
num_steps = 400

target = initialize_random_target(
    arena_width = arena_width,
    arena_height = arena_height,
    speed = 0.8
)

agent = MovingAgent(
    [0.0, y_track],
    [0.0, 0.0]
)

tracking_sheaf = build_horizontal_tracking_sheaf()

target_x_history = Float64[]
target_y_history = Float64[]

sheaf_residual_history = Float64[]

animation = @animate for _ in 1:num_steps

    update_random_target!(
        target,
        dt;
        arena_width = arena_width,
        arena_height = arena_height,
        turning_noise = 0.6
    )

    update_horizontal_agent!(
        agent,
        target,
        tracking_sheaf,
        dt;
        y_track = y_track,
        gain = 1.5,
        max_speed = 1.2,
        arena_width = arena_width
    )

    push!(target_x_history, target.position[1])
    push!(target_y_history, target.position[2])

    # Actual sheaf inconsistency δz
    z = vcat(agent.position, target.position)
    δ = coboundary_map(tracking_sheaf)
    residual = norm(δ * z)

    push!(sheaf_residual_history, residual)

    plot(
        xlims = (0, arena_width),
        ylims = (0, arena_height),
        aspect_ratio = :equal,
        xlabel = "x",
        ylabel = "y",
        title = "Horizontal tracking using CellularSheaves.jl",
        legend = :topright
    )

    # Rectangle boundary
    plot!(
        [0.0, arena_width, arena_width, 0.0, 0.0],
        [0.0, 0.0, arena_height, arena_height, 0.0];
        linewidth = 2,
        label = "Arena"
    )

    # Horizontal track
    plot!(
        [0.0, arena_width],
        [y_track, y_track];
        linestyle = :dash,
        linewidth = 2,
        label = "Agent track: y = 3"
    )

    # Target trajectory
    plot!(
        target_x_history,
        target_y_history;
        linewidth = 1.5,
        label = "Target trajectory"
    )

    # Geometric target projection onto the track
    projected_target = [
        target.position[1],
        y_track
    ]

    plot!(
        [target.position[1], projected_target[1]],
        [target.position[2], projected_target[2]];
        linestyle = :dot,
        label = false
    )

    scatter!(
        [target.position[1]],
        [target.position[2]];
        markersize = 8,
        label = "Target"
    )

    scatter!(
        [agent.position[1]],
        [agent.position[2]];
        markersize = 8,
        label = "Agent"
    )
end

gif(
    animation,
    "horizontal_sheaf_tracking.gif";
    fps = 30
)

δ = Matrix(coboundary_map(tracking_sheaf))
L = Matrix(sheaf_laplacian_matrix(tracking_sheaf))

println("Coboundary map:")
display(δ)

println("Sheaf Laplacian:")
display(L)
