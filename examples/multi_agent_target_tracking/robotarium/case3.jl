# Single-agent, two-target tracking simulation for Robotarium.
#
# The targets move in a two-dimensional rectangular arena using a
# correlated random walk. The agent is constrained to move along a
# horizontal line and uses a sheaf-Laplacian control law to track the
# target's x-coordinate.

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

# target dynamics
function initialize_random_target(;
    arena_width::Float64,
    arena_height::Float64,
    speed::Float64
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
end


function update_random_target!(
    target::Target,
    dt::Float64;
    arena_width::Float64,
    arena_height::Float64,
    turning_noise::Float64
    )

    # Randomly perturb the target's direction.
    angle_noise = turning_noise * sqrt(dt) * randn()

    rotation = [
        cos(angle_noise) -sin(angle_noise)
        sin(angle_noise)  cos(angle_noise)
    ]

    target.velocity .= rotation * target.velocity

    # Update target position.
    target.position .+= dt .* target.velocity

    # Reflect from the left and right boundaries.
    if target.position[1] < 0.0
        target.position[1] = -target.position[1]
        target.velocity[1] *= -1
    elseif target.position[1] > arena_width
        target.position[1] = 2 * arena_width - target.position[1]
        target.velocity[1] *= -1
    end

    # Reflect from the lower and upper boundaries.
    if target.position[2] < 0.0
        target.position[2] = -target.position[2]
        target.velocity[2] *= -1
    elseif target.position[2] > arena_height
        target.position[2] = 2 * arena_height - target.position[2]
        target.velocity[2] *= -1
    end

    return nothing
end


"""
Construct a sheaf for one horizontal agent and two targets.

Vertex 1: Agent position in ℝ².
Vertex 2: Target 1 position in ℝ².
Vertex 3: Target 2 position in ℝ².

Each edge compares the agent's x-coordinate with one target's
x-coordinate.
"""
function build_horizontal_tracking_sheaf()
    # Three vertices, each with a two-dimensional stalk.
    sheaf = EuclideanSheaf{Float64}([2, 2, 2])

    # Restriction map ℝ² → ℝ that extracts the x-coordinate.
    R_x = [1.0 0.0]

    # Agent and target 1.
    add_sheaf_edge!(
        sheaf,
        1,
        2,
        R_x,
        R_x
    )

    # Agent and target 2.
    add_sheaf_edge!(
        sheaf,
        1,
        3,
        R_x,
        R_x
    )

    return sheaf
end


function horizontal_sheaf_control(
    δ::AbstractMatrix{Float64},
    agent::MovingAgent,
    target_1::Target,
    target_2::Target;
    gain::Float64
  )
    z = vcat(agent.position, target_1.position, target_2.position)
    Lz = δ' * (δ * z)

    # only control over the agent's position
    agent_laplacian = Lz[1:2]

    # Negative gradient of the sheaf disagreement energy.
    control = -gain .* agent_laplacian

    return collect(control)
end


# agent dynamics
function update_horizontal_agent!(
    agent::MovingAgent,
    target_1::Target,
    target_2::Target,
    δ::AbstractMatrix{Float64},
    dt::Float64;
    gain::Float64,
    max_speed::Float64,
    arena_width::Float64
  )
  constrained_control = horizontal_sheaf_control(
        δ,
        agent,
        target_1,
        target_2;
        gain = gain
    )

    # Limit the magnitude of the agent's velocity.
    speed = norm(constrained_control)

    if speed > max_speed
        constrained_control .*= max_speed / speed
    end

    agent.velocity .= constrained_control
    agent.position .+= dt .* agent.velocity

    # Enforce the left and right arena boundaries.
    if agent.position[1] <= 0.0
        agent.position[1] = 0.0

        if agent.velocity[1] < 0.0
            agent.velocity[1] = 0.0
        end
    elseif agent.position[1] >= arena_width
        agent.position[1] = arena_width

        if agent.velocity[1] > 0.0
            agent.velocity[1] = 0.0
        end
    end

    return nothing
end

# Simulation parameters
arena_width = 10.0
arena_height = 10.0
y_track = 3.0

dt = 0.1
num_steps = 400

target_speed = 0.8
turning_noise = 0.6

control_gain = 1.5
agent_max_speed = 1.2

target_1 = initialize_random_target(
    arena_width = arena_width,
    arena_height = arena_height,
    speed = target_speed
)

target_2 = initialize_random_target(
    arena_width = arena_width,
    arena_height = arena_height,
    speed = target_speed
)


agent = MovingAgent(
    [0.0, y_track],
    [0.0, 0.0]
)

tracking_sheaf = build_horizontal_tracking_sheaf()

# The sheaf does not change during the simulation
δ = coboundary_map(tracking_sheaf)
L = sheaf_laplacian_matrix(tracking_sheaf)

target_1_x_history = Float64[]
target_1_y_history = Float64[]

target_2_x_history = Float64[]
target_2_y_history = Float64[]

residual_history = Float64[]

animation = @animate for _ in 1:num_steps
    update_random_target!(
        target_1,
        dt;
        arena_width = arena_width,
        arena_height = arena_height,
        turning_noise = turning_noise
    )

    update_random_target!(
        target_2,
        dt;
        arena_width = arena_width,
        arena_height = arena_height,
        turning_noise = turning_noise
    )

    update_horizontal_agent!(
        agent,
        target_1,
        target_2,
        δ,
        dt;
        gain = control_gain,
        max_speed = agent_max_speed,
        arena_width = arena_width
    )

    push!(target_1_x_history, target_1.position[1])
    push!(target_1_y_history, target_1.position[2])

    push!(target_2_x_history, target_2.position[1])
    push!(target_2_y_history, target_2.position[2])

    z = vcat(
    agent.position,
    target_1.position,
    target_2.position
    )

    residual = norm(δ * z)
    push!(residual_history, residual)

    plot(
        xlims = (0.0, arena_width),
        ylims = (0.0, arena_height),
        aspect_ratio = :equal,
        xlabel = "x",
        ylabel = "y",
        title = "Target Tracking with 2 Targets and 1 Agent",
        legend = :outertopright    
        )

    # Arena boundary
    plot!(
        [0.0, arena_width, arena_width, 0.0, 0.0],
        [0.0, 0.0, arena_height, arena_height, 0.0];
        linewidth = 2,
        label = false
    )

    # Horizontal agent track
    plot!(
        [0.0, arena_width],
        [y_track, y_track];
        linestyle = :dash,
        linewidth = 2,
        label = false
    )

    # Target trajectory
    plot!(
        target_1_x_history,
        target_1_y_history;
        linewidth = 1.5,
        color = "#D4A017",
        label = "Target 1 trajectory"
    )

    plot!(
        target_2_x_history,
        target_2_y_history;
        linewidth = 1.5,
        color = "#3A6EA5",
        label = "Target 2 trajectory"
    )

    # Vertical projection of the target onto the agent's track
    plot!(
        [target_1.position[1], target_1.position[1]],
        [target_1.position[2], y_track];
        linestyle = :dot,
        color = "#D4A017",
        label = false
    )

    plot!(
        [target_2.position[1], target_2.position[1]],
        [target_2.position[2], y_track];
        linestyle = :dot,
        color = "#3A6EA5",
        label = false
    )

    scatter!(
        [target_1.position[1]],
        [target_1.position[2]];
        markersize = 8,
        color = "#D4A017",
        label = "Target 1"
    )

    scatter!(
        [target_2.position[1]],
        [target_2.position[2]];
        markersize = 8,
        color = "#3A6EA5",
        label = "Target 2"
    )

    scatter!(
        [agent.position[1]],
        [agent.position[2]];
        markersize = 8,
        color = "#FF5733",
        label = "Agent"
    )
end

mp4(
    animation,
    "case3.mp4";
    fps = 30
)

println("Coboundary map:")
display(Matrix(δ))

println("Sheaf Laplacian:")
display(Matrix(L))
