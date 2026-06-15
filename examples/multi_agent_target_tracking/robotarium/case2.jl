# two-agent, single-target tracking simulation for eventual Robotarium use.
#
# The target moves in a two-dimensional rectangular arena using a
# correlated random walk. The agent is constrained to move along a
# horizontal line and uses a sheaf-Laplacian control law to track the
# target's x-coordinate. the other agent is constrained to move along a vertical line and uses a
# sheaf-Laplacian control law to track the target's y-coordinate.


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
Construct a sheaf with two vertices.
Vertex 1: Agent position in ℝ².
Vertex 2: Target position in ℝ².

The edge stalk is ℝ. Both restriction maps extract only the x-coordinate, so the sheaf measures 
horizontal disagreement.
"""
function build_horizontal_tracking_sheaf()
    # Two vertices, each with a two-dimensional stalk.
    sheaf = EuclideanSheaf{Float64}([2, 2])

    # Restriction maps ℝ² → ℝ.
    R_agent = [1.0 0.0]
    R_target = [1.0 0.0]

    add_sheaf_edge!(
        sheaf,
        1,
        2,
        R_agent,
        R_target
    )

    return sheaf
end


function horizontal_sheaf_control(
    δ_horizontal::Matrix{Float64},
    agent::MovingAgent,
    target::Target;
    gain::Float64
)
    z = vcat(agent.position, target.position)
    Lz = δ_horizontal' * (δ_horizontal * z)

    agent_laplacian = Lz[1:2]

    # Negative gradient of the sheaf disagreement energy.
    control = -gain .* agent_laplacian

    return collect(control)
end


# agent dynamics
function update_horizontal_agent!(
    agent::MovingAgent,
    target::Target,
    δ_horizontal::Matrix{Float64},
    dt::Float64;
    gain::Float64,
    max_speed::Float64,
    arena_width::Float64
)
    constrained_control = horizontal_sheaf_control(
        δ_horizontal,
        agent,
        target;
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

# vertical sheaf 
function build_vertical_tracking_sheaf()
    # Two vertices, each with a two-dimensional stalk.
    sheaf = EuclideanSheaf{Float64}([2, 2])

    # Restriction maps ℝ² → ℝ.
    R_agent = [0.0 1.0]
    R_target = [0.0 1.0]

    add_sheaf_edge!(
        sheaf,
        1,
        2,
        R_agent,
        R_target
    )

    return sheaf
end

function vertical_sheaf_control(
    δ_vertical::Matrix{Float64},
    agent::MovingAgent,
    target::Target;
    gain::Float64
)
    z = vcat(agent.position, target.position)
    Lz = δ_vertical' * (δ_vertical * z)

    agent_laplacian = Lz[1:2]

    # Negative gradient of the sheaf disagreement energy.
    control = -gain .* agent_laplacian

    return collect(control)
end

function update_vertical_agent!(
    agent::MovingAgent,
    target::Target,
    δ_vertical::Matrix{Float64},
    dt::Float64;
    gain::Float64,
    max_speed::Float64,
    arena_height::Float64
  )
    constrained_control = vertical_sheaf_control(
        δ_vertical,
        agent,
        target;
        gain = gain
    )

    # Limit the magnitude of the agent's velocity.
    speed = norm(constrained_control)

    if speed > max_speed
        constrained_control .*= max_speed / speed
    end

    agent.velocity .= constrained_control
    agent.position .+= dt .* agent.velocity

    # Enforce the lower and upper arena boundaries.
    if agent.position[2] <= 0.0
        agent.position[2] = 0.0

        if agent.velocity[2] < 0.0
            agent.velocity[2] = 0.0
        end
    elseif agent.position[2] >= arena_height
        agent.position[2] = arena_height

        if agent.velocity[2] > 0.0
            agent.velocity[2] = 0.0
        end
    end

    return nothing
end


# Simulation parameters
arena_width = 10.0
arena_height = 10.0
y_track = 3.0
x_track = 1.0

dt = 0.1
num_steps = 400

target_speed = 0.8
turning_noise = 0.6

control_gain = 1.5
agent_max_speed = 1.2

target = initialize_random_target(
    arena_width = arena_width,
    arena_height = arena_height,
    speed = target_speed
)

agent_horizontal = MovingAgent(
    [0.0, y_track],
    [0.0, 0.0]
)

agent_vertical = MovingAgent(
    [x_track, 0.0],
    [0.0, 0.0]
)

tracking_sheaf_horizontal = build_horizontal_tracking_sheaf()
tracking_sheaf_vertical = build_vertical_tracking_sheaf()

# The sheaf does not change during the simulation
δ_horizontal = coboundary_map(tracking_sheaf_horizontal)
L_horizontal = sheaf_laplacian_matrix(tracking_sheaf_horizontal)

δ_vertical = coboundary_map(tracking_sheaf_vertical)
L_vertical = sheaf_laplacian_matrix(tracking_sheaf_vertical)

target_x_history = Float64[]
target_y_history = Float64[]

animation = @animate for _ in 1:num_steps
    update_random_target!(
        target,
        dt;
        arena_width = arena_width,
        arena_height = arena_height,
        turning_noise = turning_noise
    )

    update_horizontal_agent!(
        agent_horizontal,
        target,
        tracking_sheaf_horizontal,
        dt;
        y_track = y_track,
        gain = control_gain,
        max_speed = agent_max_speed,
        arena_width = arena_width
    )

    update_vertical_agent!(
        agent_vertical,
        target,
        tracking_sheaf_vertical,
        dt;
        x_track = x_track, # track the target's x-coordinate
        gain = control_gain,
        max_speed = agent_max_speed,
        arena_height = arena_height
    )

    push!(target_x_history, target.position[1])
    push!(target_y_history, target.position[2])

    # # Compute the current sheaf disagreement for horizontal sheaf
    # z = vcat(agent_horizontal.position, target.position)
    # residual = norm(δ_horizontal * z)

    # # Compute the current sheaf disagreement for vertical sheaf
    # z = vcat(agent_vertical.position, target.position)
    # residual = norm(δ_vertical * z)

    plot(
        xlims = (0.0, arena_width),
        ylims = (0.0, arena_height),
        aspect_ratio = :equal,
        xlabel = "x",
        ylabel = "y",
        title = "Target Tracking with 1 Target and 2 Agents",
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

    # Vertical agent track
    plot!(
        [x_track, x_track],
        [0.0, arena_height];
        linestyle = :dash,
        linewidth = 2,
        color = "#3A6EA5",
        label = false
    )

    # Target trajectory
    plot!(
        target_x_history,
        target_y_history;
        linewidth = 1.5,
        color = "#D4A017",
        label = "Target trajectory"
    )

    # Vertical projection of the target onto the agent's track
    plot!(
        [target.position[1], target.position[1]],
        [target.position[2], y_track];
        linestyle = :dot,
        color = "#D4A017",
        label = false
    )

    # Horizontal projection of the target onto the agent's track
    plot!(
        [target.position[1], x_track],
        [target.position[2], target.position[2]];
        linestyle = :dot,
        color = "#D4A017",
        label = false
    ) 

    scatter!(
        [target.position[1]],
        [target.position[2]];
        markersize = 8,
        color = "#D4A017",
        label = "Target"
    )

    scatter!(
        [agent_horizontal.position[1]],
        [agent_horizontal.position[2]];
        markersize = 8,
        color = "#FF5733",
        label = "Horizontal Agent"
    )

    scatter!(
        [agent_vertical.position[1]],
        [agent_vertical.position[2]];
        markersize = 8,
        color = "#3A6EA5",
        label = "Vertical Agent"
    )
end

mp4(
    animation,
    "case2.mp4";
    fps = 30
)

println("Coboundary map:")
display(Matrix(δ_horizontal))

println("Sheaf Laplacian:")
display(Matrix(L_horizontal))

println("Coboundary map:")
display(Matrix(δ_vertical)) 

println("Sheaf Laplacian:")
display(Matrix(L_vertical))

