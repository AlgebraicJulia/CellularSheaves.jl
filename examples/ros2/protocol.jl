# protocol.jl — Julia side of the sheaf/ROS 2 wire protocol (v2).
#
# Mirrors examples/ros2/ros2_ws/src/sheaf_ros2/sheaf_ros2/protocol.py EXACTLY.
# Change the message format here AND in protocol.py together — nowhere else.
#
# v2 vs the SheafSim v1 protocol (SheafSim/julia/ground/protocol.jl):
#   v1 is unicycle-specific: state (x, y, psi, v, omega), control (v, omega).
#   v2 is dynamics-agnostic: state is a 3-D pose/twist, control is a world-frame
#   ENU velocity. The conversion to (v, omega) for a differential-drive base, or
#   to a PX4 velocity setpoint, happens in the ROS 2 node on the robot — the
#   sheaf controller never learns about wheels.
#
# Frames: everything is ENU, in the mocap / Gazebo world frame. Metres, rad, s.

using JSON3

const PROTOCOL_VERSION = 2

struct AgentState
    id   :: Int
    kind :: Symbol            # :px4 | :diffdrive
    p    :: Vector{Float64}   # position   [x, y, z]  (ENU, world)
    v    :: Vector{Float64}   # velocity   [vx, vy, vz]
    yaw  :: Float64           # heading about +z
    ok   :: Bool              # false => state is stale/invalid, hold this agent
end

struct SystemState
    t      :: Float64
    agents :: Vector{AgentState}
    stop   :: Bool            # client is going away; shut the server down cleanly
end

SystemState(t, agents) = SystemState(t, agents, false)

struct AgentControl
    id       :: Int
    u        :: Vector{Float64}   # commanded world-frame ENU velocity [ux, uy, uz]
    yaw_rate :: Float64
end

struct TargetPose
    id :: Int
    p  :: Vector{Float64}
end

struct ControlOutput
    agents  :: Vector{AgentControl}
    targets :: Vector{TargetPose}
    diag    :: Dict{String, Any}
end

ControlOutput(agents, targets) = ControlOutput(agents, targets, Dict{String,Any}())

function parse_state(raw::AbstractString)::SystemState
    d = JSON3.read(raw)
    haskey(d, :version) && d.version != PROTOCOL_VERSION &&
        error("protocol mismatch: got v$(d.version), expected v$PROTOCOL_VERSION")
    agents = [AgentState(
        a.id,
        Symbol(a.kind),
        Float64[a.p[1], a.p[2], a.p[3]],
        Float64[a.v[1], a.v[2], a.v[3]],
        Float64(a.yaw),
        Bool(a.ok),
    ) for a in d.agents]
    return SystemState(Float64(d.t), agents, get(d, :stop, false))
end

function format_control(ctrl::ControlOutput)::String
    JSON3.write(Dict(
        "version" => PROTOCOL_VERSION,
        "agents"  => [Dict(
            "id"       => a.id,
            "u"        => a.u,
            "yaw_rate" => a.yaw_rate,
        ) for a in ctrl.agents],
        "targets" => [Dict("id" => t.id, "p" => t.p) for t in ctrl.targets],
        "diag"    => ctrl.diag,
    ))
end
