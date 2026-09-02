# server.jl -- sheaf control server for the ROS 2 / Gazebo Autonomy Park demo.
#
#   julia --project=examples/ros2 examples/ros2/server.jl
#   julia --project=examples/ros2 examples/ros2/server.jl --config path/to/swarm.json
#
# Speaks protocol v2 (examples/ros2/protocol.jl) over a ZMQ REP socket. The ROS 2
# bridge sends the swarm state, this replies with one world-frame ENU velocity per
# agent. Nothing here knows about wheels, rotors, MAVLink or DDS -- that is the
# bridge's job.
#
# The control law is the sheaf Laplacian gradient flow, split into two coboundaries
# so consensus and tracking carry independent gains:
#
#     u_i = -k_c (delta_c' delta_c z)_i  -  k_t (delta_p' delta_p z)_i
#
# where z stacks the agent stalks and the target stalks, delta_c is the coboundary
# of the agent-agent (consensus) edges and delta_p that of the agent-target
# (pinning) edges. Heterogeneity is carried entirely by the restriction maps: an
# air-ground edge restricts both endpoints through the x-y projection, so a ground
# robot is asked to agree with a quadrotor about where it is over the ground and
# is never asked to match its altitude.

using CellularSheaves
using LinearAlgebra
using JSON3
using ZMQ

include(joinpath(@__DIR__, "protocol.jl"))

const STALK_DIM = 3

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

function argval(flag, default)
    i = findfirst(==(flag), ARGS)
    (isnothing(i) || i == length(ARGS)) ? default : ARGS[i+1]
end

const CONFIG_PATH = argval("--config", joinpath(@__DIR__, "config", "swarm.json"))
const CFG = JSON3.read(read(CONFIG_PATH, String))

# Optional per-tick CSV of the reported world-frame state and both residuals.
# Same column layout as scripts/fake_sim.py's --csv, so the same analysis scores
# a Gazebo run and a kinematic run -- that comparison IS the accuracy check.
const LOG_PATH = argval("--log", get(ENV, "SHEAF_LOG", ""))

const AGENTS    = CFG.agents
const TARGETS   = CFG.targets
const N_AGENTS  = length(AGENTS)
const N_TARGETS = length(TARGETS)
const KINDS     = [Symbol(a.kind) for a in AGENTS]
const V_MAX     = [Float64(a.v_max) for a in AGENTS]

"Per-agent formation offset o_i. The sheaf acts on the shifted cochain z_i = p_i - o_i,
so the global section is a rigid formation about the target rather than a rendezvous.
Without this every identity restriction map drives the agents onto the same point,
which is a correct global section and a mid-air collision."
const FORMATION = [Float64[a.formation[1], a.formation[2], a.formation[3]] for a in AGENTS]

const K_CONS = Float64(CFG.gains.consensus)
const K_TRACK = Float64(CFG.gains.tracking)
const Z_MIN = Float64(CFG.safety.z_min)
const Z_MAX = Float64(CFG.safety.z_max)
const ARENA_R = Float64(CFG.safety.arena_radius)

"""
Pairwise minimum separations for the CBF safety filter, in metres, read from an
optional `safety.min_separation` block. Defaults apply when the block is absent
so an old config still runs, and the filter still enforces something.

Three radii because the platforms are not interchangeable. Two quadrotors share
one volume of air and need the rotor wash margin; a quadrotor at 1.4 m and a
rover 1.2 m tall in the worst case cannot touch however close their ground
tracks are, so an air-ground pair is scored on HORIZONTAL distance only and gets
the smaller radius; two rovers meet in the plane by construction. Scoring an
air-ground pair in R^3 would make the shipped scenario permanently unsafe on
paper -- the rover flies directly beneath the triangle's centroid, 1.2 m away
horizontally and 1.3 m below -- and a filter that is always on is a filter that
is always fighting the controller.

`alpha` is the class-K gain in hdot >= -alpha h (Ames et al., ECC 2019). Larger
is more permissive and intervenes later; 2.0 lets a pair close to about 1 m/s of
approach speed one radius out, which is inside what these platforms do.
"""
const SEP_CFG = haskey(CFG.safety, :min_separation) ? CFG.safety.min_separation : nothing
sep_field(name::Symbol, default::Float64) =
    (SEP_CFG !== nothing && haskey(SEP_CFG, name)) ? Float64(SEP_CFG[name]) : default

const D_AIR_AIR       = sep_field(:air_air, 1.0)
const D_AIR_GROUND    = sep_field(:air_ground, 0.8)     # horizontal distance
const D_GROUND_GROUND = sep_field(:ground_ground, 0.8)  # horizontal distance
const CBF_ALPHA       = sep_field(:alpha, 2.0)

"""
Static keep-out zones for the CBF safety filter, read from an optional
`safety.keep_out` list. The list is absent by default and the shipped
`config/swarm.json` declares none, so an existing config runs unchanged with no
zone rows and the filter behaves exactly as it did.

A zone is either a `cylinder`, a vertical column measured on HORIZONTAL distance
and constraining an agent only while its altitude is inside the column's z span,
or a `sphere`, measured in R^3. The column form is the one the park actually
needs: a pillar, a net stand or a light mast is an obstacle at the height it
occupies and open air above it. Scoring such a thing in R^3 would either fence
off the airspace over it or, sized to leave that airspace usable, let a rover
drive through its base.

`applies` selects which platforms a zone constrains -- "all", "air" or "ground"
-- because the things worth avoiding are not the same for both. The pilot area
and the safety net matter to a quadrotor; a pillar matters to everything; a
landing pad you do not want overflown matters only to the air, and fencing the
rover out of it would be a bug rather than caution.
"""
struct KeepOutZone
    name::String
    column::Bool            # true: horizontal distance over a z span; false: R^3
    center::Vector{Float64}
    radius::Float64
    z_lo::Float64           # padded column span; unused for a sphere
    z_hi::Float64
    air::Bool               # does it constrain the quadrotors
    ground::Bool            # does it constrain the rovers
end

"""
Vertical margin added to each end of a column's z span before the row switches
on.

A column's row is a purely HORIZONTAL correction: its gradient has no z
component, so the row can never move an agent back across the z boundary that
switched it on, and it cannot chatter the way a self-triggering constraint
would. The pad is not there to break a feedback loop, because there is none. It
is there so that the one discontinuity that remains -- the horizontal command
stepping as an agent climbs out of the span under some other command -- happens
0.25 m clear of the obstacle rather than exactly at its lip, which is the worst
place for a command to step given the inner-loop lag these platforms have. It
also makes every column slightly taller than declared, which is the safe
direction to be wrong about an obstacle.
"""
const ZONE_Z_PAD = 0.25

function parse_keep_out()
    haskey(CFG.safety, :keep_out) || return KeepOutZone[]
    zones = KeepOutZone[]
    for (n, z) in enumerate(CFG.safety.keep_out)
        name = haskey(z, :name) ? String(z.name) : "zone$n"
        kind = String(z.kind)
        kind in ("cylinder", "sphere") ||
            error("keep_out '$name': kind must be \"cylinder\" or \"sphere\", got \"$kind\"")
        c = Float64[z.center[1], z.center[2], z.center[3]]
        r = Float64(z.radius)
        r > 0.0 || error("keep_out '$name': radius must be positive, got $r")
        column = kind == "cylinder"
        h = column ? Float64(z.height) : 0.0
        if column && !(h > 0.0)
            error("keep_out '$name': a cylinder needs a positive height, got $h")
        end
        applies = haskey(z, :applies) ? String(z.applies) : "all"
        applies in ("all", "air", "ground") ||
            error("keep_out '$name': applies must be \"all\", \"air\" or \"ground\"")
        push!(zones, KeepOutZone(name, column, c, r,
                                 c[3] - ZONE_Z_PAD, c[3] + h + ZONE_Z_PAD,
                                 applies != "ground", applies != "air"))
    end
    return zones
end

const KEEP_OUT = parse_keep_out()

"""
Hold every command at zero until every configured agent has joined the active
sheaf, read from an optional `safety.synchronized_start` and ON by default.

The default is on because of what the platforms do before they are active. A
quadrotor's adapter flies its own climb to `z_ref` and ignores the command it is
given until it is up (`Px4Adapter.send`), and while it is climbing the server
has it marked inactive, so it is invisible to the barrier as well. An agent that
is already at its slot altitude therefore crosses the airspace of a climbing one
with nothing enforcing separation, because one of the two is not being commanded
by this controller at all. Fleet arming is staggered by seconds, so that window
is not a corner case, it is every takeoff.

The gate does not stop the climbs: they happen in each vehicle's own column,
from spawn points that are already separated, which is safe by construction. It
stops TRANSLATION, which is the part that needs the other robots to exist.

The gate LATCHES open: it opens once, on the first tick when the whole swarm is
active, and never closes again. Re-arming it on any later dropout would freeze
the whole formation mid-flight every time one agent blipped out of the active
set, which is a worse failure than the one being prevented -- the remaining
agents would coast with no command while the dropped one keeps moving. After the
start, a dropout is handled where it should be: the dropped agent is commanded
zero, it leaves the coboundaries, and it stays in the barrier as an obstacle.

An agent that never becomes ready must not hold the rest of the swarm forever.
Waiting indefinitely sounds like the safe reading, and it is not: observed in
Gazebo, one quadrotor passed preflight, was later refused arming for a health
failure and exited, and the gate then held every other agent and the rover at
zero command for the whole session. The one thing still moving was the adapters'
own climb, so the operator saw a frozen swarm that did not track at all, caused
by a vehicle that had removed itself. A member missing at takeoff is the normal
case on a fleet, not an exception.

So the gate also opens on a timeout, `safety.start_timeout` seconds after the
first tick, with whoever is active at that moment, and says loudly who it gave
up on. The agents that never joined are not commanded (they are still inactive)
and they remain in the barrier as obstacles, so starting without them is not the
same as ignoring them. Set the timeout to zero to wait forever, which is the old
behaviour and a deliberate operator choice rather than the default.
"""
const SYNC_START = haskey(CFG.safety, :synchronized_start) ?
                   Bool(CFG.safety.synchronized_start) : true

"Seconds to wait for the full swarm before starting with whoever is up. Zero
waits forever. Sized well past a staggered fleet arming and climb, which is tens
of seconds, so a healthy swarm always starts together and only a genuinely
missing member trips it."
const START_TIMEOUT = haskey(CFG.safety, :start_timeout) ?
                      Float64(CFG.safety.start_timeout) : 45.0

"Latch for [`SYNC_START`](@ref). Set on the first tick the whole swarm is active."
const STARTED = Ref(false)

"Timestamp of the first tick, on the same clock the state messages carry, so the
timeout means the same thing in the kinematic harness, in a Gazebo run at half
real time, and on hardware. Wall clock would make it fire at a different point
in the flight in each of those, and would be untestable in a harness that runs
faster than real time."
const FIRST_TICK = Ref(NaN)

"SHEAF_CBF=0 runs the controller with the safety filter removed. Exists so the
filter can be measured against its own absence on the same scenario (see
config/collide.json) rather than against a remembered number."
const CBF_ON = get(ENV, "SHEAF_CBF", "1") != "0"

# ---------------------------------------------------------------------------
# Restriction maps
# ---------------------------------------------------------------------------

const IDENTITY   = Matrix{Float64}(I, STALK_DIM, STALK_DIM)
const PROJECT_XY = [1.0 0.0 0.0; 0.0 1.0 0.0]      # R^3 -> R^2, the ground plane
const PROJECT_Z  = [0.0 0.0 1.0]                   # R^3 -> R^1, altitude alone

"""
Restriction map named in the config's pinning list. These are the whole content
of the scenario: an edge carrying `PROJECT_XY` asks its endpoints to agree about
position over the ground and says nothing about altitude, while one carrying
`PROJECT_Z` asks them to agree about altitude and says nothing about where they
are. The same target vertex constrains a rover and a quadrotor completely
differently depending only on which map sits on the edge between them.
"""
function named_restriction(name::AbstractString)
    name == "xy"       && return PROJECT_XY
    name == "z"        && return PROJECT_Z
    name == "identity" && return IDENTITY
    error("unknown restriction map: $name (expected xy, z or identity)")
end

"""
Restriction maps for an edge between two vertices of the given kinds.

Same kind (air-air, ground-ground) restricts through the identity: the two stalks
are asked to agree in all of R^3. Mixed kind restricts both endpoints through the
x-y projection, so the edge only constrains the shared ground-plane coordinates.
A ground robot cannot climb, and demanding it match a quadrotor's altitude would
make the sheaf have no global section at all.
"""
function edge_maps(kind1::Symbol, kind2::Symbol)
    planar = (kind1 === :diffdrive) || (kind2 === :diffdrive)
    planar ? (PROJECT_XY, PROJECT_XY) : (IDENTITY, IDENTITY)
end

"Vertex ordering: agents 1..N_AGENTS, then targets N_AGENTS+1 .. N_AGENTS+N_TARGETS."
target_vertex(k::Int) = N_AGENTS + k

"""
Per-agent velocity feedforward operators: a list of `(target index, P)` where
`P = M' M` projects a target velocity onto the subspace that agent is actually
constrained to match.

This projection is the whole reason the scenario works. An agent pinned through
`PROJECT_Z` must rise and fall with the target but is free horizontally, so it
gets `diag(0,0,1) v` and nothing else -- handing it the full target velocity
would drag the triangle around a circle it was never asked to fly. An agent
pinned through `PROJECT_XY` gets `diag(1,1,0) v` for the same reason in reverse.
Feed the wrong one and the formation still converges, just to the wrong motion.

Agents with no pinning edge of their own inherit the operators of the nearest
pinned agent along consensus edges: they are held in formation relative to it,
so the section they belong to translates the same way.
"""
function build_feedforward()
    direct = [Tuple{Int, Matrix{Float64}}[] for _ in 1:N_AGENTS]
    for pr in CFG.pinning
        i, k = Int(pr.agent) + 1, Int(pr.target) + 1
        M = named_restriction(String(pr.map))
        push!(direct[i], (k, Matrix{Float64}(M' * M)))
    end

    adj = [Int[] for _ in 1:N_AGENTS]
    for e in CFG.edges
        i, j = Int(e[1]) + 1, Int(e[2]) + 1
        push!(adj[i], j); push!(adj[j], i)
    end

    ops = copy(direct)
    for i in 1:N_AGENTS
        isempty(ops[i]) || continue
        seen = falses(N_AGENTS); seen[i] = true
        queue = [i]
        while !isempty(queue)
            v = popfirst!(queue)
            if !isempty(direct[v])
                ops[i] = direct[v]
                break
            end
            for w in adj[v]
                seen[w] && continue
                seen[w] = true; push!(queue, w)
            end
        end
    end
    return ops
end

const FEEDFORWARD = build_feedforward()

"""
Combining several pins' feedforward. Summing `P_k v_k` is right only when the
maps are complementary (xy from one target, z from another). Two identity pins
would sum to twice the velocity and the agent overshoots every tick, which shows
up as a tracking residual that never closes. The least-squares combination
`(sum P_k)^+ sum P_k v_k` is the sum for disjoint subspaces and the average for
overlapping ones, so it is right in both cases. Precomputed per agent.
"""
const FF_SOLVE = [isempty(ops) ? zeros(3, 3) : pinv(sum(P for (_, P) in ops))
                  for ops in FEEDFORWARD]

"""
One pinning edge, resolved. A pinning entry may carry its own `offset` o_ik, in
which case the edge enforces `p_i - o_ik = T_k` instead of the default
`p_i - o_i = T_k` with the agent's formation offset.

Why this exists: consensus and pinning both act on the shifted cochain
z_i = p_i - o_i, so an agent whose consensus neighbours track a different target
than it does has no configuration satisfying both, and the flow settles at a
compromise. A per-pin offset removes the conflict for the common geometric case
of a formation tracking several targets that move rigidly together: choose
o_ik = o_i - (T_k - T_bar) and every pin agrees on z_i = T_bar, which is exactly
the section consensus wants. The scenario editor computes these.

Implemented without touching the sheaf machinery: such a pin gets its own
*virtual* target vertex whose value is T_k + (o_ik - o_i), because
p_i - o_ik = T_k rearranges to z_i = p_i - o_i = T_k + (o_ik - o_i).
"""
struct PinSpec
    agent::Int
    target::Int
    map::String
    vertex::Int
    shift::Vector{Float64}
end

const PINS = let pins = PinSpec[], nv = N_AGENTS + N_TARGETS
    for pr in CFG.pinning
        i, k = Int(pr.agent) + 1, Int(pr.target) + 1
        if haskey(pr, :offset)
            o = Float64[pr.offset[1], pr.offset[2], pr.offset[3]]
            nv += 1
            push!(pins, PinSpec(i, k, String(pr.map), nv, o .- FORMATION[i]))
        else
            push!(pins, PinSpec(i, k, String(pr.map), target_vertex(k), zeros(3)))
        end
    end
    pins
end

const NV_TOTAL = N_AGENTS + N_TARGETS + count(p -> p.vertex > N_AGENTS + N_TARGETS, PINS)

"""
Build the two sheaves whose coboundaries drive the controller, restricted to the
subgraph of *active* agents: `s_cons` keeps a consensus edge only when both
endpoints are active, `s_pin` keeps a pinning edge only when its agent is. Both
live on the full vertex set so their coboundaries act on the same stacked cochain.

The restriction is the point, not an optimization. An agent that is not flying
(failed preflight, mid-takeoff, lost mocap) must exert no force on the agents
that are: consensus with a grounded vehicle's drifting state estimate was
observed dragging the whole formation's altitude down. The active sheaf is the
subcomplex of live agents; vertices rejoin, edges and all, the tick they report
ready.
"""
function build_sheaves(mask::AbstractVector{Bool} = trues(N_AGENTS))
    s_cons = EuclideanSheaf{Float64}(fill(STALK_DIM, NV_TOTAL))
    s_pin  = EuclideanSheaf{Float64}(fill(STALK_DIM, NV_TOTAL))

    for e in CFG.edges
        i, j = Int(e[1]) + 1, Int(e[2]) + 1          # config is 0-based, Julia is 1-based
        (mask[i] && mask[j]) || continue
        rm1, rm2 = edge_maps(KINDS[i], KINDS[j])
        add_sheaf_edge!(s_cons, i, j, rm1, rm2)
    end

    for p in PINS
        mask[p.agent] || continue
        M = named_restriction(p.map)
        add_sheaf_edge!(s_pin, p.agent, p.vertex, M, M)
    end

    return s_cons, s_pin
end

"A sheaf with no edges has no coboundary worth taking; `nothing` marks a term
that contributes zero gradient."
cob_or_nothing(s) = isempty(edge_stalks(s)) ? nothing : coboundary_map(s)

"Coboundary pair for a readiness mask, cached by bitmask. Masks change a handful
of times per flight (agents coming up, an agent dropping out), so the cache stays
tiny while the steady state pays one dictionary lookup."
const MASKED_COBOUNDARIES = Dict{UInt64, Tuple{Any, Any}}()

function masked_coboundaries(mask::AbstractVector{Bool})
    key = UInt64(0)
    for b in mask
        key = (key << 1) | UInt64(b)
    end
    get!(MASKED_COBOUNDARIES, key) do
        s_cons, s_pin = build_sheaves(mask)
        (cob_or_nothing(s_cons), cob_or_nothing(s_pin))
    end
end

# ---------------------------------------------------------------------------
# Targets
# ---------------------------------------------------------------------------

"""
Position and velocity along a `waypoints` target: a piecewise-linear closed tour
of `spec.points` traversed at constant `spec.speed` [m/s]. The velocity is the
current segment's unit direction times the speed, so the feedforward is exact
everywhere except the corners, where it jumps; at the speeds this arena allows
the inner loops absorb the jump. Drawn paths from `tools/swarm_editor.html`
land here.
"""
function waypoint_state(spec, t::Float64)
    pts = [Float64[p[1], p[2], p[3]] for p in spec.points]
    n = length(pts)
    n == 1 && return pts[1], zeros(3)
    speed = Float64(spec.speed)
    segs = [pts[mod1(i + 1, n)] .- pts[i] for i in 1:n]
    lens = [norm(seg) for seg in segs]
    total = sum(lens)
    (total <= 1e-9 || speed <= 0) && return pts[1], zeros(3)
    s = mod(speed * t, total)
    for i in 1:n
        if s <= lens[i] || i == n
            d = segs[i] ./ max(lens[i], 1e-9)
            return pts[i] .+ d .* s, d .* speed
        end
        s -= lens[i]
    end
end

"Position of target `k` at time `t`. Add new `kind`s here (and, if they move,
in [`target_velocity`](@ref)) and they are immediately available to the sheaf."
function target_position(k::Int, t::Float64)
    spec = TARGETS[k]
    kind = Symbol(spec.kind)
    # `center` is read lazily: a drawn path carries its geometry in `points` and
    # has no centre at all.
    kind === :waypoints && return waypoint_state(spec, t)[1]
    c = Float64[spec.center[1], spec.center[2], spec.center[3]]
    if kind === :static
        return c
    elseif kind === :circle
        omega = 2pi / Float64(spec.period)
        r = Float64(spec.radius)
        return c .+ [r * cos(omega * t), r * sin(omega * t), 0.0]
    elseif kind === :orbit3d
        # Horizontal orbit and a vertical bob on an unrelated period, so the two
        # motions never lock: the rover's job and the quadrotors' job are visibly
        # independent even though both come from this one point.
        omega = 2pi / Float64(spec.period)
        r = Float64(spec.radius)
        wz = 2pi / Float64(spec.z_period)
        az = Float64(spec.z_amp)
        return c .+ [r * cos(omega * t), r * sin(omega * t), az * sin(wz * t)]
    elseif kind === :lissajous
        # One parameterisation covering circle, ellipse, figure-eight and the
        # full 2-D/3-D Lissajous family:
        #   p_i(t) = c_i + A_i sin(n_i w t + phi_i),   w = 2pi/period
        # A circle is A=(r,r,0), n=(1,1,0), phi=(0,90,0) degrees; a figure eight
        # is n=(1,2,0). Phases are in degrees because a human types them.
        w = 2pi / Float64(spec.period)
        return c .+ [Float64(spec.amp[i]) *
                     sin(Float64(spec.freq[i]) * w * t + deg2rad(Float64(spec.phase[i])))
                     for i in 1:3]
    else
        error("unknown target kind: $kind")
    end
end

"""
Velocity of target `k` at time `t` -- the analytic derivative of
[`target_position`](@ref). Fed forward to every active agent: the sheaf's
equilibrium section translates rigidly with the target, so adding the section's
own velocity turns the gradient flow's steady-state chase into regulation about
a fixed point of the moving frame. Without it a proportional law trails the
target by an offset proportional to speed, and -- worse, in full physics --
each platform trails by a *different* amount set by its inner-loop lag, which
shears the formation shape. Observed: feedforward off, Gazebo formation error
~0.54 m and consensus residual ~1.0; the kinematic harness shows the same
structure at ~0.16 m.
"""
function target_velocity(k::Int, t::Float64)
    spec = TARGETS[k]
    kind = Symbol(spec.kind)
    if kind === :circle || kind === :orbit3d
        omega = 2pi / Float64(spec.period)
        r = Float64(spec.radius)
        vz = 0.0
        if kind === :orbit3d
            wz = 2pi / Float64(spec.z_period)
            vz = Float64(spec.z_amp) * wz * cos(wz * t)
        end
        return [-r * omega * sin(omega * t), r * omega * cos(omega * t), vz]
    elseif kind === :lissajous
        w = 2pi / Float64(spec.period)
        return [Float64(spec.amp[i]) * Float64(spec.freq[i]) * w *
                cos(Float64(spec.freq[i]) * w * t + deg2rad(Float64(spec.phase[i])))
                for i in 1:3]
    elseif kind === :waypoints
        return waypoint_state(spec, t)[2]
    end
    return zeros(3)
end

# ---------------------------------------------------------------------------
# Control law
# ---------------------------------------------------------------------------

block(z, v) = @view z[(v-1)*STALK_DIM+1 : v*STALK_DIM]

"""
    sheaf_control(delta_c, delta_p, state) -> (controls, diagnostics)

One control tick. Stacks the reported agent positions and the current target
positions into a cochain `z`, applies the two sheaf Laplacians, and returns the
saturated per-agent velocity command.

Agents reporting `ok = false` (stale mocap, disarmed, dropped link) are commanded
to zero and their reported position is still used for their neighbours, which is
the safe reading: a robot we cannot see is a robot we do not move.
"""
const LAST_GOOD = [zeros(STALK_DIM) for _ in 1:N_AGENTS]

"Consecutive trusted ticks per agent. An agent joins the active sheaf only after
JOIN_TICKS consecutive good reports but leaves it on the first bad one --
asymmetric on purpose. Every join/leave changes the edge set, so a trust flap
from a single late odometry message would otherwise make the Laplacian (and the
residual) jump every time; distrusting instantly stays cheap because a dropped
agent is only ever commanded zero."
const JOIN_TICKS = 10
const GOOD_STREAK = zeros(Int, N_AGENTS)

"An estimate is trustworthy if it is finite and plausible for a live member of
this arena: horizontally near the geofence, vertically inside a band a robot of
that kind can actually occupy (fence height plus overshoot margin for air,
ground level for a rover). Anything outside is a broken state estimate -- a
crashed vehicle dead-reckoning at thousands of metres, or a grounded one whose
EKF drifted metres up while it waited to arm. Feeding either to the Laplacian
sends every neighbour chasing the ghost, so a distrusted agent is excluded from
the active sheaf entirely and commanded to zero. Mocap glitches on hardware hit
this same gate."
function trusted(p::Vector{Float64}, kind::Symbol)
    all(isfinite, p) || return false
    hypot(p[1], p[2]) <= ARENA_R + 2.0 || return false
    kind === :diffdrive && return abs(p[3]) <= 0.5
    return -0.5 <= p[3] <= Z_MAX + 1.0
end

function sheaf_control(state::SystemState)
    z = zeros(STALK_DIM * NV_TOTAL)

    active = falses(N_AGENTS)   # an agent absent from the message stays inactive
    seen = falses(N_AGENTS)
    # Reported this tick at a position the arena admits, whether or not the agent
    # is active. This is the OBSTACLE mask: a quadrotor still climbing, one that
    # failed preflight, one whose mocap just dropped -- none of them is being
    # commanded by this controller, and every one of them is still a solid object
    # in the airspace. `a.ok` is deliberately not required here. It says whether
    # we may steer the thing, not whether it is there.
    present = falses(N_AGENTS)
    for a in state.agents
        v = a.id + 1                                  # protocol ids are 0-based
        seen[v] = true
        present[v] = trusted(a.p, KINDS[v])
        ok = a.ok && trusted(a.p, KINDS[v])
        if ok
            LAST_GOOD[v] .= a.p
            GOOD_STREAK[v] += 1
        else
            GOOD_STREAK[v] = 0
        end
        active[v] = GOOD_STREAK[v] >= JOIN_TICKS
        z[(v-1)*STALK_DIM+1 : v*STALK_DIM] = LAST_GOOD[v] .- FORMATION[v]
    end
    for v in 1:N_AGENTS
        seen[v] || (GOOD_STREAK[v] = 0)
    end

    # Synchronized start. Latched, so this is the only place the gate can open
    # and nothing can shut it again; see [`SYNC_START`](@ref) for why.
    if SYNC_START && !STARTED[]
        isnan(FIRST_TICK[]) && (FIRST_TICK[] = state.t)
        if all(active)
            STARTED[] = true
        elseif START_TIMEOUT > 0 && state.t - FIRST_TICK[] > START_TIMEOUT
            STARTED[] = true
            missing = [String(AGENTS[i].name) for i in findall(.!active)]
            @warn "starting without $(join(missing, ", ")): not ready after $(round(Int, START_TIMEOUT)) s. They stay uncommanded and are avoided as obstacles."
        end
    end
    started = !SYNC_START || STARTED[]
    waiting = started ? Int[] : findall(.!active)
    delta_c, delta_p = masked_coboundaries(active)
    targets = [target_position(k, state.t) for k in 1:N_TARGETS]
    for k in 1:N_TARGETS
        v = target_vertex(k)
        z[(v-1)*STALK_DIM+1 : v*STALK_DIM] = targets[k]
    end
    for p in PINS
        p.vertex > N_AGENTS + N_TARGETS || continue
        z[(p.vertex-1)*STALK_DIM+1 : p.vertex*STALK_DIM] = targets[p.target] .+ p.shift
    end

    grad = zeros(length(z))
    cons_norm = 0.0
    track_norm = 0.0
    if delta_c !== nothing
        cons_residual = delta_c * z
        cons_norm = norm(cons_residual)
        grad .+= K_CONS .* (delta_c' * cons_residual)
    end
    if delta_p !== nothing
        pin_residual = delta_p * z
        track_norm = norm(pin_residual)
        grad .+= K_TRACK .* (delta_p' * pin_residual)
    end

    tvel = [target_velocity(k, state.t) for k in 1:N_TARGETS]

    # The nominal commands for every agent first, because the safety filter is a
    # statement about the swarm and cannot be assembled one agent at a time.
    u_nom = [zeros(3) for _ in 1:N_AGENTS]
    pos   = [copy(LAST_GOOD[i]) for i in 1:N_AGENTS]
    for a in state.agents
        i = a.id + 1
        pos[i] = a.p
        ff = zeros(3)
        for (k, P) in FEEDFORWARD[i]
            ff .+= P * tvel[k]
        end
        ff = FF_SOLVE[i] * ff
        u_nom[i] = ff .- Vector(block(grad, i))
    end

    # Safety filter between the sheaf gradient and the platform limits. The
    # fences stay AFTER it, unchanged: a separation correction must never be
    # able to talk a vehicle through the geofence or the ceiling, and the
    # residual risk of the reverse (the direction-preserving saturation scaling
    # a correction down) is bounded because the filter is minimally invasive and
    # the nominal command it corrects is already inside the speed limit.
    cbf_diag = cbf_filter!(u_nom, pos, active, present)

    controls = AgentControl[]
    for a in state.agents
        i = a.id + 1
        u = u_nom[i]

        u = arena_fence(u, a.p)
        u = saturate(u, V_MAX[i])
        if KINDS[i] === :diffdrive
            u[3] = 0.0                                 # ground robots do not climb
        else
            # The fence is applied AFTER saturation and clamped separately.
            # Direction-preserving saturation scales the whole vector, so a huge
            # horizontal command (e.g. chasing a distant neighbour) would scale a
            # pre-saturation fence correction down to nothing and let the vehicle
            # sink through the floor while "respecting" its speed limit.
            u[3] = clamp(altitude_fence(u[3], a.p[3]), -V_MAX[i], V_MAX[i])
        end
        # Nothing translates before the gate opens, and an agent that is not
        # active is never commanded. The climb the adapter flies on its own is
        # untouched by either: it does not come from here.
        (active[i] && started) || (u = zeros(3))

        push!(controls, AgentControl(a.id, u, 0.0))
    end

    # Per-edge residuals for the RViz overlay. Computed here, from the same
    # restriction maps the coboundaries were built from, because the sheaf math
    # has exactly one home: a visualizer that re-derives which map sits on which
    # edge will eventually disagree with the controller and then lie with
    # authority. Edges of inactive agents report active=false and res=0 -- they
    # are drawn as absent, which is what they are.
    edge_diag = Vector{Dict{String,Any}}()
    for e in CFG.edges
        i, j = Int(e[1]) + 1, Int(e[2]) + 1
        act = active[i] && active[j]
        res = 0.0
        if act
            rm1, rm2 = edge_maps(KINDS[i], KINDS[j])
            res = norm(rm1 * block(z, i) .- rm2 * block(z, j))
        end
        push!(edge_diag, Dict("type" => "cons", "a" => i - 1, "b" => j - 1,
                              "active" => act, "res" => round(res, digits = 4)))
    end
    for p in PINS
        act = active[p.agent]
        res = 0.0
        if act
            M = named_restriction(p.map)
            res = norm(M * block(z, p.agent) .- M * block(z, p.vertex))
        end
        push!(edge_diag, Dict("type" => "pin", "a" => p.agent - 1, "b" => p.target - 1,
                              "map" => p.map,
                              "active" => act, "res" => round(res, digits = 4)))
    end

    diag = Dict{String,Any}(
        "consensus_residual" => cons_norm,
        "tracking_residual"  => track_norm,
        "active_agents"      => count(active),
        "active"             => [active[i] for i in 1:N_AGENTS],
        "started"            => started,
        "waiting_for"        => [String(AGENTS[i].name) for i in waiting],
        "edges"              => edge_diag,
    )
    merge!(diag, cbf_diag)
    return controls, [TargetPose(k - 1, targets[k]) for k in 1:N_TARGETS], diag
end

"Scale a command down (never up) to respect the platform speed limit, preserving
direction so the command stays parallel to the sheaf gradient."
function saturate(u::Vector{Float64}, limit::Float64)
    n = norm(u)
    (n <= limit || n == 0.0) && return u
    return u .* (limit / n)
end

"Soft ceiling/floor. The park's Astro SOP caps altitude at 2.5 m; `z_max` in the
config should stay under that."
function altitude_fence(uz::Float64, z::Float64; k = 1.0)
    z > Z_MAX && return min(uz, 0.0) - k * (z - Z_MAX)
    z < Z_MIN && return max(uz, 0.0) + k * (Z_MIN - z)
    return uz
end

"Soft cylindrical geofence about the world origin."
function arena_fence(u::Vector{Float64}, p::Vector{Float64}; k = 1.5)
    r = hypot(p[1], p[2])
    r <= ARENA_R && return u
    n = [p[1], p[2], 0.0] ./ max(r, eps())
    outward = dot(u, n)
    return u .- max(outward, 0.0) .* n .- k * (r - ARENA_R) .* n
end

# ---------------------------------------------------------------------------
# Safety filter: control barrier functions
# ---------------------------------------------------------------------------
#
# The sheaf gradient produces a command that is correct about the FORMATION and
# says nothing about collisions: a global section is a geometric agreement, and
# two agents converging to their slots from opposite sides pass through the same
# air on the way. This filter sits between the gradient and the platform limits
# and makes separation a constraint instead of a hope.
#
# For every pair of active agents, the barrier
#
#     h_ij(p) = ||S_ij (p_i - p_j)||^2 - d_ij^2   >= 0
#
# where S_ij is the identity for an air-air pair and the x-y projection for any
# pair involving a rover (the same PROJECT_XY that sits on the mixed-kind
# consensus edges: a rover and a quadrotor only ever share the ground plane).
# The safe set stays forward invariant if hdot >= -alpha h along trajectories
# (Ames, Coogan, Egerstedt, Notomista, Sreenath, Tabuada, ECC 2019), and the
# controller is a velocity-level one -- from here every agent is a single
# integrator pdot = u -- so that condition is LINEAR in the command:
#
#     2 (S_ij dp)' S_ij (u_i - u_j)  >=  -alpha h_ij .
#
# Unlike examples/RL/lib/cbf_ipm.jl, which splits each pair's responsibility in
# half so each agent can solve its own filter onboard, this server already has
# every agent's command in one place. It therefore solves the JOINT QP
#
#     min_u  sum_i ||u_i - u_nom_i||^2   s.t. the constraints above,
#
# which is the same safe set with no conservatism from the split and no
# deadlock-breaking bias to tune. Splitting it per agent is what the
# decentralized phase will do; the constraint rows are already in that form.
#
# A row does not need both of its agents to be ours to command. An agent that is
# present but not active -- a quadrotor still climbing under its adapter's own
# control, one that failed preflight, one whose mocap dropped a beat ago -- is
# not steerable by this controller and is very much still a solid object in the
# airspace. Those get a ONE-SIDED row: the same barrier, but only the active
# agent's command appears in it, so the active agent does all the avoiding. That
# is the only honest form, because we have no say in what the other one does, and
# it is the form that matters at takeoff, when a staggered fleet always has
# someone climbing while someone else is already flying to a slot.
#
# Keep-out zones enter the SAME QP. A zone is static geometry, so only the one
# agent's velocity appears in its barrier
#
#     h_iz(p) = ||S_z (p_i - c_z)||^2 - r_z^2   >= 0,
#     2 (S_z (p_i - c_z))' S_z u_i  >=  -alpha h_iz ,
#
# with S_z the identity for a sphere and PROJECT_XY for a vertical column, and
# the same alpha as the pairwise rows. They are rows in the same problem rather
# than a second filter run afterwards because a robot squeezed between a
# neighbour and a wall must get ONE answer: two filters in sequence each undo
# part of the other's correction and the vehicle ends up doing whichever ran
# last, which is exactly the geometry where being wrong is expensive.

"Minimum separation for a pair and the map it is measured through: `true` means
the pair is scored horizontally (any pair involving a rover), `false` in R^3."
function pair_separation(i::Int, j::Int)
    air_i = KINDS[i] !== :diffdrive
    air_j = KINDS[j] !== :diffdrive
    air_i && air_j && return (D_AIR_AIR, false)
    (air_i || air_j) && return (D_AIR_GROUND, true)
    return (D_GROUND_GROUND, true)
end

"Separation vector for a pair, already restricted through the pair's map, so
`norm` of it is the distance the barrier is written in."
function pair_delta(pi_::Vector{Float64}, pj::Vector{Float64}, planar::Bool)
    dp = pi_ .- pj
    planar && (dp = Float64[dp[1], dp[2], 0.0])
    return dp
end

"""
Does zone `z` constrain an agent of kind `kind` at `p`, and if so the vector
from the zone centre to the agent, already restricted through the zone's own map
so that `norm` of it is the distance the barrier is written in. `nothing` when
the zone does not apply.

For a column that includes an agent above or below its (padded) z span: the row
is then simply absent, which is the whole point of the column form. A quadrotor
at 2.9 m is not obstructed by a 2.5 m pillar and must not be steered around it.
"""
function zone_delta(z::KeepOutZone, p::Vector{Float64}, kind::Symbol)
    air = kind !== :diffdrive
    (air ? z.air : z.ground) || return nothing
    dp = p .- z.center
    if z.column
        (z.z_lo <= p[3] <= z.z_hi) || return nothing
        dp = Float64[dp[1], dp[2], 0.0]
    end
    return dp
end

const CBF_SWEEPS = 200      # Hildreth sweeps; the shipped scenario needs 0
const CBF_TOL    = 1e-10    # dual step below which the sweep is done
const CBF_FEAS   = 1e-6     # constraint slack tolerated before falling back
const CBF_REPULSE = 1.5     # fallback repulsion gain [1/s]

"""
    cbf_filter!(u, p, active) -> Dict

Minimally-invasive correction of the nominal commands `u` (indexed by agent, and
overwritten in place) for agents at `p`, enforcing pairwise separation between
every pair of ACTIVE agents, clearance from every agent that is present but not
active, and clearance from every configured keep-out zone.

Every row needs at least one active agent, because an inactive agent is
commanded zero and a row that changes nobody's command is arithmetic with no
effect. `obstacle` marks the agents whose reported position this tick is fresh
and [`trusted`](@ref); an active agent gets a row against each of them, and only
the active agent's own command appears in it.

What is NOT constrained is a position we do not believe: an estimate outside the
trust band is a crashed vehicle dead-reckoning at thousands of metres or a
grounded one whose EKF drifted, and steering live agents around a ghost is how a
bad estimate takes down the swarm. Zone rows are built for active agents only,
for the mirror-image reason: a zone row exists to bend a command, and an
inactive agent has no command to bend.

The QP is solved in the dual. With rows `C_k` collecting the pair and zone
constraints as `C u >= b`, the KKT stationarity of `min 1/2||u-u_nom||^2 s.t.
Cu >= b` gives `u = u_nom + C' lambda`, and lambda solves the box-constrained
problem

    min_{lambda >= 0}  1/2 lambda' G lambda + lambda' w,
    G = C C',  w = C u_nom - b.

Hildreth's cyclic coordinate descent (Hildreth, "A quadratic programming
procedure", NRLQ 1957) solves that: each coordinate has a closed-form clamped
update, every iterate is dual feasible, and the objective decreases monotonically,
so the loop terminates at the sweep cap with a usable answer whatever happens.
G is at most (n choose 2) + n*|zones| square -- 6 rows for the shipped four
agents and no zones, 12 for the four agents with two zones that apply to
everything -- so the whole solve is a few hundred flops.

Returns the diagnostics the caller merges into its dict.
"""
function cbf_filter!(u::Vector{Vector{Float64}}, p::Vector{Vector{Float64}},
                     active::AbstractVector{Bool},
                     obstacle::AbstractVector{Bool} = active;
                     enabled::Bool = CBF_ON)
    ids = findall(active)
    # Present, trusted, and not one of ours to command: an obstacle only.
    obs = findall(i -> obstacle[i] && !active[i], 1:length(active))
    n = length(ids)
    min_sep = Inf
    min_margin = Inf
    min_zone = Inf

    # Constraint rows, stored by what they act on rather than densely: row k is
    # +g_k on agent A[k]'s block and, when B[k] > 0, -g_k on agent B[k]'s. A
    # keep-out zone is static and an obstacle agent is not ours to command, so
    # both carry B[k] = 0 and act on one command only; agent indices start at 1,
    # so that sentinel can never collide with a real one. This keeps the closed
    # form for G = C C' below and never allocates anything 3n-wide.
    A = Int[]; B = Int[]; g = Vector{Float64}[]; b = Float64[]
    dist = Float64[]        # distance the row is written in, for the fallback
    dmin = Float64[]        # the radius that distance is measured against
    is_zone = Bool[]        # zone row (vs. a pair or obstacle row)
    coincident = Tuple{Int,Int,Bool}[]   # degenerate pairs; Bool = one-sided
    coincident_zone = Tuple{Int,Int}[]   # (agent, zone) sitting on a zone's centre

    for ka in 1:n-1, kb in ka+1:n
        i, j = ids[ka], ids[kb]
        d_min, planar = pair_separation(i, j)
        dp = pair_delta(p[i], p[j], planar)
        d = norm(dp)
        min_sep = min(min_sep, d)
        min_margin = min(min_margin, d - d_min)
        if d < 1e-6
            # Two agents on top of each other: the barrier gradient vanishes, so
            # there is no command the constraint admits and the QP is infeasible
            # by construction. Handled by the repulsion fallback, which needs a
            # direction and takes it from the formation offsets (the geometry
            # that was pulling them apart in the first place).
            push!(coincident, (i, j, false))
            continue
        end
        h = d * d - d_min * d_min
        push!(A, i); push!(B, j); push!(g, 2.0 .* dp); push!(b, -CBF_ALPHA * h)
        push!(dist, d); push!(dmin, d_min); push!(is_zone, false)
    end

    # One-sided rows against the agents we can see but do not command. Same
    # barrier and same radii as a commanded pair, so an agent that drops out
    # mid-flight is avoided exactly as it was a tick earlier; the difference is
    # that all of the avoiding now falls on the agent still under command. Its
    # own motion is unmodelled -- we do not trust the velocity of something we
    # have declined to trust the state of -- which is part of why these radii
    # carry margin over the vehicle footprints.
    for i in ids, j in obs
        d_min, planar = pair_separation(i, j)
        dp = pair_delta(p[i], p[j], planar)
        d = norm(dp)
        min_sep = min(min_sep, d)
        min_margin = min(min_margin, d - d_min)
        if d < 1e-6
            push!(coincident, (i, j, true))
            continue
        end
        h = d * d - d_min * d_min
        push!(A, i); push!(B, 0); push!(g, 2.0 .* dp); push!(b, -CBF_ALPHA * h)
        push!(dist, d); push!(dmin, d_min); push!(is_zone, false)
    end

    n_pair_rows = length(b)

    for i in ids, (zi, z) in enumerate(KEEP_OUT)
        dp = zone_delta(z, p[i], KINDS[i])
        dp === nothing && continue
        d = norm(dp)
        min_zone = min(min_zone, d - z.radius)
        if d < 1e-6
            push!(coincident_zone, (i, zi))
            continue
        end
        h = d * d - z.radius * z.radius
        push!(A, i); push!(B, 0); push!(g, 2.0 .* dp); push!(b, -CBF_ALPHA * h)
        push!(dist, d); push!(dmin, z.radius); push!(is_zone, true)
    end

    m = length(b)
    n_active = 0
    n_zone_active = 0
    fallback = 0

    # Relative command a row is written in: the pair's difference, or the single
    # agent's own command for a static zone.
    urel(k) = B[k] > 0 ? (u[A[k]] .- u[B[k]]) : u[A[k]]

    if enabled && m > 0
        G = zeros(m, m)
        for k in 1:m, l in k:m
            s = 0.0
            A[k] == A[l] && (s += 1.0)
            (B[k] > 0 && B[k] == B[l]) && (s += 1.0)
            (B[l] > 0 && A[k] == B[l]) && (s -= 1.0)
            (B[k] > 0 && B[k] == A[l]) && (s -= 1.0)
            G[k, l] = G[l, k] = s * dot(g[k], g[l])
        end
        w = [dot(g[k], urel(k)) - b[k] for k in 1:m]

        # Nothing to do when the nominal command already satisfies every row --
        # the common case, and the one that must cost nothing: the filter is
        # inert unless agents are actually closing on each other or on a zone,
        # and then it costs one min over a handful of dot products.
        engaged = minimum(w) < -CBF_FEAS

        lambda = zeros(m)
        r = copy(w)                      # r = G*lambda + w, kept incrementally
        for _ in 1:(engaged ? CBF_SWEEPS : 0)
            step = 0.0
            for k in 1:m
                G[k, k] > 1e-12 || continue
                new = max(0.0, lambda[k] - r[k] / G[k, k])
                dl = new - lambda[k]
                if dl != 0.0
                    lambda[k] = new
                    @inbounds for l in 1:m
                        r[l] += G[l, k] * dl
                    end
                    step = max(step, abs(dl))
                end
            end
            step <= CBF_TOL && break
        end

        for k in 1:m
            lambda[k] <= 0.0 && continue
            n_active += 1
            is_zone[k] && (n_zone_active += 1)
            u[A[k]] .+= lambda[k] .* g[k]
            B[k] > 0 && (u[B[k]] .-= lambda[k] .* g[k])
        end

        # Degrade safely. A dual iterate is always feasible for the dual but the
        # primal is only safe if every row now holds; it may not, either because
        # the rows are jointly infeasible (agents already well inside d_min, or
        # already inside a zone, with conflicting escape directions) or because
        # the sweep cap was hit. Either way the answer is the same and it is not
        # an exception: push along the row's own separating axis, at a rate
        # proportional to how far inside it is. For a zone row that axis IS the
        # outward normal (g_k points from the zone centre to the agent, through
        # the zone's own map), so an agent that starts inside a zone is driven
        # straight out of it. Not minimally invasive and not optimal, but
        # unambiguously in the safe direction, and the alternative in a control
        # loop is a swarm with no commands at all.
        for k in 1:m
            (!engaged || dot(g[k], urel(k)) >= b[k] - CBF_FEAS) && continue
            fallback += 1
            dir = g[k] ./ norm(g[k])
            push_ = CBF_REPULSE * max(dmin[k] - dist[k], 0.1)
            u[A[k]] .+= push_ .* dir
            B[k] > 0 && (u[B[k]] .-= push_ .* dir)
        end
    end

    for (i, j, one_sided) in coincident
        enabled || break
        fallback += 1
        axis = FORMATION[i] .- FORMATION[j]
        d_min, planar = pair_separation(i, j)
        planar && (axis = Float64[axis[1], axis[2], 0.0])
        norm(axis) < 1e-9 && (axis = Float64[1.0, 0.0, 0.0])
        dir = axis ./ norm(axis)
        # A one-sided pair gets the whole separation rather than half of it: the
        # other body is not going to move out of the way on our account.
        u[i] .+= CBF_REPULSE * d_min * (one_sided ? 2.0 : 1.0) .* dir
        one_sided || (u[j] .-= CBF_REPULSE * d_min .* dir)
    end

    for (i, zi) in coincident_zone
        enabled || break
        fallback += 1
        z = KEEP_OUT[zi]
        # Dead on a sphere's centre, or on a column's axis: the barrier gradient
        # vanishes and the geometry offers no outward normal to read off. Leave
        # along the horizontal direction from the arena origin to the zone
        # centre. It is deterministic (no tie for two agents to break the same
        # way), it points out of the middle of the field rather than across it,
        # and it exists for every zone that is not itself centred on the origin;
        # one that is gets +x.
        axis = Float64[z.center[1], z.center[2], 0.0]
        norm(axis) < 1e-9 && (axis = Float64[1.0, 0.0, 0.0])
        u[i] .+= CBF_REPULSE * z.radius .* (axis ./ norm(axis))
    end

    n_zone_rows = m - n_pair_rows + length(coincident_zone)
    n_obs_rows = count(k -> !is_zone[k] && B[k] == 0, 1:m) +
                 count(c -> c[3], coincident)

    return Dict{String,Any}(
        "cbf_enabled"          => enabled,
        "cbf_constraints"      => enabled ? m + length(coincident) + length(coincident_zone) : 0,
        "cbf_active"           => n_active,
        "cbf_fallback"         => fallback,
        "cbf_obstacle_rows"    => enabled ? n_obs_rows : 0,
        "cbf_zone_constraints" => enabled ? n_zone_rows : 0,
        "cbf_zone_active"      => n_zone_active,
        "min_separation"       => isfinite(min_sep) ? round(min_sep, digits = 4) : -1.0,
        "min_margin"           => isfinite(min_margin) ? round(min_margin, digits = 4) : -1.0,
        # Smallest signed clearance to any zone that applied this tick, negative
        # inside. -1.0 is also the "no zone was evaluated" sentinel, matching
        # min_margin above; read it together with cbf_zone_constraints, which is
        # 0 exactly when no zone row existed.
        "min_zone_margin"      => isfinite(min_zone) ? round(min_zone, digits = 4) : -1.0,
    )
end

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------

"JIT-warm every function on the serve path before binding the socket. Without
this the bridge's first request lands mid-compilation: the reply takes seconds,
the bridge trips its watchdog, and the compilation CPU spike coincides with the
quadrotors' takeoff -- the least forgiving moment of the whole flight. Observed
as `timesync RTT 1992 ms` and an attitude failure at takeoff before this existed."
function warmup()
    for all_ok in (true, false), _ in 1:(JOIN_TICKS + 1)   # past the join debounce
        agents = [Dict("id" => i - 1, "kind" => String(KINDS[i]),
                       "p" => [0.4i, -0.3i, 1.0], "v" => zeros(3),
                       "yaw" => 0.0, "ok" => all_ok) for i in 1:N_AGENTS]
        req = JSON3.write(Dict("version" => PROTOCOL_VERSION, "t" => 0.0, "agents" => agents))
        state = parse_state(req)
        controls, targets, diag = sheaf_control(state)
        format_control(ControlOutput(controls, targets, diag))
    end
    foreach(lg -> lg .= 0.0, LAST_GOOD)   # do not leak warmup positions into the run
    fill!(GOOD_STREAK, 0)
    STARTED[] = false                     # nor a warmup-opened start gate
    return nothing
end

function run_server(endpoint::AbstractString)
    # Scripts default to exit-on-SIGINT, which skips the `finally` below and dumps
    # a backtrace on every ordinary Ctrl-C. Turn it into a catchable exception.
    Base.exit_on_sigint(false)

    warmup()

    ctx  = Context()
    sock = Socket(ctx, REP)
    ZMQ.bind(sock, endpoint)

    println("sheaf server  endpoint=$endpoint")
    println("  agents=$N_AGENTS  kinds=$(KINDS)  targets=$N_TARGETS")
    println("  consensus edges=$(length(CFG.edges))  pinning=$(length(CFG.pinning))")
    println("  k_cons=$K_CONS  k_track=$K_TRACK  z in [$Z_MIN, $Z_MAX]  arena_r=$ARENA_R")
    println("  synchronized_start=$(SYNC_START ? "on" : "OFF")")
    println("  cbf=$(CBF_ON ? "on" : "OFF")  d_min air-air=$D_AIR_AIR  air-ground=$D_AIR_GROUND " *
            "(horizontal)  ground-ground=$D_GROUND_GROUND  alpha=$CBF_ALPHA")
    if isempty(KEEP_OUT)
        println("  keep-out zones: none")
    else
        for z in KEEP_OUT
            shape = z.column ?
                "cylinder r=$(z.radius) z in [$(round(z.z_lo, digits=2)), $(round(z.z_hi, digits=2))] (padded, horizontal)" :
                "sphere r=$(z.radius) (R^3)"
            who = z.air && z.ground ? "all" : (z.air ? "air" : "ground")
            println("  keep-out '$(z.name)'  $shape  c=$(z.center)  applies=$who")
        end
    end
    flush(stdout)

    logio = isempty(LOG_PATH) ? nothing : open(LOG_PATH, "w")
    if logio !== nothing
        cols = ["t"; [string(AGENTS[i].name, "_", c) for i in 1:N_AGENTS for c in ("x", "y", "z")];
                     [string(TARGETS[k].name, "_", c) for k in 1:N_TARGETS for c in ("x", "y", "z")]]
        println(logio, join(vcat(cols, ["cons_residual", "track_residual"]), ","))
    end

    n = 0
    t_start = time()
    try
        while true
            state = parse_state(String(ZMQ.recv(sock)))
            if state.stop
                ZMQ.send(sock, format_control(ControlOutput(AgentControl[], TargetPose[])))
                println("Client requested shutdown after $n steps.")
                break
            end
            controls, targets, diag = sheaf_control(state)
            ZMQ.send(sock, format_control(ControlOutput(controls, targets, diag)))
            n += 1
            if logio !== nothing
                ps = [zeros(3) for _ in 1:N_AGENTS]
                for a in state.agents; ps[a.id+1] = a.p; end
                print(logio, round(state.t, digits = 3))
                for q in ps, c in q; print(logio, ",", round(c, digits = 4)); end
                for tp in targets, c in tp.p; print(logio, ",", round(c, digits = 4)); end
                println(logio, ",", round(diag["consensus_residual"], digits = 6),
                               ",", round(diag["tracking_residual"], digits = 6))
                n % 40 == 0 && flush(logio)
            end
            if n % 100 == 0
                rate = n / (time() - t_start)
                println("  step $n  t=$(round(state.t, digits=2))s  " *
                        "cons=$(round(diag["consensus_residual"], digits=4))  " *
                        "track=$(round(diag["tracking_residual"], digits=4))  " *
                        "rate=$(round(Int, rate)) Hz")
                flush(stdout)
            end
        end
    catch e
        e isa InterruptException || rethrow(e)
        println("Server stopped.")
    finally
        logio === nothing || close(logio)
        ZMQ.close(sock)
        ZMQ.close(ctx)
    end
end

if abspath(PROGRAM_FILE) == (@__FILE__)
    run_server(String(CFG.zmq_endpoint))
end
