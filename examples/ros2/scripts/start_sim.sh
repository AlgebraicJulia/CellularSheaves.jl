#!/usr/bin/env bash
# start_sim.sh : bring up Gazebo Harmonic, the ground robot and the PX4 SITL
# instances for the Autonomy Park sheaf demo.
#
# Gazebo and PX4 are started here as plain processes rather than from the ROS
# launch file, so that the controller can be restarted without restarting the
# simulator. The ROS 2 side (mavros, ros_gz_bridge, sheaf_bridge) comes up
# separately via scripts/start_bridge.sh.
#
#   ./start_sim.sh            # headless server plus GUI
#   HEADLESS=1 ./start_sim.sh # no GUI, for CI or a remote box
#
# Requires: gz-harmonic, a built PX4-Autopilot (PX4_DIR), ros_gz_sim.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$HERE")"
CONFIG="${SHEAF_CONFIG:-$ROOT/config/swarm.json}"
PX4_DIR="${PX4_DIR:-$HOME/PX4-Autopilot}"
WORLD_NAME="autonomy_park"
RUNDIR="${RUNDIR:-/tmp/sheaf_sim}"
HEADLESS="${HEADLESS:-0}"

mkdir -p "$RUNDIR"
: > "$RUNDIR/pids"

log() { printf '[start_sim] %s\n' "$*"; }
record() { echo "$1" >> "$RUNDIR/pids"; }

[[ -f "$CONFIG" ]]        || { echo "no config at $CONFIG" >&2; exit 1; }
[[ -d "$PX4_DIR" ]]       || { echo "no PX4-Autopilot at $PX4_DIR (set PX4_DIR)" >&2; exit 1; }
command -v gz >/dev/null  || { echo "gz not on PATH; install gz-harmonic" >&2; exit 1; }

# PX4 resolves PX4_GZ_WORLD against its own worlds directory, so the world has to
# live there. Copy rather than symlink: PX4's build tree is often a container
# layer that outlives the mount the repo was on.
PX4_WORLDS="$PX4_DIR/Tools/simulation/gz/worlds"
mkdir -p "$PX4_WORLDS"
cp "$ROOT/worlds/$WORLD_NAME.sdf" "$PX4_WORLDS/$WORLD_NAME.sdf"

# Real-time factor knob. PX4 lockstep follows the gz clock but gz does not wait
# for PX4, so on a loaded machine a PX4 process that loses the CPU for a second
# leaves its quadrotor flying open-loop on the last latched motor command --
# observed as a mid-flight attitude failure with a "timesync RTT too high"
# warning in the same breath. Running the world slower than real time buys the
# same margin lockstep would: every consumer downstream (bridge, brain, MAVROS)
# is sim-time driven, so behavior is identical, just slower on a wall clock.
RTF="${RTF:-1.0}"
sed -i "s|<real_time_factor>[0-9.]*</real_time_factor>|<real_time_factor>$RTF</real_time_factor>|" \
    "$PX4_WORLDS/$WORLD_NAME.sdf"
log "real-time factor $RTF"
log "cameras ${CAMERAS:-0}"

export GZ_SIM_RESOURCE_PATH="$ROOT/models:$PX4_DIR/Tools/simulation/gz/models:${GZ_SIM_RESOURCE_PATH:-}"
export GZ_SIM_SYSTEM_PLUGIN_PATH="$PX4_DIR/build/px4_sitl_default/build_gz:${GZ_SIM_SYSTEM_PLUGIN_PATH:-}"

log "world  $PX4_WORLDS/$WORLD_NAME.sdf"
log "config $CONFIG"

# RENDER_ENGINE: ogre2 is the default and wants a real GPU. Under software
# rendering (llvmpipe -- which is what you get when the container's Mesa is
# older than the host's GPU) ogre2 is unusably slow, while ogre1 stays
# interactive. Irrelevant headless: the server does not render at all.
GZ_ARGS=(-r -v 2)
[[ -n "${RENDER_ENGINE:-}" ]] && GZ_ARGS+=(--render-engine "$RENDER_ENGINE")
GZ_ARGS+=("$PX4_WORLDS/$WORLD_NAME.sdf")
[[ "$HEADLESS" == "1" ]] && GZ_ARGS=(-s "${GZ_ARGS[@]}")
gz sim "${GZ_ARGS[@]}" > "$RUNDIR/gz.log" 2>&1 &
record $!
log "gazebo pid $! (log: $RUNDIR/gz.log)"

# Wait for the world to actually exist before anything tries to spawn into it.
for _ in $(seq 1 60); do
  if gz topic -l 2>/dev/null | grep -q "/world/$WORLD_NAME"; then break; fi
  sleep 1
done
gz topic -l 2>/dev/null | grep -q "/world/$WORLD_NAME" \
  || { echo "gazebo never advertised /world/$WORLD_NAME; see $RUNDIR/gz.log" >&2; exit 1; }
log "world up"

# Ground robots: spawn the model, one process per agent, named so that the
# DiffDrive plugin's auto-scoped topics match what the adapters subscribe to.
python3 - "$CONFIG" <<'PY' > "$RUNDIR/ground.txt"
import json, sys
cfg = json.load(open(sys.argv[1]))
for a in cfg["agents"]:
    if a["kind"] == "diffdrive":
        x, y, z = a["spawn"]
        print(a.get("gz_model", a["name"]), x, y, z, a.get("platform", "jackal"))
PY

# The model is chosen by the agent's `platform`. A platform without a model in
# this tree still runs, with the Jackal's body and the config's own
# v_max/omega_max, and says so rather than failing at spawn.
while read -r name x y z platform; do
  [[ -z "${name:-}" ]] && continue
  sdf="$ROOT/models/apark_${platform:-jackal}/model.sdf"
  if [[ ! -f "$sdf" ]]; then
    [[ "${platform:-jackal}" != "jackal" ]] &&
      log "no model for platform '${platform}': spawning $name with the Jackal body (kinematics still come from the config)"
    sdf="$ROOT/models/apark_jackal/model.sdf"
  fi
  log "spawning ground robot $name at ($x, $y, $z)"
  ros2 run ros_gz_sim create -world "$WORLD_NAME" \
      -file "$sdf" -name "$name" -x "$x" -y "$y" -z "$z" >> "$RUNDIR/spawn.log" 2>&1
done < "$RUNDIR/ground.txt"

# Air robots: one PX4 SITL instance each, attaching to the running world.
#
# The variable names below are read straight out of this PX4 version's
# px4-rc.simulator, not from a tutorial: v1.15.4 selects the airframe with
# PX4_SIM_MODEL (the "gz_" prefix is stripped to get the model name), and
# PX4_GZ_MODEL_NAME means something else entirely (attach to an *existing*
# model rather than spawn one). There is no PX4_GZ_MODEL in this release.
#
# PX4_GZ_STANDALONE makes PX4 wait for the Gazebo we already started. Note that
# in standalone mode PX4 does not auto-detect the running world, so PX4_GZ_WORLD
# has to be set explicitly or gz_bridge attaches to nothing.
python3 - "$CONFIG" <<'PY' > "$RUNDIR/air.txt"
import json, sys
cfg = json.load(open(sys.argv[1]))
for i, a in enumerate(cfg["agents"]):
    if a["kind"] == "px4":
        x, y, z = a["spawn"]
        inst = a.get("instance", i)
        print(inst, a["name"], x, y, z, a.get("platform", "x500"),
              a.get("gz_model", f"x500_{inst}"))
PY

# Start one PX4 instance and wait for its own preflight to pass. gz-transport
# discovery races when several processes register against one world: an unlucky
# instance comes up with some sensor topics never connected ("Compass Sensor 0
# missing", "Battery unhealthy"), its EKF drifts on the ground for the whole
# session, and if it ever passes preflight late it joins the formation with a
# meters-wrong state estimate. PX4's "Ready for takeoff!" line is the ground
# truth that every sensor arrived and the EKF aligned, so an instance that does
# not print it in time is killed and started again rather than trusted.
#
# Killing a failed instance does not remove the model it spawned: Gazebo then
# refuses the retry's spawn ("Entity named [x500_2] already exists ... Entity not
# spawned"), the new PX4 attaches to nothing, and every retry fails the same way
# the first attempt did. Observed as three identical "never became ready" lines
# in a row. The stale model has to go before the next attempt.
remove_gz_model() {
  gz service -s "/world/$WORLD_NAME/remove" --reqtype gz.msgs.Entity \
     --reptype gz.msgs.Boolean --timeout 3000 \
     --req "name: \"$1\" type: MODEL" >/dev/null 2>&1 || true
}

# A platform with its own body in models/ is spawned here rather than by PX4,
# because PX4 names whatever it spawns "${PX4_SIM_MODEL#gz_}_${instance}" and the
# config's gz_model (and every topic derived from it) says x500_N. Spawning it
# ourselves is the only way to get the park's own body under the name the rest of
# the stack already agreed on: PX4_GZ_MODEL_NAME then makes gz_bridge attach to
# that model instead of spawning a second one.
#
# What is NOT changed is how it flies. Each apark_* air model is x500_base's
# links, inertias, rotor placement and motor plugins with only base_link's
# visuals replaced, and PX4_SYS_AUTOSTART stays 4001, so the airframe parameters
# still match the geometry. A camera instance keeps the x500_camN wrapper: the
# camera is the point of that run, and only the x500 body carries one.
start_px4_instance() {
  local inst="$1" name="$2" x="$3" y="$4" z="$5" platform="$6" gz_model="$7"
  local attempt pid sim_model body_sdf attach_name=""
  sim_model="$([ "$inst" -lt "${CAMERAS:-0}" ] && echo "gz_x500_cam$inst" || echo gz_x500)"
  body_sdf="$ROOT/models/apark_${platform:-x500}/model.sdf"
  if [[ "$inst" -ge "${CAMERAS:-0}" && -f "$body_sdf" ]]; then
    attach_name="$gz_model"
  elif [[ "${platform:-x500}" != "x500" && "$inst" -ge "${CAMERAS:-0}" ]]; then
    log "no model for platform '${platform}': flying $name with the x500 body (limits still come from the config)"
  fi
  for attempt in 1 2 3; do
    log "starting PX4 instance $inst ($name) at ($x, $y, $z), attempt $attempt"
    # clears a leftover from a previous attempt, or from a crashed earlier run
    remove_gz_model "${attach_name:-${sim_model#gz_}_$inst}"
    if [[ -n "$attach_name" ]]; then
      log "spawning $platform body for $name as $attach_name"
      ros2 run ros_gz_sim create -world "$WORLD_NAME" \
          -file "$body_sdf" -name "$attach_name" -x "$x" -y "$y" -z "$z" \
          >> "$RUNDIR/spawn.log" 2>&1
    fi
    (
      cd "$PX4_DIR"
      PX4_GZ_STANDALONE=1 \
      PX4_GZ_WORLD="$WORLD_NAME" \
      PX4_GZ_WORLDS="$PX4_WORLDS" \
      PX4_SYS_AUTOSTART=4001 \
      PX4_SIM_MODEL="$sim_model" \
      PX4_GZ_MODEL_NAME="$attach_name" \
      PX4_GZ_MODEL_POSE="$x,$y,$z,0,0,0" \
      ./build/px4_sitl_default/bin/px4 -i "$inst" -d \
        > "$RUNDIR/px4_$inst.log" 2>&1 &
      echo $! > "$RUNDIR/px4_$inst.pid"
    )
    pid="$(cat "$RUNDIR/px4_$inst.pid")"
    for _ in $(seq 1 45); do
      grep -q "Ready for takeoff" "$RUNDIR/px4_$inst.log" 2>/dev/null && break
      kill -0 "$pid" 2>/dev/null || break
      sleep 1
    done
    if grep -q "Ready for takeoff" "$RUNDIR/px4_$inst.log" 2>/dev/null; then
      echo "$pid" >> "$RUNDIR/pids"
      log "instance $inst ($name) preflight passed"
      return 0
    fi
    log "instance $inst ($name) never became ready; restarting it"
    kill "$pid" 2>/dev/null || true
    sleep 2
    kill -9 "$pid" 2>/dev/null || true
    sleep 1
  done
  echo "PX4 instance $inst failed preflight after 3 attempts; see $RUNDIR/px4_$inst.log" >&2
  return 1
}

while read -r inst name x y z platform gz_model; do
  [[ -z "${inst:-}" ]] && continue
  start_px4_instance "$inst" "$name" "$x" "$y" "$z" "$platform" "$gz_model"
done < "$RUNDIR/air.txt"

log "up. pids in $RUNDIR/pids, logs in $RUNDIR/"
log "next: ./start_bridge.sh   (and, separately, the Julia server)"
