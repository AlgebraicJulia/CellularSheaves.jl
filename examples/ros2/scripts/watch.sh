#!/usr/bin/env bash
# watch.sh : run the whole demo with the Gazebo GUI open, and hold until Ctrl-C.
#
# Same components as smoke_test.sh in the same order and for the same reasons
# (brain first so its JIT storm cannot starve a PX4 instance aligning its EKF),
# but it renders, it does not assert, and it does not stop. Use it to watch;
# use `make run` to get a verdict.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$HERE")"
RUNDIR="${RUNDIR:-/tmp/sheaf_sim}"
export SHEAF_CONFIG="${SHEAF_CONFIG:-$ROOT/config/swarm.json}"

# shellcheck source=scripts/lib.sh
source "$HERE/lib.sh"
export RTF="${RTF:-0.5}"
export RENDER_ENGINE="${RENDER_ENGINE:-ogre}"
ENDPOINT="$(python3 -c 'import json,sys;print(json.load(open(sys.argv[1]))["zmq_endpoint"])' "$SHEAF_CONFIG")"
require_free_zmq_port "$ENDPOINT" || exit 1
reset_rundir "$RUNDIR"

cleanup() {
  echo
  echo "[watch] shutting down"
  [[ -n "${BRIDGE_PID:-}" ]] && kill "$BRIDGE_PID" 2>/dev/null
  [[ -n "${BRAIN_PID:-}"  ]] && kill -INT "$BRAIN_PID" 2>/dev/null
  "$HERE/stop_sim.sh" >/dev/null 2>&1
}
trap cleanup EXIT INT TERM

# ROS 2's setup.bash reads unset variables, so -u stands down while sourcing.
set +u
source /opt/ros/humble/setup.bash
[[ -f /ws/ros2_ws/install/setup.bash ]] && source /ws/ros2_ws/install/setup.bash
set -u

echo "[watch] preparing the Julia environment"
"$HERE/setup_julia.sh" || exit 1

echo "[watch] starting the sheaf server"
export SHEAF_LOG="$RUNDIR/state.csv"
julia --heap-size-hint=1G --project="${JULIA_ENV_DIR:-/ws/julia_env}" "$ROOT/server.jl" \
    > "$RUNDIR/brain.log" 2>&1 &
BRAIN_PID=$!
for _ in $(seq 1 120); do
  grep -a -q "sheaf server" "$RUNDIR/brain.log" 2>/dev/null && break
  sleep 1
done
grep -a -q "sheaf server" "$RUNDIR/brain.log" || { echo "brain failed; see $RUNDIR/brain.log"; exit 1; }

echo "[watch] starting Gazebo (render engine: $RENDER_ENGINE, real-time factor: $RTF)"
"$HERE/start_sim.sh" || exit 1

echo "[watch] starting mavros + bridges"
"$HERE/start_bridge.sh" > "$RUNDIR/bridge.log" 2>&1 &
BRIDGE_PID=$!

cat <<'BANNER'

  ------------------------------------------------------------------
  Running. The Gazebo window should be open.

  Three quadrotors arm, climb to staggered altitudes and assemble a
  triangle around an orbiting target; the ground robot holds its slot
  beneath them, coupled through the x-y-projected restriction maps.
  Give it about a minute to arm and climb.

  Live numbers:  ros2 topic echo /sheaf_bridge/residuals
  Logs:          $RUNDIR/{brain,bridge,px4_*}.log
  Ctrl-C to stop.
  ------------------------------------------------------------------

BANNER

# Report convergence every 15 s so there is something to read while watching.
while sleep 15; do
  line="$(grep -a "step " "$RUNDIR/brain.log" | tail -1)"
  [[ -n "$line" ]] && echo "[watch] ${line#  }"
done
