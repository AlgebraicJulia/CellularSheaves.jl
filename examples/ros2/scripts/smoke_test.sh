#!/usr/bin/env bash
# smoke_test.sh : run the whole stack headless and assert it actually flew.
#
# Run inside the container:  ./scripts/smoke_test.sh
#
# Brings up Gazebo (no GUI), PX4, the Julia brain and the ROS 2 bridge, lets the
# swarm settle, then checks three things that a GUI would only let you check by
# squinting at it:
#
#   1. every quadrotor reached ARMED + OFFBOARD
#   2. the ground robot is publishing odometry
#   3. the consensus residual fell
#
# Tears everything down on any exit path, including failure.

set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$HERE")"
REPO="$(cd "$ROOT/../.." && pwd)"
RUNDIR="${RUNDIR:-/tmp/sheaf_sim}"
SETTLE="${SETTLE:-75}"          # wall seconds; sim seconds = SETTLE * RTF
export RTF="${RTF:-0.5}"        # see start_sim.sh: margin against CPU stalls
export SHEAF_CONFIG="${SHEAF_CONFIG:-$ROOT/config/swarm.json}"

# shellcheck source=scripts/lib.sh
source "$HERE/lib.sh"

ENDPOINT="$(python3 -c 'import json,sys;print(json.load(open(sys.argv[1]))["zmq_endpoint"])' "$SHEAF_CONFIG")"
require_free_zmq_port "$ENDPOINT" || exit 1
reset_rundir "$RUNDIR"
fails=0

# Memory pressure is a first-class failure mode here, not an inconvenience: a
# page-fault stall of a PX4 process leaves its quadrotor flying open-loop on the
# last latched motor command until the process wakes -- observed as random
# mid-flight attitude failures with "timesync RTT too high" in the same breath,
# one victim per run on a swapping host. Name it up front.
avail_kb="$(awk '/MemAvailable/{print $2}' /proc/meminfo)"
swap_used_kb="$(awk '/SwapTotal/{t=$2} /SwapFree/{f=$2} END{print t-f}' /proc/meminfo)"
printf '  info  MemAvailable %d MiB, swap in use %d MiB\n' "$((avail_kb/1024))" "$((swap_used_kb/1024))"
if (( avail_kb < 4000000 )); then
  printf '  warn  host is memory-pressured; PX4 stalls and random attitude failures are likely.\n'
  printf '        Run the containers with --memory/--memory-swap pinning, or free RAM first.\n'
fi
ok()  { printf '  \033[32mok\033[0m    %s\n' "$*"; }
bad() { printf '  \033[31mFAIL\033[0m  %s\n' "$*"; fails=$((fails+1)); }

cleanup() {
  echo "== teardown =="
  [[ -n "${BRIDGE_PID:-}" ]] && kill "$BRIDGE_PID" 2>/dev/null
  [[ -n "${BRAIN_PID:-}"  ]] && kill -INT "$BRAIN_PID" 2>/dev/null
  "$HERE/stop_sim.sh" >/dev/null 2>&1
}
trap cleanup EXIT

# ROS 2's setup.bash reads unset variables (AMENT_TRACE_SETUP_FILES and friends),
# so `set -u` has to stand down while it is sourced. Sourcing is the only place
# this applies; every check below still runs under -u.
set +u
source /opt/ros/humble/setup.bash
[[ -f /ws/ros2_ws/install/setup.bash ]] && source /ws/ros2_ws/install/setup.bash
set -u

# The brain comes up BEFORE the simulator, deliberately. Julia's JIT is a
# multi-core CPU storm for tens of seconds; when it overlapped a running sim it
# starved whichever PX4 instance was aligning its EKF at that moment, and that
# instance never passed preflight again for the whole run. The server's own
# warmup() runs before it binds the socket, so once "sheaf server" appears the
# storm is over and the simulator boots on a quiet machine.
echo "== julia environment =="
"$HERE/setup_julia.sh" || { echo "setup_julia.sh failed"; exit 1; }

echo "== brain =="
export SHEAF_LOG="$RUNDIR/state.csv"
( julia --heap-size-hint=1G --project="${JULIA_ENV_DIR:-/ws/julia_env}" "$ROOT/server.jl" > "$RUNDIR/brain.log" 2>&1 ) &
BRAIN_PID=$!
for _ in $(seq 1 120); do
  grep -a -q "sheaf server" "$RUNDIR/brain.log" 2>/dev/null && break
  sleep 1
done
grep -a -q "sheaf server" "$RUNDIR/brain.log" || { echo "brain never started; see $RUNDIR/brain.log"; exit 1; }
ok "julia server listening (JIT warm)"

echo "== simulator =="
HEADLESS=1 "$HERE/start_sim.sh" || { echo "start_sim.sh failed"; exit 1; }

echo "== bridge =="
"$HERE/start_bridge.sh" > "$RUNDIR/bridge.log" 2>&1 &
BRIDGE_PID=$!

echo "== settling for ${SETTLE}s =="
sleep "$SETTLE"

echo "== checks =="
python3 - "$SHEAF_CONFIG" <<'PY' > "$RUNDIR/names.txt"
import json, sys
cfg = json.load(open(sys.argv[1]))
for a in cfg["agents"]:
    print(a["kind"], a["name"], a.get("gz_model", a["name"]))
PY

while read -r kind name model; do
  case "$kind" in
    px4)
      st="$(timeout 10 ros2 topic echo --once "/$name/mavros/state" 2>/dev/null)"
      armed="$(grep -m1 '^armed:' <<<"$st" | awk '{print $2}')"
      mode="$(grep -m1 '^mode:'  <<<"$st" | awk '{print $2}')"
      if [[ "$armed" == "true" && "$mode" == "OFFBOARD" ]]; then
        ok "$name armed in OFFBOARD"
      else
        bad "$name armed=${armed:-?} mode=${mode:-?} (wanted true/OFFBOARD)"
      fi
      ;;
    diffdrive)
      if timeout 10 ros2 topic echo --once "/model/$model/odometry" >/dev/null 2>&1; then
        ok "$model publishing odometry"
      else
        bad "$model published no odometry"
      fi
      ;;
  esac
done < "$RUNDIR/names.txt"

# Convergence is judged by scripts/analyze_run.py over the trailing 20% of the
# state log -- residuals oscillate a little through the physical inner loop, so a
# single 20 Hz sample is noise; the tail mean is the honest measure. The peak ->
# last line from the brain log stays as context. These gates pass vacuously if an
# agent never joins (its edges are masked out of the residual); that failure mode
# belongs to the armed/odometry checks above, which is why both kinds exist.
peak="$(grep -a "step " "$RUNDIR/brain.log" | sed -n 's/.*cons=\([0-9.]*\).*/\1/p' | sort -g | tail -1)"
last="$(grep -a "step " "$RUNDIR/brain.log" | tail -1 | sed -n 's/.*cons=\([0-9.]*\).*/\1/p')"
echo "  info  consensus residual peak ${peak:-?} -> last ${last:-?}"
if python3 "$HERE/analyze_run.py" "$RUNDIR/state.csv" --config "$SHEAF_CONFIG" --max-cons 1.0 --max-track 1.0 --max-edge 0.6 --max-pin 0.6; then
  ok "tail-mean residuals converged (analyze_run.py)"
else
  bad "tail-mean residuals did not converge (see table above)"
fi

echo
if (( fails == 0 )); then echo "SMOKE TEST PASSED"; else echo "SMOKE TEST: $fails failure(s). Logs in $RUNDIR/"; fi
exit $(( fails > 0 ))
