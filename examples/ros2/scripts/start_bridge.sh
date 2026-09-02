#!/usr/bin/env bash
# start_bridge.sh : start the ROS 2 side (mavros, ros_gz bridges, sheaf_bridge).
# Safe to restart independently of the simulator.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$HERE")"
export SHEAF_CONFIG="${SHEAF_CONFIG:-$ROOT/config/swarm.json}"
exec ros2 launch sheaf_ros2 sim.launch.py config:="$SHEAF_CONFIG" "$@"
