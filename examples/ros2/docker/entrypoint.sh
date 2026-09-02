#!/usr/bin/env bash
# Sources ROS and the demo workspace (building it on first run), then execs.
set -e

source /opt/ros/humble/setup.bash

WS=/ws/ros2_ws
REPO_PKG=/ws/repo/examples/ros2/ros2_ws/src/sheaf_ros2

if [[ -d "$REPO_PKG" ]]; then
  mkdir -p "$WS/src"
  # Symlink rather than copy so edits in the mounted repo take effect on rebuild.
  ln -sfn "$REPO_PKG" "$WS/src/sheaf_ros2"
  # Always rebuild: it takes about a second, and a build gated on "first run"
  # means a persistent workspace volume silently pins the console scripts and
  # data files of whatever revision existed when the volume was created. A new
  # node added to setup.py entry_points simply would not exist at launch time.
  echo "[entrypoint] building sheaf_ros2"
  ( cd "$WS" && colcon build --symlink-install --packages-select sheaf_ros2 )
  # shellcheck disable=SC1091
  source "$WS/install/setup.bash"
fi

exec "$@"
