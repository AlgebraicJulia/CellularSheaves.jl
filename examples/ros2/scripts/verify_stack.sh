#!/usr/bin/env bash
# verify_stack.sh : confirm the container really has the park's stack before you
# spend time debugging a simulation that was never going to start.
#
# Run inside the container:  ./scripts/verify_stack.sh
#
# Checks each piece independently and reports every failure rather than stopping
# at the first, because "gz is missing" and "PX4 never built" have completely
# different fixes and you want to learn about both in one pass.

set -uo pipefail

PX4_DIR="${PX4_DIR:-/opt/PX4-Autopilot}"
fails=0

ok()   { printf '  \033[32mok\033[0m    %s\n' "$*"; }
bad()  { printf '  \033[31mFAIL\033[0m  %s\n' "$*"; fails=$((fails+1)); }
note() { printf '        %s\n' "$*"; }

echo "== base =="
if [[ -f /etc/os-release ]] && grep -q 'VERSION_CODENAME=jammy' /etc/os-release; then
  ok "Ubuntu 22.04 jammy"
else
  bad "not jammy (the park pins Ubuntu 22.04.5)"
fi

echo "== ROS 2 =="
if [[ -f /opt/ros/humble/setup.bash ]]; then
# ROS 2's setup.bash reads unset variables (AMENT_TRACE_SETUP_FILES and friends),
# so `set -u` has to stand down while it is sourced. Sourcing is the only place
# this applies; every check below still runs under -u.
  set +u; source /opt/ros/humble/setup.bash; set -u
  ok "ROS 2 humble at /opt/ros/humble"
else
  bad "no /opt/ros/humble/setup.bash"
fi
for pkg in ros_gz_sim ros_gz_bridge mavros mavros_msgs; do
  if ros2 pkg prefix "$pkg" >/dev/null 2>&1; then ok "ros2 package $pkg"; else bad "missing ros2 package $pkg"; fi
done

echo "== Gazebo =="
if command -v gz >/dev/null 2>&1; then
  ver="$(gz sim --versions 2>/dev/null | head -1)"
  if [[ "$ver" == 8.* ]]; then ok "Gazebo Harmonic (gz-sim $ver)"; else bad "gz-sim $ver is not Harmonic (expected 8.x)"; fi
else
  bad "gz not on PATH"
fi

echo "== PX4 =="
if [[ -x "$PX4_DIR/build/px4_sitl_default/bin/px4" ]]; then
  ok "px4 sitl binary built"
else
  bad "no px4 binary at $PX4_DIR/build/px4_sitl_default/bin/px4"
fi
if compgen -G "$PX4_DIR/Tools/simulation/gz/models/x500*" >/dev/null; then
  ok "x500 airframe model present"
else
  bad "x500 model missing from $PX4_DIR/Tools/simulation/gz/models"
fi
if compgen -G "$PX4_DIR/build/px4_sitl_default/build_gz/*.so" >/dev/null 2>&1; then
  ok "gz system plugin built"
else
  note "no build_gz/*.so (PX4 may build it lazily on first run; not fatal)"
fi

echo "== Julia =="
if command -v julia >/dev/null 2>&1; then
  ok "julia $(julia --version | awk '{print $3}')"
  if julia -e 'using ZMQ' >/dev/null 2>&1; then
    ok "ZMQ.jl loads"
  else
    note "ZMQ.jl not instantiated yet (run the brain service once; it instantiates on start)"
  fi
else
  bad "julia not on PATH"
fi

echo "== python =="
if python3 -c 'import zmq' 2>/dev/null; then ok "pyzmq"; else bad "pyzmq missing"; fi

echo "== demo package =="
set +u; [[ -f /ws/ros2_ws/install/setup.bash ]] && source /ws/ros2_ws/install/setup.bash; set -u
if ros2 pkg prefix sheaf_ros2 >/dev/null 2>&1; then
  ok "sheaf_ros2 built and sourced"
else
  note "sheaf_ros2 not sourced; the entrypoint builds it on first run into /ws/ros2_ws"
fi

echo
if (( fails == 0 )); then
  echo "stack ok."
else
  echo "$fails check(s) failed."
fi
exit $(( fails > 0 ))
