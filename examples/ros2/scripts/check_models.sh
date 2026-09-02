#!/usr/bin/env bash
# check_models.sh : prove every models/apark_* body parses, spawns and renders.
#
# Run INSIDE the container. Three separate questions get three separate answers,
# because passing one of them says nothing about the others:
#   1. does libsdformat parse the file at all
#   2. does every <uri> it references exist on this image
#   3. does gz actually create the entity, and does a rendering server load the
#      meshes without falling back to a missing-mesh warning
# Step 3 renders headlessly with a camera in the world: a server started with -s
# and no rendering sensor never touches a visual mesh, so a broken mesh URI would
# pass unnoticed.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$HERE")"
PX4_DIR="${PX4_DIR:-/opt/PX4-Autopilot}"
WORLD_NAME="model_check"
RUNDIR="${RUNDIR:-/tmp/model_check}"
MODELS=("${@:-}")
[[ -z "${MODELS[0]:-}" ]] && MODELS=(apark_jackal apark_husky apark_b1 apark_go1 \
                                    apark_freefly apark_homebrew apark_bebop)

mkdir -p "$RUNDIR"
export GZ_SIM_RESOURCE_PATH="$ROOT/models:$PX4_DIR/Tools/simulation/gz/models:${GZ_SIM_RESOURCE_PATH:-}"
# gz-transport discovers every gz on the same network, so without a partition of
# its own this check happily lists the models of whatever demo is already running
# and reports them as its own.
export GZ_PARTITION="model_check_$$"
fail=0

echo "== 1. sdf parse =="
for m in "${MODELS[@]}"; do
  if gz sdf -k "$ROOT/models/$m/model.sdf" 2>"$RUNDIR/$m.sdferr"; then
    echo "  ok    $m"
  else
    echo "  FAIL  $m"; cat "$RUNDIR/$m.sdferr"; fail=1
  fi
done

echo "== 2. mesh uris resolve =="
for m in "${MODELS[@]}"; do
  missing=0
  while read -r uri; do
    p="${uri#file://}"
    if [[ "$uri" == model://* ]]; then
      rel="${uri#model://}"
      p=""
      IFS=: read -ra dirs <<< "$GZ_SIM_RESOURCE_PATH"
      for d in "${dirs[@]}"; do [[ -f "$d/$rel" ]] && p="$d/$rel" && break; done
    fi
    if [[ -z "$p" || ! -f "$p" ]]; then echo "  MISSING $m -> $uri"; missing=1; fi
    # the gazebo-classic material scripts inherited from PX4's x500_base are not
    # files on this image and gz resolves them internally; they are not meshes
  done < <(grep -o '<uri>[^<]*</uri>' "$ROOT/models/$m/model.sdf" | sed 's|</\?uri>||g' \
           | grep -v 'materials/scripts')
  [[ $missing -eq 0 ]] && echo "  ok    $m" || fail=1
done

echo "== 3. spawn into a rendering server =="
cat > "$RUNDIR/$WORLD_NAME.sdf" <<'SDF'
<?xml version="1.0" ?>
<sdf version="1.9">
  <world name="model_check">
    <plugin filename="gz-sim-physics-system" name="gz::sim::systems::Physics"/>
    <plugin filename="gz-sim-user-commands-system" name="gz::sim::systems::UserCommands"/>
    <plugin filename="gz-sim-scene-broadcaster-system" name="gz::sim::systems::SceneBroadcaster"/>
    <plugin filename="gz-sim-sensors-system" name="gz::sim::systems::Sensors"/>
    <light type="directional" name="sun">
      <diffuse>1 1 1 1</diffuse><direction>-0.5 0.1 -0.9</direction>
    </light>
    <model name="ground_plane"><static>true</static><link name="link">
      <collision name="collision"><geometry><plane><normal>0 0 1</normal></plane></geometry></collision>
      <visual name="visual"><geometry><plane><normal>0 0 1</normal><size>100 100</size></plane></geometry></visual>
    </link></model>
    <model name="observer"><static>true</static>
      <link name="link"><pose>-8 0 3 0 0.3 0</pose>
        <sensor name="cam" type="camera">
          <camera><horizontal_fov>1.2</horizontal_fov>
            <image><width>320</width><height>240</height></image>
            <clip><near>0.1</near><far>100</far></clip></camera>
          <update_rate>2</update_rate><always_on>1</always_on><topic>check_cam</topic>
        </sensor>
      </link>
    </model>
  </world>
</sdf>
SDF

# ogre2, not ogre1: ogre1 opens an X display even under --headless-rendering,
# while ogre2 goes through EGL and runs on the container's software Mesa.
gz sim -s -r -v 3 --headless-rendering --render-engine "${RENDER_ENGINE:-ogre2}" \
    "$RUNDIR/$WORLD_NAME.sdf" > "$RUNDIR/gz.log" 2>&1 &
gzpid=$!
for _ in $(seq 1 60); do
  gz topic -l 2>/dev/null | grep -q "/world/$WORLD_NAME" && break
  sleep 1
done

i=0
for m in "${MODELS[@]}"; do
  ros2 run ros_gz_sim create -world "$WORLD_NAME" -file "$ROOT/models/$m/model.sdf" \
      -name "test_$m" -x $((i * 3)) -y 0 -z 0.5 > "$RUNDIR/$m.spawn" 2>&1
  echo "  create $m: $(grep -aiE 'success|error|fail' "$RUNDIR/$m.spawn" | tail -1)"
  i=$((i + 1))
done
sleep 10

echo "-- gz model --list"
gz model --list 2>&1 | sed 's/^/  /'
for m in "${MODELS[@]}"; do
  if gz model --list 2>/dev/null | grep -q "test_$m"; then
    echo "  spawned      $m"
  else
    echo "  NOT SPAWNED  $m"; fail=1
  fi
done

echo "-- the camera actually rendered a frame with all this in front of it"
# Minutes, not seconds, and that is not a hang: the renderer parses every mesh on
# one thread, the Go1 trunk alone is 64 MB of Collada, and with no GPU in the
# container this all goes through software Mesa. The first frame is the proof
# that all of it got through the loader.
#
# Written to a file rather than piped into head: the pipeline runs under
# pipefail, head closing early SIGPIPEs gz topic, and the whole check then
# reports a failure it did not have.
timeout "${RENDER_TIMEOUT:-900}" gz topic -e -t /check_cam -n 1 > "$RUNDIR/frame.txt" 2>/dev/null
if grep -q '^width:' "$RUNDIR/frame.txt"; then
  echo "  ok: /check_cam published"
else
  echo "  NO FRAME on /check_cam: the meshes above were never handed to a renderer"
  fail=1
fi

echo "-- mesh and sdf complaints in the server log"
if grep -Ei "unable to find file|failed to load|missing|error|\[Err\]" "$RUNDIR/gz.log" \
     | grep -v "OGRE EXCEPTION.*RenderSystem_GL" | sed 's/^/  /' | grep .; then
  fail=1
else
  echo "  none"
fi

echo "-- resting pose after settling. Every model puts its own links at the right"
echo "   height above its origin, so a ground robot dropped from z=0.5 comes to"
echo "   rest with its origin at 0: below that it has sunk, above it, it floats."
for m in "${MODELS[@]}"; do
  printf '  %-16s %s\n' "$m" \
    "$(gz model -m "test_$m" -p 2>/dev/null | grep -A3 -m1 'Pose' | tr '\n' ' ' | tr -s ' ')"
done

kill "$gzpid" 2>/dev/null
wait "$gzpid" 2>/dev/null
echo "== verdict: $([[ $fail -eq 0 ]] && echo PASS || echo FAIL) (log: $RUNDIR/gz.log)"
exit $fail
