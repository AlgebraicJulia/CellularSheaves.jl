# Autonomy Park ROS 2 / Gazebo demo.
#
# Run from the HOST:
#   make test        unit tests, no ROS/Julia/Docker needed          (seconds)
#   make loop        the real control loop, kinematic, no ROS        (~2 min)
#   make image       build the Humble + Harmonic + PX4 + Julia image (~30 min, once)
#   make run         full Gazebo simulation, headless, with a verdict (~8 min)
#   make shell       drop into the container to run it by hand / with a GUI
#
# Run INSIDE the container (make shell first):
#   make verify      check the image really has the park's stack
#   make sim         Gazebo + ground robot + PX4 instances
#   make bridge      mavros + ros_gz + sheaf_bridge      (second shell)
#   make brain       the Julia sheaf server              (third shell)
#   make smoke       headless end-to-end run with assertions
#   make down        stop whatever start_sim.sh started
#
# Knobs (see README): RTF=0.5 halves sim speed for CPU margin, SETTLE=N wall
# seconds to fly before scoring, MEM=8g the container RAM pin.

ROOT    := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
REPO    := $(abspath $(ROOT)../..)
PKG     := $(ROOT)ros2_ws/src/sheaf_ros2
CSV     ?= /tmp/sheaf_loop.csv
DURATION?= 90

IMAGE   ?= sheaf-apark:humble
RTF     ?= 0.5
SETTLE  ?= 300
MEM     ?= 8g

# Everything the container needs, in one place so `run` and `shell` cannot drift.
# The empty-manifest mount masks the repo's untracked root Manifest.toml, which
# dev-depends on sibling checkouts that do not exist in an image; see
# scripts/setup_julia.sh. Memory is pinned rather than merely limited because a
# swapped-out PX4 flies its quadrotor open-loop until it wakes.
# --ipc host and CycloneDDS are load-bearing, not preferences. FastDDS ships
# same-host data over shared memory, and each container gets a PRIVATE /dev/shm:
# two containers then discover each other's topics but exchange no data -- RViz
# shows an empty world with no error anywhere. Sharing the host IPC namespace
# fixes FastDDS; pinning CycloneDDS (pure UDP loopback, no SHM) removes the
# failure mode entirely and matches what the compose file already declares.
DOCKER_ARGS = --rm --network host --ipc host -e RMW_IMPLEMENTATION=rmw_cyclonedds_cpp \
  --memory $(MEM) --memory-swap $(MEM) \
  -v $(REPO):/ws/repo:rw \
  -v $(ROOT)docker/empty-manifest.toml:/ws/repo/Manifest.toml:ro \
  -v sheaf-julia-depot:/ws/julia_depot -v sheaf-julia-env:/ws/julia_env \
  -e JULIA_DEPOT_PATH=/ws/julia_depot -e JULIA_ENV_DIR=/ws/julia_env \
  -e RUNDIR=/ws/repo/examples/ros2/.run -e RTF=$(RTF) \
  -w /ws/repo/examples/ros2

.PHONY: test loop image run shell rviz verify smoke sim bridge brain down clean

test:
	cd $(PKG) && python3 -m pytest test -q

# Waits for the server's banner rather than sleeping a fixed time. A blind sleep
# hides the one failure that matters here: if another run already holds the ZMQ
# port, this server dies on bind and the harness cheerfully drives the OTHER
# run's controller, producing plausible numbers for the wrong experiment.
loop:
	@cd $(REPO) && julia --project=examples/ros2 examples/ros2/server.jl > /tmp/sheaf_loop.log 2>&1 & \
	 echo $$! > /tmp/sheaf_server.pid; \
	 for i in $$(seq 1 90); do grep -aq "sheaf server" /tmp/sheaf_loop.log 2>/dev/null && break; sleep 1; done; \
	 if ! grep -aq "sheaf server" /tmp/sheaf_loop.log; then \
	   echo "server.jl never started:"; tail -5 /tmp/sheaf_loop.log; \
	   kill $$(cat /tmp/sheaf_server.pid) 2>/dev/null; exit 1; \
	 fi; \
	 python3 $(ROOT)scripts/fake_sim.py --duration $(DURATION) --csv $(CSV); \
	 status=$$?; \
	 wait $$(cat /tmp/sheaf_server.pid) 2>/dev/null || true; \
	 exit $$status

image:
	cd $(ROOT) && docker build -f docker/Dockerfile -t $(IMAGE) .

# The headless one-shot: brings up Gazebo, PX4, the brain and the bridge inside
# the container, flies the formation, scores it and exits non-zero if it did not
# converge. Logs and the state CSV land in examples/ros2/.run/ on the host.
run:
	docker run $(DOCKER_ARGS) -e HEADLESS=1 -e SETTLE=$(SETTLE) $(IMAGE) ./scripts/smoke_test.sh

# X11 + GPU passthrough shared by the interactive targets. Run
# `xhost +local:docker` on the host once per login first.
#
# GPU=1 uses the NVIDIA runtime, which injects the driver's own GL userspace and
# makes the container's Mesa version irrelevant. It needs the driver loaded
# (`modprobe nvidia`, or prime-select on-demand) and nvidia-container-toolkit
# installed. Without it we pass /dev/dri and fall back to whatever the container's
# Mesa can do -- on a GPU newer than that Mesa, that means software llvmpipe, so
# RENDER_ENGINE defaults to ogre1 which stays interactive where ogre2 does not.
ifeq ($(GPU),1)
  GUI_ARGS = --gpus all -e NVIDIA_DRIVER_CAPABILITIES=all \
             -e __NV_PRIME_RENDER_OFFLOAD=1 -e __GLX_VENDOR_LIBRARY_NAME=nvidia
  RENDER_ENGINE ?= ogre2
  # CAMERAS is a count: the first N drones carry one. Default is a single
  # camera -- one feed tells the story, and each additional render competes
  # with the GUI and the physics on a machine that is already swapping.
  # CAMERAS=3 for all, CAMERAS=0 for none.
  CAMERAS ?= 1
else
  GUI_ARGS = --device /dev/dri --group-add $(shell getent group render | cut -d: -f3)
  RENDER_ENGINE ?= ogre
endif
CAMERAS ?= 0
X11_ARGS = -e DISPLAY=$(DISPLAY) -e QT_X11_NO_MITSHM=1 \
           -v /tmp/.X11-unix:/tmp/.X11-unix:rw -e RENDER_ENGINE=$(RENDER_ENGINE) \
           -e CAMERAS=$(CAMERAS)

shell:
	docker run $(DOCKER_ARGS) -it $(X11_ARGS) $(GUI_ARGS) $(IMAGE) bash

# One command to watch it fly: Gazebo with its GUI, the ground robot, the PX4
# fleet, the Julia brain and the ROS 2 bridge, all in one container. Ctrl-C to
# stop. Add GPU=1 once the NVIDIA driver and container toolkit are in place.
watch:
	docker run $(DOCKER_ARGS) -it $(X11_ARGS) $(GUI_ARGS) $(IMAGE) ./scripts/watch.sh

# The algorithm view: attaches to a RUNNING demo (make watch / make run in
# another shell) over the shared host-network DDS and draws the sheaf itself --
# active edges colored by residual, slots, targets, commanded velocities.
# Gazebo shows the robots; this shows what the controller thinks of them.
rviz:
	docker run $(DOCKER_ARGS) -it $(X11_ARGS) $(GUI_ARGS) $(IMAGE) bash -lc '\
	  set +u && source /opt/ros/humble/setup.bash && set -u && \
	  rviz2 -d /ws/repo/examples/ros2/config/sheaf.rviz --ros-args -p use_sim_time:=true'

verify:
	$(ROOT)scripts/verify_stack.sh

smoke:
	$(ROOT)scripts/smoke_test.sh

sim:
	$(ROOT)scripts/start_sim.sh

bridge:
	$(ROOT)scripts/start_bridge.sh

brain:
	cd $(REPO) && julia --project=examples/ros2 examples/ros2/server.jl

down:
	$(ROOT)scripts/stop_sim.sh

clean:
	rm -rf $(PKG)/.pytest_cache $(PKG)/**/__pycache__ /tmp/sheaf_sim
