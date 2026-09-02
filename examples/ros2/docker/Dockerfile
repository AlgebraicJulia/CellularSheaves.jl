# Autonomy Park sheaf demo: ROS 2 Humble + Gazebo Harmonic + PX4 SITL + Julia.
#
# The versions here are not preferences, they are the park's. Its wiki pins
# Ubuntu 22.04 (Jammy) and ROS 2 Humble, installs Gazebo via `gz-harmonic`, and
# runs PX4 SITL as `make px4_sitl gz_x500` with MAVROS and QGroundControl. Humble
# with Harmonic is a non-default pairing, so ros_gz comes from
# packages.osrfoundation.org rather than the ROS index; that is expected, and it
# is why this image exists instead of a one-line apt install.
#
#   cd examples/ros2 && docker build -f docker/Dockerfile -t sheaf-apark:humble .
#   (or just `make -C examples/ros2 image`)
#
# Expect a long first build: PX4 is compiled from source.

FROM ros:humble-ros-base-jammy

SHELL ["/bin/bash", "-lc"]
ENV DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 \
    PX4_DIR=/opt/PX4-Autopilot \
    GZ_VERSION=harmonic

ARG PX4_TAG=v1.15.4
ARG JULIA_VERSION=1.11.5

RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates curl wget gnupg lsb-release git sudo \
        python3-colcon-common-extensions \
        build-essential cmake ninja-build ccache \
        python3-pip python3-venv python3-jinja2 python3-empy \
        libgstreamer1.0-dev gstreamer1.0-plugins-bad gstreamer1.0-libav \
        x11-apps mesa-utils libgl1-mesa-dri \
    && rm -rf /var/lib/apt/lists/*

# Gazebo Harmonic and the Humble/Harmonic ros_gz build, both from OSRF.
RUN curl -fsSL https://packages.osrfoundation.org/gazebo.gpg \
        -o /usr/share/keyrings/pkgs-osrf-archive-keyring.gpg \
    && echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/pkgs-osrf-archive-keyring.gpg] http://packages.osrfoundation.org/gazebo/ubuntu-stable $(lsb_release -cs) main" \
        > /etc/apt/sources.list.d/gazebo-stable.list \
    && apt-get update && apt-get install -y --no-install-recommends \
        gz-harmonic libgz-sim8-dev \
        ros-humble-ros-gzharmonic \
    && rm -rf /var/lib/apt/lists/*

# MAVROS, plus the geographiclib geoid datasets it refuses to start without.
RUN apt-get update && apt-get install -y --no-install-recommends \
        ros-humble-mavros ros-humble-mavros-extras ros-humble-mavros-msgs \
        ros-humble-rmw-cyclonedds-cpp \
    && rm -rf /var/lib/apt/lists/* \
    && wget -q https://raw.githubusercontent.com/mavlink/mavros/master/mavros/scripts/install_geographiclib_datasets.sh \
        -O /tmp/install_geographiclib_datasets.sh \
    && bash /tmp/install_geographiclib_datasets.sh \
    && rm /tmp/install_geographiclib_datasets.sh

# PX4 SITL. --no-nuttx skips the cross toolchain we will never use; --no-sim-tools
# skips the simulators PX4 would otherwise pull in alongside the one we want.
RUN git clone --depth 1 --branch ${PX4_TAG} --recursive \
        https://github.com/PX4/PX4-Autopilot.git ${PX4_DIR} \
    && bash ${PX4_DIR}/Tools/setup/ubuntu.sh --no-nuttx --no-sim-tools \
    && rm -rf /var/lib/apt/lists/*

# Build only. `make px4_sitl gz_x500` selects an airframe *and launches it*, and
# DONT_RUN=1 is not honoured on that path, so it wedges the build at the pxh>
# prompt forever. `px4_sitl_default` is the plain build target; the gz_bridge
# module is part of that board config, and the airframe is chosen at run time via
# PX4_GZ_MODEL instead.
RUN cd ${PX4_DIR} && make px4_sitl_default && ccache -C

# Julia. Architecture is resolved at build time so the same Dockerfile produces a
# working image on an x86 workstation and on a Jetson Orin Nano, which is where
# this has to run once the controller moves onboard.
RUN set -e; \
    arch="$(dpkg --print-architecture)"; \
    case "$arch" in \
      amd64) jarch=x86_64; jdir=x64 ;; \
      arm64) jarch=aarch64; jdir=aarch64 ;; \
      *) echo "unsupported architecture: $arch" >&2; exit 1 ;; \
    esac; \
    short="${JULIA_VERSION%.*}"; \
    wget -q "https://julialang-s3.julialang.org/bin/linux/${jdir}/${short}/julia-${JULIA_VERSION}-linux-${jarch}.tar.gz" -O /tmp/julia.tgz; \
    mkdir -p /opt/julia && tar -xzf /tmp/julia.tgz -C /opt/julia --strip-components=1; \
    ln -sf /opt/julia/bin/julia /usr/local/bin/julia; \
    rm /tmp/julia.tgz

RUN pip3 install --no-cache-dir pyzmq==26.2.0

# RViz for the sheaf overlay (ros-base deliberately excludes GUI tools). A
# separate layer at the end so adding it never invalidates the PX4 build cache.
RUN apt-get update && apt-get install -y --no-install-recommends \
        ros-humble-rviz2 \
    && rm -rf /var/lib/apt/lists/*

ENV ROS_DOMAIN_ID=42 \
    GZ_SIM_RESOURCE_PATH=/ws/repo/examples/ros2/models \
    SHEAF_CONFIG=/ws/repo/examples/ros2/config/swarm.json \
    JULIA_DEPOT_PATH=/ws/julia_depot

WORKDIR /ws
COPY docker/entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["bash"]

# Real robot bodies for the ground platforms. Both layers sit at the very end so
# adding a platform never invalidates the PX4 build above.
#
# Clearpath ships the A200 (Husky) and J100 (Jackal) meshes in this one package,
# which is a plain apt install and needs no workspace build.
RUN apt-get update && apt-get install -y --no-install-recommends \
        ros-humble-clearpath-description \
    && rm -rf /var/lib/apt/lists/*

# Unitree publishes the B1 and Go1 meshes only as part of unitree_ros, which is
# ~1 GB of history and a dozen robots. A blobless sparse clone of just the two
# description packages we reference keeps this under 150 MB.
RUN git clone --depth 1 --filter=blob:none --sparse \
        https://github.com/unitreerobotics/unitree_ros.git /opt/unitree_ros \
    && cd /opt/unitree_ros \
    && git sparse-checkout set robots/b1_description/meshes robots/go1_description/meshes \
    && rm -rf /opt/unitree_ros/.git

# Airframe meshes. PX4 ships only the x500 body, so the other three park quads
# borrow real bodies from two long-standing open descriptions. Both clones are
# blobless and sparse down to the individual meshes referenced by the models,
# because rotors_description alone is over 100 MB of geometry we do not use.
RUN git clone --depth 1 --filter=blob:none --sparse \
        https://github.com/ethz-asl/rotors_simulator.git /tmp/rotors \
    && cd /tmp/rotors \
    && git sparse-checkout set --no-cone \
        /rotors_description/meshes/neo11.dae \
        /rotors_description/meshes/iris.dae \
    && mkdir -p /opt/rotors_description \
    && cp -r /tmp/rotors/rotors_description/meshes /opt/rotors_description/ \
    && rm -rf /tmp/rotors

RUN git clone --depth 1 --filter=blob:none --sparse \
        https://github.com/gsilano/BebopS.git /tmp/bebops \
    && cd /tmp/bebops \
    && git sparse-checkout set --no-cone /meshes/bebop.dae \
    && mkdir -p /opt/bebop_description \
    && cp -r /tmp/bebops/meshes /opt/bebop_description/ \
    && rm -rf /tmp/bebops
