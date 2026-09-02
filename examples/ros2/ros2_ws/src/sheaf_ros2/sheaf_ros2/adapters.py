"""Per-robot ROS 2 adapters.

Each adapter owns exactly one robot: it reads that robot's state into the
protocol's world-frame ENU representation, and writes a commanded world-frame
ENU velocity back out in whatever dialect the platform speaks. The sheaf server
sees only the protocol, so adding a platform means adding an adapter here and
nothing else.

Two are implemented:

  Px4Adapter        MAVROS offboard velocity setpoints. Runs the arm/OFFBOARD
                    handshake and an initial climb to `z_ref`, and only reports
                    `ready` once the robot is actually accepting commands.
  DiffDriveAdapter  Gazebo `DiffDrive` (and, unchanged, a real Jackal): converts
                    the planar velocity to (v, omega) through the near-identity
                    diffeomorphism in transforms.py.

Both accept an optional `pose_topic`. When set, position and yaw come from that
PoseStamped instead of odometry -- that is the hook for the park's OptiTrack
stream, and it is the only change needed to move from Gazebo odometry to mocap.
"""

from __future__ import annotations

from typing import Optional

from geometry_msgs.msg import PoseStamped, Twist
from nav_msgs.msg import Odometry
from rclpy.node import Node
from rclpy.qos import qos_profile_sensor_data

from .protocol import AgentState
from .transforms import quat_to_yaw, saturate_norm, world_velocity_to_unicycle


class BaseAdapter:
    """Common state bookkeeping. Subclasses fill in `send` and `ready`."""

    def __init__(self, node: Node, index: int, spec: dict):
        self.node = node
        self.index = index
        self.spec = spec
        self.name = spec["name"]
        self.kind = spec["kind"]
        self.v_max = float(spec["v_max"])

        self.p = [0.0, 0.0, 0.0]
        self.v = [0.0, 0.0, 0.0]
        self.yaw = 0.0
        self.last_state_stamp: Optional[float] = None
        self._pose_override = False

        # Simulator odometry is spawn-relative, not world-frame: PX4's EKF2 puts
        # the local origin wherever the vehicle booted, and Gazebo's DiffDrive
        # odometry starts at zero at the spawn pose. Consensus over raw odometry
        # is therefore consensus across N different frames -- the residual still
        # falls, but toward a world-frame formation distorted by the spawn
        # offsets. Adding the known spawn back recovers the world frame exactly
        # (all spawns have yaw 0; a nonzero spawn yaw would need the rotation
        # too, so refuse it rather than silently mis-rotate). A mocap pose_topic
        # is world-frame by construction and skips this correction entirely --
        # on hardware this offset is dead code, which is the point.
        self.origin = [float(x) for x in spec.get("spawn", (0.0, 0.0, 0.0))]
        if abs(float(spec.get("spawn_yaw", 0.0))) > 1e-9:
            raise ValueError(
                f"agent {self.name!r}: nonzero spawn yaw is not supported by the "
                "odometry frame correction; use a mocap pose_topic instead"
            )

        topic = spec.get("pose_topic", "")
        if topic:
            self._pose_override = True
            node.create_subscription(PoseStamped, topic, self._on_pose, qos_profile_sensor_data)

    # -- state in ----------------------------------------------------------
    def _on_odom(self, msg: Odometry) -> None:
        tw = msg.twist.twist.linear
        self.v = [tw.x, tw.y, tw.z]
        if self._pose_override:
            # Position comes from mocap; odometry only contributes velocity and
            # must NOT refresh the state stamp. Otherwise a live odometry stream
            # vouches for a dead mocap feed and the controller flies the whole
            # swarm on positions frozen at zero -- observed, not hypothetical.
            return
        pos = msg.pose.pose.position
        q = msg.pose.pose.orientation
        self.p = [
            pos.x + self.origin[0],
            pos.y + self.origin[1],
            pos.z + self.origin[2] - getattr(self, "_z_corr", 0.0),
        ]
        self.yaw = quat_to_yaw(q.x, q.y, q.z, q.w)
        self.last_state_stamp = self.node.get_clock().now().nanoseconds * 1e-9

    def _on_pose(self, msg: PoseStamped) -> None:
        pos = msg.pose.position
        q = msg.pose.orientation
        self.p = [pos.x, pos.y, pos.z]
        self.yaw = quat_to_yaw(q.x, q.y, q.z, q.w)
        self.last_state_stamp = self.node.get_clock().now().nanoseconds * 1e-9

    def state(self, stale_after: float) -> AgentState:
        now = self.node.get_clock().now().nanoseconds * 1e-9
        fresh = self.last_state_stamp is not None and (now - self.last_state_stamp) < stale_after
        return AgentState(
            id=self.index,
            kind=self.kind,
            p=list(self.p),
            v=list(self.v),
            yaw=self.yaw,
            ok=bool(fresh and self.ready),
        )

    # -- control out -------------------------------------------------------
    @property
    def ready(self) -> bool:
        raise NotImplementedError

    def send(self, u, yaw_rate: float) -> None:
        raise NotImplementedError

    def stop(self) -> None:
        self.send([0.0, 0.0, 0.0], 0.0)


class DiffDriveAdapter(BaseAdapter):
    """Gazebo DiffDrive / Clearpath Jackal."""

    def __init__(self, node: Node, index: int, spec: dict):
        super().__init__(node, index, spec)
        model = spec.get("gz_model", spec["name"])
        ns = spec.get("ns", f"/model/{model}")
        self.omega_max = float(spec.get("omega_max", 2.0))
        self.offset = float(spec.get("offset", 0.2))

        self.cmd_pub = node.create_publisher(Twist, f"{ns}/cmd_vel", 10)
        node.create_subscription(Odometry, f"{ns}/odometry", self._on_odom, qos_profile_sensor_data)

    @property
    def ready(self) -> bool:
        return self.last_state_stamp is not None

    def send(self, u, yaw_rate: float) -> None:
        v, w = world_velocity_to_unicycle(
            u[0], u[1], self.yaw, self.offset, self.v_max, self.omega_max
        )
        msg = Twist()
        msg.linear.x = v
        msg.angular.z = w
        self.cmd_pub.publish(msg)


class Px4Adapter(BaseAdapter):
    """PX4 multirotor over MAVROS.

    PX4 refuses OFFBOARD until setpoints are already streaming, so this adapter
    always publishes -- a zero velocity when it has nothing better to say -- and
    retries the mode change and the arming command on a slow cadence until they
    take. `ready` stays False until the vehicle is armed, in OFFBOARD, and has
    climbed to within `climb_tol` of `z_ref`, which keeps the sheaf from trying
    to run a formation through a takeoff.
    """

    MODE_OFFBOARD = "OFFBOARD"
    RETRY_PERIOD = 1.0
    # Pad calibration and sanity bounds, applied at the arm moment when the
    # vehicle is known to sit at its spawn pose. The SITL EKF's grounded
    # altitude wanders by metres under host scheduling jitter (observed -1.8 to
    # +3.4 across runs); since the sheaf commands velocity, not position, PX4's
    # internal bias matters only through the feedback WE read -- so zero it on
    # the pad like a pilot zeroes an altimeter, and refuse to fly only when the
    # estimate is wrong in ways a pad calibration cannot fix: a horizontal
    # position far from the spawn, or a vertical bias beyond any plausible
    # drift. A commanded "climb" computed against an uncorrected +1 m bias is a
    # commanded descent into the floor; that is the bug this replaces.
    GROUND_XY_TOL = 1.0
    GROUND_Z_INSANE = 3.0
    # After this many arm/disarm cycles the vehicle is left on the ground for
    # good; endlessly re-arming a robot that keeps failing is noise, not help.
    MAX_ARM_CYCLES = 5

    def __init__(self, node: Node, index: int, spec: dict):
        super().__init__(node, index, spec)
        from mavros_msgs.msg import PositionTarget, State
        from mavros_msgs.srv import CommandBool, SetMode

        self._PositionTarget = PositionTarget
        ns = spec.get("ns", f"/{spec['name']}")
        self.z_ref = float(spec.get("z_ref", 1.5))
        self.climb_tol = float(spec.get("climb_tol", 0.15))
        self.climb_speed = float(spec.get("climb_speed", 0.6))

        self.armed = False
        self.mode = ""
        self._climbed = False
        self._last_retry = 0.0
        self._ground_checked = False
        self._unfit = False
        self._arm_cycles = 0
        self._z_corr = 0.0

        self.sp_pub = node.create_publisher(PositionTarget, f"{ns}/mavros/setpoint_raw/local", 10)
        node.create_subscription(State, f"{ns}/mavros/state", self._on_mavros_state, qos_profile_sensor_data)
        node.create_subscription(
            Odometry, f"{ns}/mavros/local_position/odom", self._on_odom, qos_profile_sensor_data
        )
        self.arm_cli = node.create_client(CommandBool, f"{ns}/mavros/cmd/arming")
        self.mode_cli = node.create_client(SetMode, f"{ns}/mavros/set_mode")
        self._CommandBool = CommandBool
        self._SetMode = SetMode

    def _on_mavros_state(self, msg) -> None:
        was_armed = self.armed
        self.armed = bool(msg.armed)
        self.mode = str(msg.mode)
        if was_armed and not self.armed:
            # Any disarm (landing detector, failsafe, RC kill) invalidates the
            # climb: on re-arm the vehicle is back on the ground, and feeding it
            # lateral sheaf commands down there is how quadrotors tip over.
            # The ground check and fitness reset too -- the estimator re-inits
            # on the ground, so it gets a fresh chance, bounded by the cycle cap.
            self._climbed = False
            self._ground_checked = False
            self._unfit = False
            self._z_corr = 0.0
            self._arm_cycles += 1

    @property
    def offboard(self) -> bool:
        return self.mode == self.MODE_OFFBOARD

    @property
    def ready(self) -> bool:
        return self.armed and self.offboard and self._climbed and not self._unfit

    def _publish_velocity(self, u, yaw_rate: float) -> None:
        PT = self._PositionTarget
        msg = PT()
        msg.header.stamp = self.node.get_clock().now().to_msg()
        msg.coordinate_frame = PT.FRAME_LOCAL_NED  # MAVROS converts our ENU input
        msg.type_mask = (
            PT.IGNORE_PX | PT.IGNORE_PY | PT.IGNORE_PZ
            | PT.IGNORE_AFX | PT.IGNORE_AFY | PT.IGNORE_AFZ
            | PT.IGNORE_YAW
        )
        msg.velocity.x, msg.velocity.y, msg.velocity.z = (float(c) for c in u)
        msg.yaw_rate = float(yaw_rate)
        self.sp_pub.publish(msg)

    def _retry_handshake(self) -> None:
        now = self.node.get_clock().now().nanoseconds * 1e-9
        if now - self._last_retry < self.RETRY_PERIOD:
            return
        self._last_retry = now
        if not self.offboard and self.mode_cli.service_is_ready():
            req = self._SetMode.Request()
            req.custom_mode = self.MODE_OFFBOARD
            self.mode_cli.call_async(req)
        elif not self.armed and self.arm_cli.service_is_ready():
            req = self._CommandBool.Request()
            req.value = True
            self.arm_cli.call_async(req)

    def send(self, u, yaw_rate: float) -> None:
        """Stream a setpoint every tick regardless of readiness -- that stream is
        itself the precondition PX4 checks before allowing OFFBOARD."""
        if not (self.armed and self.offboard):
            if self._arm_cycles < self.MAX_ARM_CYCLES:
                self._retry_handshake()
            elif self._arm_cycles == self.MAX_ARM_CYCLES:
                self._arm_cycles += 1  # log the giving-up exactly once
                self.node.get_logger().error(
                    f"{self.name}: {self.MAX_ARM_CYCLES} arm cycles failed; leaving it grounded"
                )
            self._publish_velocity([0.0, 0.0, 0.0], 0.0)
            return

        if not self._ground_checked:
            self._ground_checked = True
            import math as _math
            xy_err = _math.hypot(self.p[0] - self.origin[0], self.p[1] - self.origin[1])
            z_bias = self.p[2] - self.origin[2]
            if xy_err > self.GROUND_XY_TOL or abs(z_bias) > self.GROUND_Z_INSANE:
                self._unfit = True
                self.node.get_logger().error(
                    f"{self.name}: grounded pose implausible (xy off {xy_err:.2f} m, "
                    f"z off {z_bias:.2f} m); not flying it"
                )
            elif self._pose_override:
                # Position comes from mocap, which is world-frame truth. There is
                # no estimator bias to null out, and nulling one would inject the
                # pad's own mounting height as a permanent altitude error.
                self._z_corr = 0.0
            else:
                self._z_corr = z_bias
                if abs(z_bias) > 0.2:
                    self.node.get_logger().warn(
                        f"{self.name}: zeroed altimeter on the pad (bias {z_bias:+.2f} m)"
                    )
        if self._unfit:
            self._publish_velocity([0.0, 0.0, 0.0], 0.0)
            return

        if not self._climbed:
            err = self.z_ref - self.p[2]
            if abs(err) < self.climb_tol:
                self._climbed = True
            else:
                vz = max(-self.climb_speed, min(self.climb_speed, err))
                self._publish_velocity([0.0, 0.0, vz], 0.0)
                return

        self._publish_velocity(saturate_norm(u, self.v_max), yaw_rate)


ADAPTERS = {"px4": Px4Adapter, "diffdrive": DiffDriveAdapter}


def make_adapter(node: Node, index: int, spec: dict) -> BaseAdapter:
    kind = spec["kind"]
    if kind not in ADAPTERS:
        raise ValueError(f"unknown agent kind {kind!r}; expected one of {sorted(ADAPTERS)}")
    return ADAPTERS[kind](node, index, spec)
