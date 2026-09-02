"""Trot gait for the Unitree bodies, driven by how fast they are actually moving.

The quadrupeds are driven as unicycles: a DiffDrive plugin on two hidden wheels
carries the body and publishes the odometry the sheaf controller closes on. That
is deliberate and this node does not touch it. What it adds is the missing half
of the picture, which is that a legged robot sliding along on frozen legs reads
as a bug to anyone watching. So the legs are animated to match the motion the
wheels are already producing.

The one thing that has to be right is that the gait is a function of *distance*,
not of wall time. A fixed-frequency leg cycle looks like moonwalking as soon as
the body speed stops matching the stride the animation assumed. Here the phase
advances at f = duty * speed / stride, so one full cycle happens per stride
travelled, and a stopped robot simply stops with its feet down. Outside the
cadence band the stride gives instead of the frequency: above f_max the steps
lengthen, below f_min they shorten, and either way the foot still travels at the
body's own speed while it is on the ground.

Within a cycle each foot is placed rather than each joint being posed: the stance
foot is held at a fixed point on the ground and swept backwards through the body
frame at exactly the body's own speed (this is what keeps the feet from skating),
the swing foot is lifted and carried forward, and the thigh and calf angles come
out of the two-link inverse kinematics of the leg. Yaw rate enters as a per-leg
speed, v - w * y, so the outside legs of a turn take longer steps and the robot
can step in place while spinning.

Targets go out as gz.msgs.Double on each joint's own JointPositionController
topic, over gz-transport directly, the same route mocap_node.py takes. There is
no ROS message for "one joint position" that ros_gz would bridge without a
parameter_bridge process per joint.
"""

from __future__ import annotations

import json
import math
import os

import rclpy
from nav_msgs.msg import Odometry
from rclpy.node import Node
from rclpy.qos import qos_profile_sensor_data

# Leg name -> (phase offset, sign of the leg's y offset). A trot moves the
# diagonal pairs together: FR with RL, FL with RR, half a cycle apart.
LEGS = {
    "FR": (0.0, -1.0),
    "FL": (0.5, +1.0),
    "RR": (0.5, -1.0),
    "RL": (0.0, +1.0),
}

# Everything here is read off the same xacro the model geometry came from, except
# the three gait numbers at the end, which are a look choice.
#
#   l_thigh, l_calf  link lengths, so the inverse kinematics matches the model
#   hip_height       hip joint to foot centre in the standing pose; the body
#                    rides at a fixed height on the wheels, so this is constant
#   y_leg            lateral offset of the leg's sagittal plane from the centre
#                    line (leg_offset_y + thigh_offset), which is the moment arm
#                    that turns yaw rate into a per-leg speed
#   stride           distance covered per step at the nominal cadence
#   lift             how high the swing foot clears the ground
#   f_max            cadence ceiling; above it the stride lengthens instead, so
#                    the feet still do not skate
#   f_min            cadence floor; below it the stride shortens instead, which
#                    is what makes a robot rolling to a halt take smaller and
#                    smaller steps rather than one enormous slow one that has to
#                    be abandoned mid-stride when it finally stops
PLATFORMS = {
    "go1": dict(l_thigh=0.213, l_calf=0.213, hip_height=0.3202, y_leg=0.12675,
                stride=0.22, lift=0.07, f_max=3.0, f_min=0.5),
    "b1": dict(l_thigh=0.35, l_calf=0.35, hip_height=0.5263, y_leg=0.19875,
               stride=0.36, lift=0.11, f_max=2.2, f_min=0.4),
}

# Fraction of the cycle a foot spends on the ground. 0.5 is what makes a trot a
# trot: the two diagonal pairs hand off with neither pair ever both airborne.
DUTY = 0.5

# The joint controllers hold the last target they were sent, so this rate is the
# resolution of the animation and not a sampling of it, and it has to be set
# against the fastest the legs ever move rather than against a comfortable
# average. At 50 Hz a swinging thigh travels so far between targets that the leg
# rings: measured on the Go1 at 0.3 m/s, the foot swung 23 cm up instead of 9 and
# then punched 4.5 cm through the floor. 100 Hz is clean at a walk but still puts
# the foot a centimetre and a half under the ground at the Go1's 1.5 m/s ceiling;
# 400 Hz holds the planted foot to a tenth of a millimetre at every speed, for
# about a third of a core with twelve doubles per robot per tick.
RATE_HZ = 400.0

# Below this the robot counts as parked and the legs hold the neutral stance,
# rather than creeping through a cycle on odometry noise.
SPEED_DEADBAND = 0.02


def leg_ik(x: float, u: float, l1: float, l2: float) -> tuple[float, float]:
    """Thigh and calf angles that put the foot at (x forward, u below the hip).

    Sign convention is the xacro's: both joints turn about +y, so a positive
    thigh angle swings the leg backwards and the calf angle is negative (knees
    fold rearwards). Zero on both is the leg pointing straight down, which is
    also how the links are laid out in the model file."""
    reach = l1 + l2
    r = math.hypot(x, u)
    # Never ask for a straight or an inverted knee: at full extension the
    # inverse kinematics is singular and the leg snaps between solutions.
    r = min(max(r, 0.25 * reach), 0.98 * reach)
    scale = r / max(math.hypot(x, u), 1e-9)
    x, u = x * scale, u * scale

    cos_knee = (r * r - l1 * l1 - l2 * l2) / (2.0 * l1 * l2)
    knee = math.acos(min(max(cos_knee, -1.0), 1.0))
    # Angle of the hip-to-foot line off vertical, minus the angle the thigh sits
    # off that line once the knee is folded by `knee`.
    gamma = math.atan2(x, u)
    delta = math.atan2(l2 * math.sin(knee), l1 + l2 * math.cos(knee))
    return delta - gamma, -knee


def foot_target(phase: float, stride: float, hip_height: float, lift: float) -> tuple[float, float]:
    """Where this leg's foot should be, in the body frame, at cycle position `phase`.

    Stance is a straight sweep from the front of the step to the back at constant
    speed: that is the half that has to be exact, because any error in it is a
    foot sliding across the ground. Swing is a raised cosine forward, which only
    has to look unhurried."""
    phase %= 1.0
    if phase < DUTY:
        s = phase / DUTY
        return stride * (0.5 - s), hip_height
    s = (phase - DUTY) / (1.0 - DUTY)
    x = stride * (-0.5 + 0.5 * (1.0 - math.cos(math.pi * s)))
    return x, hip_height - lift * math.sin(math.pi * s)


class GaitDriver(Node):
    def __init__(self):
        super().__init__("gait_driver")
        self.declare_parameter("config", "")

        config_path = self.get_parameter("config").value or os.environ.get("SHEAF_CONFIG", "")
        if not config_path:
            raise RuntimeError("set the `config` parameter (or SHEAF_CONFIG)")
        with open(config_path) as fh:
            cfg = json.load(fh)

        self.robots = []
        for spec in cfg["agents"]:
            geom = PLATFORMS.get(spec.get("platform"))
            model = spec.get("gz_model") or spec.get("name")
            if geom is None or not model:
                continue
            self.robots.append(dict(model=model, geom=geom, phase=0.0, v=0.0, w=0.0,
                                    parked=False))

        if not self.robots:
            # The launch file already skips this node for a quadruped-free
            # config; this covers being run by hand against one.
            self.get_logger().info("gait: no b1/go1 agents in this config, nothing to animate")
            return

        from gz.transport13 import Node as GzNode
        from gz.msgs10.double_pb2 import Double

        self._Double = Double
        self._gz = GzNode()   # must outlive the publishers

        for r in self.robots:
            r["pubs"] = {}
            for leg in LEGS:
                for joint in ("hip", "thigh", "calf"):
                    name = f"{leg}_{joint}_joint"
                    topic = f"/model/{r['model']}/joint/{name}/0/cmd_pos"
                    r["pubs"][name] = self._gz.advertise(topic, Double)
            self.create_subscription(
                Odometry, f"/model/{r['model']}/odometry",
                lambda msg, rr=r: self._on_odom(msg, rr), qos_profile_sensor_data)
            self.get_logger().info(f"gait: animating {r['model']} from /model/{r['model']}/odometry")

        self._last = self.get_clock().now()
        self.create_timer(1.0 / RATE_HZ, self._tick)

    def _on_odom(self, msg: Odometry, robot: dict) -> None:
        # DiffDrive reports its twist in the body frame, which is already the
        # frame the gait needs: forward speed and yaw rate, no rotation applied.
        robot["v"] = msg.twist.twist.linear.x
        robot["w"] = msg.twist.twist.angular.z

    def _tick(self) -> None:
        now = self.get_clock().now()
        dt = (now - self._last).nanoseconds * 1e-9
        self._last = now
        # A paused or rewound sim clock must not teleport the gait.
        if not 0.0 < dt < 0.2:
            return

        for r in self.robots:
            g = r["geom"]
            v, w = r["v"], r["w"]
            # Cadence follows ground speed, with the widest leg's contribution
            # from yaw folded in so a robot spinning on the spot still steps.
            speed = abs(v) + abs(w) * g["y_leg"]
            if speed < SPEED_DEADBAND:
                self._park(r)
                continue
            r["parked"] = False
            freq = min(max(DUTY * speed / g["stride"], g["f_min"]), g["f_max"])
            r["phase"] = (r["phase"] + freq * dt) % 1.0

            for leg, (offset, y_sign) in LEGS.items():
                # Ground speed under this particular hip. The stride each leg
                # takes is whatever its own foot has to cover during stance, so
                # the inside legs of a turn take short steps and the outside
                # legs long ones, both without sliding.
                v_leg = v - w * (y_sign * g["y_leg"])
                stride = v_leg * DUTY / freq
                # Swing height in proportion to step length. Inside the cadence
                # band the stride is the nominal one and this is simply the full
                # lift; it only bites at a crawl, where a short step paired with
                # a full-height lift would read as prancing.
                lift = g["lift"] * min(1.0, abs(stride) / g["stride"])
                x, u = foot_target(r["phase"] + offset, stride, g["hip_height"], lift)
                thigh, calf = leg_ik(x, u, g["l_thigh"], g["l_calf"])
                self._send(r, leg, 0.0, thigh, calf)

    def _park(self, robot: dict) -> None:
        # The controllers latch the last target, so the neutral stance is sent
        # once and then left alone rather than restated four hundred times a
        # second at a robot that is standing still.
        if robot["parked"]:
            return
        robot["parked"] = True
        g = robot["geom"]
        thigh, calf = leg_ik(0.0, g["hip_height"], g["l_thigh"], g["l_calf"])
        for leg in LEGS:
            self._send(robot, leg, 0.0, thigh, calf)

    def _send(self, robot: dict, leg: str, hip: float, thigh: float, calf: float) -> None:
        for joint, value in (("hip", hip), ("thigh", thigh), ("calf", calf)):
            msg = self._Double()
            msg.data = value
            robot["pubs"][f"{leg}_{joint}_joint"].publish(msg)


def main(args=None) -> None:
    rclpy.init(args=args)
    node = GaitDriver()
    try:
        rclpy.spin(node)
    except KeyboardInterrupt:
        pass
    finally:
        node.destroy_node()
        rclpy.shutdown()


if __name__ == "__main__":
    main()
