"""mocap_natnet -- OptiTrack/Motive as the world-frame truth for the real park.

This is the hardware twin of `mocap_node.py`. That node republishes Gazebo's
ground truth onto each agent's `pose_topic`; this one puts Motive's NatNet
stream on exactly the same topics, in the same frame, message type and units, so
`BaseAdapter` cannot tell which one is upstream. That interchangeability is the
whole design: the adapters, the bridge, the protocol and server.jl are identical
in the lab and at the park, and the only thing that changes is which of these two
nodes the launch file starts.

It publishes to two places, and both are necessary:

  `pose_topic`                     what the sheaf bridge reads. Covers every
                                   agent, quadrotors and ground robots alike.
  `<ns>/mavros/vision_pose/pose`   PX4 agents only. MAVROS turns this into
                                   VISION_POSITION_ESTIMATE and EKF2 fuses it.

Skipping the second is the mistake worth spelling out. A quadrotor's attitude and
velocity loops run on EKF2, not on anything we publish; without external vision
that estimate is dead reckoning from an IMU, and it walks. The sheaf would then
be servoing a perfectly observed position through an inner loop that believes it
is somewhere else, and the formation error at the park would be the estimator's
drift, not the controller's. See the failure ledger's "spawn-relative frames"
entry: this is the same bug one layer lower down.

MAVROS offers two sinks and this uses the vision one. `mocap_pose_estimate`
sends ATT_POS_MOCAP, which PX4 still accepts but treats as a legacy path;
`vision_pose` sends VISION_POSITION_ESTIMATE, which is what the current EKF2
external-vision parameters are written against. Both plugins load by default (a
launch of this stack prints "Plugin mocap_pose_estimate initialized" either way),
so picking one is a matter of which topic gets published, not which plugin runs.

Airframe side, which no code here can enforce: EKF2_EV_CTRL must enable
horizontal and vertical position fusion, EKF2_HGT_REF must select vision (not
baro) if the altitude is to come from mocap, and the mocap yaw only reaches the
estimator if yaw fusion is enabled too. Bring one vehicle up and confirm PX4's
local position tracks the Motive numbers before arming anything else.

Run it against the loopback Motive stand-in with no hardware at all:

    ros2 run sheaf_ros2 fake_natnet --config config/swarm.json &
    ros2 run sheaf_ros2 mocap_natnet --ros-args -p config:=config/swarm.json
"""

from __future__ import annotations

import json
import os

import rclpy
from geometry_msgs.msg import PoseStamped
from rclpy.node import Node
from rclpy.qos import qos_profile_sensor_data

from .mocap_routing import is_live, mocap_settings, plan_routes, resolve_ids, unclaimed
from .natnet import NatNetClient
from .transforms import natnet_to_enu


class MocapNatNet(Node):
    def __init__(self):
        super().__init__("mocap_natnet")
        self.declare_parameter("config", "")
        self.declare_parameter("frame_id", "world")

        config_path = self.get_parameter("config").value or os.environ.get("SHEAF_CONFIG", "")
        if not config_path:
            raise RuntimeError("set the `config` parameter (or SHEAF_CONFIG)")
        with open(config_path) as fh:
            cfg = json.load(fh)

        self.frame_id = str(self.get_parameter("frame_id").value)
        self.settings = mocap_settings(cfg)
        self.routes = plan_routes(cfg)
        self.client = None
        self._pubs = {}
        self._vision_pubs = {}
        self._last_vision = {}
        self._seen = set()
        self._untracked = set()
        self._by_id = {}
        self._last_frame_wall = None
        self._decode_errors = 0

        if not is_live(cfg):
            # Inert by construction. The launch file does not start this node for
            # a simulation config, but running it by hand against one must be a
            # no-op rather than a second publisher fighting mocap_node for the
            # same topics.
            self.get_logger().warn(
                "mocap.source is 'sim': this config wants Gazebo ground truth, so "
                "mocap_natnet is standing down. Set mocap.source to 'natnet' in "
                "swarm.json to talk to Motive."
            )
            return

        if not self.routes:
            raise RuntimeError(
                "mocap.source is 'natnet' but no agent has a `pose_topic`, so there is "
                "nowhere to put the poses"
            )

        for r in self.routes:
            self._pubs[r.index] = self.create_publisher(
                PoseStamped, r.pose_topic, qos_profile_sensor_data
            )
            if r.vision_topic:
                self._vision_pubs[r.index] = self.create_publisher(
                    PoseStamped, r.vision_topic, qos_profile_sensor_data
                )

        s = self.settings
        self.client = NatNetClient(
            server=s["server"],
            local=s["local"],
            command_port=s["command_port"],
            data_port=s["data_port"],
            multicast=s["multicast"],
            on_frame=self._on_frame,
            on_error=self._on_decode_error,
        )
        info = self.client.connect()
        self.get_logger().info(
            f"mocap: {info.name} at {s['server']}, NatNet "
            f"{'.'.join(str(c) for c in info.natnet_version)}, {s['up_axis']}-up"
        )
        descriptions = self.client.request_descriptions()
        self._by_id = resolve_ids(self.routes, descriptions)
        for body_id, r in sorted(self._by_id.items()):
            sink = f" (+ {r.vision_topic})" if r.vision_topic else ""
            self.get_logger().info(
                f"mocap: rigid body {r.body!r} (id {body_id}) -> {r.pose_topic}{sink}"
            )
        extra = unclaimed(self.routes, descriptions)
        if extra:
            self.get_logger().info(f"mocap: ignoring unclaimed rigid bodies {extra}")
        self.client.start()

        self._vision_period = 1.0 / max(s["vision_rate_hz"], 1e-3)
        self.create_timer(1.0, self._watch_stream)

    # -- stream in ---------------------------------------------------------
    def _on_frame(self, frame) -> None:
        """Called on the NatNet reader thread. rclpy publishers are safe here,
        which is the same arrangement mocap_node uses for the gz-transport
        thread, so the two nodes stay structurally identical."""
        now = self.get_clock().now()
        stamp = now.to_msg()
        wall = now.nanoseconds * 1e-9
        self._last_frame_wall = wall
        for body in frame.bodies:
            route = self._by_id.get(body.body_id)
            if route is None:
                continue
            if not body.tracking_valid:
                # Motive keeps streaming the last known pose of a body it has
                # lost. Publishing that would hand the sheaf a robot frozen in
                # mid-air that looks perfectly fresh. Staying silent instead lets
                # the bridge's staleness check mask the agent out, which is the
                # behaviour the ghost-agent ledger entry asks for.
                if body.body_id not in self._untracked:
                    self._untracked.add(body.body_id)
                    self.get_logger().warn(
                        f"mocap: lost tracking on {route.body!r} ({route.agent}); "
                        "holding its pose topic silent"
                    )
                continue
            if body.body_id in self._untracked:
                self._untracked.discard(body.body_id)
                self.get_logger().info(f"mocap: {route.body!r} tracking again")
            if body.body_id not in self._seen:
                self._seen.add(body.body_id)
                self.get_logger().info(f"mocap: first fix for {route.body!r}")

            p, q = natnet_to_enu(body.position, body.quaternion, self.settings["up_axis"])
            msg = PoseStamped()
            msg.header.stamp = stamp
            msg.header.frame_id = self.frame_id
            msg.pose.position.x, msg.pose.position.y, msg.pose.position.z = p
            (
                msg.pose.orientation.x,
                msg.pose.orientation.y,
                msg.pose.orientation.z,
                msg.pose.orientation.w,
            ) = q
            self._pubs[route.index].publish(msg)

            vision = self._vision_pubs.get(route.index)
            if vision is not None and (
                wall - self._last_vision.get(route.index, -1e9) >= self._vision_period
            ):
                self._last_vision[route.index] = wall
                vision.publish(msg)

    def _on_decode_error(self, exc: Exception) -> None:
        self._decode_errors += 1
        if self._decode_errors in (1, 10, 100) or self._decode_errors % 1000 == 0:
            self.get_logger().error(
                f"mocap: undecodable NatNet datagram (#{self._decode_errors}): {exc}"
            )

    def _watch_stream(self) -> None:
        """Motive going away is silence, not an error, and silence is exactly
        what a healthy but stationary volume looks like to a naive check. Report
        the gap so the operator knows which system stopped."""
        if self._last_frame_wall is None:
            self.get_logger().warn("mocap: no frames yet from Motive", throttle_duration_sec=5.0)
            return
        age = self.get_clock().now().nanoseconds * 1e-9 - self._last_frame_wall
        if age > self.settings["max_age_s"]:
            self.get_logger().warn(
                f"mocap: no frame for {age:.2f} s", throttle_duration_sec=5.0
            )

    def destroy_node(self) -> None:
        if self.client is not None:
            self.client.stop()
        super().destroy_node()


def main(args=None) -> None:
    rclpy.init(args=args)
    node = MocapNatNet()
    try:
        rclpy.spin(node)
    except KeyboardInterrupt:
        pass
    finally:
        node.destroy_node()
        rclpy.shutdown()


if __name__ == "__main__":
    main()
