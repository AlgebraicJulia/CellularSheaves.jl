"""Gazebo ground truth as a stand-in for the Autonomy Park's OptiTrack.

Subscribes to Gazebo's `/world/<world>/dynamic_pose/info` over **gz-transport
directly** and republishes one `PoseStamped` per robot on the topic its adapter
watches. Do not "simplify" this to go through ros_gz_bridge: its Pose_V ->
TFMessage conversion drops every entity name (24 transforms, all with
`child_frame_id: ''`, verified against this container's ros-humble-ros-gzharmonic),
which silently unmatches every robot and leaves the swarm flying on stale zeros.
The raw gz message keeps the names.

Why mocap at all: without it each robot's position comes from its own estimator --
PX4's EKF2 for the quadrotors, wheel dead-reckoning for the rover. The rover
cannot even see itself being dragged; each quadrotor's frame is anchored wherever
it booted and drifts. Consensus over drifting, differently-anchored estimates is
consensus over a lie each robot tells itself, which is precisely why the park
runs OptiTrack. Swapping this node for the real mocap driver is a launch change:
point each agent's `pose_topic` at the driver's output instead of ours.
"""

from __future__ import annotations

import json
import os

import rclpy
from geometry_msgs.msg import PoseStamped
from rclpy.node import Node
from rclpy.qos import qos_profile_sensor_data


class MocapBridge(Node):
    def __init__(self):
        super().__init__("mocap_bridge")
        self.declare_parameter("config", "")
        self.declare_parameter("world", os.environ.get("GZ_WORLD", "autonomy_park"))

        config_path = self.get_parameter("config").value or os.environ.get("SHEAF_CONFIG", "")
        if not config_path:
            raise RuntimeError("set the `config` parameter (or SHEAF_CONFIG)")
        with open(config_path) as fh:
            cfg = json.load(fh)
        world = self.get_parameter("world").value

        # Gazebo names a PX4-spawned quadrotor after its airframe and instance
        # (x500_0 ...), not after our agent name, so the mapping is explicit in
        # the config. The pose message also carries every child link (base_link,
        # wheels); exact-matching on model names filters those out for free.
        self.pubs = {}
        self.spawns = {}
        for a in cfg["agents"]:
            model = a.get("gz_model")
            if not model:
                continue
            topic = a.get("pose_topic") or f"/mocap/{model}/pose"
            pub = self.create_publisher(PoseStamped, topic, qos_profile_sensor_data)
            # Aliases cover alternate physical models for the same agent (the
            # camera-equipped x500 spawns under a different gz name); whichever
            # name shows up in the world feeds the same pose topic.
            for name in [model] + list(a.get("gz_aliases", [])):
                self.pubs[name] = pub
                self.spawns[name] = list(a.get("spawn", [0.0, 0.0, 0.0]))
            self.get_logger().info(f"mocap: gazebo model {model!r} (+{len(a.get('gz_aliases', []))} aliases) -> {topic}")

        from gz.transport13 import Node as GzNode
        from gz.msgs10.boolean_pb2 import Boolean
        from gz.msgs10.pose_pb2 import Pose
        from gz.msgs10.pose_v_pb2 import Pose_V

        self._GzPose, self._GzBool = Pose, Boolean
        self._world = world
        self._rescued_at = {}

        self.seen = set()
        self._gz = GzNode()   # must outlive the subscription
        topic = f"/world/{world}/dynamic_pose/info"
        if not self._gz.subscribe(Pose_V, topic, self._on_poses):
            raise RuntimeError(f"could not subscribe to {topic} over gz-transport")
        self.get_logger().info(f"mocap: subscribed to {topic} (gz-transport, names intact)")

    def _on_poses(self, msg) -> None:
        # Runs on the gz-transport thread; rclpy publishers are safe to call here.
        stamp = self.get_clock().now().to_msg()
        for pose in msg.pose:
            pub = self.pubs.get(pose.name)
            if pub is None:
                continue
            if pose.name not in self.seen:
                self.seen.add(pose.name)
                self.get_logger().info(f"mocap: first fix for {pose.name}")
            if pose.position.z < -0.7:
                self._rescue(pose.name)
                continue   # do not publish a subterranean pose as truth
            out = PoseStamped()
            out.header.stamp = stamp
            out.header.frame_id = "world"
            out.pose.position.x = pose.position.x
            out.pose.position.y = pose.position.y
            out.pose.position.z = pose.position.z
            out.pose.orientation.x = pose.orientation.x
            out.pose.orientation.y = pose.orientation.y
            out.pose.orientation.z = pose.orientation.z
            out.pose.orientation.w = pose.orientation.w
            pub.publish(out)


    def _rescue(self, name: str) -> None:
        """Sim referee: teleport a fallen model back to its pad.

        The GUI translate tool can drop a model below the floor, after which it
        falls forever while its agent gets masked and the demo quietly loses a
        robot (observed at z = -2800). Real mocap would report the same garbage
        for a robot carried out of the volume; the sim can do what a human at
        the park does -- put the robot back on its pad. Cooldown guards against
        re-teleporting every tick while it settles."""
        import time
        now = time.monotonic()
        if now - self._rescued_at.get(name, -10.0) < 5.0:
            return
        self._rescued_at[name] = now
        spawn = self.spawns[name]
        self.get_logger().warn(f"referee: {name} fell out of the world; returning it to its pad")
        msg = self._GzPose()
        msg.name = name
        msg.position.x, msg.position.y = spawn[0], spawn[1]
        msg.position.z = max(spawn[2], 0.15)
        msg.orientation.w = 1.0
        try:
            self._gz.request(f"/world/{self._world}/set_pose", msg,
                             self._GzPose, self._GzBool, 500)
        except Exception:  # noqa: BLE001 - the referee must never kill mocap
            pass


def main(args=None) -> None:
    rclpy.init(args=args)
    node = MocapBridge()
    try:
        rclpy.spin(node)
    except KeyboardInterrupt:
        pass
    finally:
        node.destroy_node()
        rclpy.shutdown()


if __name__ == "__main__":
    main()
