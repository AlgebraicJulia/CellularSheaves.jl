"""sheaf_bridge -- ROS 2 node that runs the swarm against the Julia sheaf server.

One timer tick is one control step:

    gather every adapter's state  ->  ZMQ REQ to the Julia server  ->  dispatch
    the returned world-frame velocities back through the adapters.

The node is deliberately thin. It holds no model of the sheaf, no gains, and no
knowledge of which robots are neighbours -- all of that lives in server.jl and in
swarm.json. What it does own is safety: if the server does not answer within
`watchdog_ms`, every robot is commanded to zero and stays there until a reply
arrives. A brain that stops thinking must not leave the body coasting.

Deployment note: for the centralized phase this runs once, on the ground station.
For the decentralized phase the same node runs on each robot with a single-agent
config, and the agents exchange restricted stalk values over ROS 2 topics instead
of being gathered here. The adapter and protocol code is unchanged either way.
"""

from __future__ import annotations

import json
import os

import rclpy
import zmq
from rclpy.node import Node
from std_msgs.msg import Float32MultiArray

from .adapters import make_adapter
from .protocol import format_state, parse_control
from .rviz_markers import SheafMarkers
from .viz import TargetViz


class SheafBridge(Node):
    def __init__(self):
        super().__init__("sheaf_bridge")

        self.declare_parameter("config", "")
        self.declare_parameter("endpoint", "")
        self.declare_parameter("stale_after", 0.5)

        config_path = self.get_parameter("config").value or os.environ.get("SHEAF_CONFIG", "")
        if not config_path:
            raise RuntimeError("set the `config` parameter (or SHEAF_CONFIG) to swarm.json")
        with open(config_path) as fh:
            self.cfg = json.load(fh)

        self.stale_after = float(self.get_parameter("stale_after").value)
        self.watchdog_s = float(self.cfg.get("watchdog_ms", 250)) * 1e-3
        endpoint = self.get_parameter("endpoint").value or self.cfg["zmq_endpoint"]

        self.adapters = [make_adapter(self, i, spec) for i, spec in enumerate(self.cfg["agents"])]

        self.ctx = zmq.Context()
        self.sock = self._connect(endpoint)
        self.endpoint = endpoint

        self.diag_pub = self.create_publisher(Float32MultiArray, "~/residuals", 10)

        # Targets are cochain values, not bodies; nothing in Gazebo knows they
        # exist until this draws them.
        self.viz = None
        if self.cfg.get("visualize", False):
            self.viz = TargetViz(
                self,
                self.cfg.get("targets", []),
                os.environ.get("GZ_WORLD", "autonomy_park"),
            )

        # The RViz overlay is always on: publishing an unsubscribed MarkerArray
        # costs nothing, and it is the only view that shows the sheaf itself.
        self.sheaf_markers = SheafMarkers(self, self.cfg)

        self.t0 = self.get_clock().now().nanoseconds * 1e-9
        self.faults = 0
        period = 1.0 / float(self.cfg["rate_hz"])
        self.timer = self.create_timer(period, self.tick)

        self.get_logger().info(
            f"bridge up: {len(self.adapters)} agents -> {endpoint} at {self.cfg['rate_hz']} Hz"
        )

    def _connect(self, endpoint: str):
        """A REQ socket that has timed out is stuck in the wrong protocol state; the
        only clean recovery is to throw it away and open a new one."""
        sock = self.ctx.socket(zmq.REQ)
        sock.setsockopt(zmq.RCVTIMEO, int(self.watchdog_s * 1000))
        sock.setsockopt(zmq.SNDTIMEO, int(self.watchdog_s * 1000))
        sock.setsockopt(zmq.LINGER, 0)
        sock.connect(endpoint)
        return sock

    def _all_stop(self, reason: str) -> None:
        self.faults += 1
        if self.faults == 1 or self.faults % 20 == 0:
            self.get_logger().warn(f"all-stop ({self.faults}): {reason}")
        for a in self.adapters:
            a.stop()
        self.sock.close()
        self.sock = self._connect(self.endpoint)

    def tick(self) -> None:
        t = self.get_clock().now().nanoseconds * 1e-9 - self.t0
        states = [a.state(self.stale_after) for a in self.adapters]

        try:
            self.sock.send(format_state(t, states))
            out = parse_control(self.sock.recv())
        except zmq.Again:
            self._all_stop("no reply from sheaf server within watchdog")
            return
        except (ValueError, KeyError) as exc:
            self._all_stop(f"malformed reply: {exc}")
            return

        if self.faults:
            self.get_logger().info("sheaf server responding again")
            self.faults = 0

        by_id = {c.id: c for c in out.agents}
        for i, adapter in enumerate(self.adapters):
            cmd = by_id.get(i)
            if cmd is None:
                adapter.stop()
                continue
            adapter.send(cmd.u, cmd.yaw_rate)

        if self.viz is not None and out.targets:
            self.viz.update([t.p for t in out.targets])

        self.sheaf_markers.publish(
            positions=[a.p for a in states],
            targets=[t.p for t in out.targets],
            controls={c.id: c.u for c in out.agents},
            diag=out.diag,
        )

        msg = Float32MultiArray()
        msg.data = [
            float(out.diag.get("consensus_residual", 0.0)),
            float(out.diag.get("tracking_residual", 0.0)),
        ]
        self.diag_pub.publish(msg)

    def destroy_node(self) -> bool:
        if self.viz is not None:
            self.viz.stop()
        for a in self.adapters:
            a.stop()
        self.sock.close()
        self.ctx.term()
        return super().destroy_node()


def main(args=None) -> None:
    rclpy.init(args=args)
    node = SheafBridge()
    try:
        rclpy.spin(node)
    except KeyboardInterrupt:
        pass
    finally:
        node.destroy_node()
        rclpy.shutdown()


if __name__ == "__main__":
    main()
