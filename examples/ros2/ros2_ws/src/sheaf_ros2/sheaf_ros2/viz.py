"""Render the sheaf's targets in Gazebo.

The targets live only in the controller -- they are points in a cochain, not
bodies -- so without this the demo shows robots orbiting nothing. This spawns one
emissive sphere per target and moves it to the position the controller reports.

Implementation notes worth keeping:

* Spheres are *models* created through the server's `/world/<w>/create`, not GUI
  markers. The marker service is provided by the GUI process, so markers vanish
  in headless runs and cannot be checked by the smoke test; models exist on the
  server, show up in `/world/<w>/pose/info`, and are therefore verifiable.
* They are `<static>` and have no `<collision>`. A target is a mathematical
  object; it must not fall, and a quadrotor must not be able to hit it.
* Updates run on a daemon thread, never on the control tick. `set_pose` is a
  blocking service call, and a stalled Gazebo would otherwise stall the control
  loop -- visualization must not be able to hurt flight.
"""

from __future__ import annotations

import threading
import time
from typing import Dict, List, Sequence

_SDF = """<?xml version="1.0" ?>
<sdf version="1.9">
  <model name="{name}">
    <static>true</static>
    <link name="link">
      <visual name="visual">
        <geometry><sphere><radius>{radius}</radius></sphere></geometry>
        <material>
          <ambient>{r} {g} {b} 1</ambient>
          <diffuse>{r} {g} {b} 1</diffuse>
          <emissive>{er} {eg} {eb} 1</emissive>
        </material>
      </visual>
    </link>
  </model>
</sdf>"""


class TargetViz:
    """Spawns a sphere per target and keeps it on the reported position."""

    def __init__(self, node, targets: Sequence[dict], world: str, rate_hz: float = 15.0):
        self.node = node
        self.world = world
        self.rate = rate_hz
        self.names: List[str] = []
        self._latest: Dict[str, Sequence[float]] = {}
        self._lock = threading.Lock()
        self._alive = False

        try:
            from gz.transport13 import Node as GzNode
            from gz.msgs10.boolean_pb2 import Boolean
            from gz.msgs10.entity_factory_pb2 import EntityFactory
            from gz.msgs10.pose_pb2 import Pose
        except ImportError as exc:
            node.get_logger().warn(f"target visualization off (no gz python bindings: {exc})")
            return

        self._Pose, self._Boolean = Pose, Boolean
        self._gz = GzNode()

        for i, spec in enumerate(targets):
            name = spec.get("name", f"target{i}")
            colour = spec.get("color", [0.95, 0.12, 0.10])
            req = EntityFactory()
            req.sdf = _SDF.format(
                name=name, radius=spec.get("radius_viz", 0.26),
                r=colour[0], g=colour[1], b=colour[2],
                # Emissive at half intensity so the sphere reads as lit from
                # within against the dark ground without blowing out to white.
                er=colour[0] * 0.5, eg=colour[1] * 0.5, eb=colour[2] * 0.5,
            )
            # A drawn path has no centre; start its marker on the first waypoint
            # so it does not visibly jump in from the origin on the first tick.
            pts = spec.get("points") or []
            c = spec.get("center") or (list(pts[0]) if pts else [0.0, 0.0, 1.0])
            req.pose.position.x, req.pose.position.y, req.pose.position.z = c[0], c[1], c[2]
            req.pose.orientation.w = 1.0
            ok, rep = self._gz.request(f"/world/{world}/create", req, EntityFactory, Boolean, 5000)
            if ok and rep.data:
                node.get_logger().info(f"spawned target marker {name!r}")
            else:
                # Almost always "already exists" from a previous run in the same
                # world; the pose updates below will drive it either way.
                node.get_logger().warn(f"could not spawn {name!r} (may already exist)")
            self.names.append(name)
            self._latest[name] = list(c)

        self._alive = True
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()

    def update(self, positions: Sequence[Sequence[float]]) -> None:
        """Called from the control tick. Stores only; never blocks."""
        if not self._alive:
            return
        with self._lock:
            for name, p in zip(self.names, positions):
                self._latest[name] = list(p)

    def _run(self) -> None:
        period = 1.0 / self.rate
        while self._alive:
            with self._lock:
                snapshot = dict(self._latest)
            for name, p in snapshot.items():
                msg = self._Pose()
                msg.name = name
                msg.position.x, msg.position.y, msg.position.z = p[0], p[1], p[2]
                msg.orientation.w = 1.0
                try:
                    self._gz.request(f"/world/{self.world}/set_pose", msg,
                                     self._Pose, self._Boolean, 300)
                except Exception:  # noqa: BLE001 - visualization must never kill the node
                    pass
            time.sleep(period)

    def stop(self) -> None:
        self._alive = False
