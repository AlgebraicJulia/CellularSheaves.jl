"""RViz overlay for the sheaf: the algorithm's state, drawn.

Gazebo renders bodies; nothing in it shows the mathematics being run against
them. This publishes a MarkerArray on /sheaf/viz with:

  * consensus edges  -- thick lines between agents, green -> red by residual
  * pinning edges    -- thinner lines from agent to target, same colormap
  * slots            -- translucent spheres where each agent is asked to be
  * targets          -- solid spheres at the controller's target positions
  * velocity arrows  -- the commanded u per agent
  * name labels      -- above each agent, gray when the agent is masked out

Everything comes from data the bridge already holds each tick: adapter
positions, the server reply (targets, per-edge residuals, active flags) and the
config (formation offsets, edge lists). Nothing here re-derives sheaf math; the
server's diag is the single source of truth, so when an agent is cut from the
active subcomplex its edges *visibly vanish* -- which is the point of drawing it.

All markers are stamped in the fixed frame "world" (the mocap frame), so no TF
tree is needed and the same overlay works in simulation and at the park.
"""

from __future__ import annotations

from typing import Dict, List, Sequence

from geometry_msgs.msg import Point
from std_msgs.msg import ColorRGBA
from visualization_msgs.msg import Marker, MarkerArray

# One residual worth of disagreement maps green -> red across this many metres.
RES_FULL_SCALE = 0.5


def residual_color(res: float, alpha: float = 1.0) -> ColorRGBA:
    """Green at zero disagreement, red at RES_FULL_SCALE metres and beyond."""
    f = max(0.0, min(1.0, res / RES_FULL_SCALE))
    return ColorRGBA(r=f, g=1.0 - f, b=0.08, a=alpha)


GRAY = ColorRGBA(r=0.45, g=0.45, b=0.45, a=0.9)


class SheafMarkers:
    def __init__(self, node, cfg: dict):
        self.node = node
        self.cfg = cfg
        self.names = [a["name"] for a in cfg["agents"]]
        self.offsets = [a["formation"] for a in cfg["agents"]]
        self.pub = node.create_publisher(MarkerArray, "/sheaf/viz", 10)

    def _marker(self, ns: str, mid: int, mtype: int) -> Marker:
        m = Marker()
        m.header.frame_id = "world"
        m.header.stamp = self.node.get_clock().now().to_msg()
        m.ns = ns
        m.id = mid
        m.type = mtype
        m.action = Marker.ADD
        m.pose.orientation.w = 1.0
        # A short lifetime instead of explicit DELETEs: if the bridge dies, the
        # overlay disappears within a second rather than freezing a stale lie.
        m.lifetime.nanosec = 800_000_000
        return m

    def publish(self, positions: Sequence[Sequence[float]],
                targets: Sequence[Sequence[float]],
                controls: Dict[int, Sequence[float]],
                diag: dict) -> None:
        arr = MarkerArray()
        active: List[bool] = list(diag.get("active", [True] * len(self.names)))

        def pt(p) -> Point:
            return Point(x=float(p[0]), y=float(p[1]), z=float(p[2]))

        # --- edges, straight from the server's diag -------------------------
        for n, e in enumerate(diag.get("edges", [])):
            kind = e.get("type", "cons")
            m = self._marker(f"edge_{kind}", n, Marker.LINE_LIST)
            m.scale.x = 0.035 if kind == "cons" else 0.018
            if not e.get("active", False):
                continue  # absent edge: not drawn, exactly like the Laplacian sees it
            a = int(e["a"])
            pa = positions[a]
            pb = targets[int(e["b"])] if kind == "pin" else positions[int(e["b"])]
            m.points = [pt(pa), pt(pb)]
            c = residual_color(float(e.get("res", 0.0)), alpha=1.0 if kind == "cons" else 0.75)
            m.colors = [c, c]
            arr.markers.append(m)

        # --- targets and slots ----------------------------------------------
        for k, tp in enumerate(targets):
            m = self._marker("target", k, Marker.SPHERE)
            m.pose.position = pt(tp)
            m.scale.x = m.scale.y = m.scale.z = 0.3
            m.color = ColorRGBA(r=0.95, g=0.12, b=0.10, a=1.0)
            arr.markers.append(m)

            for i, off in enumerate(self.offsets):
                s = self._marker("slot", k * len(self.names) + i, Marker.SPHERE)
                kinds = self.cfg["agents"][i]["kind"]
                if kinds == "diffdrive":
                    s.pose.position = Point(x=tp[0] + off[0], y=tp[1] + off[1], z=0.12)
                else:
                    s.pose.position = pt([tp[0] + off[0], tp[1] + off[1], tp[2] + off[2]])
                s.scale.x = s.scale.y = s.scale.z = 0.18
                s.color = ColorRGBA(r=0.25, g=0.55, b=0.95, a=0.35)
                arr.markers.append(s)

        # --- commanded velocities and labels --------------------------------
        for i, p in enumerate(positions):
            u = controls.get(i)
            if u is not None and any(abs(c) > 1e-3 for c in u):
                m = self._marker("cmd_vel", i, Marker.ARROW)
                m.points = [pt(p), pt([p[0] + u[0], p[1] + u[1], p[2] + u[2]])]
                m.scale.x, m.scale.y = 0.03, 0.07
                m.color = ColorRGBA(r=0.95, g=0.75, b=0.1, a=0.9) if (i < len(active) and active[i]) else GRAY
                arr.markers.append(m)

            t = self._marker("label", i, Marker.TEXT_VIEW_FACING)
            t.pose.position = pt([p[0], p[1], p[2] + 0.42])
            t.scale.z = 0.16
            t.text = self.names[i] if (i < len(active) and active[i]) else f"{self.names[i]} [masked]"
            t.color = ColorRGBA(r=1.0, g=1.0, b=1.0, a=0.95) if (i < len(active) and active[i]) else GRAY
            arr.markers.append(t)

        self.pub.publish(arr)
