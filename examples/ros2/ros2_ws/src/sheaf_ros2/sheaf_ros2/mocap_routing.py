"""swarm.json -> "which Motive rigid body feeds which ROS topics". Pure, no ROS.

Kept apart from the node so the part that is easy to get wrong -- and that no
amount of desk testing against real hardware can cover -- is the part that runs
under pytest. Two mistakes this file exists to make loud rather than silent:

  * a rigid body in Motive that no agent claims (harmless, but worth logging),
  * an agent whose body is not in the take at all (NOT harmless: that agent
    would sit with a never-updating pose, be masked as stale by the bridge, and
    the operator would be left guessing which of the two systems is at fault).

The second case raises with the full list of both sides, because at the park the
person reading it is standing in a net enclosure holding a controller.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence

# Motive streams at 100-240 Hz. MAVROS forwards vision_pose straight into PX4's
# EKF2, which fuses external vision at 30-50 Hz and gains nothing above that;
# pushing 240 Hz at it just burns the telemetry link the setpoints share. The
# node drops whole mocap frames rather than resampling, so the delivered rate
# quantizes DOWN to a divisor of the mocap rate: 30 Hz asked of a 100 Hz stream
# measures 25 Hz. That is the intended behaviour, not undershoot -- every
# forwarded pose is a real measurement and none is interpolated.
DEFAULT_VISION_RATE_HZ = 30.0

# A body Motive has lost keeps streaming its last pose with the tracking flag
# cleared. We drop those, so an occlusion shows up as silence -- and the bridge's
# own staleness check (`stale_after`) then masks the agent out of the sheaf. This
# is the age past which the node warns that it is happening.
DEFAULT_MAX_AGE_S = 0.25

DEFAULTS = {
    "source": "sim",
    "server": "127.0.0.1",
    "local": "0.0.0.0",
    "multicast": "239.255.42.99",
    "command_port": 1510,
    "data_port": 1511,
    "up_axis": "y",
    "vision_rate_hz": DEFAULT_VISION_RATE_HZ,
    "max_age_s": DEFAULT_MAX_AGE_S,
}


@dataclass(frozen=True)
class MocapRoute:
    """One agent's path from a Motive asset to the topics it must appear on."""

    index: int
    agent: str
    kind: str
    body: str                     # rigid body name as typed into Motive
    body_id: Optional[int]        # set only when the config pins a streaming id
    pose_topic: str               # what BaseAdapter subscribes to
    vision_topic: Optional[str]   # MAVROS vision_pose sink, PX4 agents only
    spawn: Sequence[float]


def mocap_settings(cfg: dict) -> dict:
    """The `mocap` block with defaults filled in. Absent block means simulation."""
    out = dict(DEFAULTS)
    out.update({k: v for k, v in cfg.get("mocap", {}).items() if not k.startswith("_")})
    unknown = set(out) - set(DEFAULTS)
    if unknown:
        raise ValueError(f"unknown keys in the config's `mocap` block: {sorted(unknown)}")
    out["command_port"] = int(out["command_port"])
    out["data_port"] = int(out["data_port"])
    out["vision_rate_hz"] = float(out["vision_rate_hz"])
    out["max_age_s"] = float(out["max_age_s"])
    if out["up_axis"] not in ("y", "z"):
        raise ValueError(f"mocap.up_axis must be 'y' or 'z', got {out['up_axis']!r}")
    if out["source"] not in ("sim", "natnet"):
        raise ValueError(f"mocap.source must be 'sim' or 'natnet', got {out['source']!r}")
    return out


def is_live(cfg: dict) -> bool:
    """True when this config asks for a real mocap system rather than Gazebo."""
    return mocap_settings(cfg)["source"] == "natnet"


def plan_routes(cfg: dict) -> List[MocapRoute]:
    """One route per agent that has a `pose_topic`.

    An agent without a `pose_topic` has opted out of mocap and runs on its own
    odometry, so it gets no route rather than a route nobody listens to.

    The Motive asset name comes from `mocap_body`, defaulting to the agent name.
    The default is the useful one: naming the rigid body after the agent is the
    thing an operator does anyway, and `mocap_body` exists for the take that was
    calibrated before anyone agreed on names.
    """
    routes: List[MocapRoute] = []
    seen: Dict[str, str] = {}
    for i, a in enumerate(cfg["agents"]):
        topic = a.get("pose_topic", "")
        if not topic:
            continue
        body = a.get("mocap_body") or a["name"]
        if body in seen:
            raise ValueError(
                f"agents {seen[body]!r} and {a['name']!r} both claim Motive rigid body "
                f"{body!r}; set `mocap_body` on one of them"
            )
        seen[body] = a["name"]
        vision = None
        if a["kind"] == "px4":
            # PX4 has to be TOLD where it is. The sheaf reading mocap does not
            # help the vehicle: EKF2 keeps integrating its own accelerometers and
            # flies the attitude and velocity loops off that estimate, so a
            # quadrotor with no vision fusion drifts underneath a controller that
            # can see it perfectly. MAVROS republishes this topic as
            # VISION_POSITION_ESTIMATE; EKF2_HGT_REF and EKF2_EV_CTRL must select
            # it on the airframe side or the messages arrive and are ignored.
            ns = a.get("ns", f"/{a['name']}")
            vision = f"{ns}/mavros/vision_pose/pose"
        routes.append(
            MocapRoute(
                index=i,
                agent=a["name"],
                kind=a["kind"],
                body=body,
                body_id=int(a["mocap_id"]) if a.get("mocap_id") is not None else None,
                pose_topic=topic,
                vision_topic=vision,
                spawn=[float(c) for c in a.get("spawn", (0.0, 0.0, 0.0))],
            )
        )
    return routes


def resolve_ids(
    routes: Sequence[MocapRoute], descriptions: Dict[int, str]
) -> Dict[int, MocapRoute]:
    """Streaming id -> route, matching on the names Motive reported.

    `mocap_id` in the config wins over discovery and does not have to appear in
    the descriptions at all, which is the escape hatch for a Motive setup whose
    model definitions this codec declines to parse.
    """
    by_name = {name: body_id for body_id, name in descriptions.items()}
    out: Dict[int, MocapRoute] = {}
    missing: List[str] = []
    for r in routes:
        body_id = r.body_id
        if body_id is None:
            body_id = by_name.get(r.body)
        if body_id is None:
            missing.append(f"{r.agent} -> {r.body!r}")
            continue
        if body_id in out:
            raise ValueError(
                f"agents {out[body_id].agent!r} and {r.agent!r} both resolve to Motive "
                f"streaming id {body_id}; check `mocap_id` in the config"
            )
        out[body_id] = r
    if missing:
        raise ValueError(
            "Motive is not streaming a rigid body for: "
            + ", ".join(missing)
            + ". Assets in the take: "
            + (", ".join(f"{v!r} (id {k})" for k, v in sorted(descriptions.items())) or "none")
            + ". Rename the asset in Motive, set `mocap_body` on the agent, or pin "
            "`mocap_id`."
        )
    return out


def unclaimed(routes: Sequence[MocapRoute], descriptions: Dict[int, str]) -> List[str]:
    """Names Motive streams that no agent asked for. Informational only: a take
    normally carries props and calibration squares nobody flies."""
    claimed = {r.body for r in routes}
    return sorted(name for name in descriptions.values() if name not in claimed)
