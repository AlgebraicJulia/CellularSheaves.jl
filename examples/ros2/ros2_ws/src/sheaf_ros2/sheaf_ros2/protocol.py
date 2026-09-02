"""Python side of the sheaf/ROS 2 wire protocol (v2).

Mirrors examples/ros2/protocol.jl EXACTLY. Change the message format here AND in
protocol.jl together -- nowhere else.

Frames: everything is ENU, in the mocap / Gazebo world frame. Metres, rad, s.

This module is deliberately free of any ROS or numpy import so it can be unit
tested (and reused by scripts/fake_sim.py) without a ROS 2 installation.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from typing import Any, Dict, List, Sequence

PROTOCOL_VERSION = 2


@dataclass
class AgentState:
    id: int
    kind: str  # "px4" | "diffdrive"
    p: Sequence[float]  # position [x, y, z]
    v: Sequence[float]  # velocity [vx, vy, vz]
    yaw: float
    ok: bool  # False => state is stale/invalid, controller should hold this agent


@dataclass
class AgentControl:
    id: int
    u: List[float]  # commanded world-frame ENU velocity
    yaw_rate: float


@dataclass
class TargetPose:
    id: int
    p: List[float]


@dataclass
class ControlOutput:
    agents: List[AgentControl]
    targets: List[TargetPose] = field(default_factory=list)
    diag: Dict[str, Any] = field(default_factory=dict)


def format_state(t: float, agents: Sequence[AgentState], stop: bool = False) -> bytes:
    """Serialize the request sent from the ROS 2 bridge to the Julia server.

    `stop=True` asks the server to finish and exit; it replies once with an empty
    control set. Used so a client can shut the brain down without a signal, which
    keeps Ctrl-C semantics out of the normal teardown path.
    """
    return json.dumps(
        {
            "version": PROTOCOL_VERSION,
            "t": float(t),
            "stop": bool(stop),
            "agents": [
                {
                    "id": int(a.id),
                    "kind": str(a.kind),
                    "p": [float(x) for x in a.p],
                    "v": [float(x) for x in a.v],
                    "yaw": float(a.yaw),
                    "ok": bool(a.ok),
                }
                for a in agents
            ],
        }
    ).encode()


def parse_control(raw: bytes | str) -> ControlOutput:
    """Deserialize the reply from the Julia server."""
    d = json.loads(raw)
    version = d.get("version", PROTOCOL_VERSION)
    if version != PROTOCOL_VERSION:
        raise ValueError(f"protocol mismatch: got v{version}, expected v{PROTOCOL_VERSION}")
    return ControlOutput(
        agents=[
            AgentControl(
                id=int(a["id"]),
                u=[float(x) for x in a["u"]],
                yaw_rate=float(a.get("yaw_rate", 0.0)),
            )
            for a in d["agents"]
        ],
        targets=[TargetPose(id=int(t["id"]), p=[float(x) for x in t["p"]]) for t in d.get("targets", [])],
        diag=dict(d.get("diag", {})),
    )
