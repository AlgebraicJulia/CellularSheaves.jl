import json
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from sheaf_ros2.protocol import (  # noqa: E402
    PROTOCOL_VERSION,
    AgentState,
    format_state,
    parse_control,
)


def test_format_state_shape():
    raw = format_state(0.5, [AgentState(0, "px4", [1.0, 2.0, 3.0], [0.1, 0.0, 0.0], 0.25, True)])
    d = json.loads(raw)
    assert d["version"] == PROTOCOL_VERSION
    assert d["t"] == 0.5
    assert d["agents"][0] == {
        "id": 0,
        "kind": "px4",
        "p": [1.0, 2.0, 3.0],
        "v": [0.1, 0.0, 0.0],
        "yaw": 0.25,
        "ok": True,
    }


def test_parse_control_reads_a_server_reply():
    """Byte-for-byte a reply emitted by examples/ros2/server.jl."""
    raw = (
        '{"version":2,'
        '"agents":[{"id":0,"u":[0.1,-0.2,0.05],"yaw_rate":0.0},'
        '{"id":1,"u":[0.0,0.0,0.0],"yaw_rate":0.0}],'
        '"targets":[{"id":0,"p":[1.4,0.0,1.5]}],'
        '"diag":{"consensus_residual":0.31,"tracking_residual":0.12}}'
    )
    out = parse_control(raw)
    assert [a.id for a in out.agents] == [0, 1]
    assert out.agents[0].u == [0.1, -0.2, 0.05]
    assert out.targets[0].p == [1.4, 0.0, 1.5]
    assert out.diag["consensus_residual"] == pytest.approx(0.31)


def test_parse_control_rejects_a_version_mismatch():
    with pytest.raises(ValueError, match="protocol mismatch"):
        parse_control('{"version":1,"agents":[]}')


def test_missing_optional_fields_default():
    out = parse_control('{"version":2,"agents":[{"id":3,"u":[0,0,0]}]}')
    assert out.agents[0].yaw_rate == 0.0
    assert out.targets == []
    assert out.diag == {}
