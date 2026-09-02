"""NatNet codec, frame conversion and Motive-to-agent routing.

No ROS and no OptiTrack: everything here is bytes and dictionaries, plus one
end-to-end pass over real UDP sockets against the loopback Motive stand-in. What
is deliberately NOT covered is whether a real Motive agrees with this codec --
that is the one claim only the park can settle.
"""

import math
import os
import socket
import sys
import time

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from sheaf_ros2 import natnet  # noqa: E402
from sheaf_ros2.mocap_routing import (  # noqa: E402
    is_live,
    mocap_settings,
    plan_routes,
    resolve_ids,
    unclaimed,
)
from sheaf_ros2.transforms import natnet_to_enu, quat_to_yaw  # noqa: E402


# -- frame conversion ------------------------------------------------------


def test_yup_position_maps_east_north_up():
    """Motive Y-up (x right, y up, z toward front) -> ENU is (x, -z, y)."""
    p, _ = natnet_to_enu((1.0, 2.0, 3.0), (0.0, 0.0, 0.0, 1.0), "y")
    assert p == [1.0, -3.0, 2.0]


def test_zup_is_the_identity():
    p, q = natnet_to_enu((1.0, 2.0, 3.0), (0.1, 0.2, 0.3, 0.9), "z")
    assert p == [1.0, 2.0, 3.0]
    assert q == [0.1, 0.2, 0.3, 0.9]


def test_yup_rotation_about_up_becomes_enu_yaw():
    """A heading change in Motive must read as the same yaw in ENU, or the
    unicycle transform steers the rover into the netting."""
    theta = 0.7
    h = theta / 2.0
    _, q = natnet_to_enu((0.0, 0.0, 0.0), (0.0, math.sin(h), 0.0, math.cos(h)), "y")
    assert quat_to_yaw(*q) == pytest.approx(theta, abs=1e-9)


def test_yup_conversion_preserves_handedness():
    """The basis change must have determinant +1: a reflection would leave
    positions right and silently mirror every rotation."""
    e = [natnet_to_enu(v, (0.0, 0.0, 0.0, 1.0), "y")[0] for v in
         ((1, 0, 0), (0, 1, 0), (0, 0, 1))]
    cx = [e[0][1] * e[1][2] - e[0][2] * e[1][1],
          e[0][2] * e[1][0] - e[0][0] * e[1][2],
          e[0][0] * e[1][1] - e[0][1] * e[1][0]]
    assert cx == e[2]


def test_yup_conversion_is_isometric():
    p, _ = natnet_to_enu((0.3, -1.2, 2.5), (0.0, 0.0, 0.0, 1.0), "y")
    assert math.sqrt(sum(c * c for c in p)) == pytest.approx(
        math.sqrt(0.3 ** 2 + 1.2 ** 2 + 2.5 ** 2)
    )


def test_unknown_up_axis_is_refused():
    with pytest.raises(ValueError, match="up_axis"):
        natnet_to_enu((0, 0, 0), (0, 0, 0, 1), "x")


# -- codec -----------------------------------------------------------------


def test_message_framing_round_trip():
    raw = natnet.pack_message(natnet.NAT_FRAMEOFDATA, b"payload")
    msg_id, payload = natnet.unpack_message(raw)
    assert msg_id == natnet.NAT_FRAMEOFDATA
    assert payload == b"payload"


def test_connect_command_carries_the_ping_payload():
    """Motive matches on the literal "Ping"; an empty Connect is ignored."""
    assert natnet.pack_command(natnet.NAT_CONNECT) == b"\x00\x00\x05\x00Ping\x00"


def test_server_info_round_trip():
    info = natnet.ServerInfo("Motive", (3, 0, 1, 0), (4, 0, 0, 0))
    out = natnet.unpack_server_info(natnet.pack_server_info(info))
    assert out.name == "Motive"
    assert out.natnet_version == (4, 0, 0, 0)
    assert out.major == 4


def test_short_server_info_is_refused():
    with pytest.raises(natnet.NatNetError):
        natnet.unpack_server_info(b"\x00" * 100)


def test_model_def_round_trip_gives_id_to_name():
    payload = natnet.pack_model_def([(1, "uav0"), (7, "jackal0")], major=3)
    assert natnet.unpack_model_def(payload, 3) == {1: "uav0", 7: "jackal0"}


def test_model_def_refuses_descriptor_types_it_cannot_measure():
    """A force-plate descriptor cannot be length-skipped safely, and guessing
    resynchronises the parser onto garbage rather than failing."""
    import struct

    payload = struct.pack("<ii", 1, 3)  # one descriptor, type 3 (force plate)
    with pytest.raises(natnet.NatNetError, match="mocap_id"):
        natnet.unpack_model_def(payload, 3)


def test_frame_round_trip_keeps_pose_and_tracking_flag():
    frame = natnet.MocapFrame(
        42,
        (
            natnet.RigidBody(1, (0.1, 0.2, 0.3), (0.0, 0.0, 0.0, 1.0), 0.0005, True),
            natnet.RigidBody(2, (1.0, 2.0, 3.0), (0.0, 0.5, 0.0, 0.866), 0.002, False),
        ),
    )
    out = natnet.unpack_frame(natnet.pack_frame(frame), 3)
    assert out.frame_number == 42
    assert [b.body_id for b in out.bodies] == [1, 2]
    assert out.bodies[0].position == pytest.approx((0.1, 0.2, 0.3), abs=1e-6)
    assert out.bodies[0].tracking_valid is True
    assert out.bodies[1].tracking_valid is False


def test_frame_parsing_walks_past_marker_sections():
    """Marker sets and legacy markers sit BEFORE the rigid bodies, so getting
    their lengths wrong shifts every position that follows."""
    import struct

    payload = struct.pack("<i", 9)
    payload += struct.pack("<i", 1) + b"all\0" + struct.pack("<i", 2) + struct.pack("<6f", *range(6))
    payload += struct.pack("<i", 3) + struct.pack("<9f", *range(9))
    payload += struct.pack("<i", 1)
    payload += struct.pack("<i", 5) + struct.pack("<3f", 7.0, 8.0, 9.0)
    payload += struct.pack("<4f", 0.0, 0.0, 0.0, 1.0) + struct.pack("<f", 0.0)
    payload += struct.pack("<h", 1)
    out = natnet.unpack_frame(payload, 3)
    assert out.bodies[0].body_id == 5
    assert out.bodies[0].position == pytest.approx((7.0, 8.0, 9.0))


def test_truncated_frame_raises_instead_of_reading_past_the_end():
    good = natnet.pack_frame(
        natnet.MocapFrame(1, (natnet.RigidBody(1, (0, 0, 0), (0, 0, 0, 1), 0.0, True),))
    )
    with pytest.raises(natnet.NatNetError, match="truncated"):
        natnet.unpack_frame(good[:20], 3)


def test_natnet_2_is_refused_rather_than_misparsed():
    with pytest.raises(natnet.NatNetError, match="not supported"):
        natnet.unpack_frame(b"\x00" * 64, 2)


# -- config routing --------------------------------------------------------


def _cfg(**overrides):
    cfg = {
        "agents": [
            {
                "name": "uav0",
                "kind": "px4",
                "ns": "/uav0",
                "spawn": [0.6, 1.4, 0.0],
                "pose_topic": "/mocap/x500_0/pose",
            },
            {
                "name": "jackal0",
                "kind": "diffdrive",
                "spawn": [-2.6, 2.6, 0.1],
                "pose_topic": "/mocap/jackal0/pose",
            },
        ]
    }
    cfg.update(overrides)
    return cfg


def test_settings_default_to_simulation():
    cfg = _cfg()
    assert mocap_settings(cfg)["source"] == "sim"
    assert is_live(cfg) is False


def test_settings_reject_a_typo_rather_than_ignoring_it():
    """A misspelled key that is silently dropped is a mocap host nobody notices
    is still pointing at localhost."""
    with pytest.raises(ValueError, match="unknown keys"):
        mocap_settings(_cfg(mocap={"source": "natnet", "srever": "10.0.0.4"}))


def test_settings_reject_an_impossible_up_axis():
    with pytest.raises(ValueError, match="up_axis"):
        mocap_settings(_cfg(mocap={"up_axis": "x"}))


def test_underscore_keys_are_notes_not_settings():
    s = mocap_settings(_cfg(mocap={"source": "natnet", "_note": "written for humans"}))
    assert s["source"] == "natnet"


def test_routes_default_the_body_name_to_the_agent_name():
    routes = plan_routes(_cfg())
    assert [r.body for r in routes] == ["uav0", "jackal0"]


def test_mocap_body_overrides_the_default():
    cfg = _cfg()
    cfg["agents"][0]["mocap_body"] = "Rigid Body 1"
    assert plan_routes(cfg)[0].body == "Rigid Body 1"


def test_only_px4_agents_get_a_vision_pose_sink():
    """The ground robot has no EKF to feed; the quadrotor's EKF2 drifts without
    one, whatever the sheaf can see."""
    routes = plan_routes(_cfg())
    assert routes[0].vision_topic == "/uav0/mavros/vision_pose/pose"
    assert routes[1].vision_topic is None


def test_an_agent_without_a_pose_topic_gets_no_route():
    cfg = _cfg()
    del cfg["agents"][1]["pose_topic"]
    assert [r.agent for r in plan_routes(cfg)] == ["uav0"]


def test_two_agents_cannot_claim_one_rigid_body():
    cfg = _cfg()
    cfg["agents"][1]["mocap_body"] = "uav0"
    with pytest.raises(ValueError, match="both claim"):
        plan_routes(cfg)


def test_ids_resolve_by_name():
    routes = plan_routes(_cfg())
    by_id = resolve_ids(routes, {3: "jackal0", 8: "uav0"})
    assert by_id[8].agent == "uav0"
    assert by_id[3].pose_topic == "/mocap/jackal0/pose"


def test_a_missing_asset_names_both_sides():
    routes = plan_routes(_cfg())
    with pytest.raises(ValueError) as exc:
        resolve_ids(routes, {3: "jackal0"})
    assert "uav0" in str(exc.value) and "jackal0" in str(exc.value)


def test_mocap_id_wins_over_discovery():
    """The escape hatch for a take whose descriptors this codec will not parse."""
    cfg = _cfg()
    cfg["agents"][0]["mocap_id"] = 11
    by_id = resolve_ids(plan_routes(cfg), {3: "jackal0"})
    assert by_id[11].agent == "uav0"


def test_props_in_the_take_are_reported_not_fatal():
    routes = plan_routes(_cfg())
    descriptions = {1: "uav0", 2: "jackal0", 9: "calib_square"}
    assert sorted(r.agent for r in resolve_ids(routes, descriptions).values()) == [
        "jackal0",
        "uav0",
    ]
    assert unclaimed(routes, descriptions) == ["calib_square"]


# -- sockets ---------------------------------------------------------------


def test_client_handshakes_and_receives_frames_from_the_loopback_server():
    """Everything except the far end of the wire: real UDP, real handshake, real
    model-definition discovery, real reader thread."""
    free = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    free.bind(("127.0.0.1", 0))
    data_port = free.getsockname()[1]
    free.close()

    server = natnet.FakeNatNetServer(
        bodies={4: ("uav0", (1.0, 2.0, 1.4), 0.0), 5: ("jackal0", (-2.6, 2.6, 0.1), 0.0)},
        command_port=0,
        client_data_port=data_port,
        rate_hz=200.0,
    )
    server.start()
    frames = []
    client = natnet.NatNetClient(
        server="127.0.0.1",
        local="127.0.0.1",
        command_port=server.command_port,
        data_port=data_port,
        multicast="",
        on_frame=frames.append,
    )
    try:
        info = client.connect(timeout=3.0)
        assert info.major >= 3
        assert client.request_descriptions(timeout=3.0) == {4: "uav0", 5: "jackal0"}
        client.start()
        deadline = time.monotonic() + 5.0
        while not frames and time.monotonic() < deadline:
            time.sleep(0.02)
        assert frames, "no NatNet frames arrived from the loopback server"
        bodies = {b.body_id: b for b in frames[-1].bodies}
        p, _ = natnet_to_enu(bodies[4].position, bodies[4].quaternion, "y")
        assert p == pytest.approx([1.0, 2.0, 1.4], abs=1e-5)
        p, _ = natnet_to_enu(bodies[5].position, bodies[5].quaternion, "y")
        assert p == pytest.approx([-2.6, 2.6, 0.1], abs=1e-5)
    finally:
        client.stop()
        server.stop()


def test_connect_to_nothing_fails_fast_with_an_actionable_message():
    free = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    free.bind(("127.0.0.1", 0))
    dead_port = free.getsockname()[1]
    free.close()
    client = natnet.NatNetClient(server="127.0.0.1", local="127.0.0.1", command_port=dead_port)
    try:
        with pytest.raises(natnet.NatNetError, match="Streaming pane"):
            client.connect(timeout=0.4)
    finally:
        client.stop()


def test_shipped_config_is_still_a_simulation_config():
    """Guard rail: the demo config must never ship pointed at real hardware."""
    import json

    here = os.path.dirname(__file__)
    path = os.path.join(here, "..", "..", "..", "..", "config", "swarm.json")
    with open(os.path.normpath(path)) as fh:
        cfg = json.load(fh)
    assert is_live(cfg) is False
    # Every agent must be routable, but NOT which agents those are: swarm.json is
    # the file the scenario editor writes to, so pinning the exact roster here
    # fails the suite every time someone saves a scenario, which is the intended
    # workflow rather than a defect. The guard that matters is the one above.
    routes = plan_routes(cfg)
    assert len(routes) == len(cfg["agents"])
    assert {r.agent for r in routes} == {a["name"] for a in cfg["agents"]}
