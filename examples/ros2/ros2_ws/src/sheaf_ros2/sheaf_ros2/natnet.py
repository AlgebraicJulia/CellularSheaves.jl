"""NatNet 3/4 codec, client and loopback server. No ROS, no numpy, no SDK.

Why a NatNet client and not `mocap4r2` or MAVROS `vision_pose`:

  * The park's mocap is OptiTrack driven by Motive, and NatNet is the protocol
    Motive already speaks on the wire. Everything else on the shortlist ends up
    speaking it too, one layer down.
  * `mocap4r2_optitrack_driver` is a thin wrapper over Motive's NatNet SDK: it
    pulls in the mocap4r2 stack plus a vendored, version-pinned closed-source
    SDK binary, and then emits `mocap4r2_msgs/RigidBodies`, which is not the
    `PoseStamped` the adapters read. So that route is this file plus two extra
    packages in the image plus a translating node. The wrapper buys nothing we
    do not have to write anyway, and it adds a binary we cannot run in CI.
  * MAVROS `vision_pose` is not a source at all, it is a *sink*. It only exists
    for the PX4 vehicles and it carries pose *into* an EKF; it cannot tell the
    ground robot where it is. It is complementary, not an alternative, which is
    why `mocap_natnet_node` publishes to it in addition to the pose topics.

The parsing here is deliberately partial: a frame is read as far as the rigid
bodies and then abandoned. Labeled markers, skeletons, force plates and the
timing block all sit *after* the rigid bodies in the packet, and we need none of
them, so refusing to parse them removes most of the version-to-version churn in
this protocol. The one thing lost is Motive's own transmit timestamp, so latency
is measured at arrival instead; on the park's wired network that costs a couple
of milliseconds against a 50 ms control tick.

NatNet major < 3 is rejected rather than guessed at: the 2.x frame packs marker
data inside each rigid body, so a 2.x stream parsed as 3.x silently yields
plausible-looking garbage positions. Motive 2.0 and newer stream NatNet 3+.
"""

from __future__ import annotations

import math
import socket
import struct
import threading
import time
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence, Tuple

NAT_CONNECT = 0
NAT_SERVERINFO = 1
NAT_REQUEST = 2
NAT_RESPONSE = 3
NAT_REQUEST_MODELDEF = 4
NAT_MODELDEF = 5
NAT_REQUEST_FRAMEOFDATA = 6
NAT_FRAMEOFDATA = 7
NAT_MESSAGESTRING = 8
NAT_DISCONNECT = 9
NAT_KEEPALIVE = 10
NAT_UNRECOGNIZED_REQUEST = 100

DESC_MARKERSET = 0
DESC_RIGIDBODY = 1
DESC_SKELETON = 2

DEFAULT_COMMAND_PORT = 1510
DEFAULT_DATA_PORT = 1511
DEFAULT_MULTICAST = "239.255.42.99"

# Motive marks a rigid body it has lost as "not tracking" but keeps streaming its
# last known pose. Bit 0 of the per-body params word is that flag.
TRACKING_VALID = 0x01


@dataclass(frozen=True)
class RigidBody:
    """One Motive rigid body, still in Motive's own frame and units."""

    body_id: int
    position: Tuple[float, float, float]
    quaternion: Tuple[float, float, float, float]  # x, y, z, w
    mean_error: float
    tracking_valid: bool


@dataclass(frozen=True)
class MocapFrame:
    frame_number: int
    bodies: Tuple[RigidBody, ...]


@dataclass(frozen=True)
class ServerInfo:
    name: str
    app_version: Tuple[int, int, int, int]
    natnet_version: Tuple[int, int, int, int]

    @property
    def major(self) -> int:
        return self.natnet_version[0]


class NatNetError(RuntimeError):
    """Raised for anything the codec refuses to guess at."""


# -- byte-level helpers ----------------------------------------------------


class _Reader:
    """Cursor over a payload. Every read is bounds checked, because a truncated
    UDP datagram parsed with bare struct.unpack_from reads whatever follows it in
    the buffer and reports it as a robot position."""

    def __init__(self, data: bytes):
        self.data = data
        self.pos = 0

    def take(self, n: int) -> bytes:
        if n < 0 or self.pos + n > len(self.data):
            raise NatNetError(
                f"truncated NatNet payload: wanted {n} bytes at offset {self.pos} "
                f"of {len(self.data)}"
            )
        out = self.data[self.pos : self.pos + n]
        self.pos += n
        return out

    def int32(self) -> int:
        return struct.unpack("<i", self.take(4))[0]

    def int16(self) -> int:
        return struct.unpack("<h", self.take(2))[0]

    def float32(self) -> float:
        return struct.unpack("<f", self.take(4))[0]

    def floats(self, n: int) -> Tuple[float, ...]:
        return struct.unpack(f"<{n}f", self.take(4 * n))

    def cstring(self) -> str:
        end = self.data.find(b"\0", self.pos)
        if end < 0:
            raise NatNetError("unterminated string in NatNet payload")
        out = self.data[self.pos : end].decode("utf-8", "replace")
        self.pos = end + 1
        return out

    def skip(self, n: int) -> None:
        self.take(n)


def _cstring(text: str) -> bytes:
    return text.encode("utf-8") + b"\0"


def pack_message(msg_id: int, payload: bytes) -> bytes:
    """NatNet framing: little-endian id, little-endian payload size, payload."""
    return struct.pack("<HH", msg_id, len(payload)) + payload


def unpack_message(data: bytes) -> Tuple[int, bytes]:
    """Split a datagram into (message id, payload).

    The declared size is trusted only as far as the datagram actually goes;
    Motive has been seen to under-report it on the model-definition packet, and
    truncating to a short header would drop rigid bodies with no error at all.
    """
    if len(data) < 4:
        raise NatNetError(f"NatNet packet shorter than its header ({len(data)} bytes)")
    msg_id, size = struct.unpack("<HH", data[:4])
    payload = data[4:]
    if size and size <= len(payload):
        payload = payload[:size]
    return msg_id, payload


def pack_command(command: int, command_str: str = "") -> bytes:
    """Build one of the short command packets a client sends to port 1510.

    `Connect` carries the literal payload "Ping", which is what Motive matches
    on; the request commands carry an empty payload. This mirrors the reference
    SDK's `send_command` exactly, because Motive ignores anything else.
    """
    if command == NAT_CONNECT:
        return pack_message(command, _cstring("Ping"))
    if command in (NAT_REQUEST_MODELDEF, NAT_REQUEST_FRAMEOFDATA, NAT_KEEPALIVE, NAT_DISCONNECT):
        return pack_message(command, b"")
    if command == NAT_REQUEST:
        return pack_message(command, _cstring(command_str))
    raise NatNetError(f"unsupported NatNet command {command}")


# -- server info -----------------------------------------------------------


def pack_server_info(info: ServerInfo) -> bytes:
    name = info.name.encode("utf-8")[:255]
    return (
        name.ljust(256, b"\0")
        + bytes(info.app_version)
        + bytes(info.natnet_version)
    )


def unpack_server_info(payload: bytes) -> ServerInfo:
    """`char name[256]; uint8 appVersion[4]; uint8 natNetVersion[4]`."""
    if len(payload) < 264:
        raise NatNetError(f"ServerInfo payload is {len(payload)} bytes, expected >= 264")
    name = payload[:256].split(b"\0", 1)[0].decode("utf-8", "replace")
    app = tuple(payload[256:260])
    net = tuple(payload[260:264])
    return ServerInfo(name=name, app_version=app, natnet_version=net)


# -- model definitions -----------------------------------------------------


def pack_model_def(bodies: Sequence[Tuple[int, str]], major: int = 3) -> bytes:
    """Encode a rigid-body-only dataset, the shape Motive sends when the only
    assets in the take are rigid bodies. Used by the loopback server and the
    tests; the client never needs to encode this."""
    out = struct.pack("<i", len(bodies))
    for body_id, name in bodies:
        out += struct.pack("<i", DESC_RIGIDBODY)
        out += _cstring(name)
        out += struct.pack("<ii", body_id, -1)
        out += struct.pack("<3f", 0.0, 0.0, 0.0)
        if major >= 3:
            out += struct.pack("<i", 0)  # marker count
    return out


def _skip_rigid_body_desc(r: _Reader, major: int) -> Tuple[int, str]:
    name = r.cstring() if major >= 2 else ""
    body_id = r.int32()
    r.int32()      # parent id
    r.floats(3)    # offset from parent
    if major >= 3:
        n_markers = r.int32()
        r.skip(12 * n_markers)  # marker positions
        r.skip(4 * n_markers)   # active labels
        if major >= 4:
            for _ in range(n_markers):
                r.cstring()
    return body_id, name


def unpack_model_def(payload: bytes, major: int = 3) -> Dict[int, str]:
    """Rigid-body streaming id -> name, from a NAT_MODELDEF packet.

    Marker-set and skeleton descriptors are walked past because their lengths
    are computable. Force plate, device and camera descriptors are not walked:
    their layouts changed across 3.x and a mis-skip would resynchronise onto
    garbage rather than fail. If the park streams any of those, pin the ids in
    the config with `mocap_id` and skip discovery altogether -- which is what
    the error message says.
    """
    r = _Reader(payload)
    count = r.int32()
    out: Dict[int, str] = {}
    for _ in range(count):
        kind = r.int32()
        if kind == DESC_RIGIDBODY:
            body_id, name = _skip_rigid_body_desc(r, major)
            out[body_id] = name
        elif kind == DESC_MARKERSET:
            r.cstring()
            for _ in range(r.int32()):
                r.cstring()
        elif kind == DESC_SKELETON:
            r.cstring()
            r.int32()  # skeleton id
            for _ in range(r.int32()):
                _skip_rigid_body_desc(r, major)
        else:
            raise NatNetError(
                f"model definition contains descriptor type {kind}, whose length this "
                "codec does not compute; set `mocap_id` per agent in swarm.json to "
                "name rigid bodies by streaming id instead of by discovery"
            )
    return out


# -- frames ----------------------------------------------------------------


def pack_frame(frame: MocapFrame, major: int = 3) -> bytes:
    """Encode the prefix of a data frame this codec reads back: frame number,
    an empty marker-set section, an empty legacy-marker section, then the rigid
    bodies. A trailing filler block stands in for everything real Motive puts
    after the rigid bodies, so the loopback server exercises the "stop parsing
    early" path rather than hiding it."""
    if major < 3:
        raise NatNetError("pack_frame only emits NatNet 3+ frames")
    out = struct.pack("<i", frame.frame_number)
    out += struct.pack("<i", 0)  # marker set count
    out += struct.pack("<i", 0)  # legacy other-marker count
    out += struct.pack("<i", len(frame.bodies))
    for b in frame.bodies:
        out += struct.pack("<i", b.body_id)
        out += struct.pack("<3f", *b.position)
        out += struct.pack("<4f", *b.quaternion)
        out += struct.pack("<f", b.mean_error)
        out += struct.pack("<h", TRACKING_VALID if b.tracking_valid else 0)
    out += b"\x00" * 64  # skeletons, labeled markers, timing: never parsed
    return out


def unpack_frame(payload: bytes, major: int = 3) -> MocapFrame:
    """Read a NAT_FRAMEOFDATA packet as far as the rigid bodies."""
    if major < 3:
        raise NatNetError(
            f"NatNet major version {major} is not supported; its frame packs marker "
            "data inside each rigid body and misreading that yields plausible but "
            "wrong positions. Motive 2.0 and newer stream NatNet 3 or later."
        )
    r = _Reader(payload)
    frame_number = r.int32()
    for _ in range(r.int32()):        # marker sets
        r.cstring()
        r.skip(12 * r.int32())
    r.skip(12 * r.int32())            # legacy unlabeled markers
    bodies: List[RigidBody] = []
    for _ in range(r.int32()):
        body_id = r.int32()
        px, py, pz = r.floats(3)
        qx, qy, qz, qw = r.floats(4)
        mean_error = r.float32()
        params = r.int16()
        bodies.append(
            RigidBody(
                body_id=body_id,
                position=(px, py, pz),
                quaternion=(qx, qy, qz, qw),
                mean_error=mean_error,
                tracking_valid=bool(params & TRACKING_VALID),
            )
        )
    return MocapFrame(frame_number=frame_number, bodies=tuple(bodies))


# -- client ----------------------------------------------------------------


class NatNetClient:
    """Minimal Motive client: handshake on the command port, frames on the data
    port, one callback per frame off a reader thread.

    Multicast is the Motive default and is what `multicast` selects; setting it
    to "" puts the client in unicast mode, where it simply binds the data port
    and waits for whatever Motive (or the loopback server) sends there. Unicast
    is what the tests use, because a container without a multicast route would
    otherwise make them environment dependent.
    """

    def __init__(
        self,
        server: str,
        local: str = "0.0.0.0",
        command_port: int = DEFAULT_COMMAND_PORT,
        data_port: int = DEFAULT_DATA_PORT,
        multicast: str = DEFAULT_MULTICAST,
        on_frame: Optional[Callable[[MocapFrame], None]] = None,
        on_error: Optional[Callable[[Exception], None]] = None,
    ):
        self.server = server
        self.local = local
        self.command_port = int(command_port)
        self.data_port = int(data_port)
        self.multicast = multicast or ""
        self.on_frame = on_frame
        self.on_error = on_error

        self.info: Optional[ServerInfo] = None
        self.descriptions: Dict[int, str] = {}
        self._major = 3
        self._cmd_sock: Optional[socket.socket] = None
        self._data_sock: Optional[socket.socket] = None
        self._thread: Optional[threading.Thread] = None
        self._stop = threading.Event()

    # The handshake is separate from `start` so a caller can fail fast on a
    # misconfigured address instead of discovering it as silence at flight time.
    def connect(self, timeout: float = 3.0) -> ServerInfo:
        sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        sock.bind((self.local, 0))
        sock.settimeout(timeout)
        self._cmd_sock = sock
        addr = (self.server, self.command_port)

        sock.sendto(pack_command(NAT_CONNECT), addr)
        info = None
        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            try:
                data, _ = sock.recvfrom(65535)
            except socket.timeout:
                break
            msg_id, payload = unpack_message(data)
            if msg_id == NAT_SERVERINFO:
                info = unpack_server_info(payload)
                break
        if info is None:
            raise NatNetError(
                f"no ServerInfo from {self.server}:{self.command_port} within {timeout} s; "
                "check that Motive's Streaming pane is enabled and that this host can "
                "reach it (Motive answers only the interface it is bound to)"
            )
        self.info = info
        self._major = info.major
        if self._major < 3:
            raise NatNetError(
                f"Motive is streaming NatNet {'.'.join(map(str, info.natnet_version))}; "
                "this client needs 3 or later"
            )
        return info

    def request_descriptions(self, timeout: float = 3.0) -> Dict[int, str]:
        """Ask for the asset list so rigid bodies can be matched by NAME.

        Names are the thing an operator actually sets in Motive; streaming ids
        get renumbered whenever assets are recreated, so a config that pinned
        ids would silently swap two robots after a re-calibration.
        """
        if self._cmd_sock is None:
            raise NatNetError("connect() first")
        self._cmd_sock.sendto(pack_command(NAT_REQUEST_MODELDEF), (self.server, self.command_port))
        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            try:
                data, _ = self._cmd_sock.recvfrom(65535)
            except socket.timeout:
                break
            msg_id, payload = unpack_message(data)
            if msg_id == NAT_MODELDEF:
                self.descriptions = unpack_model_def(payload, self._major)
                return self.descriptions
        raise NatNetError("Motive did not answer the model-definition request")

    def start(self) -> None:
        sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        if self.multicast:
            sock.bind(("", self.data_port))
            mreq = socket.inet_aton(self.multicast) + socket.inet_aton(self.local)
            sock.setsockopt(socket.IPPROTO_IP, socket.IP_ADD_MEMBERSHIP, mreq)
        else:
            sock.bind((self.local, self.data_port))
        sock.settimeout(0.5)
        self._data_sock = sock
        self._stop.clear()
        self._thread = threading.Thread(target=self._read_loop, name="natnet", daemon=True)
        self._thread.start()

    @property
    def data_port_bound(self) -> int:
        """The port actually bound, which differs from `data_port` when 0 was
        requested. Tests bind 0 so they can run concurrently."""
        if self._data_sock is None:
            raise NatNetError("start() first")
        return self._data_sock.getsockname()[1]

    def _read_loop(self) -> None:
        while not self._stop.is_set():
            try:
                data, _ = self._data_sock.recvfrom(65535)
            except socket.timeout:
                continue
            except OSError:
                return
            try:
                msg_id, payload = unpack_message(data)
                if msg_id != NAT_FRAMEOFDATA:
                    continue
                frame = unpack_frame(payload, self._major)
            except Exception as exc:  # noqa: BLE001 - one bad datagram must not end the stream
                if self.on_error is not None:
                    self.on_error(exc)
                continue
            if self.on_frame is not None:
                self.on_frame(frame)

    def stop(self) -> None:
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=2.0)
            self._thread = None
        for sock in (self._data_sock, self._cmd_sock):
            if sock is not None:
                try:
                    sock.close()
                except OSError:
                    pass
        self._data_sock = self._cmd_sock = None


# -- loopback server -------------------------------------------------------


@dataclass
class FakeNatNetServer:
    """A Motive stand-in that speaks the real protocol on real sockets.

    This exists because the one thing that cannot be tested at a desk is the
    OptiTrack system, and a node tested only through mocked-out internals is a
    node whose sockets, threading and handshake have never run. Pointing the
    real client at this server exercises all of it; only the bytes on the far
    side are synthetic.

    Bodies are given in ENU (the frame the node is expected to publish) and
    encoded back into Motive's Y-up convention on the way out, so a test can
    assert that what comes out of the node equals what went into the server.
    """

    bodies: Dict[int, Tuple[str, Tuple[float, float, float], float]]
    host: str = "127.0.0.1"
    command_port: int = 0
    data_port: int = 0
    client_data_port: int = DEFAULT_DATA_PORT
    rate_hz: float = 100.0
    up_axis: str = "y"
    motion: str = "static"
    natnet_version: Tuple[int, int, int, int] = (3, 1, 0, 0)
    _cmd: Optional[socket.socket] = field(default=None, init=False)
    _out: Optional[socket.socket] = field(default=None, init=False)
    _threads: List[threading.Thread] = field(default_factory=list, init=False)
    _stop: threading.Event = field(default_factory=threading.Event, init=False)
    frames_sent: int = field(default=0, init=False)

    def start(self) -> None:
        self._cmd = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        self._cmd.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        self._cmd.bind((self.host, self.command_port))
        self._cmd.settimeout(0.2)
        self.command_port = self._cmd.getsockname()[1]
        self._out = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        self._out.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        self._out.setsockopt(socket.IPPROTO_IP, socket.IP_MULTICAST_TTL, 2)
        self._stop.clear()
        for target in (self._command_loop, self._stream_loop):
            t = threading.Thread(target=target, daemon=True)
            t.start()
            self._threads.append(t)

    def _command_loop(self) -> None:
        info = ServerInfo("FakeMotive", (0, 0, 0, 0), self.natnet_version)
        while not self._stop.is_set():
            try:
                data, peer = self._cmd.recvfrom(65535)
            except socket.timeout:
                continue
            except OSError:
                return
            try:
                msg_id, _ = unpack_message(data)
            except NatNetError:
                continue
            if msg_id == NAT_CONNECT:
                self._cmd.sendto(pack_message(NAT_SERVERINFO, pack_server_info(info)), peer)
            elif msg_id == NAT_REQUEST_MODELDEF:
                defs = [(bid, name) for bid, (name, _, _) in sorted(self.bodies.items())]
                self._cmd.sendto(
                    pack_message(NAT_MODELDEF, pack_model_def(defs, self.natnet_version[0])), peer
                )

    def _encode(self, enu: Sequence[float], yaw: float) -> RigidBody:
        x, y, z = enu
        half = yaw / 2.0
        if self.up_axis == "z":
            pos = (x, y, z)
            quat = (0.0, 0.0, math.sin(half), math.cos(half))
        else:
            # Inverse of transforms.natnet_to_enu: ENU (x, y, z) came from Motive
            # (x, z, -y), and the yaw axis is Motive's +Y.
            pos = (x, z, -y)
            quat = (0.0, math.sin(half), 0.0, math.cos(half))
        return RigidBody(0, pos, quat, 0.0002, True)

    def _stream_loop(self) -> None:
        period = 1.0 / max(self.rate_hz, 1e-3)
        n = 0
        t0 = time.monotonic()
        while not self._stop.is_set():
            t = time.monotonic() - t0
            bodies = []
            for body_id, (_, enu, yaw) in sorted(self.bodies.items()):
                p = list(enu)
                if self.motion == "circle":
                    p[0] += 0.5 * math.cos(0.4 * t)
                    p[1] += 0.5 * math.sin(0.4 * t)
                rb = self._encode(p, yaw)
                bodies.append(
                    RigidBody(body_id, rb.position, rb.quaternion, rb.mean_error, True)
                )
            packet = pack_message(
                NAT_FRAMEOFDATA,
                pack_frame(MocapFrame(n, tuple(bodies)), self.natnet_version[0]),
            )
            try:
                self._out.sendto(packet, (self.host, self.client_data_port))
            except OSError:
                return
            n += 1
            self.frames_sent = n
            time.sleep(period)

    def stop(self) -> None:
        self._stop.set()
        for t in self._threads:
            t.join(timeout=2.0)
        self._threads = []
        for sock in (self._cmd, self._out):
            if sock is not None:
                try:
                    sock.close()
                except OSError:
                    pass
        self._cmd = self._out = None


def fake_main(argv=None) -> None:
    """`fake_natnet` console script: stream a swarm.json's agents as rigid bodies.

    Poses come straight from each agent's `spawn`, so whatever the node
    republishes can be checked against the config by eye or by script. This is
    the desk substitute for the park's Motive machine.
    """
    import argparse
    import json

    from .mocap_routing import mocap_settings, plan_routes

    ap = argparse.ArgumentParser(description="loopback OptiTrack/Motive stand-in")
    ap.add_argument("--config", required=True, help="path to swarm.json")
    ap.add_argument("--command-port", type=int, default=None)
    ap.add_argument("--data-port", type=int, default=None, help="port to stream TO")
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--rate", type=float, default=100.0)
    ap.add_argument("--motion", choices=("static", "circle"), default="static")
    ap.add_argument("--seconds", type=float, default=0.0, help="0 runs until interrupted")
    args = ap.parse_args(argv)

    with open(args.config) as fh:
        cfg = json.load(fh)
    settings = mocap_settings(cfg)
    routes = plan_routes(cfg)
    bodies = {
        i + 1: (r.body, tuple(r.spawn), 0.0)
        for i, r in enumerate(routes)
    }
    server = FakeNatNetServer(
        bodies=bodies,
        host=args.host,
        command_port=args.command_port if args.command_port is not None else settings["command_port"],
        client_data_port=args.data_port if args.data_port is not None else settings["data_port"],
        rate_hz=args.rate,
        up_axis=settings["up_axis"],
        motion=args.motion,
    )
    server.start()
    print(
        f"fake Motive: command {args.host}:{server.command_port}, streaming "
        f"{len(bodies)} rigid bodies to {args.host}:{server.client_data_port} "
        f"at {args.rate:g} Hz, {settings['up_axis']}-up, motion {args.motion}",
        flush=True,
    )
    for body_id, (name, enu, _) in sorted(bodies.items()):
        print(f"  id {body_id}: {name} at ENU {tuple(round(c, 3) for c in enu)}", flush=True)
    try:
        if args.seconds > 0:
            time.sleep(args.seconds)
        else:
            while True:
                time.sleep(0.5)
    except KeyboardInterrupt:
        pass
    finally:
        server.stop()


if __name__ == "__main__":
    fake_main()
