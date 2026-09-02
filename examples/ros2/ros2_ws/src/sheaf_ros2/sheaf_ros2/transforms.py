"""Frame and dynamics conversions for the sheaf bridge.

The sheaf controller emits a world-frame ENU velocity per agent. Turning that
into something a particular robot accepts is a *local* concern and lives here:

  * PX4 multirotor  -- takes the ENU velocity directly (MAVROS converts to NED).
  * Differential drive -- is nonholonomic, so we use the standard near-identity
    diffeomorphism: control the point a distance `offset` ahead of the wheel
    axle instead of the axle centre.

Pure math, no ROS and no numpy, so it is unit testable anywhere.
"""

from __future__ import annotations

import math
from typing import Sequence, Tuple


def quat_to_yaw(x: float, y: float, z: float, w: float) -> float:
    """Yaw (rotation about +z) of a quaternion, in radians, wrapped to (-pi, pi]."""
    siny_cosp = 2.0 * (w * z + x * y)
    cosy_cosp = 1.0 - 2.0 * (y * y + z * z)
    return math.atan2(siny_cosp, cosy_cosp)


def saturate_norm(vec: Sequence[float], limit: float) -> list:
    """Scale `vec` down (never up) so its Euclidean norm is at most `limit`.

    Scaling the whole vector rather than clipping components preserves the
    commanded direction, which matters: componentwise clipping bends the
    velocity away from the sheaf gradient and can stall consensus.
    """
    if limit <= 0.0:
        return [0.0] * len(vec)
    n = math.sqrt(sum(c * c for c in vec))
    if n <= limit or n == 0.0:
        return [float(c) for c in vec]
    s = limit / n
    return [float(c) * s for c in vec]


def world_velocity_to_unicycle(
    ux: float,
    uy: float,
    yaw: float,
    offset: float,
    v_max: float,
    omega_max: float,
) -> Tuple[float, float]:
    """Near-identity diffeomorphism: planar world velocity -> (v, omega).

    With the controlled point at p = [x + l*cos(yaw), y + l*sin(yaw)],

        [v; omega] = [ cos(yaw)      sin(yaw)     ] [ux]
                     [-sin(yaw)/l    cos(yaw)/l   ] [uy]

    which is exactly invertible for l > 0. `offset` is l. Outputs are clamped
    to the platform limits *after* the transform, so the robot degrades to
    "turn first, then drive" rather than driving off in the wrong direction.
    """
    if offset <= 0.0:
        raise ValueError("offset (l) must be > 0 for the unicycle transform")
    c, s = math.cos(yaw), math.sin(yaw)
    v = c * ux + s * uy
    omega = (-s * ux + c * uy) / offset
    v = max(-v_max, min(v_max, v))
    omega = max(-omega_max, min(omega_max, omega))
    return v, omega


def clamp_altitude(uz: float, z: float, z_min: float, z_max: float, k: float = 1.0) -> float:
    """Soft altitude fence: cancel any climb above z_max / descent below z_min and
    push back toward the band. Keeps the Astro SOP ceiling enforceable in sim."""
    if z > z_max:
        return min(uz, 0.0) - k * (z - z_max)
    if z < z_min:
        return max(uz, 0.0) + k * (z_min - z)
    return uz


def natnet_to_enu(
    position: Sequence[float],
    quaternion: Sequence[float],
    up_axis: str = "y",
) -> Tuple[list, list]:
    """Motive's streaming frame -> the ENU world frame everything else here uses.

    Motive streams in whichever up axis the Streaming pane is set to, and the
    default is **Y-up**: +X right, +Y up, +Z toward the front of the volume,
    right-handed. Our world frame is ENU: +X east, +Y north, +Z up. Taking
    East = X, North = -Z, Up = Y is the unique axis permutation that lines the
    up axes up while keeping the handedness, so

        p_enu = (x, -z, y)

    and, because a change of basis by a proper rotation P rotates a quaternion's
    vector part and leaves its scalar alone,

        q_enu = (qx, -qz, qy, qw).

    Set Motive's Streaming pane to Z-up and this is the identity: Z-up Motive is
    already ENU. Both are supported because which one a park is configured for
    is an operator setting we cannot check from here, and getting it wrong is
    not a crash -- it is a swarm that flies a formation rotated 90 degrees into
    the netting. `up_axis` therefore has to be stated in the config, not guessed.

    Units are metres and the quaternion is (x, y, z, w), matching both NatNet
    and geometry_msgs, so no scaling or reordering happens here.
    """
    x, y, z = (float(c) for c in position)
    qx, qy, qz, qw = (float(c) for c in quaternion)
    axis = up_axis.lower()
    if axis == "z":
        return [x, y, z], [qx, qy, qz, qw]
    if axis == "y":
        return [x, -z, y], [qx, -qz, qy, qw]
    raise ValueError(f"up_axis must be 'y' or 'z', got {up_axis!r}")
