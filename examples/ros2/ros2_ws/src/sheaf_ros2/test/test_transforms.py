import math
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from sheaf_ros2.transforms import (
    clamp_altitude,
    quat_to_yaw,
    saturate_norm,
    world_velocity_to_unicycle,
)


def test_quat_to_yaw_identity():
    assert abs(quat_to_yaw(0.0, 0.0, 0.0, 1.0)) < 1e-12


def test_quat_to_yaw_half_pi():
    h = math.pi / 4
    assert abs(quat_to_yaw(0.0, 0.0, math.sin(h), math.cos(h)) - math.pi / 2) < 1e-9


def test_saturate_norm_passthrough():
    assert saturate_norm([0.3, 0.4, 0.0], 1.0) == [0.3, 0.4, 0.0]


def test_saturate_norm_scales_and_keeps_direction():
    out = saturate_norm([3.0, 4.0, 0.0], 1.0)
    assert abs(math.sqrt(sum(c * c for c in out)) - 1.0) < 1e-12
    assert abs(out[0] / out[1] - 3.0 / 4.0) < 1e-12


def test_unicycle_aligned_is_pure_forward():
    v, w = world_velocity_to_unicycle(0.5, 0.0, 0.0, 0.2, 10.0, 10.0)
    assert abs(v - 0.5) < 1e-12
    assert abs(w) < 1e-12


def test_unicycle_perpendicular_is_pure_turn():
    v, w = world_velocity_to_unicycle(0.0, 0.5, 0.0, 0.25, 10.0, 10.0)
    assert abs(v) < 1e-12
    assert abs(w - 2.0) < 1e-12


def test_unicycle_is_inverse_of_offset_point_kinematics():
    """Round trip: (v, omega) -> offset-point velocity -> (v, omega)."""
    yaw, l = 0.7, 0.18
    ux, uy = -0.31, 0.44
    v, w = world_velocity_to_unicycle(ux, uy, yaw, l, 10.0, 10.0)
    # forward kinematics of the offset point
    fx = v * math.cos(yaw) - l * w * math.sin(yaw)
    fy = v * math.sin(yaw) + l * w * math.cos(yaw)
    assert abs(fx - ux) < 1e-12
    assert abs(fy - uy) < 1e-12


def test_unicycle_respects_limits():
    v, w = world_velocity_to_unicycle(9.0, 9.0, 0.0, 0.2, 0.8, 2.0)
    assert abs(v) <= 0.8 + 1e-12
    assert abs(w) <= 2.0 + 1e-12


def test_clamp_altitude_band():
    assert clamp_altitude(0.5, 1.5, 0.3, 2.0) == 0.5
    assert clamp_altitude(0.5, 2.4, 0.3, 2.0) < 0.0   # above ceiling: pushed down
    assert clamp_altitude(-0.5, 0.1, 0.3, 2.0) > 0.0  # below floor: pushed up
