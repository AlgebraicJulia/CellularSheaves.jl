"""Bring up the ROS 2 half of the Autonomy Park sheaf demo.

Everything here is generated from config/swarm.json so that the agent list, the
graph and the gains have exactly one home. Per agent kind this launches:

  px4        one mavros_node, namespaced /<name>, talking to PX4 SITL instance i
             on the conventional UDP port pair (14540+i in, 14580+i out)
  diffdrive  one ros_gz_bridge parameter_bridge for that model's cmd_vel/odometry

plus the /clock bridge, the sheaf_bridge node itself, one mocap source picked by
the config's `mocap.source` (Gazebo ground truth or the park's OptiTrack), and --
when the config has a Unitree body in it -- the gait driver that animates its
legs.

Gazebo and the PX4 SITL instances are *not* started here -- they are plain
processes with their own environment, launched by scripts/start_sim.sh. Keeping
them out of the ROS launch graph means you can restart the controller without
restarting the simulator, which is most of what you do while tuning.

    ros2 launch sheaf_ros2 sim.launch.py config:=/path/to/swarm.json
"""

import json
import os

from ament_index_python.packages import get_package_share_directory
from launch import LaunchDescription
from launch.actions import DeclareLaunchArgument, OpaqueFunction
from launch.substitutions import LaunchConfiguration
from launch_ros.actions import Node

from sheaf_ros2.mocap_routing import is_live


def mavros_launch(filename: str) -> str:
    """Path to one of the yaml files mavros ships. PX4 SITL works without them,
    but the denylist in px4_pluginlists.yaml turns off plugins that spew errors
    against a simulated airframe with no camera or rangefinder."""
    return os.path.join(get_package_share_directory("mavros"), "launch", filename)


def _default_config() -> str:
    env = os.environ.get("SHEAF_CONFIG")
    if env:
        return env
    return os.path.join(get_package_share_directory("sheaf_ros2"), "config", "swarm.json")


def _spawn(context, *args, **kwargs):
    config_path = LaunchConfiguration("config").perform(context)
    use_sim_time = LaunchConfiguration("use_sim_time").perform(context).lower() in ("true", "1")

    with open(config_path) as fh:
        cfg = json.load(fh)

    world = os.environ.get("GZ_WORLD", "autonomy_park")

    nodes = [
        Node(
            package="ros_gz_bridge",
            executable="parameter_bridge",
            name="clock_bridge",
            arguments=["/clock@rosgraph_msgs/msg/Clock[gz.msgs.Clock"],
            output="screen",
        ),
    ]

    # Exactly one thing may own the pose topics. `mocap.source` in swarm.json
    # picks it: "sim" gets the Gazebo stand-in (which reads gz-transport
    # directly, because the ros_gz TFMessage route drops entity names -- see
    # mocap_node.py), "natnet" gets the real OptiTrack client. Running both
    # would put two publishers on every `pose_topic` and the adapters would
    # interleave simulated and measured poses without complaining once.
    if is_live(cfg):
        nodes.append(
            Node(
                package="sheaf_ros2",
                executable="mocap_natnet",
                name="mocap_natnet",
                output="screen",
                # No use_sim_time: Motive frames are timestamped off the wall
                # clock, and there is no Gazebo publishing /clock at the park.
                parameters=[{"config": config_path}],
            )
        )
    else:
        nodes.append(
            Node(
                package="sheaf_ros2",
                executable="mocap_bridge",
                name="mocap_bridge",
                output="screen",
                parameters=[{"config": config_path, "world": world, "use_sim_time": use_sim_time}],
            )
        )

    for i, spec in enumerate(cfg["agents"]):
        kind = spec["kind"]
        name = spec["name"]

        if kind == "px4":
            instance = int(spec.get("instance", i))
            if instance < int(os.environ.get("CAMERAS", "0") or 0) and spec.get("camera_topic"):
                cam = spec["camera_topic"]
                nodes.append(
                    Node(
                        package="ros_gz_image",
                        executable="image_bridge",
                        name=f"{name}_camera_bridge",
                        arguments=[cam],
                        remappings=[(cam, f"/{name}/camera")],
                        parameters=[{"use_sim_time": use_sim_time}],
                        output="screen",
                    )
                )
            nodes.append(
                Node(
                    package="mavros",
                    executable="mavros_node",
                    # No `name=`. mavros_node is a composed container: it runs a
                    # router node plus one node per plugin. Setting `name` remaps
                    # __node for all of them at once, so two plugins end up
                    # creating the same topic with different types and the process
                    # aborts with "create_subscription() called for existing topic
                    # name .../status with incompatible type".
                    #
                    # The namespace carries the whole layout instead: mavros 2 does
                    # not add a `mavros` level of its own, so /uav0 alone would give
                    # /uav0/state, not /uav0/mavros/state.
                    namespace=f"/{name}/mavros",
                    output="screen",
                    # A MAVROS that dies (observed: abort after its FCU went
                    # insane in free fall) takes its whole vehicle out of the
                    # demo permanently unless someone restarts it.
                    respawn=True,
                    respawn_delay=3.0,
                    parameters=[
                        mavros_launch("px4_pluginlists.yaml"),
                        mavros_launch("px4_config.yaml"),
                        {
                            "fcu_url": f"udp://:{14540 + instance}@127.0.0.1:{14580 + instance}",
                            "gcs_url": "",
                            # These are the names mavros_node actually declares
                            # (see mavros/launch/node.launch). ROS 2 rejects
                            # undeclared parameters, so target_system_id and
                            # friends do not merely get ignored -- the node
                            # refuses to start.
                            "tgt_system": 1 + instance,
                            "tgt_component": 1,
                            "fcu_protocol": "v2.0",
                            "use_sim_time": use_sim_time,
                        },
                    ],
                )
            )

        elif kind == "diffdrive":
            model = spec.get("gz_model", name)
            nodes.append(
                Node(
                    package="ros_gz_bridge",
                    executable="parameter_bridge",
                    name=f"{model}_bridge",
                    arguments=[
                        f"/model/{model}/cmd_vel@geometry_msgs/msg/Twist]gz.msgs.Twist",
                        f"/model/{model}/odometry@nav_msgs/msg/Odometry[gz.msgs.Odometry",
                    ],
                    parameters=[{"use_sim_time": use_sim_time}],
                    output="screen",
                )
            )

        else:
            raise ValueError(f"agent {name!r} has unknown kind {kind!r}")

    # The Unitree bodies are driven by the same DiffDrive plugin as the rovers,
    # so nothing above knows their legs exist. This animates them off the
    # odometry that plugin already publishes; a config with no quadruped in it
    # does not get the node at all.
    if any(spec.get("platform") in ("b1", "go1") for spec in cfg["agents"]):
        nodes.append(
            Node(
                package="sheaf_ros2",
                executable="gait_driver",
                name="gait_driver",
                output="screen",
                parameters=[{"config": config_path, "use_sim_time": use_sim_time}],
            )
        )

    nodes.append(
        Node(
            package="sheaf_ros2",
            executable="sheaf_bridge",
            name="sheaf_bridge",
            output="screen",
            parameters=[{"config": config_path, "use_sim_time": use_sim_time}],
        )
    )
    return nodes


def generate_launch_description() -> LaunchDescription:
    return LaunchDescription(
        [
            DeclareLaunchArgument("config", default_value=_default_config()),
            DeclareLaunchArgument("use_sim_time", default_value="true"),
            OpaqueFunction(function=_spawn),
        ]
    )
