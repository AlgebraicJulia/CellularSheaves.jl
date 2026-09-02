import os
from glob import glob

from setuptools import setup

package_name = "sheaf_ros2"

setup(
    name=package_name,
    version="0.1.0",
    packages=[package_name],
    data_files=[
        ("share/ament_index/resource_index/packages", ["resource/" + package_name]),
        (os.path.join("share", package_name), ["package.xml"]),
        (os.path.join("share", package_name, "launch"), glob("launch/*.launch.py")),
    ],
    install_requires=["setuptools"],
    zip_safe=True,
    maintainer="Itay Kadosh",
    maintainer_email="adamkadosh@gmail.com",
    description="ROS 2 bridge between a Julia cellular-sheaf controller and PX4 / differential-drive robots.",
    license="MIT",
    tests_require=["pytest"],
    entry_points={
        "console_scripts": [
            "sheaf_bridge = sheaf_ros2.bridge_node:main",
            "mocap_bridge = sheaf_ros2.mocap_node:main",
            "mocap_natnet = sheaf_ros2.mocap_natnet_node:main",
            "fake_natnet = sheaf_ros2.natnet:fake_main",
            "gait_driver = sheaf_ros2.gait_node:main",
        ],
    },
)
