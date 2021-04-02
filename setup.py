"""Installation of angle_visualization

Author: Romain Fayat, April 2021
"""
import os
from setuptools import setup


def read(fname):
    "Read a file in the current directory."
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(
    name="angle_visualization",
    version="0.1",
    author="Romain Fayat",
    author_email="r.fayat@gmail.com",
    description="Code for visualizing 3D angles with matplotlib",
    install_requires=["numpy",
                      "matplotlib",
                      "scipy",
                      "cartopy",
                      "numpy-quaternion",
                      ],
    packages=["angle_visualization"],
    long_description=read('README.md'),
)
