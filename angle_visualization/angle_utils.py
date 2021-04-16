"""Utils functions when dealing with 3-dimensional angles.

Author: Romain Fayat, March 2021
"""
import numpy as np
from .utils import unpack_first_arg, list_args_to_array
RAD2DEG = 180. / np.pi
DEG2RAD = np.pi / 180.


@list_args_to_array
@unpack_first_arg
def cartesian_to_polar(x, y, z, deg=False):
    "Convert 3D cartesian coordinates to polar (radius, inclination, azimuth)"
    r = np.sqrt(x**2 + y**2 + z**2)
    inclination = np.arccos(z / r)

    azimuth = np.arctan2(y, x)  # result in [-pi, pi]

    # Conversion to degrees if needed
    if deg:
        return r, RAD2DEG * inclination, RAD2DEG * azimuth
    else:
        return r, inclination, azimuth


@unpack_first_arg
def cartesian_to_latitude_longitude(x, y, z, deg=False):
    "Convert 3D cartesian coordinates to longitude and lattitude"
    r, inclination, azimuth = cartesian_to_polar(x, y, z)
    latitude = np.pi / 2 - inclination

    if deg:
        return RAD2DEG * latitude, RAD2DEG * azimuth
    else:
        return latitude, azimuth
