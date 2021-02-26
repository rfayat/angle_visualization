"""
Objects used for plotting 3D shapes with matplotlib.

Author: Romain Fayat, February 2021
"""
import numpy as np
import quaternion
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa W0611
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
from copy import copy
from .utils import grab_current_axis, set_3Dlim


class Arrow3D(FancyArrowPatch):
    "3D arrow"

    def __init__(self, x, y, z, origin=[0, 0, 0], *args, **kwargs):
        "Create an Arrows3D from an origin and xyz coordinates"
        FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
        self._origin = origin
        self._vect = x, y, z

    @property
    def coord(self):
        "Return the coordinates of the arrow's extremities for each axis"
        x, y, z = self._vect
        xc = [self._origin[0], self._origin[0] + x]
        yc = [self._origin[1], self._origin[1] + y]
        zc = [self._origin[2], self._origin[2] + z]
        return xc, yc, zc

    @property
    def coord_range(self):
        "Return the range of values covered by the arrow on the different axis"
        xr, yr, zr = copy(self.coord)
        xr.sort()
        yr.sort()
        zr.sort()

        return xr, yr, zr

    def draw(self, renderer):
        xs, ys, zs = proj3d.proj_transform(*self.coord, renderer.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        FancyArrowPatch.draw(self, renderer)

    @grab_current_axis
    def extend_ax_lim(self, ax=None, margin=0):
        "Extend the axis limits so that the arrow is fully visible"
        # Grab the current axis limits and the range of values of the coord
        xlim, ylim, zlim = ax.get_xlim(), ax.get_ylim(), ax.get_zlim()
        xr, yr, zr = self.coord_range

        # Compute the new axis limits
        xlim_new = [min(xr[0] - margin, xlim[0]), max(xr[1] + margin, xlim[1])]
        ylim_new = [min(yr[0] - margin, ylim[0]), max(yr[1] + margin, ylim[1])]
        zlim_new = [min(zr[0] - margin, zlim[0]), max(zr[1] + margin, zlim[1])]

        # Set the axis limits
        set_3Dlim(xlim_new, ylim_new, zlim_new, ax=ax)

    @grab_current_axis
    def fit_ax_lim(self, ax=None, margin=0):
        "Fit the axis limits so that they are centered on the arrow"
        # Grab the range of values of the coordinates
        xr, yr, zr = self.coord_range

        # Add a margin around the range of values
        xlim_new = [xr[0] - margin, xr[1] + margin]
        ylim_new = [yr[0] - margin, yr[1] + margin]
        zlim_new = [zr[0] - margin, zr[1] + margin]

        # Set the axis limits
        set_3Dlim(xlim_new, ylim_new, zlim_new, ax=ax)

    @classmethod
    @grab_current_axis
    def add_arrow(cls, *args, ax=None, adjust_ax_lim=True, margin=0, **kwargs):
        "Create an Arrow3D object and add it to the provided axis"
        arrow = cls(*args, **kwargs)
        ax.add_artist(arrow)

        if adjust_ax_lim:
            arrow.extend_ax_lim(ax=ax, margin=margin)

        return arrow

    @classmethod
    def add_xyz(cls, V=np.eye(3), ax=None, colors="rgb", **kwargs):
        "Plot the columns of V with a common origin"
        arrow_all = []

        for i, v in enumerate(V.T):
            arrow = cls.add_arrow(*v,
                                  color=colors[i],
                                  ax=ax, **kwargs)
            arrow_all.append(arrow)

        return arrow_all


class Sphere3D():
    "Class for plotting spheres"

    def __init__(self, radius=1., origin=[0, 0, 0]):
        "Create a 3D sphere"
        self._radius = radius
        self._origin = origin

    def coord(self, steps=10):
        "Return the xyz coordinates of a sphere as 3 (2*steps, steps) arrays."
        u, v = np.mgrid[0:2 * np.pi:2 * steps * 1j, 0:np.pi:steps * 1j]
        x = self._radius * np.cos(u) * np.sin(v) + self._origin[0]
        y = self._radius * np.sin(u) * np.sin(v) + self._origin[1]
        z = self._radius * np.cos(v) + self._origin[2]
        return x, y, z

    @classmethod
    @grab_current_axis
    def add_sphere(cls, radius=1., origin=[0, 0, 0], steps=10, ax=None, **kwargs):  # noqa E501
        "Plot a sphere as a surface and return it."
        sphere = cls(radius, origin)
        x, y, z = sphere.coord(steps=steps)
        ax.plot_surface(x, y, z, **kwargs)
        return sphere


class Quaternion3D():
    "Class for plotting rotations / vectors represented by quaternions"

    def __init__(self, q, tol=1e-7):
        "Initialize the object from a np.quaternion object"
        self._q = q
        self.is_vect = np.abs(self._q.real) < tol
        self.is_rot = not self.is_vect and np.abs(self._q.norm() - 1) < tol
        if not self.is_vect and not self.is_rot:
            raise ValueError("Expected quaternion representing a vector or a rotation.")  # noqa E501

        if self.is_rot:
            self._angle = self._q.angle

    def __repr__(self):
        "Print the quaternion"
        return self._q.__repr__()

    @classmethod
    @grab_current_axis
    def add_quaternion(cls, q, radius=1., ax=None, tol=1e-7, origin=[0, 0, 0]):
        quat = cls(q, tol=tol)
        if quat.is_vect:
            quat._arrow = Arrow3D.add_arrow(*quat._q.vec, ax=ax, color="k",
                                            lw=3, arrowstyle="-|>",
                                            mutation_scale=20)
        elif quat.is_rot:
            vec_normalized = 2 * radius * quat._q.vec / np.linalg.norm(quat._q.vec)  # noqa E501
            origin_arrow = [o - v / 2 for o, v in zip(origin, vec_normalized)]
            quat._arrow = Arrow3D.add_arrow(*vec_normalized, ax=ax, color="k",
                                            lw=3, arrowstyle="-|>",
                                            mutation_scale=20,
                                            origin=origin_arrow)
            quat._sphere = Sphere3D.add_sphere(origin=origin, radius=radius,
                                               ax=ax, color="grey", alpha=.5)

        return quat


if __name__ == "__main__":
    # Create a rotation quaternion
    angle = np.pi / 6
    angle_cos = np.cos(angle / 2)
    orientation = np.array([1., .25, .75])
    normalization = np.sqrt(((1 - angle_cos**2) / (orientation**2).sum()))
    orientation_normalized = normalization * orientation

    r = np.quaternion(angle_cos, *orientation_normalized)
    M = quaternion.as_rotation_matrix(r)

    fig = plt.figure(figsize=(15, 15))
    ax = fig.add_subplot(111, projection='3d')
    Quaternion3D.add_quaternion(q=r, origin=[.5, .5, 0])
    ax.set_xlabel('x_values')
    ax.set_ylabel('y_values')
    ax.set_zlabel('z_values')
    plt.draw()
    fig.show()

    str(input("Press enter to stop."))
