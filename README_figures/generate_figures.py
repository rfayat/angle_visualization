"""Script for generating the figures for angle_visualization' readme.

Author: Romain Fayat, April 2021
"""
import matplotlib.pyplot as plt
import numpy as np
import scipy.spatial
import angle_visualization
import angle_visualization.utils
import angle_visualization.triangulation
from angle_visualization.angle_utils import cartesian_to_latitude_longitude
from pathlib import Path
import seaborn as sns
import pandas as pd

OUTPUT_DIR = Path(__file__).parent.absolute()


def plot_triangulated_cube(*args, **kwargs):
    "3D plot of a cube triangulated using Delaunay_Complete"
    # Generate the cube's coordinates
    c = [0, 1]  # Verticies boundaries
    coord = np.array([[i, j, k] for i in c for j in c for k in c])
    # Triangulation using scipy
    cube = scipy.spatial.Delaunay(coord)
    cube_faces = np.array([cube._points[e] for e in cube.simplices[:, 1:]])
    # Triangulation using Delaunay_Complete
    cube_complete = angle_visualization.triangulation.Delaunay_Complete(coord)

    fig = plt.figure(figsize=(5, 3))
    ax = fig.add_subplot(121, projection='3d')
    angle_visualization.plot_faces(cube_faces, ax=ax, edge_colors="k",
                                   linewidths=1, alpha=.8)
    ax.axis('off')
    ax.set_title("scipy.spatial's\nDelaunay", size=8)

    ax = fig.add_subplot(122, projection='3d')
    angle_visualization.plot_faces(cube_complete.faces, ax=ax,
                                   edge_colors="k", linewidths=1, alpha=.8)
    ax.axis('off')
    ax.set_title("angle_visualization's\nDelaunay_complete", size=8)

    fig.suptitle("Triangulated cube using Delaunay triangulation", size=12)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.savefig(*args, **kwargs)
    plt.close(fig)


def scatterplot_fibonacci_sphere(*args, **kwargs):
    "3D scatter plot of the euclidean coordinates of a fibonacci sphere"
    # Generate the coordinates of a fibonacci sphere with n_points points
    n_points = 1000
    x, y, z = angle_visualization.triangulation.fibonacci_sphere(n_points)
    # Display points in a scatter plot
    fig = plt.figure(figsize=(4, 4))
    ax = fig.add_subplot(111, projection='3d')
    ax.view_init(0, 90)  # orientation of the sphere, y axis toward us
    # Scatter plot of the sphere
    colors = np.cos(np.pi * x) * np.sin(2 * np.pi * z)  # dummy colors
    size = 50 * (1 + y)  # dummy size of the markers
    ax.scatter(x, y, z, c=colors, s=size, alpha=1)
    # Change the axes limits
    angle_visualization.utils.set_3Dlim(*[[-.8, .8] for _ in range(3)], ax=ax)
    ax.axis('off')
    ax.set_title(f"Scatterplot of a Fibonacci\nsphere with {n_points} points")
    fig.tight_layout()
    fig.savefig(*args, **kwargs)
    plt.close(fig)


def plot_triangulated_sphere(*args, **kwargs):
    "Plot the faces of a triangulated fibonacci sphere"
    # Delaunay triangulation of the sphere
    n_points = 3000
    d = angle_visualization.triangulation.Delaunay_Sphere(n_points)

    # Compute dummy colors from the triangles' centroids
    # Normalization to get the projection on the sphere
    center_all = d.face_centroids
    cmap = sns.mpl_palette("viridis", as_cmap=True)
    colors = cmap(np.sin(center_all[:, 1] * np.pi) *
                  np.sin(center_all[:, 2] * np.pi * 2))

    fig = plt.figure(figsize=(4, 4))
    ax = fig.add_subplot(111, projection='3d')

    angle_visualization.plot_faces(d.faces, ax=ax, linewidths=1,
                                   face_colors=colors, edge_colors=colors)

    angle_visualization.utils.set_3Dlim(*[[-.8, .8] for _ in range(3)], ax=ax)
    # Final tweaking and save the figure
    ax.axis('off')
    ax.set_title("Triangulated Fibonacci sphere")
    fig.tight_layout()
    fig.savefig(*args, **kwargs)
    plt.close(fig)


def plot_3D_histogram(data, *args, **kwargs):
    "Plot a 3D histogram for input 3D euclidean coordinates"
    # Create a Delaunay sphere and compute the 3D histogram
    n_points = 5000
    d = angle_visualization.triangulation.Delaunay_Sphere(n_points)
    triangle_idx_all, counts = d.spherical_histogram(data.values,
                                                     block_size=50000)
    # Create the colormap for each face of the sphere (here log scale)
    cmap = sns.mpl_palette("viridis", as_cmap=True)
    log_count = np.full(len(counts), -1.)
    log_count[counts != 0] = np.log(counts[counts != 0]) / np.log(counts[counts != 0]).max()  # noqa
    colors = cmap(log_count)
    # Create the figure and rotate the axis
    fig = plt.figure(figsize=(4, 4))
    ax = fig.add_subplot(111, projection='3d')
    ax.view_init(90, 0)
    # Plot the faces of the sphere
    angle_visualization.plot_faces(d.faces, ax=ax, linewidths=1,
                                   face_colors=colors, edge_colors=colors)
    angle_visualization.utils.set_3Dlim(*[[-.8, .8] for _ in range(3)], ax=ax)
    # Final tweaking and save the figure
    ax.axis("off")
    ax.set_title("3D histogram for an example time series")
    fig.tight_layout()
    fig.savefig(*args, **kwargs)
    plt.close(fig)


def plot_projected_3D_histogram(data, *args, **kwargs):
    "Plot the 2D projection of a 3D histogram"
    # Create a Delaunay sphere and compute the 3D histogram
    n_points = 5000
    d = angle_visualization.triangulation.Delaunay_Sphere(n_points)
    triangle_idx_all, counts = d.spherical_histogram(data.values,
                                                     block_size=50000)
    # Create the colormap for each face of the sphere (here log scale)
    cmap = sns.mpl_palette("viridis", as_cmap=True)
    log_count = np.full(len(counts), -1.)
    log_count[counts != 0] = np.log(counts[counts != 0]) / np.log(counts[counts != 0]).max()  # noqa
    colors = cmap(log_count)

    # Create the projected axis using cartopy
    fig = plt.figure(figsize=(4, 4))
    projection = angle_visualization.LAEA  # Lambert Azimuthal Equal Area projection  # noqa
    ax = fig.add_subplot(1, 1, 1, projection=projection)
    ax.gridlines()

    # Plot the projected faces
    latlon_face_all = np.array([
        cartesian_to_latitude_longitude(f, deg=True) for f in d.faces
    ])
    angle_visualization.fill_projected_faces(latlon_face_all[:, 0, :],
                                             latlon_face_all[:, 1, :],
                                             face_colors=colors, alpha=1)

    # Final tweaking and save the figure
    ax.axis("off")
    ax.set_title("3D histogram projected using LAEA projection")
    fig.tight_layout()
    fig.savefig(*args, **kwargs)
    plt.close(fig)


if __name__ == "__main__":
    sample_data = pd.read_csv("README_figures/sample_data.csv")
    # Plot of a triangulated cube
    p = OUTPUT_DIR / "triangulated_cube.png"
    plot_triangulated_cube(p, transparent=True)
    # Scatterplot of the euclidean coordinates of a fibonacci sphere
    p = OUTPUT_DIR / "fibonacci_sphere.png"
    scatterplot_fibonacci_sphere(p, transparent=True)
    # Plot the faces of a triangulated Fibonacci sphere
    p = OUTPUT_DIR / "triangulated_fibonacci_sphere.png"
    plot_triangulated_sphere(p, transparent=True)
    # Plot a 3D histogram on a triangulated sphere
    p = OUTPUT_DIR / "histogram_3D.png"
    plot_3D_histogram(sample_data, p, transparent=True)
    # Plot a projected 3D histogram
    p = OUTPUT_DIR / "projected_histogram_3D.png"
    plot_projected_3D_histogram(sample_data, p, transparent=True)
