"Useful functions for 3D plots"
import numpy as np
import mpl_toolkits.mplot3d.axes3d as ax3d
from .utils import grab_current_axis
import cartopy.crs as ccrs
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
LAEA = ccrs.LambertAzimuthalEqualArea(central_latitude=90)


@grab_current_axis
def plot_faces(faces, face_colors=None, edge_colors=None, ax=None, **kwargs):
    """Plot a list of faces using Poly3DCollection.

    Parameters
    ----------
    faces: list
        List of the coordinates of the faces to plot, each element of `faces`
        is a numpy array with shape (n_points, 3).

    face_colors: color or sequence of colors (default: None)
        Set the facecolor(s) of the collection. face_colors can be a color (all
        patches have same color), or a sequence of colors; if it is a sequence
        the patches will cycle through the sequence.

    edge_colors: color or sequence of colors or 'face'
        The collection edgecolor(s).  If a sequence, the patches cycle
        through it.  If 'face', match the facecolor.

    ax: mplot3d.axes3d.Axes3D
        Axis to use fort the plot, by default grab the current axis.

    **kwargs:
        Additional key-word arguments that will be passed to Poly3DCollection.

    Returns
    -------
    tri: mplot3d.art3d.Poly3DCollection
        The resulting matplotlib Poly3DCollection.

    """
    tri = ax3d.art3d.Poly3DCollection(faces, **kwargs)
    if face_colors is not None:
        tri.set_facecolor(face_colors)
    if edge_colors is not None:
        tri.set_edgecolor(edge_colors)
    ax.add_collection3d(tri)
    return tri


@grab_current_axis
def fill_projected_faces(lat_face_all, lon_face_all, face_colors=None,
                         ax=None, **kwargs):
    """Plot a list of faces using a Lambert azimuthal equal area projection

    Example:
    -------
    import matplotlib.pyplot as plt
    import numpy as np
    import angle_visualization
    from angle_visualization.angle_utils import cartesian_to_latitude_longitude
    from angle_visualization.triangulation import Delaunay_Sphere
    import cartopy.crs as ccrs

    # Generate the sphere faces
    n_points = 300
    d = Delaunay_Sphere(n_points)
    latlon_face_all = np.array([cartesian_to_latitude_longitude(f, deg=True) for f in d.faces])  # noqa E501
    random_colors = [np.random.choice(["r","g","b"]) for _ in d.faces]

    # Create a figure and axis with Lambert Azimuthal Equal Area projection
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=angle_visualization.LAEA)
    ax.gridlines()

    # Fill the faces of the triangulated sphere
    angle_visualization.fill_projected_faces(latlon_face_all[:, 0, :],
                                             latlon_face_all[:, 1, :],
                                             face_colors = random_colors,
                                             alpha=.3)

    # Add a scatterplot of all points
    lat, lon = cartesian_to_latitude_longitude(d._points, deg=True)
    ax.scatter(lon, lat, transform=ccrs.PlateCarree())

    """
    # Replicate the input face_color if needed
    if face_colors is not None and len(face_colors) != len(lat_face_all):
        face_colors = [face_colors for _ in lat_face_all]

    # Loop over the coordinates of the faces
    patches = []  # list of patches that will be converted to a PatchCollection
    for face_idx, (lat, lon) in enumerate(zip(lat_face_all, lon_face_all)):
        # Remove the points with longitude 0 at the south pole
        if not (np.any(lat <= -85) and np.any(np.abs(lon) <= 5)):
            filling_transform = ccrs.Geodetic()
            # Duplicate the 1st element of the coordinates to close the shapes
            lat, lon = np.append(lat, lat[0]), np.append(lon, lon[0])
            # Add the resulting patch to the list of patches
            patches += plt.fill(lon, lat, transform=filling_transform,
                                color=face_colors[face_idx], **kwargs)
    ax.add_collection(PatchCollection(patches))
