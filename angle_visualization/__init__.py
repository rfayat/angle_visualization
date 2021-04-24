"Useful functions for 3D plots"
import numpy as np
import mpl_toolkits.mplot3d.axes3d as ax3d
from .utils import grab_current_axis
from .angle_utils import cartesian_to_latitude_longitude
import cartopy.crs as ccrs
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
from collections.abc import Iterable
LAEA = ccrs.LambertAzimuthalEqualArea(central_latitude=90)
TRANSFORM = ccrs.PlateCarree()


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
def fill_projected_faces(lat_face_all, lon_face_all,
                         face_colors=None, edge_colors=None,
                         ax=None, **kwargs):
    """Plot a list of faces from geographic coordinates using a LAEA projection

    The faces are plotted as a matplotlib PatchCollection, each one of the
    patches being genereted using matplotlib's `fill` function.

    Parameters
    ----------
    lat_face_all: array, shape = (n_faces, n_vertices)
        The latitudes of the `n_faces` vertices for the `n_faces` faces in deg

    lon_face_all: array, shape = (n_faces, n_vertices)
        The longitudes of the `n_faces` vertices for the `n_faces` faces in deg

    face_colors: string, iterable of len n_faces or None (default: None)
        The color for the faces, if face_colors is an iterable of length
        `n_faces`, each element of face_colors is passed to `plt.fill` color.
        Else, face_colors is used for generating each patch of the
        PatchCollection.

    ax: matplotlib axis (default: None)
        Matplotlib axis with a cartopy LambertAzimuthalEqualArea projection. If
        None is provided (default), the current axis is grabed.

    **kwargs
        Additional key-word arguments passed to matplotlib's fill function for
        each facet.

    Returns
    -------
    patch_collection: matplotlib PatchCollection
        The resulting PatchCollection

    Notes
    -----
    The patches overlapping the south pole and the prime meridian are ignored
    due to rendering issue. The criterion for ignoring such facets is having
    one vertex with latitude inferior to -85 deg and one vertex whose longitude
    has an absolute value inferior to 5 deg.

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
    n_faces = len(lat_face_all)

    if kwargs.get("facecolor") is not None:
        face_colors = [kwargs.get("facecolor")] * n_faces
    elif kwargs.get("fc") is not None:
        face_colors = [kwargs.get("fc")] * n_faces
    elif face_colors is None or (isinstance(face_colors, Iterable) and len(face_colors)) != n_faces:  # noqa 501
        face_colors = [face_colors] * n_faces

    # Replicate the edge color if needed or use the face color
    if kwargs.get("ec") is not None:
        edge_colors = [kwargs.get("ec")] * n_faces
        kwargs.pop("ec")
    elif kwargs.get("edgecolor") is not None:
        edge_colors = [kwargs.get("edgecolor")] * n_faces
        kwargs.pop("edgecolor")
    elif isinstance(edge_colors, Iterable) and len(edge_colors) != n_faces:
        edge_colors = [edge_colors] * n_faces
    elif edge_colors is None:
        edge_colors = face_colors


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
                                fc=face_colors[face_idx],
                                ec=edge_colors[face_idx],
                                **kwargs)
    patch_collection = PatchCollection(patches)
    ax.add_collection(patch_collection)
    return patch_collection


def fill_projected_faces_euclidean(x, y, z, **kwargs):
    """Plot a list of faces from euclidean coordinates using a LAEA projection

    The x, y, z euclidean coordinates of the faces are converted to geographic
    coordinates before using `fill_projected_faces` to plot them. See this
    function's documentation for more information about the function's use.

    Parameters
    ----------
    x, y, z, arrays, shape = (n_points, n_vertices)
        The 3-dimensional euclidean coordinates of the vertices of the faces
        that will be plotted.

    **kwargs
        Additional key-word arguments passed to `fill_projected_faces`.

    Returns
    -------
    patch_collection: matplotlib PatchCollection
        The resulting PatchCollection

    Example:
    -------
    import matplotlib.pyplot as plt
    import numpy as np
    import angle_visualization
    from angle_visualization.triangulation import Delaunay_Sphere

    # Generate the sphere faces
    n_points = 300
    d = Delaunay_Sphere(n_points)
    random_colors = [np.random.choice(["r","g","b"]) for _ in d.faces]

    # Create a figure and axis with Lambert Azimuthal Equal Area projection
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=angle_visualization.LAEA)
    ax.gridlines()

    # Fill the faces of the triangulated sphere
    angle_visualization.fill_projected_faces_euclidean(
        *d.faces.transpose(2, 0, 1), face_colors = random_colors, alpha=.3
    )

    """
    # Convert the euclidean coordinates to geographic coordinates
    latlon_face_all = np.array(
        [cartesian_to_latitude_longitude(xf, yf, zf, deg=True) for (xf, yf, zf) in zip(x, y, z)]  # noqa E501
    )
    # Plot the result using fill_projected_faces
    return fill_projected_faces(latlon_face_all[:, 0, :],
                                latlon_face_all[:, 1, :],
                                **kwargs)


@grab_current_axis
def plot_projected(lon, lat, ax=None, **kwargs):
    "Plot using a crs projection."
    projected = LAEA.transform_points(TRANSFORM, lon, lat)
    return ax.plot(*projected.T, transform=LAEA, **kwargs)


def plot_projected_euclidean(x, y, z, **kwargs):
    "Plot from euclidean coordinates using a crs projection."
    lat, lon = cartesian_to_latitude_longitude(x, y, z, deg=True)
    return plot_projected(lon, lat, **kwargs)

@grab_current_axis
def scatter_projected(lon, lat, ax=None, **kwargs):
    "Scatterplot using a crs projection."
    return ax.scatter(lon, lat, transform=TRANSFORM, **kwargs)


def scatter_projected_euclidean(x, y, z, **kwargs):
    "Scatterplot from euclidean coordinates using a crs projection."
    lat, lon = cartesian_to_latitude_longitude(x, y, z, deg=True)
    return scatter_projected(lon, lat, **kwargs)
