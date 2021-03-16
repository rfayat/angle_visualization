"Useful functions for 3D plots"
import mpl_toolkits.mplot3d.axes3d as ax3d
from .utils import grab_current_axis


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
