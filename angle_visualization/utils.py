"Utils functions."
import matplotlib.pyplot as plt


def grab_current_axis(f):
    "Grab the current axis if None is provided"

    def g(*args, **kwargs):
        "Call f after grabbing the current axis if None is provided"
        if kwargs.get("ax") is None:
            kwargs.update({"ax": plt.gca()})
        return f(*args, **kwargs)

    return g


@grab_current_axis
def set_3Dlim(xlim=None, ylim=None, zlim=None, ax=None):
    "Set the x, y and z axis limits of a 3D axis"
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    if zlim is not None:
        ax.set_zlim(zlim)
    return ax
