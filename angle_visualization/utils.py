"Utils functions."
import matplotlib.pyplot as plt
import numpy as np
from functools import wraps


def unpack_first_arg(f):
    "Treat the second dimension of the first arg as independent inputs"

    @wraps(f)
    def g(*args, **kwargs):
        if len(args) == 1 and isinstance(args[0], np.ndarray):
            return f(*(args[0].T), **kwargs)
        else:
            return f(*args, **kwargs)
    return g


def get_integer_dtype(max_value, min_value=0):
    "Return the best numpy dtype for integers with input max value"
    if min_value >= 0:
        if max_value < 2**8:
            return np.uint8
        elif max_value < 2**16:
            return np.uint16
        elif max_value < 2**32:
            return np.uint32
        else:
            return np.uint64
    else:
        max_value = max(max_value, abs(min_value))
        if max_value < 2**7:
            return np.int8
        elif max_value < 2**15:
            return np.int16
        elif max_value < 2**31:
            return np.int32
        else:
            return np.int64


def store_value_on_first_computation(f):
    "Compute and store the result of a method its first call"

    @wraps(f)
    def g(self, *args, **kwargs):
        "Store the value of f on the first call and return it"
        method_name = f.__name__
        stored_result_name = "_" + method_name

        if getattr(self, stored_result_name, None) is None:
            setattr(self, stored_result_name, f(self, *args, **kwargs))

        return getattr(self, stored_result_name)

    return g


def grab_current_axis(f):
    "Grab the current axis if None is provided"

    @wraps(f)
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
