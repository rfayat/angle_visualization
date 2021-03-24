"Triangulation for points of a closed shape in 3D"
import scipy
import scipy.spatial
import numpy as np
from . import utils
from .utils import store_value_on_first_computation


def fibonacci_sphere(n_points: int):
    "Return the coordinates of the points of a Fibonacci sphere"
    ga = (3 - np.sqrt(5)) * np.pi  # golden angle

    # Create a list of golden angle increments along the range of n_points
    theta = ga * np.arange(n_points)

    # Z is a split into a range of -1 to 1 in order to create a unit circle
    z = np.linspace(1 / n_points - 1, 1 - 1 / n_points, n_points)

    # a list of the radii at each height step of the unit circle
    radius = np.sqrt(1 - z * z)

    # Determine where xy fall on the sphere from the azimuthal and polar angles
    x = radius * np.cos(theta)
    y = radius * np.sin(theta)

    return np.vstack((x, y, z))


def allocate_spherical_triangle(inv, order_all, V):
    """Allocate 3D-vectors to a spherical triangle from precomputed variables

    Parameters
    ----------
    inv: array, shape = (n_triangles, 3, 3)
        For each triangle, the inverse of (3, 3) matrix corresponding to the
        verticies coordinates as columns.

    order_all: array of integers, shape = (n_triangles, n_triangles)
        For each triangle (rows), the indexes for the triangles sorted by
        centroid distances.

    V: array, shape = (n_points, 3)
        The range of 3-dimensional arrays that need to be mapped on the sphere.

    Returns
    -------
    triangle_idx_all: array of integers, shape = (n_points)
        For each element of V, the index of the spherical triangle it belongs
        to.

    """
    triangle_idx_all = np.zeros(len(V), dtype=np.uint)
    triangle_idx_t_minus_1 = 0  # preallocation
    for i, v in enumerate(V):
        # The triangles will be tested starting at the closest of the
        # triangle allocated a t - 1
        order = order_all[triangle_idx_t_minus_1]

        # Loop over all spherical triangles according to order
        for triangle_idx, inverted_triangle_coordinates in zip(order, inv[order]):  # noqa E501
            if np.all(np.dot(inverted_triangle_coordinates, v) >= 0):
                triangle_idx_t_minus_1 = triangle_idx
                triangle_idx_all[i] = triangle_idx
                break  # end the loop when the matching triangle has been found
    return triangle_idx_all


def allocate_spherical_triangle_block(inv, order_all, V, block_size=10000):
    """Allocate 3D-vectors to a spherical triangle from precomputed variables

    Parameters
    ----------
    inv: array, shape = (n_triangles, 3, 3)
        For each triangle, the inverse of (3, 3) matrix corresponding to the
        verticies coordinates as columns.

    order_all: array of integers, shape = (n_triangles, n_triangles)
        For each triangle (rows), the indexes for the triangles sorted by
        centroid distances.

    V: array, shape = (n_points, 3)
        The range of 3-dimensional arrays that need to be mapped on the sphere.

    block_size: int, default = 10000
        Number of arrays of V to proces simultaneously. Dot products between
        a (3, 3) and a (3, block_size) arrays will be performed when looking
        the matching triangle of the elements of V, taking advantage of
        numpy's parallelization.

    Returns
    -------
    triangle_idx_all: array of integers, shape = (n_points)
        For each element of V, the index of the spherical triangle it belongs
        to.

    """
    triangle_idx_all = np.zeros(len(V), dtype=np.uint)

    # Preallocation using the first data point
    idx = allocate_spherical_triangle(inv, order_all, V[[0]])
    triangle_idx_t_minus_1 = idx[0]

    block_start = np.arange(0, len(V), block_size)
    block_end = np.append(block_start[1:], len(V))

    for start, stop in zip(block_start, block_end):
        # The triangles will be tested starting at the closest of the
        # triangle allocated a t - 1
        order = order_all[triangle_idx_t_minus_1]
        allocated = np.zeros(stop - start, dtype=bool)

        # Loop over all spherical triangles according to order
        for triangle_idx, inverted_triangle_coordinates in zip(order, inv[order]):  # noqa E501
            # Check if the vectors to allocate belong to the triangle
            to_allocate = np.argwhere(~allocated).flatten() + start
            prod = np.dot(inverted_triangle_coordinates, V[to_allocate].T)
            good_triangle = np.all(prod >= 0., axis=0)
            # Change allocated to True for the vectors that have been allocated
            allocated[~allocated] = good_triangle
            # Store the triangle's index for the vectors that belong to it
            triangle_idx_all[to_allocate[good_triangle]] = triangle_idx
            # Break the loop once all vectors have been assigned to a triangle
            if np.all(allocated):
                break

        # Update the triangle for sorting using the last value of the block
        triangle_idx_t_minus_1 = triangle_idx_all[stop - 1]
    return triangle_idx_all


class Delaunay_Complete(scipy.spatial.Delaunay):
    "Extends Delaunay triangulation for closed 3D objects"

    def __init__(self, *args, **kwargs):
        "Delaunay tessellation in 3 dimensions for closed objects."
        super().__init__(*args, **kwargs)
        self.faces_indexes = self.get_faces_indexes()

    @property
    @store_value_on_first_computation
    def edges_indexes(self):
        """Compute and store and return the edge indexes

        This computation might be a bit long which is why we don't directly
        perform it when initializing the object.
        """
        return self.get_connected_pairs(self.faces_indexes)

    def get_faces_indexes(self):
        "Return an array of the indexes of the points to use for each face"
        faces_indexes = self.simplices[:, 1:]

        # find the points with missing triangles around the left-out point
        unique, counts = np.unique(self.simplices[:, 1:], return_counts=True)
        reference_point_neighbors = unique[counts == counts.min()]

        # add faces for the points already connected
        for i, p1 in enumerate(reference_point_neighbors):
            for p2 in reference_point_neighbors[i + 1:]:
                if self.points_are_connected(p1, p2, self.vertices[:, 1:]):
                    new_face_idx = [self.reference_point_idx, p1, p2]
                    faces_indexes = np.vstack((faces_indexes, new_face_idx))

        return faces_indexes

    @property
    @store_value_on_first_computation
    def faces(self):
        "Coordinate of the triangles including the ones for the left-out point"
        return np.array([self._points[e] for e in self.faces_indexes])

    @property
    @store_value_on_first_computation
    def edges(self):
        return np.array([self._points[e] for e in self.edges_indexes])

    @property
    @store_value_on_first_computation
    def face_centroids(self):
        "Return the average of the triangles' points coordinates"
        return np.vstack(tuple(e.mean(axis=0) for e in self.faces))

    @property
    @store_value_on_first_computation
    def face_centroids_sorted(self):
        "Return the indexes of the sorted pairwise distances of the centroids"
        # Compute and sort the pairwise distance matrix
        pairwise_dist = scipy.spatial.distance.pdist(self.face_centroids)
        pairwise_dist = scipy.spatial.distance.squareform(pairwise_dist)
        pairwise_dist_sorted = np.argsort(pairwise_dist, axis=1)

        # Change the dtype for the sake of memory
        dtype = utils.get_integer_dtype(len(pairwise_dist_sorted))
        pairwise_dist_sorted = pairwise_dist_sorted.astype(dtype)
        return pairwise_dist_sorted

    @property
    def reference_point_idx(self):
        "Return the index of the point which has no triangle"
        return self.simplices[0, 0]

    def get_connected_pairs(self, vertices):
        "Return all pairs of points from triplets"
        pairs_all = np.vstack((vertices[:, (0, 1)],
                               vertices[:, (0, 2)],
                               vertices[:, (1, 2)]))
        pairs_all = np.sort(pairs_all, axis=1)
        return np.unique(pairs_all, axis=0)

    def points_are_connected(self, idx1, idx2, vertices):
        "Indicate whether points with input indexes are connected"
        has_idx1 = np.any(vertices == idx1, axis=1)
        has_idx2 = np.any(vertices == idx2, axis=1)
        has_both = has_idx1 & has_idx2
        return np.any(has_both)


class Delaunay_Sphere(Delaunay_Complete):
    "A triangulated unit sphere computed using Fibonacci's formula"

    def __init__(self, npoints=300, *args, **kwargs):
        "Create a triangulated unit sphere with the given number of vertices"
        fib = fibonacci_sphere(npoints).T
        # Triangulation
        super().__init__(fib, *args, **kwargs)

    @property
    @store_value_on_first_computation
    def inverted_face_coordinates(self):
        """Return the invert of each triangle's verticies coordinates.

        More precisely, for each face, we compute the invert of the (3, 3)
        matrix containing the verticies as columns.

        """
        return np.linalg.inv(self.faces.transpose(0, 2, 1))

    def spherical_histogram(self, V, block_size=None):
        """Allocate each 3-dimensional vector in V to its face of the sphere

        Notes
        -----
        For the vector with timestamp t stored in V, we test all triangles
        starting by the neighbors of the triangle allocated to vector with
        timestamp t-1. The pipeline is therefore relatively efficient for
        continous timeseries but can become quite heavy if V[t-1] is not
        related to V[t].

        """
        inv = self.inverted_face_coordinates

        # Allocation of each 3D vector of V to its spherical triangle
        order_all = self.face_centroids_sorted
        if len(V) == 0:
            return np.array([]), np.zeros(len(self.faces), dtype=int)
        elif block_size is None:
            triangle_idx_all = allocate_spherical_triangle(inv, order_all, V)
        else:
            block_size = min(block_size, len(V))
            triangle_idx_all = allocate_spherical_triangle_block(
                inv, order_all, V, block_size=block_size
            )

        # Compute the histogram of values
        histogram = np.zeros(len(self.faces), dtype=np.int)
        unique, counts = np.unique(triangle_idx_all, return_counts=True)
        histogram[unique] = counts
        # Return the spherical triangle index and histogram
        return triangle_idx_all, histogram
