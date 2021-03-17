"Triangulation for points of a closed shape in 3D"
import scipy
import numpy as np


class Delaunay_Complete(scipy.spatial.Delaunay):
    "Extends Delaunay triangulation for closed 3D objects"

    def __init__(self, *args, **kwargs):
        "Delaunay tessellation in 3 dimensions for closed objects."
        super().__init__(*args, **kwargs)
        self.faces_indexes = self.get_faces_indexes()

    def get_faces_indexes(self):
        "Return an array of the indexes of the points to use for each face"
        faces_indexes = self.simplices[:, 1:]

        # find the points with missing triangles around the left-out point
        unique, counts = np.unique(self.simplices[:, 1:], return_counts=True)
        reference_point_neighbors = unique[counts == counts.min()]

        # add faces for the points already connected
        for i, p1 in enumerate(reference_point_neighbors):
            for p2 in reference_point_neighbors[i + 1:]:
                if self.points_are_connected(p1, p2):
                    new_face_idx = [self.reference_point_idx, p1, p2]
                    faces_indexes = np.vstack((faces_indexes, new_face_idx))

        return faces_indexes

    @property
    def faces(self):
        "Coordinate of the triangles including the ones for the left-out point"
        return [self._points[e] for e in self.faces_indexes]

    @property
    def face_centroids(self):
        "Return the average of the triangles' points coordinates"
        return np.vstack(tuple(e.mean(axis=0) for e in self.faces))

    @property
    def reference_point_idx(self):
        "Return the index of the point which has no triangle"
        return self.simplices[0, 0]

    def points_are_connected(self, idx1, idx2):
        "Indicate whether points with input indexes are connected"
        has_idx1 = np.any(self.vertices[:, 1:] == idx1, axis=1)
        has_idx2 = np.any(self.vertices[:, 1:] == idx2, axis=1)
        has_both = has_idx1 & has_idx2
        return np.any(has_both)
