# Angle visualization
A few python tools for visualizing 3D angles, shapes and histograms on the unit sphere with matplotlib.

## Installation
### Installation with pip
Using a python>=3.6 environment, simply clone the repository and install it with pip:

```bash
$ git clone https://github.com/rfayat/angle_visualization.git
$ cd angle_visualization
$ pip install .
```

### Troubleshooting
The installation pipeline was tested on ubuntu 20.04 and should be working fine on linux. The module [cartopy](https://scitools.org.uk/cartopy/docs/latest/index.html) is however a dependency for plotting 2D projections of the sphere and might cause hurdles during installation especially on windows and MacOS. Try following [cartopy's installation guidelines](https://scitools.org.uk/cartopy/docs/latest/installing.html) if it is the case. One potential source of issues is the dependency of cartopy on [GEOS](https://trac.osgeo.org/geos), an open source C++ geometry engine which will probably need to be installed independently.

### Figure generation
To generate the figures shown below, install the additional requirement (seaborn), uncomment the desired function at the end of [this script](README_figures/generate_figures.py) (after `if __name__=="__main__":`) and run:
```bash
$ python -m README_figures.generate_figures
```
## Delaunay triangulation of a Fibonacci sphere
### Complete Delaunay triangulation for 3D shapes

An extension of [scipy's Delaunay triangulation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.Delaunay.html) routine, (`angle_visualization.triangulation.Delaunay_Complete`), allows to obtain a triangulation of closed 3D objects. The triangulation can be performed using a 2D array containing the x, y and z coordinates of the vertices with shape `(n_vertices, 3)`:

```python
# Generate the coordinates of a Fibonacci sphere with n_points points
from angle_visualization.triangulation import Delaunay_Complete
# vertices is an array of shape (n_vertices, 3)
triangulated = angle_visualization.triangulation.Delaunay_Complete(vertices)
```

A comparison of the triangulation resulting from scipy's routine and its extension is shown below. The corresponding code can be found in [generate_figures.py](README_figures/generate_figures.py) (`plot_triangulated_cube` function):

![Delaunay_cube](README_figures/triangulated_cube.png)


### Triangulated Fibonacci sphere
#### Fibonacci sphere

The euclidean coordinates of a Fibonacci sphere can be generated using `angle_visualization.triangulation.fibonacci_sphere`,

```python
# Generate the coordinates of a Fibonacci sphere with n_points points
from angle_visualization.triangulation import fibonacci_sphere
n_points = 1000
x, y, z = fibonacci_sphere(n_points)
```

![Fibonacci_Sphere](README_figures/fibonacci_sphere.png)


#### Generating a triangulated Fibonacci sphere
A triangulated Fibonacci sphere with `n_points=3000` vertices can be generated as follows:  
```python
from angle_visualization.triangulation import Delaunay_Sphere
# Delaunay triangulation of the sphere
n_points = 3000
d = Delaunay_Sphere(n_points)
```
It has among others properties storing the edges (`Delaunay_Sphere.edges`), faces (`Delaunay_Sphere.faces`) and face centroids (`Delaunay_Sphere.face_centroids`) coordinates as numpy arrays.

#### Visualization

The facets of the Delaunay sphere can be visualized using `angle_visualization.plot_faces` as follows:
```python
import angle_visualization
import angle_visualization.triangulation
# Create a triangulated Fibonacci sphere with 3000 vertices
n_points = 3000
d = angle_visualization.triangulation.Delaunay_Sphere(n_points)

# Plot the result
fig = plt.figure(figsize=(4, 4))
ax = fig.add_subplot(111, projection='3d')
angle_visualization.plot_faces(d.faces,
                               ax=ax,
                               linewidths=...,
                               face_colors=...,
                               edge_colors=...)
```
An example for plotting faces with colors depending on the face centroids' coordinates can be found in [generate_figures.py](README_figures/generate_figures.py) (`plot_triangulated_sphere` function):

![Triangulated_Fibonacci_Sphere](README_figures/triangulated_fibonacci_sphere.png)



## Angle histogram on a sphere
### 3D angle histogram on a Fibonacci sphere
The `Delaunay_Sphere` can allocate each element of a 2D array with shape `(n_arrays, 3)` where axis 1 corresponds to the x, y and z coordinates of vectors on the unit sphere as follows:

```python
import angle_visualization
import angle_visualization.triangulation
# Create a triangulated Fibonacci sphere with 3000 vertices
n_points = 3000
d = angle_visualization.triangulation.Delaunay_Sphere(n_points)
# Compute the histogram on the sphere for the input array data
# data is a 2D array of unit 3D vectors with shape (n_arrays, 3)
triangle_idx_all, counts = d.spherical_histogram(data)
```

After converting the counts to a colorscale, the resulting histogram on the unit sphere can be visualized using the `angle_visualization.plot_faces` function. An example can be found in [generate_figures.py](README_figures/generate_figures.py) (`plot_histogram_3D` function):


![histogram_3D](README_figures/histogram_3D.png)
### 2D projection of the sphere
Leveraging [cartopy](https://scitools.org.uk/cartopy/docs/latest/)'s functionalities, the resulting histogram on the unit sphere can be projected in 2D using the `angle_visualization.fill_projected_faces_euclidean` function. For instance, using the outputs of `Delaunay_Sphere.spherical_histogram`, a [Lambert Azimuthal Equal-Area (LAEA)](https://en.wikipedia.org/wiki/Lambert_azimuthal_equal-area_projection) projection of the histogram can be visualized as follows:

```python
# Create the projected axis using cartopy
fig = plt.figure(figsize=(4, 4))
projection = angle_visualization.LAEA  # Lambert Azimuthal Equal-Area projection  # noqa
ax = fig.add_subplot(1, 1, 1, projection=projection)
ax.gridlines()

# Plot the projected faces
angle_visualization.fill_projected_faces_euclidean(
    *d.faces.transpose(2, 0, 1), face_colors=..., alpha=1
)
```

An example can be found in [generate_figures.py](README_figures/generate_figures.py) (`plot_projected_histogram_3D` function) and results in:

![projected_histogram_3D](README_figures/projected_histogram_3D.png)

## Other 3D-shapes
### Sphere
```python
# Plot a sphere
sphere = angle_visualization.shapes.Sphere3D.add_sphere(
    radius=1.,
    origin=[0, 0, 0],
    steps=30,
    ax=ax, color="blue", alpha=.7,
)
```


![sphere](README_figures/sphere.png)
### Arrows

`angle_visualization.shapes.Arrow3D` inherits from matplotlib's `FancyArrowPatch`. Plotting a single arrow an xyz origin can be done as follows:
```python
# Plot a single 3-dimensional arrow
arrow = angle_visualization.shapes.Arrow3D.add_arrow(
    2, .3, 3,
    origin=[-.5, 0, -1],
    color="purple",
    mutation_scale=20,
    arrowstyle="-|>",
    lw=5,
    adjust_ax_lim=False,
    ax=ax
)

# Plot an xyz origin in 3D
xyz = angle_visualization.shapes.Arrow3D.add_xyz(
    V=np.eye(3),
    origin=[0, 0, 0],
    mutation_scale=10,
    arrowstyle="-|>",
    lw=3,
    adjust_ax_lim=False,
    ax=ax
)
```

![arrows](README_figures/arrows.png)

### Rotations / Quaternions
**Not implemented yet.** The idea here would be to have a simple way to visualize unit quaternions / rotations (axis and magnitude) using matplotlib similarly to [WolframAlpha's representation of rotations](https://www.wolframalpha.com/input/?i=quaternion%3A+0%2B2i-j-3k&lk=3), taking as input a [scipy Rotation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.transform.Rotation.html) object.

If you feel like developing this functionality, don't hesitate opening a pull request or having a look at the (very) preliminary `angle_visualization.shapes.Quaternion3D` class, I'd be glad to help.
