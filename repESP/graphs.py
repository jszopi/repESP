import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_pdf import PdfPages
import os
import numpy as np
import random
import math

# This was necessary to prevent y-axis label from being cut off when plotting
# http://stackoverflow.com/a/17390833
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

from . import field_comparison

DIR_LABELS = ['x', 'y', 'z']


def _plot_common(dimension, title, guideline=False):
    """Set up plot of correct dimensionality and return related objects"""
    fig = plt.figure()
    if dimension == 3:
        ax = fig.add_subplot(111, projection='3d')
    elif dimension == 2:
        ax = fig.add_subplot(111)
    else:
        raise NotImplementedError("Plotting of dimension {0} not implemented"
                                  .format(dimension))

    # Add a horizontal line at 0 for 2D plots
    if guideline and dimension == 2:
        ax.axhline(color='k', linestyle='--')

    if title is not None:
        plt.title(title)

    return fig, ax


def plot(*fields, color=None, color_span=None, dist_field_filter=None,
         exclusion_dist=0, rand_skim=0.01, extra_filter=None, save_to=None,
         axes_limits=None, title=None, get_limits=None):

    assert 2 <= len(fields) <= 3

    # Pack fields and color together
    if color is not None:
        fields_and_color = list(fields) + [color]
    else:
        fields_and_color = fields

    field_comparison._check_grids(*fields_and_color)
    field_comparison._check_fields_for_nans(*fields_and_color)
    # Necessary, as the original Field will be overwritten when filtering
    dist_field_filter_type = dist_field_filter.field_type

    fig, ax = _plot_common(len(fields), title, guideline=True)
    _set_axis_labels(ax, *fields)

    # This function got really fat due to all that filtering and it can still
    # handle only one additional filter. Some refactoring is due. TODO
    if extra_filter is not None:
        fields_and_color = extra_filter(*fields_and_color)
        # This filtering step changes all Fields to np.arrays. As a result, in
        # the next filtering step, by dist_field_filter, a mixture of np.arrays
        # and Fields is passed, which is not handled by the filters. While that
        # deficiency was intentional, I don't think there's a reason it should
        # not be handled (TODO). But for now, a kludge:
        if dist_field_filter is not None:
            dist_field_filter = dist_field_filter.values
    if dist_field_filter is not None:
        if dist_field_filter_type != 'dist':
            print("WARNING: The field selected for filtering is not of type "
                  "'dist' but ", dist_field_filter.field_type)
        dist_field_filter, *fields_and_color = field_comparison.filter_by_dist(
            exclusion_dist, *([dist_field_filter] + fields_and_color))
    elif exclusion_dist:
        print("WARNING: exclusion distance specified but no Field passed to "
              "filter by.")
    fields_and_color = field_comparison.skim(rand_skim, *fields_and_color)
    fields_and_color = list(map(field_comparison._flatten_no_nans,
                            fields_and_color))

    if color is not None:
        cmap = _get_cmap(len(fields), color.field_type)
        cmap_name = color.lookup_name()
        *fields, color = fields_and_color
        # ax.scatter has to be inside of the 'color is not None' conditional
        # because an error occurs when the kwarg ``c`` is explicitly set to
        # None, even though it's the default value.
        vmin, vmax = color_span if color_span is not None else None, None
        image = ax.scatter(*fields, c=color, cmap=cmap, vmin=vmin, vmax=vmax)
        cbar = fig.colorbar(image, label=cmap_name)
    else:
        fields = fields_and_color
        ax.scatter(*fields)

    _set_axes_limits(len(fields), ax, axes_limits)
    _save_or_display(save_to)

    # Save limits to get_limits. This is useful when they are to be reused in
    # other plots. Saving the limits to an argument was more intuitive than
    # returning them.
    if get_limits is not None:
        # Code copied from _set_axes_limits (TODO: DRY)
        limits = []
        for dir_label in DIR_LABELS[:len(fields)]:
            # Get current limits
            limits.append(getattr(ax, "get_" + dir_label + "lim")())
        get_limits[:] = limits


def _set_axes_limits(dimension, ax, axes_limits):
    """Set axes limits"""
    if axes_limits is None:
        return
    # Smaller lengths are allowed, will be interpreted as the first few axes.
    # This should be an Exception not assertion though.
    assert len(axes_limits) <= dimension

    for axis_limits, dir_label in zip(axes_limits, DIR_LABELS):
        # Get current limits
        limits = list(getattr(ax, "get_" + dir_label + "lim")())
        for i, axis_limit in enumerate(axis_limits):
            if axis_limit is not None:
                limits[i] = axis_limit
        getattr(ax, "set_" + dir_label + "lim")(limits)

    # Although **not for my purposes at the moment** (I only want to set limits
    # so that different plots can be easily compared, so both axes will be
    # getting set), it would be nice to rescale the axes which were not
    # modified. However, when autoscaling, matplotlib always uses all the data.
    # ax.relim() with ax.autoscale_view() seemed to be relevant but they do not
    # easily operate on datapoints I think.


def _set_axis_labels(ax, *fields):
    """Set axis labels based on free-form names of Fields being plotted"""
    for field, dir_label in zip(fields, DIR_LABELS):
        getattr(ax, "set_" + dir_label + "label")(field.lookup_name())


def _get_cmap(dimension, field_type):
    """Return a color map based on plot dimensionality and field type"""
    if field_type == 'dist':
        if dimension != 3:
            # Shading by distance is more intuitive
            return plt.get_cmap('Blues_r')
        else:
            print("WARNING: Shading by distance doesn't look good on a 3D "
                  "plot. Colouring instead.")
    # This colour map is fairly intuitive, looks good in greysale and
    # to coloublind people (except for those with the very rare tritanomaly),
    # according to a simulation on:
    # http://www.color-blindness.com/coblis-color-blindness-simulator/
    return plt.get_cmap('plasma')
    # 'gnuplot2' was also considered but in greyscale it yields black and white
    # at the ends and IMO contains unnecessarily many hues. 'plasma' is simpler
    # and uses similar colours.


def _save_or_display(save_to=None):
    """Save the plot or display it if save_to is None"""
    if save_to is None:
        plt.show()
    else:
        if type(save_to) is PdfPages:
            # Need to check the type first, because it may be a file object if
            # a pdf is to be created, see:
            # http://matplotlib.org/faq/howto_faq.html#save-multiple-plots-to-one-pdf-file
            plt.savefig(save_to, format="pdf")
        elif os.path.isfile(save_to):
            raise FileExistsError("File exists: " + save_to)
        else:
            # DPI may need to be increased
            plt.savefig(save_to)
    plt.close()


def plot_points(points_field, dimension, title=None, color_span=None,
                axes_limits=None, save_to=None, rand_skim=1, plane_eqn=None,
                dist_thresh=None, molecule=None, atom_dist_threshs=None,
                atom_format=None, show_all_atoms=False):

    project_onto_plane = _check_args(dimension, plane_eqn, dist_thresh)
    field_comparison._check_fields_for_nans(points_field)
    fig, ax = _plot_common(dimension, title)

    # Skimming, filtering and projecting
    points, values = _points_dist_filter(
        points_field.points, points_field.values, plane_eqn, dist_thresh)
    points, values = _points_rand_skim(points, values, rand_skim)
    points = _project_points(points, project_onto_plane, dimension, plane_eqn)

    _plot_atoms(molecule, ax, dimension, plane_eqn, project_onto_plane,
                atom_dist_threshs, atom_format, show_all_atoms)

    cmap_name = points_field.lookup_name()
    cmap = plt.get_cmap('plasma')
    vmin, vmax = color_span if color_span is not None else None, None

    image = ax.scatter(*list(zip(*points))[:dimension], c=values,
                       cmap=cmap, vmin=vmin, vmax=vmax)
    cbar = fig.colorbar(image, label=cmap_name)

    _set_axis_labels2(ax, dimension, project_onto_plane, plane_eqn)
    _set_axes_limits(dimension, ax, axes_limits)
    if dimension == 2:
        plt.axes().set_aspect('equal')
    _save_or_display(save_to)


def _check_args(dimension, plane_eqn, dist_thresh):
    """Checks arguments and decides whether to project points"""
    if dimension == 3:
        project_onto_plane = False
    elif dimension == 2:
        if plane_eqn is None:
            project_onto_plane = False
        else:
            project_onto_plane = True
    else:
        raise ValueError("Parameter `dimension` needs to be either 2 or 3 but "
                         "{0} was given.".format(dimension))

    if dist_thresh is not None and plane_eqn is None:
        raise ValueError("`dist_thresh` was specified but no `plane_eqn` was "
                         "given.")

    if dist_thresh is None and dimension == 2:
        print("WARNING: A 2D plot will look cluttered without cut-off value "
              "for the distance from the specified plane (`dist_thresh`).")

    return project_onto_plane


def _set_axis_labels2(ax, dimension, project_onto_plane, plane_eqn):
    if project_onto_plane:
        ax.set_xlabel(r'Coordinates mapped onto plane ${0:.2f}x {1:+.2f}y '
                      '{2:+.2f}z {3:+.2f} = 0$'.format(*plane_eqn))
    else:
        # Zip with dimension to stop early if it's less than 3 dimensions
        for dir_label, dim in zip(DIR_LABELS, range(dimension)):
            getattr(ax, "set_" + dir_label + "label")(dir_label)


def _plot_atoms(molecule, ax, dimension, plane_eqn, project_onto_plane,
                atom_dist_threshs, atom_format, show_all_atoms):
    if molecule is None:
        return

    # Default values for formatting
    if atom_format is None:
        atom_format = [
            {'color': 'red', 'bbox': dict(facecolor='none', edgecolor='red')},
            {'color': 'red', 'bbox': dict(facecolor='none', edgecolor='red',
             linestyle='dashed')}, {'color': 'grey'}, {'color': 'grey'}]
    if atom_dist_threshs is None:
        atom_dist_threshs = [0, 0.5, 1]

    # This is outside of the loop to take advantage of projecting all atoms at
    # once with _project_points
    coords = [atom.coords for atom in molecule]
    coords = _project_points(coords, project_onto_plane, dimension, plane_eqn)

    for atom, coord in zip(molecule, coords):
        assert 0 <= len(atom_format) - len(atom_dist_threshs) <= 1
        atom_string = '{0}{1}'.format(atom.identity, atom.label)
        if plane_eqn is None:
            _plot_atom(ax, coord, atom_string, dimension, {'color': 'red'},
                       marker_fill='k')
        else:
            # This big for-else loop checks into which threshold range fits the
            # atom's distance
            for curr_thresh, curr_format in zip(atom_dist_threshs,
                                                atom_format):
                dist = _plane_point_dist(plane_eqn, atom.coords)
                # Prevent retyping _plot_atom arguments by creating a lambda
                plot_atom = lambda curr_format, marker_fill: _plot_atom(
                    ax, coord, atom_string, dimension, curr_format,
                    marker_fill=marker_fill)
                if _check_dist(dist, curr_thresh):
                    plot_atom(curr_format, 'k')
                    break
            else:
                # If it doesn't fit into any threshold, check if such atoms
                # should be plotted and if their plotting arguments have been
                # supplied as the additional, hanging element of `atom_format`
                if (len(atom_format) == len(atom_dist_threshs) + 1 and
                        show_all_atoms):
                    plot_atom(atom_format[-1], 'grey')


def _plot_atom(ax, coords, atom_string, dimension, curr_format, marker='D',
               marker_fill='b'):
    """Plot atom as text and optionally marker"""
    ax.text(*coords[:dimension], atom_string, **curr_format)
    if marker is not None:
        ax.scatter(*coords[:dimension], marker='D', c=marker_fill)


def _plane_point_dist(equation, point):
    """Calculate the distance between a point and a plane given by equation

    Parameters
    ----------
    equation : List[float]
        A list of coefficients of the equation describing the plane :math:
        `ax + by + cz + d = 0`. The length should hence be 4. For example, the
        plane :math: `z = 0` corresponds to the argument `[0, 0, 1, 0]`.

    point : List[float]
        The coordinates of the point `[x, y, z]`. A list of length 3.

    Returns
    -------
    float
        The calculated distance according to the equation:

        .. math:: \frac{ax + by + cz + d}{\root{a^2 + b^2 + c^2}}

        Returning the signed value of this expression allows to distinguish
        between points lying on the opposite sides of the plane.
    """
    result = equation[3]
    denom = 0
    for coef, coord in zip(equation, point):
        result += coef*coord
        denom += coef**2
    return result/np.sqrt(denom)


def _plane_through_points(point1, point2, point3):
    point1 = np.array(point1)
    point2 = np.array(point2)
    point3 = np.array(point3)

    u = point2 - point1
    v = point3 - point1
    cross = np.cross(u, v)

    if not np.count_nonzero(cross):
        raise ValueError("The supplied points appear to be colinear.")

    a, b, c = cross[:3]
    d = - (a*point1[0] + b*point1[1] + c*point1[2])
    return a, b, c, d


def plane_through_atoms(molecule, label1, label2, label3):
    points = [molecule[label - 1].coords for label in [label1, label2, label3]]
    return _plane_through_points(*points)


def _project_point_onto_plane(plane, point):
    """Calculate coordinates of a point perpendicularly projected onto plane"""
    point = np.array(point)
    plane = np.array(plane)
    normal = plane[:3]
    return point - np.dot(point, normal)*normal/np.linalg.norm(normal)**2


def _get_alt_coords(plane_eqn):
    """Create new coordinate system with z-axis orthogonal to given plane"""
    # Normal to the plane
    normal = np.array(plane_eqn[:3])
    # Normalize (set magnitude to 1)
    normal = normal/np.linalg.norm(normal)
    # Set suggested direction of i, here it is the old x-direction.
    i = np.array([1, 0, 0])
    # Check if the normal coincides with the x-direction. If so, the i vector
    # needs to be initially pointed in a different direction, e.g. that of y.
    if np.dot(i, normal) == 1:
        i = np.array([0, 1, 0])
    # Select direction as close to the suggested one but orthogonal to the
    # normal vector. This is done by taking the *rejected* vector when
    # projecting i onto normal (the subtrahend).
    i_prime = i - np.dot(i, normal)*normal
    # Normalize
    i_prime = i_prime/np.linalg.norm(i_prime)
    # Find vector orthogonal to both i and the normal vector by taking their
    # cross product. The order there is significant and was chosen to obtain a
    # right-handed coordinate system (i, j, normal), just like (x, y, z)
    j_prime = np.cross(normal, i_prime)
    # No need to normalize
    return i_prime, j_prime, normal


def _new_coord_matrix(new_coord_system):
    """Calculate matrix of transformation from old to new coordinate system"""
    old = np.identity(3)
    # The templates are used to create the elements of the desired 3x3 matrix.
    # They are populated in a meshgrid fashion, but I couldn't get it to work
    # due to the nesting, so I settled on a list comprehension kludge.
    # To simplify the code, the 9 elements of the matrix are kept as a
    # contiguous array of 9 vectors, hence the reshaping.
    old_template = np.array([old]*3).reshape(9, 3)
    new_template = np.array([[elem]*3 for elem in new_coord_system])
    new_template = new_template.reshape(9, 3)
    # The desired matrix is calculated as an element-wise dot product
    matrix = np.array([np.dot(old_elem, new_elem) for old_elem, new_elem in
                       zip(old_template, new_template)])
    return matrix.reshape(3, 3)


def _project_points(points, project_onto_plane, dimension, plane_eqn):
    """Project points onto the given plane (3D) or its new coordinate system"""
    if project_onto_plane:
        if dimension == 3:
            # Simple perpendicular projection onto 3D plane (this is expected
            # to be rarely used and is not accessible through the 'public'
            # `plot_points` as it switches projection off in 3D
            points = [_project_point_onto_plane(plane_eqn, point) for point in
                      points]
        elif dimension == 2:
            # This is actually more than a projection, as 'looking' at the
            # plane in a perpendicular manner requires a change of coordinate
            # system. Otherwise the points would then be projected onto the
            # (x, y) plane when flattening for plotting.
            matrix = _new_coord_matrix(_get_alt_coords(plane_eqn))
            points = [np.dot(matrix, point) for point in points]
    return points


def _check_dist(dist, thresh):
    """Check if a distance is below the given threshold value"""
    # The second condition ensures that floats are rounded correctly. With some
    # of the grids some points may lie on the threshold value but would not be
    # caught by the first condition due to float precision.
    # Absolute tolerance was selected as one decimal place fewer than what
    # seems to be the precision of Gaussian .esp coordinates.
    return abs(dist) <= thresh or math.isclose(abs(dist), thresh, abs_tol=1e-4)


def _points_dist_filter(points, values, plane_eqn, dist_thresh):
    if dist_thresh is None or plane_eqn is None:
        return points, values
    _points, _values = [], []
    for point, value in zip(points, values):
        dist = _plane_point_dist(plane_eqn, point)
        if _check_dist(dist, dist_thresh):
            _points.append(point)
            _values.append(value)
    return _points, _values


def _points_rand_skim(points, values, rand_skim):
    if rand_skim == 1:
        return points, values
    _points, _values = [], []
    for point, value in zip(points, values):
        if random.random() <= rand_skim:
            _points.append(point)
            _values.append(value)
    return _points, _values
