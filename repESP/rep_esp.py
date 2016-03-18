from scipy.spatial.distance import euclidean
import numpy as np

from .resp_helpers import Points, NonGridField
from .cube_helpers import GridField, Grid, angstrom_per_bohr


def _calc_field(molecule, points, field_func, *field_func_args):
    """Calculate values at points to a function

    Note that this method returns as many `NonGridField` objects as the
    calculation function returns values. The function is selected from this
    module based on a name description (through `globals`).
    """
    # PERFORMANCE CONSIDERATIONS
    # This will likely be a bottleneck and, while a mathematical trick to
    # reduce the complexity from O(n*g^3) may exist, I don't know one.
    # Still, this method should be optimized but **only once its proven a
    # bottleneck** by profiling.
    # (1) Some straightforward optimization has been performed. The reason
    # why this method expects a field_func returning more than one value is
    # that calling field_func only once prevents the reevaluation of some
    # common parts of its body. For _rep_esp_func and _dist_func, it is the
    # `euclidean' distance, and for the latter also the logic behind
    # choosing the closest atom.
    # (2) However, `euclidean' still gets reevaluated for different calls
    # of this method. It could be memoized as a distance field resulting
    # from a new Molecule method. TODO
    # (3) np.array should be preferred over the intermediate `results` list
    # (4) Then it may be worth mapping field_func onto array elements
    # instead of iterating over it but that's a disputable topic:
    # https://wiki.python.org/moin/PythonSpeed/PerformanceTips#Loops
    # https://www.python.org/doc/essays/list2str/
    # (5) Finally, iterating the points has a good potential for
    # parallelization.

    # This is important to prevent the coordinates being given as a list and
    # hence new Points objects being created upon initialization of
    # NonGridField
    if type(points) is not Points:
        raise TypeError("A Points object was expected, given {0}.".format(
            type(points)))

    for point in points:
        values = field_func(molecule, point[0], point[1], point[2],
                            *field_func_args)
        # A convoluted hack to create the list of `results` only when its
        # length is known i.e. the first set of `values` is evaluated.
        while True:
            try:
                for value, result in zip(values, results):
                    result.append(value)
                break
            except NameError:
                results = [[] for i in range(len(values))]
    return results


def calc_grid_field(molecule, grid, field_func, *field_func_args):
    # Some code overlap between this function and the non-grid one, but IMO
    # trying to DRY would overcomplicate it due to slight differences.
    if type(grid) is not Grid:
        raise TypeError("A Grid object was expected, given {0}.".format(
            type(grid)))

    points = []
    for ix in range(grid.axes[0].point_count):
        x = grid.origin_coords[0] + ix*grid.dir_intervals[0]
        for iy in range(grid.axes[1].point_count):
            y = grid.origin_coords[1] + iy*grid.dir_intervals[1]
            for iz in range(grid.axes[2].point_count):
                z = grid.origin_coords[2] + iz*grid.dir_intervals[2]
                points.append((x, y, z))
    points = Points(points)

    field_func, field_types, field_infos = _field_func_helper(
        field_func, *field_func_args)
    results = _calc_field(molecule, points, field_func, *field_func_args)

    fields = []
    for result, field_type, field_info in zip(results, field_types,
                                              field_infos):
        field = np.array(result)
        field.resize(grid.points_on_axes)
        fields.append(GridField(field, grid, field_type, field_info))

    return fields


def calc_non_grid_field(molecule, points, field_func, *field_func_args):

    field_func, field_types, field_infos = _field_func_helper(
        field_func, *field_func_args)
    results = _calc_field(molecule, points, field_func, *field_func_args)

    fields = []
    for result, field_type, field_info in zip(results, field_types,
                                              field_infos):
        field = np.array(result)
        fields.append(NonGridField(field, points, field_type, field_info))

    return fields


def _field_func_helper(field_func, *field_func_args):
    # Any new field_func must leave details in here
    if field_func == 'rep_esp':
        field_types = ['rep_esp']*len(field_func_args[0])
        field_infos = [[elem] for elem in field_func_args[0]]
        func = _rep_esp_func
    elif field_func == 'dist':
        field_types = ['parent_atom', 'dist']
        # Voronoi means closest-atom in molecular partitioning lingo
        field_infos = [['Voronoi']]*2
        func = _dist_func
    else:
        raise NotImplementedError("The requested function is not "
                                  "implemented: {0}".format(field_func))
    return func, field_types, field_infos


def _rep_esp_func(molecule, x, y, z, charge_types):
    """Calculate ESP value at given point due to charges on atoms"""
    values = [0]*len(charge_types)
    for atom in molecule:
        dist = euclidean([x, y, z], atom.coords)
        for i, charge_type in enumerate(charge_types):
            values[i] += atom.charges[charge_type]/(dist/angstrom_per_bohr)
    return values


def _dist_func(molecule, x, y, z):
    """For a given point, find the closest atom and its distance"""
    min_dist = float('inf')
    min_atom = None
    for atom in molecule:
        dist = euclidean([x, y, z], atom.coords)
        if dist < min_dist:
            min_dist = dist
            min_atom = atom.label
    return (min_atom, min_dist)
