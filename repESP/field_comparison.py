import random
import numpy as np
import inspect
from warnings import warn
from cube_helpers import GridError, Field
from operator import attrgetter
import os


def difference(field1, field2, relative=False, absolute=False):
    _check_grids(field1, field2)
    # The first element contains information about whether the difference is
    # relative or absolute, the second contains the free-form names of the
    # compared fields.
    info = [[], []]
    if absolute:
        func = lambda val1, val2: abs(val1 - val2)
        info[0].append('abs')
    else:
        func = lambda val1, val2: val1 - val2

    if relative:
        func = lambda val1, val2: func(val1, val2)/val1
        info[0].append('rel')

    info[1].append(field1.lookup_name())
    info[1].append(field2.lookup_name())

    values = np.vectorize(func)(field1.values, field2.values)

    return Field(values, field1.grid, 'diff', info)


def _check_grids(field1, *fields):
    """Check if field grids or array dimensions match."""
    if not len(fields):
        warn('No fields to be compared! As this is a helper function, the '
             'issue is likely due to the caller: ' + inspect.stack()[1][3])

    attr_dict = {Field: 'grid', np.ndarray: 'shape'}
    for field in fields:
        if type(field) != type(field1):
            raise TypeError("Expected '{0}' objects to be compared, received "
                            "'{1}'.".format(type(field1), type(field)))
        else:
            try:
                attr = attr_dict[type(field)]
            except KeyError:
                raise TypeError("Checking dimensions of type '{0}' is not "
                                "supported.".format(type(field)))
            if attrgetter(attr)(field1) != attrgetter(attr)(field):
                raise GridError('Grids of fields to be compared do not match.')


def _flatten_no_nans(ndarray_input):
    """Flatten ndarray and remove None elements."""
    return [elem for elem in ndarray_input.flat if not np.isnan(elem)]


def filter_by_dist(exclusion_dist, dist, *fields, assign_val=None):
    """Filter the data points in input fields by their distances.

    This function takes any number of field-like objects, whose data points
    correspond to each other between the fields. The `dist` object specifies
    their distance transform, e.g. from an electron density isosurface. All
    data points, whose distances are smaller or equal a threshold value, will
    be filtered out (actually, replaced with a special value `assign_val`). The
    input objects are not affected, as new ndarrays are returned.

    Parameters
    ----------
    exclusion_dist : float
        data points at distances smaller or equal this number will be filtered.
    dist : Field or np.ndarray
        The object specifying distance of each data point (distance transform).
        The data points in this object also get filtered.
    *fields : Field or np.ndarray
        The objects containing the data points to be filtered by distance.
    assign_val : Optional
        The data points within the exclusion distance will be replaced by this
        value. Defaults to None.

    Returns
    -------
    List[np.ndarray]
        The filtered values of `dist` and `fields`.

    """
    condition = lambda elems: elems[0] <= exclusion_dist
    fields = [dist] + list(fields)
    return _iterate_fields(condition, assign_val, *fields)


def filter_by_atom(molecule, atom_label, method, *fields, assign_val=None):
    # _iterate_fields will check grids anyway
    grid = fields[0].grid
    if method == 'dist':
        closest_atom = molecule.calc_field(grid, 'dist')[0]
    elif method == 'qtaim':
        # TODO: It's not a good idea to assume the location of those cubes.
        # Currently, all the paths should be specified by the caller, i.e. the
        # path should be one of the arguments of this function. However, in the
        # future the issue of paths will be resolved in some smart way, so I'll
        # leave it as it is now.
        path = os.path.dirname(molecule.parent_cube.cube_fn) + '/'
        closest_atom = molecule.extract_qtaim_basins(grid, path)
    else:
        raise NotImplementedError("Method {0} is not implemented.".format(
            method))
    condition = lambda elems: elems[0] != atom_label
    fields = [closest_atom] + list(fields)
    return _iterate_fields(condition, assign_val, *fields)


def skim(rand_skim, *fields, assign_val=None):
    """Skim the number of data points in input fields.

    This function takes any number of field-like objects, whose data points
    correspond to each other between the fields. This function will randomly
    remove some of the corresponding data points and retain a fraction of all
    data points given by the `rand_skim` argument. The data points are not
    actually removed but replaced with a special value `assign_val`. The input
    objects are not affected, as new ndarrays are returned.

    This function is useful for reducing the number of data points in a 3D
    plot, whose interactivity or clarity would otherwise be compromised.

    Parameters
    ----------
    rand_skim : float
        A number between 0 and 1. The fraction of data points to be retained.
    *fields : Field or np.ndarray
        The objects containing the data points to be skimmed.
    assign_val : Optional
        The data points to be removed will be replaced by this value. Defaults
        to None.

    Returns
    -------
    List[np.ndarray]
        The skimmed fields.

    """
    condition = lambda elems: random.random() > rand_skim
    return _iterate_fields(condition, assign_val, *fields)


def _iterate_fields(condition, assign_val, *fields):
    """Iterate and remove corresponding elements in ndarrays.

    This function takes any number of field-like objects, whose data points
    correspond to each other between the fields. It will remove the
    corresponding data points from all the fields depending on whether the list
    of corresponding data points satisfies the passed `condition` function. The
    data points are not actually removed but replaced with a special value
    `assign_val`. The input objects are not affected, as new ndarrays are
    returned.

    Parameters
    ----------
    condition : bool
        A function which takes the list of corresponding elements in the
        iterated fields and return a boolean value stating whether the list
        satisfies the condition.
    assign_val
        The data points to be removed will be replaced by this value.
    *fields : Field or np.ndarray
        The fields containing the data points to be iterated through.

    Returns
    -------
    List[np.ndarray]
        The transformed field values.

    """
    _check_grids(*fields)
    fields = [field.values.copy() if type(field) is Field else field.copy() for
              field in fields]
    # Numpy array iteration with nditer has a convenient write mode, but it
    # would affect the input arrays, so copies were made earlier.
    for field_elems in np.nditer(fields, op_flags=['readwrite']):
        if condition(field_elems):
            for field_elem in field_elems:
                field_elem[...] = assign_val
    return fields
