import random
import numpy as np
import inspect
from warnings import warn
from cube_helpers import GridError, Field, _check_for_nans
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


def _check_fields_for_nans(*fields):
    for field in fields:
        _check_for_nans(field.values)


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
    # Set assign_val to be zero for the atomic label field to prevent an error
    # about assigning None to an np.array of an integer dtype. This is expected
    # to be the only case when atomic label fields are filtered, so a local
    # workaround is fine. However, if this turns out to be more common, a
    # general solution will be more elegant.
    assign_val = [assign_val]*len(fields)
    assign_val[0] = 0
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


def _iterate_fields(condition, assign_vals, *fields):
    """Iterate and remove corresponding elements in ndarrays.

    This function takes any number of field-like objects, whose data points
    correspond to each other between the fields. It will remove the
    corresponding data points from all the fields depending on whether the list
    of corresponding data points satisfies the passed `condition` function. The
    data points are not actually removed but replaced with a special value
    `assign_vals`. The input objects are not affected, as new ndarrays are
    returned.

    Parameters
    ----------
    condition : bool
        A function which takes the list of corresponding elements in the
        iterated fields and return a boolean value stating whether the list
        satisfies the condition.
    assign_vals
        The data points to be removed will be replaced by this value. This
        argument can be a single value or a list of the same length as the
        number of fields passed to specify a different value for each field.
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

    # Check if assign_vals is a list
    try:
        # Perhaps this should raise a more informative error, but since this is
        # an internal-use function, an assertion is fine
        assert len(assign_vals) == len(fields)
    except TypeError:
        # If not a list, make it a list of identical values
        assign_vals = [assign_vals] * len(fields)

    # Numpy array iteration with nditer has a convenient write mode, but it
    # would affect the input arrays, so copies were made earlier.
    for field_elems in np.nditer(fields, op_flags=['readwrite']):
        if condition(field_elems):
            # An assertion just to make sure that the array iteration logic is
            # correct.
            assert len(field_elems) == len(assign_vals)
            for field_elem, assign_val in zip(field_elems, assign_vals):
                try:
                    field_elem[...] = assign_val
                except TypeError as e:
                    raise TypeError(
                        "Could not assign `{0}` as an element of a np.array of"
                        " dtype `{1}`. The original, tricky example of this "
                        "error is assigning `None` to an integer np.array. The"
                        " caller should have specified a suitable value to "
                        "assign, e.g. 0 for atomic label field. The handled "
                        "error message was : \nTypeError: ".format(
                            assign_val, field_elem.dtype) + str(e))
                    # Atomic label fields are only expected to come up when
                    # filtering by atomic label (filter_by_atom), so that
                    # function manually sets 0 for that field.

    return fields
