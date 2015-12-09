import ipdb
import random
import numpy as np
import inspect
from warnings import warn
from cube_helpers import GridError, Field
from operator import attrgetter


def difference(field1, field2, relative=False, absolute=False):
    _check_grids(field1, field2)
    if absolute:
        func = lambda val1, val2: abs(val1 - val2)
        name = 'abs_'
    else:
        func = lambda val1, val2: val1 - val2
        name = ''

    if relative:
        func = lambda val1, val2: func(val1, val2)/val1
        name += 'rel_'

    name += 'diff'
    values = np.vectorize(func)(field1.values, field2.values)

    return Field(values, field1.grid, name)


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


def filter_by_dist(exclusion_dist, dist, *fields, assign_val=None):
    """Filter the datapoints in input fields by their distances.

    This function takes any number of field-like objects which describe
    datapoints corresponding between the fields. The `dist` object specifies
    their distances from something, e.g. an electron density isosurface. This
    function will filter out the values of all datapoints, whose distances are
    smaller than a threshold value, will be filtered or actually replaced by a
    placeholder value `assign_val`. The input objects are not affected, as new
    np.ndarrays are returned.


    Parameters
    ----------
    exclusion_dist : float
        Datapoints at distances smaller or equal this number will be filtered.
    dist : Field or np.ndarray
        The Field or ndarrays specifying distance of each datapoint. This
        object also gets filtered.
    *fields : Field or np.ndarray
        The Fields or ndarrays containing the datapoints to be filtered by
        distance.
    assign_val : Optional
        The datapoints within the exclusion distance will be replaced by this
        value. Defaults to None.

    Returns
    -------
    List[np.ndarray]
        The filtered values of `dist` and `fields`.

    """
    condition = lambda elems: elems[0] <= exclusion_dist
    fields = [dist] + list(fields)
    return _iterate_fields(condition, assign_val, *fields)


def skim(rand_skim, *fields, assign_val=None):
    """Skim the number of datapoints in input fields.

    This function takes any number of field-like objects which describe
    datapoints corresponding between the fields. This function will randomly
    remove some of the datapoints (the corresponding points between the
    different fields) and retain a fraction of all datapoints given by the
    `rand_skim` argument. The datapoints are not actually removed but replaced
    with placeholder values `assign_val`. The input objects are not affected,
    as new np.ndarrays are returned.

    This function is useful for reducing the number of datapoints in a 3D plot,
    whose interactivity or clarity would otherwise be affected.

    Parameters
    ----------
    rand_skim : float
        A number between 0 and 1. The fraction of datapoints to be retained.
    *fields : Field or np.ndarray
        The Fields or ndarrays containing the datapoints to be skimmed.
    assign_val : Optional
        The datapoints to be removed will be replaced by this value. Defaults
        to None.

    Returns
    -------
    List[np.ndarray]
        The skimmed fields.

    """
    condition = lambda elems: random.random() > rand_skim
    return _iterate_fields(condition, assign_val, *fields)


def _iterate_fields(condition, assign_val, *fields):
    """Iterate and remove corresponding elements of ndarrays

    This function takes any number of field-like objects which describe
    datapoints corresponding between the fields. It will remove the
    corresponding datapoints from all the fields depending whether they satisfy
    the passed `condition` function. The datapoints are not actually removed
    but replaced with placeholder values `assign_val`. The input objects are
    not affected, as new np.ndarrays are returned.

    Parameters
    ----------
    condition : bool
        A function which takes the list of corresponding elements in the
        iterated fields and return a boolean value stating whether the list
        satisfies the condition.
    assign_val
        The datapoints to be removed will be replaced by this value.
    *fields : Field or np.ndarray
        The fields containing the datapoints to be iterated.

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
