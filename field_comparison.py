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
    # Note the code repetition from the `skim' function. The logic is somewhat
    # different here so I wasn't sure how I could generalize and decided that
    # modularity was more important. Note that the arrays are traversed twice
    # if both this function and `skim' are called consecutively.
    _check_grids(dist, *fields)
    if 'dist' not in dist.field_type:
        raise TypeError("The field passed was of type '{0}' but one of the "
                        "'dist' types was expected.".format(dist.field_type))

    fields = [field.values.copy() if type(field) is Field else field.copy() for
              field in fields]
    dist = dist.values.copy() if type(dist) is Field else dist.copy()
    for dist_elem, *field_elems in np.nditer([dist] + fields,
                                             op_flags=['readwrite']):
        if dist_elem <= exclusion_dist:
            dist_elem[...] = assign_val
            for field_elem in field_elems:
                field_elem[...] = assign_val
    # Return all the input fields, including distance fields, as a list
    return [dist] + fields


def skim(rand_skim, *fields, assign_val=None):
    _check_grids(*fields)
    fields = [field.values.copy() if type(field) is Field else field.copy() for
              field in fields]
    # Numpy array iteration with nditer has a convenient write mode, but it
    # would affect the input arrays, so copies were made earlier.
    for field_elems in np.nditer(fields, op_flags=['readwrite']):
        if random.random() > rand_skim:
            for field_elem in field_elems:
                field_elem[...] = assign_val
    return fields
