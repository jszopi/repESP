import ipdb
import random
import numpy as np
import inspect
from warnings import warn
from cube_helpers import GridError, Field


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
    if not len(fields):
        warn('No fields to be compared! As this is a helper function, the '
             'issue is likely due to the caller: ' + inspect.stack()[1][3])
    for field in fields:
        if field1.grid != field.grid:
            raise GridError('Grids of the fields to be compared do not match.')


def filter_by_dist(dist, *fields, exclusion_dist=-1, rand_skim=1,
                   assign_val=None):
    _check_grids(dist, *fields)
    fields = [field.values.copy() for field in fields]
    dist = dist.values.copy()
    # Numpy array iteration with nditer has a convenient write mode, but it
    # would affect the input arrays, so copies were made earlier.
    for dist_elem, *field_elems in np.nditer([dist] + fields,
                                             op_flags=['readwrite']):
        if dist_elem <= exclusion_dist or random.random() > rand_skim:
            dist_elem[...] = assign_val
            for field_elem in field_elems:
                field_elem[...] = assign_val
    # Return all the input fields as a list
    return [dist] + fields
