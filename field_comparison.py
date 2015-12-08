import ipdb
import numpy as np
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
    for field in fields:
        if field1.grid != field.grid:
            raise GridError('Grids of the fields to be compared do not match.')
