import ipdb
import numpy as np
from cube_helpers import GridError, Field


def difference(field1, field2, relative=False):
    _check_grids(field1, field2)
    if relative:
        func = lambda val1, val2: (val1 - val2)/val1
        name = 'rel_diff'
    else:
        func = lambda val1, val2: val1 - val2
        name = 'diff'

    values = np.vectorize(func)(field1.values, field2.values)

    return Field(values, field1.grid, name)


def _check_grids(field1, field2):
    if field1.grid != field2.grid:
        raise GridError('Grids of the fields to be compared do not match.')
