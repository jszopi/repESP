import ipdb
import numpy as np
from scipy.spatial.distance import euclidean


def reproduce_field(atom_list, charge_type, grid):
    field = []
    for ix in range(grid.axes[0].point_count):
        x = grid.origin_coords[0] + ix*grid.dir_intervals[0]
        for iy in range(grid.axes[1].point_count):
            y = grid.origin_coords[1] + iy*grid.dir_intervals[1]
            for iz in range(grid.axes[2].point_count):
                z = grid.origin_coords[2] + iz*grid.dir_intervals[2]
                value = 0
                for atom in atom_list:
                    value += atom.charges[charge_type]/euclidean([x, y, z],
                                                                 atom.coords)
                field.append(value)

    field = np.array(field)
    field.resize(grid.points_on_axes)

    return field
