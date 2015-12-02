import ipdb
import numpy as np
from operator import attrgetter
from scipy.ndimage.morphology import distance_transform_edt as scipy_edt

AXES = ['x', 'y', 'z']


class GridError(Exception):
    pass


class Cube(object):

    def __init__(self, cube_fn):

        with open(cube_fn, 'r') as f:

            self.gaussian_input = f.readline().rstrip('\n')
            self.cube_title = f.readline().rstrip('\n')
            self.atom_count, *origin_coords, nval = f.readline().split()
            if float(nval) != 1:
                raise GridError('NVal in the cube is different than 1. Not '
                                'sure what it means in practice.')
            self.atom_count = int(self.atom_count)

            self.grid = Grid([f.readline().split() for i in range(3)])
            self.grid.origin_coords = [float(coord) for coord in origin_coords]

            self.atoms = []
            # The atoms will be added to the list in the order of occurence in
            # the input, which is assumed to correspond to Gaussian labels.
            for label in range(self.atom_count):
                atom_temp = f.readline().split()
                for index in range(4):
                    atom_temp[index+1] = float(atom_temp[index+1])

                new_atom = Atom(int(label)+1, int(atom_temp[0]), atom_temp[2:])
                new_atom.charges['cube'] = atom_temp[1]
                self.atoms.append(new_atom)

            # This may be unfeasible for very large cubes
            self.field = f.read().split()

        self.field = np.array(list(map(float, self.field)))

        if len(self.field) != np.prod(self.grid.points_on_axes):
            raise GridError('The number of points in the cube {0} is not equal'
                            ' to the product of number of points in the XYZ '
                            'directions given in the cube header: {1}.'
                            .format(len(self.field), self.grid.points_on_axes))

        self.field.resize(self.grid.points_on_axes)

    def distance_transform(self, isovalue):
        """This should only be applied to the electron density cube."""

        if not self.grid.aligned_to_coord:
            raise GridError('Distance transform not implemented with cube not '
                            'aligned with the coordinate system.')

        # Select isosurface and its interior as a 3D solid of 0s.
        select_iso = lambda x: 1 if x < isovalue else 0
        field = np.vectorize(select_iso)(self.field)
        return scipy_edt(field, sampling=self.grid.dir_intervals)


class Atom(object):

    # http://www.science.co.il/PTelements.asp
    periodic = [('H', 'Hydrogen'),
                ('He', 'Helium'),
                ('Li', 'Lithium'),
                ('Be', 'Beryllium'),
                ('B', 'Boron'),
                ('C', 'Carbon'),
                ('N', 'Nitrogen'),
                ('O', 'Oxygen'),
                ('F', 'Fluorine'),
                ('Ne', 'Neon'),
                ('Na', 'Sodium'),
                ('Mg', 'Magnesium'),
                ('Al', 'Aluminum'),
                ('Si', 'Silicon'),
                ('P', 'Phosphorus'),
                ('S', 'Sulfur'),
                ('Cl', 'Chlorine'),
                ('Ar', 'Argon')]

    def __init__(self, label, atomic_no, coords=None):
        self.label = label
        self.atomic_no = atomic_no
        try:
            self.identity = Atom.periodic[atomic_no-1][0]
        except IndexError:
            print('WARNING: Element of atomic number {0} not implemented. '
                  'Setting its identity to atomic number'.format(atomic_no))
            self.identity = str(atomic_no)

        self.charges = {}
        self.coords = coords

    def __str__(self):
        return 'Atom ' + str(self.label) + ' (' + self.identity + ')'

    def __repr__(self):
        return str(self)


class Grid(object):

    def __init__(self, grid_input):

        self.origin_coords = None

        if np.shape(grid_input) != (3, 4):
            raise GridError('Incorrect grid formatting. Expected a list of '
                            'shape 3x4, instead got: ' + str(grid_input))

        self.axes = [GridAxis(label) for label in AXES]
        self.aligned_to_coord = True

        for axis_number, input_axis in enumerate(grid_input):
            aligned_to_axis = self._add_axis(axis_number, input_axis)
            self.aligned_to_coord = self.aligned_to_coord and aligned_to_axis

        self.dir_intervals = []
        if self.aligned_to_coord:
            for axis in range(3):
                self.dir_intervals.append(self.axes[axis].dir_interval)
        else:
            print('WARNING: The cube is not aligned with coordinate system.')

        self.points_on_axes = [axis.point_count for axis in self.axes]

    def _add_axis(self, axis_number, input_axis):
        axis_to_set = self.axes[axis_number]
        axis_to_set.set_point_count(input_axis.pop(0))
        return axis_to_set.set_intervals(input_axis)


class GridAxis(object):

    def __init__(self, label):
        self.label = label
        self.point_count = None
        self.intervals = []  # xyz
        self.dir_interval = None  # Interval in its 'own' direction

    def set_point_count(self, point_count):

        if int(point_count) != float(point_count):
            raise GridError('Number of points in direction {0} is not an '
                            'integer: {1}'.format(self.label, point_count))
        if self.label == 'x' and float(point_count) < 0:
            raise GridError('Gaussian requested distance in angstroms, which '
                            'is not currently supported.')

        self.point_count = int(point_count)

    def set_intervals(self, intervals):

        aligned_to_coord_axis = True

        for direction, interval in enumerate(intervals):
            self.intervals.append(float(interval))
            if AXES[direction] == self.label:
                self.dir_interval = float(interval)
            elif float(interval) != 0:
                aligned_to_coord_axis = False

        if not aligned_to_coord_axis:
            print('INFO: Cube axis {0} is not aligned to its coordinate'
                  ' axis: The intervals are: {1}'.format(self.label,
                                                         intervals))

        return aligned_to_coord_axis
