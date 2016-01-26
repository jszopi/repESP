import ipdb
import numpy as np
from operator import attrgetter
from scipy.ndimage.morphology import distance_transform_edt as scipy_edt
from scipy.spatial.distance import euclidean

AXES = ['x', 'y', 'z']


class GridError(Exception):
    pass


class Cube(object):

    title_to_type = {
        ' Electrostatic pot': 'esp',
        ' Electron density ': 'ed',
        }

    def __init__(self, cube_fn):
        with open(cube_fn, 'r') as f:

            self.gaussian_input = f.readline().rstrip('\n')
            self.title = f.readline().rstrip('\n')
            self.atom_count, *origin_coords, nval = f.readline().split()
            if float(nval) != 1:
                raise GridError('NVal in the cube is different than 1. Not '
                                'sure what it means in practice.')
            self.atom_count = int(self.atom_count)

            grid = Grid([f.readline().split() for i in range(3)])
            grid.origin_coords = [float(coord) for coord in origin_coords]

            self.molecule = Molecule()
            # The atoms will be added to the Molecule in the order of occurence
            # in the input, which is assumed to correspond to Gaussian labels.
            for label in range(self.atom_count):
                atom_temp = f.readline().split()
                for index in range(4):
                    atom_temp[index+1] = float(atom_temp[index+1])

                new_atom = Atom(int(label)+1, int(atom_temp[0]), atom_temp[2:])
                new_atom.charges['cube'] = atom_temp[1]
                self.molecule.append(new_atom)

            # This may be unfeasible for very large cubes
            field = f.read().split()

        try:
            self.cube_type = Cube.title_to_type[self.title[:18]]
        except KeyError:
            raise NotImplementedError("Cube title '" + self.title + "' is not "
                                      "associated with a known cube type.")

        self.field = Field(Cube.field_from_raw(field, grid), grid,
                           self.cube_type)

    @staticmethod
    def field_from_raw(raw_field, grid):
        field = np.array(list(map(float, raw_field)))
        if len(field) != np.prod(grid.points_on_axes):
            raise GridError('The number of points in the cube {0} is not equal'
                            ' to the product of number of points in the XYZ '
                            'directions given in the cube header: {1}.'
                            .format(len(field), grid.points_on_axes))

        field.resize(grid.points_on_axes)
        return field


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

    def print_with_charge(self, charge_type):
        print(self, ', charge: {0: .4f}'.format(self.charges[charge_type]),
              sep='')

    def __str__(self):
        return 'Atom {0:2}:  {1:2}'.format(self.label, self.identity)

    def __repr__(self):
        return str(self)


class Molecule(list):
    """A list of atoms with extra functionalities."""

    def __init__(self, *args):
        list.__init__(self, *args)

    def calc_field(self, grid, field_func, *field_func_args):
        """Calculate field values point-wise according to a function

        Note that this method returns as many fields as the calculation
        function returns values. The function is selected from this class based
        on a name description (through `getattr').
        """
        # PERFORMANCE CONSIDERATIONS
        # This will likely be a bottleneck and, while a mathematical trick to
        # reduce the complexity from O(n*g^3) may exist, I don't know one.
        # Still, this method be optimized but **only once its proven a
        # bottleneck** by profiling.
        # (1) Some straightforward optimization has been performed. The reason
        # why this method expects a field_func returning more than one value is
        # that calling field_func only once prevents the reevaluation of some
        # common parts of its body. For _rep_esp_func and _dist_func, it is the
        # `euclidean' distance, and for the latter also the logic behind
        # choosing the closest atom.
        # (2) However, `euclidean' still gets reevaluated for different calls
        # of this method. It could be memoized as a distance field resulting
        # from a new Molecule method. Since that would better operate on grid
        # rather than coordinates, the loop in *this* method should iterate the
        # grid indices and pass all the grid calculations to the distance
        # method. TODO
        # (3) np.array should be preferred over the intermediate list here
        # (4) Then it may be worth mapping field_func onto array elements
        # instead of iterating over it but that's a disputable topic:
        # https://wiki.python.org/moin/PythonSpeed/PerformanceTips#Loops
        # https://www.python.org/doc/essays/list2str/
        # (5) Finally, iterating the grid has a good potential for
        # parallelization.

        field_func, field_types, field_infos = self._field_func_helper(
            field_func, *field_func_args)

        for ix in range(grid.axes[0].point_count):
            x = grid.origin_coords[0] + ix*grid.dir_intervals[0]
            for iy in range(grid.axes[1].point_count):
                y = grid.origin_coords[1] + iy*grid.dir_intervals[1]
                for iz in range(grid.axes[2].point_count):
                    z = grid.origin_coords[2] + iz*grid.dir_intervals[2]
                    values = field_func(x, y, z, *field_func_args)
                    while True:
                        try:
                            for value, result in zip(values, results):
                                result.append(value)
                            break
                        except NameError:
                            results = [[] for i in range(len(values))]

        fields = []
        for result, field_type, field_info in zip(results, field_types,
                                                  field_infos):
            field = np.array(result)
            field.resize(grid.points_on_axes)
            fields.append(Field(field, grid, field_type, field_info))

        return fields

    def _field_func_helper(self, field_func, *field_func_args):
        # Any new field_func must leave details in here
        if field_func == 'rep_esp':
            field_types = ['rep_esp']*len(field_func_args[0])
            field_infos = field_func_args[0]
            func = self._rep_esp_func
        elif field_func == 'dist':
            field_types = ['closest_atom', 'closest_atom_dist']
            field_infos = ['own']*2
            # Own implementation as opposed to Henkelman's
            func = self._dist_func
        else:
            raise NotImplementedError("The requested function is not "
                                      "implemented: {0}".format(field_func))
        return func, field_types, field_infos

    def _rep_esp_func(self, x, y, z, charge_types):
        """Calculate ESP value at given point due to charges on atoms"""
        values = [0]*len(charge_types)
        for atom in self:
            dist = euclidean([x, y, z], atom.coords)
            for i, charge_type in enumerate(charge_types):
                values[i] += atom.charges[charge_type]/dist
        return values

    def _dist_func(self, x, y, z):
        """For a given point, find the closest atom and its distance"""
        min_dist = float('inf')
        min_atom = None
        for atom in self:
            dist = euclidean([x, y, z], atom.coords)
            if dist < min_dist:
                min_dist = dist
                min_atom = atom.label
        return (min_atom, min_dist)


class Field(object):

    def __init__(self, values, grid, field_type, field_info=None):
        self.values = values
        self.grid = grid
        self.field_type = field_type
        self.field_info = field_info

    def distance_transform(self, isovalue):
        """This should only be applied to the electron density cube."""

        if self.field_type != 'ed':
            print("WARNING: Distance transform should only be applied to "
                  "electron density fields, attempted on field type: '{0}'."
                  .format(self.field_type))

        if not self.grid.aligned_to_coord:
            raise GridError('Distance transform not implemented for grid not '
                            'aligned with the coordinate system.')

        # Select isosurface and its interior as a 3D solid of 0s.
        select_iso = lambda x: 1 if x < isovalue else 0
        field = np.vectorize(select_iso)(self.values)
        dist = scipy_edt(field, sampling=self.grid.dir_intervals)
        return Field(dist, self.grid, 'ed_dist')

    def write_cube(self, output_fn, molecule, charge_type=None):
        """Write the field as a Gaussian cube file.

        Raises FileExistsError when the file exists.
        """
        with open(output_fn, 'x') as f:
            f.write(' Cube file generated by repESP.\n')
            f.write(' Cube file for field of type {0}.\n'.format(
                self.field_type))
            f.write(' {0:4}   {1: .6f}   {2: .6f}   {3: .6f}    1\n'.format(
                len(molecule), *self.grid.origin_coords))
            for axis in self.grid.axes:
                f.write(' {0:4}   {1: .6f}   {2: .6f}   {3: .6f}\n'.format(
                    axis.point_count, *axis.intervals))
            for atom in molecule:
                if charge_type is None:
                    charge = atom.atomic_no
                else:
                    charge = atom.charges[charge_type]
                f.write(' {0:4}   {1: .6f}   {2: .6f}   {3: .6f}   {4: .6f}\n'
                        .format(atom.atomic_no, charge, *atom.coords))
            i = 1
            for value in self.values.flatten():
                f.write(' {0: .5E}'.format(value))
                if not i % 6:
                    f.write('\n')
                if not i % self.grid.axes[2].point_count:
                    f.write('\n')
                    i = 1
                else:
                    i += 1


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

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

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

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

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
