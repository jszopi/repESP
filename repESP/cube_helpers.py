import numpy as np
from operator import attrgetter
from scipy.ndimage.morphology import distance_transform_edt as scipy_edt
from scipy.spatial.distance import euclidean
import glob
import os
import sys

# http://www.gaussian.com/g_tech/g_ur/k_constants.htm
angstrom_per_bohr = 0.5291772086


def _check_for_nans(values):
    try:
        values = values.flat
    except AttributeError:
        pass
    # http://stackoverflow.com/a/6736970
    if np.isnan(np.sum(values)):
        raise InputFormatError("Values contain NANs!")


class Molecule(list):
    """A list of atoms with extra functionalities."""

    def __init__(self, parent_cube, *args):
        list.__init__(self, *args)
        self.parent_cube = parent_cube

    def verbose_compare(self, other):
        if self == other:
            print("The molecules are the same.")
            return
        # Otherwise:
        print("The molecules differ at the following atoms:")
        for atom, other_atom in zip(self, other):
            if atom != other_atom:
                print("{0} != {1}".format(atom, other_atom))
        if len(self) != len(other):
            which = self if len(self) > len(other) else other
            which_str = 'first' if len(self) > len(other) else 'second'

            print("The {0} molecule has {1} more atoms:".format(
                which_str, abs(len(other) - len(self))))
            for atom in which[min(len(self), len(other)):]:
                print(atom)

    def extract_qtaim_basins(self, grid, path):
        """Extract QTAIM basins from Henkelman group's ``bader`` program

        The ``bader`` command needed to generate input cube files is::

            bader -p all_atom -vac off density.cube

        Assigning low density points to vacuum needs to be switched off in
        order to allow the basins to extend to infinity.

        This method returns a field with atomic labels indicating which basin
        each point belongs to.
        """
        output_files = glob.glob(path + 'BvAt*.cube')
        expected = [path + 'BvAt{0:04}.cube'.format(i+1) for i in
                    range(len(self))]

        if sorted(output_files) != expected:
            if len(output_files) == 0:
                msg = "No ``bader`` output cube files found!"
            else:
                for output_file, expected_file in zip(output_files, expected):
                    if output_file != expected_file:
                        msg += "Missing expected ``bader`` output cube file: "
                        msg += os.path.basename(expected)
                        break

            raise InputFormatError(msg + " To generate the files use the "
                                   "command: ``bader -p all_atom -vac off "
                                   "density.cube``")

        cubes = [Cube(expected_file) for expected_file in expected]
        # Compare grids with that provided. TODO: Would be better to use the
        # function field_comparison._check_grids, but can't import that module
        # here and won't be able to pass a special message. All that requires
        # some refactoring but is a sensible thing to do.
        for i, cube in enumerate(cubes):
            if cube.field.grid != grid:
                raise GridError("The grid of `bader' cube number {0} is "
                                "different from that of the molecule "
                                "requesting extraction.".format(i+1))

        result = []
        # Iterate all the cubes element-wise and produce a field with atomic
        # labels indicating which basin each point belongs to.
        # (This probably isn't the numpy way of doing this. It operates on
        # iterators though, so should be memory efficient.)
        for point in zip(*[cube.field.values.flat for cube in cubes]):
            point_bool = [True if elem else False for elem in point]
            if sum(point_bool) == 0:
                raise InputFormatError("Found point not assigned to any atom "
                                       "by the ``bader`` program. Maybe the "
                                       "``-vac off`` option was not set?")
            elif sum(point_bool) > 1:
                raise InputFormatError("Found point assigned to many atoms "
                                       "by the ``bader`` program. Possible "
                                       "numerical inconsistency in algorithm.")
            result.append(point_bool.index(True)+1)

        result = np.array(result)
        result.resize(self.parent_cube.field.grid.points_on_axes)
        return GridField(result, self.parent_cube.field.grid, 'parent_atom',
                         ['qtaim'])


class GridField(Field):

    def __init__(self, values, grid, field_type, field_info=None,
                 check_nans=True):
        self.grid = grid
        super().__init__(values, field_type, field_info, check_nans)

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
        return GridField(dist, self.grid, 'dist', ['ed', isovalue])

    def write_cube(self, output_fn, molecule, charge_type=None,
                   write_coords_in_bohr=True):
        """Write the field as a Gaussian cube file.

        Raises FileExistsError when the file exists.
        """
        with open(output_fn, 'x') as f:
            f.write(' Cube file generated by repESP.\n')
            f.write(' Cube file for field of type {0}.\n'.format(
                self.field_type))
            origin_coords = self.grid.origin_coords
            if write_coords_in_bohr:
                origin_coords = [elem/angstrom_per_bohr for elem in
                                 origin_coords]
            f.write(' {0:4}   {1: .6f}   {2: .6f}   {3: .6f}    1\n'.format(
                len(molecule), *origin_coords))
            for axis in self.grid.axes:
                axis_intervals = axis.intervals
                if write_coords_in_bohr:
                    axis_intervals = [elem/angstrom_per_bohr for elem in
                                      axis_intervals]
                f.write(' {0:4}   {1: .6f}   {2: .6f}   {3: .6f}\n'.format(
                    axis.point_count, *axis_intervals))
            for atom in molecule:
                if charge_type is None:
                    charge = atom.atomic_no
                else:
                    charge = atom.charges[charge_type]
                atom_coords = atom.coords
                if write_coords_in_bohr:
                    atom_coords = [coord/angstrom_per_bohr for coord in
                                   atom_coords]
                f.write(' {0:4}   {1: .6f}   {2: .6f}   {3: .6f}   {4: .6f}\n'
                        .format(atom.atomic_no, charge, *atom_coords))
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
