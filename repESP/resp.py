from fortranformat import FortranRecordWriter

from .cube_helpers import InputFormatError, Atom, Molecule, Field
from .cube_helpers import angstrom_per_bohr


class DuplicateEntryError(Exception):
    pass


class InputValueError(Exception):
    pass


class G09_esp(object):

    def __init__(self, fn, coords_in_bohr=True, allow_dupes=False):
        self._read_in(fn, coords_in_bohr, allow_dupes)

    def _read_in(self, fn, coords_in_bohr, allow_dupes):
        with open(fn, 'r') as f:
            # Checks two first lines
            self._read_header(fn, f)
            self._read_atoms(f, coords_in_bohr)
            self._read_moments(f)
            self._read_esp_points(f, coords_in_bohr, allow_dupes)

    def _read_header(self, fn, f):
        line = f.readline().rstrip('\n')
        if line != " ESP FILE - ATOMIC UNITS":
            raise InputFormatError("The input file {0} does not seem to be the"
                                   " G09 .esp format. Generate by specifying "
                                   "Pop=MK/CHelp(G) with IOp(6/50=1)".format(
                                       fn))
        line = f.readline().split()
        self.charge = int(line[2])
        self.multip = int(line[-1])

    def _read_atoms(self, f, coords_in_bohr):
        line = f.readline().split()
        atom_count = int(line[-1])
        self.molecule = Molecule(self)
        for i in range(atom_count):
            line = f.readline().split()
            identity = line[0]
            atomic_no = Atom.inv_periodic[identity]
            coords = [float(coord.replace('D', 'E')) for coord in line[1:4]]
            # Neglect the ESP value at atoms, which is given by last value
            self.molecule.append(Atom(i+1, atomic_no, coords, coords_in_bohr))

    def _read_moments(self, f):
        assert f.readline().rstrip('\n') == " DIPOLE MOMENT:"
        # Currently not implemented, the lines are just skipped
        for i in range(4):
            f.readline()

    def _read_esp_points(self, f, coords_in_bohr, allow_dupes):
        line = f.readline().split()
        expected = "ESP VALUES AND GRID POINT COORDINATES. #POINTS ="
        assert ' '.join(line[:-1]) == expected
        expected_points_count = int(line[-1])

        points_coords = []
        values = []
        for line in f:
            line = [val.replace('D', 'E') for val in line.split()]
            points_coords.append(tuple(line[1:4]))
            values.append(float(line[0]))

        if len(points_coords) != expected_points_count:
            raise InputFormatError(
                "The number of ESP points {0} does not agree with that "
                "specified at the top of the input file: {1}".format(
                    len(points_coords), expected_points_count))

        try:
            self.field = NonGridField(
                values, points_coords, 'esp', field_info=['input'],
                allow_dupes=allow_dupes, coords_in_bohr=coords_in_bohr)

        except DuplicateEntryError:
            raise InputFormatError(
                "Duplicate points in the input file. This might be an artefact"
                " of the algorithm which produced the points. If these points "
                "are to be counted twice, the NonGridField needs to be called "
                "with `allow_dupes=True`")
        except InputValueError as e:
            # Translate the errors when creating a field to errors due to input
            # file format
            raise InputFormatError(e)


class Points(list):

    def __init__(self, points_coords, coords_in_bohr=False, allow_dupes=True):
        super().__init__()
        self.allow_dupes = allow_dupes
        if not self.allow_dupes:
            self.points_dict = {}

        for point_coords in points_coords:
            self.append(self._check_and_create_point(point_coords,
                                                     coords_in_bohr))

    def _check_and_create_point(self, point_coords, coords_in_bohr):
        if len(point_coords) != 3:
            raise InputValueError(
                "Encountered a point with a number of coordinates {0}, which "
                "is different from 3.".format(len(point_coords)))

        if not self.allow_dupes:
            for point_coord in point_coords:
                if type(point_coord) is not str:
                    raise InputValueError(
                        "If no duplicates are allowed, `points` must be"
                        " a list of tuples of *strings*. Encountered type {0} "
                        "instead.".format(type(point_coord)))

            if point_coords in self.points_dict:
                raise DuplicateEntryError("Encountered a duplicate point: {0}"
                                          .format(point_coords))

        try:
            result = [float(point_coord) for point_coord in point_coords]
        except ValueError:
            raise InputValueError("Couldn't convert coordinates {0} to float."
                                  .format(point_coords))

        if coords_in_bohr:
            result = [angstrom_per_bohr*point_coord for point_coord in result]

        return result


class NonGridField(Field):

    def __init__(self, values, points, field_type, field_info=None,
                 allow_dupes=False, coords_in_bohr=True):
        """Create a NonGridField from given coordinates and values

        Parameters
        ----------
        points : List[Tuple]
            The inner tuples represent coordinates and should hence have
            lengths of 3. If the coordinates are given in Angstrom, the
            `coords_in_bohr` parameter must be set to False. If the
            `allow_dupes` parameter is False (default), the values must be
            given as strings to enable checking for duplicates.

        values : List
            The list of values at *corresponding* coordinates.
        """
        super().__init__(values, field_type, field_info, check_nans=False)

        if type(points) is Points:
            self.points = points
        else:
            self.points = Points(points, coords_in_bohr, allow_dupes)

        if len(points) != len(values):
            raise ValueError("The number of points {0} is different from the "
                             "number of values {1}.".format(len(points),
                                                            len(values)))

    def write_to_file(self, output_fn, molecule, write_coords_in_bohr=True):
        # Numeric formats specified in resp input specification
        # http://upjv.q4md-forcefieldtools.org/RED/resp/#other3
        header_format = FortranRecordWriter('2I5')
        atoms_format = FortranRecordWriter('17X,3E16.7')
        esp_points_format = FortranRecordWriter('1X,4E16.7')

        with open(output_fn, 'x') as f:
            f.write(header_format.write([len(molecule),
                                         len(self.points)]) + "\n")
            for atom in molecule:
                if write_coords_in_bohr:
                    coords = [atom_coord/angstrom_per_bohr for atom_coord in
                              atom.coords]
                else:
                    coords = atom.coords
                f.write(atoms_format.write(coords) + "\n")
            for point_coords, esp_val in zip(self.points, self.values):
                if write_coords_in_bohr:
                    point_coords = [point_coord/angstrom_per_bohr for
                                    point_coord in point_coords]
                f.write(esp_points_format.write([esp_val] + point_coords)+"\n")
