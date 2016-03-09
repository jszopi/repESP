from fortranformat import FortranRecordWriter
import textwrap

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
            raise InputFormatError(
                "The input file {0} does not seem to be the G09 .esp format. "
                "Generate by specifying Pop=MK/CHelp(G) with IOp(6/50=1)"
                .format(fn))
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
            points = Points(points_coords, coords_in_bohr, allow_dupes)
        except DuplicateEntryError:
            raise InputFormatError(
                "Duplicate points in the input file. This might be an artefact"
                " of the algorithm which produced the points. If these points "
                "are to be counted twice, the NonGridField needs to be called "
                "with `allow_dupes=True`")
        except InputValueError as e:
            # Translate the errors when creating Points to errors due to input
            # file format
            raise InputFormatError(e)

        self.field = NonGridField(values, points, 'esp', field_info=['input'])


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

    def __init__(self, values, points, field_type, field_info=None):
        """Create a NonGridField from given coordinates and values

        Parameters
        ----------
        points : Points

        values : List
            The list of values at *corresponding* coordinates.
        """
        super().__init__(values, field_type, field_info, check_nans=False)
        self.points = points

        if type(points) is not Points:
            raise TypeError("Expected type Points for the points argument, "
                            "instead got {0}".format(type(points)))

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


def read_respin2(fn, ref_molecule=None):
    with open(fn, 'r') as inp:
        # Rewind to end of cntrl section
        line = inp.readline()
        while "&end" not in line:
            line = inp.readline()
        # Skip two lines ...
        for i in range(3):
            line = inp.readline()
        # ... and the third one will be `charge, iuniq`
        charge, iuniq = [int(elem) for elem in line.split()]

        # Create a molecule
        molecule = Molecule(None)
        for i, line in enumerate(inp):
            if len(line.split()) != 2:
                break
            atom = Atom(i+1, int(line.split()[0]))
            # Crucial bit: reading in ivary
            atom.ivary = int(line.split()[1])
            molecule.append(atom)

    # Check input file consistency
    if len(molecule) != iuniq:
        raise InputFormatError("The number of atoms {0} doesn't agree with"
                               " the `iuniq` value in the input file: {1}"
                               .format(len(molecule), iuniq))
    # Check the molecule against a reference molecule
    if ref_molecule is not None and molecule != ref_molecule:
        molecule.verbose_compare(ref_molecule)
        raise InputFormatError("The molecule in the .respin2 file differs "
                               "from the other molecule as shown above.")

    return molecule, charge, iuniq


def write_modified_respin2(molecule, charge, iuniq, fn_out, check_ivary=True):
    with open(fn_out, 'w') as out:
        out.write(textwrap.dedent(
            """\
            Resp charges for organic molecule

             &cntrl

             nmol = 1,
             ihfree = 1,
             ioutopt = 1,
             iqopt = 2,

             &end
                1.0
            Resp charges for organic molecule
            """))
        numbers = FortranRecordWriter('2I5')
        # `charge, iuniq` line
        print(numbers.write([charge, iuniq]), file=out)

        if check_ivary:
            print("\nPlease check if the following generated RESP input is "
                  "what you want. Note that the hydrogens to be equivalenced "
                  "were selected automatically by the program which generated "
                  "the `.respin1` file (likely `respgen`).\n")
        for atom in molecule:
            # Freeze non-hydrogens
            ivary = atom.ivary if atom.atomic_no == 1 else -1
            if check_ivary:
                _print_ivary_action(atom, ivary, molecule)
            print(numbers.write([atom.atomic_no, ivary]), file=out)

        print(file=out)


def _print_ivary_action(atom, ivary, molecule):
    print(atom, end='')
    if ivary == -1:
        print(", frozen")
    elif ivary > 0:
        print(", equivalenced to atom", molecule[ivary-1].label)
    else:
        print()
