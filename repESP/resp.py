from collections import OrderedDict
from fortranformat import FortranRecordWriter

from .cube_helpers import InputFormatError, Atom, Molecule

# http://www.gaussian.com/g_tech/g_ur/k_constants.htm
angstrom_per_bohr = 0.5291772086


class G09_esp(object):

    def __init__(self, fn):
        self._read_in(fn)

    def _read_in(self, fn):
        # Note: distances are assumed to be given (and are written by the
        # write_to_file method) in Bohr radii, but are internally manipulated
        # in Angstroms.
        with open(fn, 'r') as f:
            # Checks two first lines
            self._read_header(fn, f)
            self._read_atoms(f)
            self._read_moments(f)
            self._read_esp_points(f)

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

    def _read_atoms(self, f):
        line = f.readline().split()
        atom_count = int(line[-1])
        self.molecule = Molecule(self)
        for i in range(atom_count):
            line = f.readline().split()
            identity = line[0]
            atomic_no = Atom.inv_periodic[identity]
            coords = line[1:4]
            coords = [angstrom_per_bohr*float(coord.replace('D', 'E')) for
                      coord in line[1:4]]
            # Neglect the ESP value at atoms, which is given by last value
            self.molecule.append(Atom(i+1, atomic_no, coords))

    def _read_moments(self, f):
        assert f.readline().rstrip('\n') == " DIPOLE MOMENT:"
        # Currently not implemented, the lines are just skipped
        for i in range(4):
            f.readline()

    def _read_esp_points(self, f):
        line = f.readline().split()
        expected = "ESP VALUES AND GRID POINT COORDINATES. #POINTS ="
        assert ' '.join(line[:-1]) == expected
        self.esp_points_count = int(line[-1])
        # Use OrderedDict for easy eye-checking if the written values are
        # correct. If it causes performance issues, it can be replaced with a
        # regular dictionary.
        self.esp_points = OrderedDict()
        for line in f:
            line = [float(val.replace('D', 'E')) for val in line.split()]
            coords = (angstrom_per_bohr*val for val in line[1:4])
            esp_val = line[0]
            if coords in self.esp_points:
                # Since they're tuples of *floats*, duplicates may not be
                # spotted this way! TODO
                raise InputFormatError(
                    "Duplicate points in the input file. This might be an "
                    "artefact of the algorithm which produced the points. If "
                    "these points are to be counted twice, the program needs "
                    "to be modified.")
            else:
                self.esp_points[coords] = esp_val
        assert len(self.esp_points) == self.esp_points_count

    def write_to_file(self, fn):
        # Numeric formats specified in resp input specification
        # http://upjv.q4md-forcefieldtools.org/RED/resp/#other3
        header_format = FortranRecordWriter('2I5')
        atoms_format = FortranRecordWriter('17X,3E16.7')
        esp_points_format = FortranRecordWriter('1X,4E16.7')

        # Distances are also written in Bohr
        with open(fn, 'x') as f:
            f.write(header_format.write([len(self.molecule),
                                         len(self.esp_points)]) + "\n")
            for atom in self.molecule:
                coords = [val/angstrom_per_bohr for val in atom.coords]
                f.write(atoms_format.write(coords) + "\n")
            for esp_coords, esp_val in self.esp_points.items():
                esp_coords = [val/angstrom_per_bohr for val in esp_coords]
                f.write(esp_points_format.write([esp_val] + esp_coords) + "\n")
