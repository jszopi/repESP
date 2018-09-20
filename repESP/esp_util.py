from .types import Atom, Esp, Coords, NonGridMesh, Molecule, NumericField, NumericFieldValue
from .types import make_coords, make_esp
from .charges import Charge, Dipole, MoleculeWithCharges, Quadrupole
from .charges import make_charge, make_dipole_moment, make_quadrupole_moment
from .exceptions import InputFormatError

from dataclasses import dataclass
from fortranformat import FortranRecordWriter as FW, FortranRecordReader as FR
from typing import List, TextIO, Tuple
import re


@dataclass
class GaussianEspData:
    # Gaussian 09 .esp file format generated with IOp(6/50=1).
    charge: float
    multiplicity: int
    molecule_with_charges: MoleculeWithCharges
    dipole: Dipole
    quadrupole: Quadrupole
    field: NumericField


@dataclass
class EspData:
    atoms_coords: List[Coords]
    field: NumericField

    @classmethod
    def from_gaussian(cls, gaussian_esp_data: GaussianEspData):
        return EspData(
            [atom.coords for atom in gaussian_esp_data.molecule_with_charges.molecule.atoms],
            gaussian_esp_data.field
        )


def parse_gaussian_esp(f: TextIO) -> GaussianEspData:

    get_line = lambda: f.readline().rstrip('\n')

    charge, multiplicity, atom_count = _parse_prelude([get_line() for i in range(3)])

    atoms_and_charges = [_parse_atom(get_line()) for _ in range(atom_count)]
    molecule_with_charges = MoleculeWithCharges(
        Molecule([atom_and_charge[0] for atom_and_charge in atoms_and_charges]),
        [atom_and_charge[1] for atom_and_charge in atoms_and_charges]
    )

    if get_line() != " DIPOLE MOMENT:":
        raise InputFormatError("Expected dipole moment section header.")

    dipole = _parse_dipole(get_line())

    if get_line() != " TRACELESS QUADRUPOLE MOMENT:":
        raise InputFormatError("Expected quadrupole moment section header.")

    quadrupole = _parse_quadrupole([get_line(), get_line()])

    points_header_re = re.compile(" ESP VALUES AND GRID POINT COORDINATES. #POINTS =\s+([0-9]+)")
    points_header_match = points_header_re.match(get_line())

    if points_header_match is None:
        raise InputFormatError("Expected ESP points section header.")

    point_count = int(points_header_match.group(1))
    field = _parse_esp_points(f)

    if len(field.mesh) != point_count:
        raise InputFormatError(
            f"The number of ESP points ({len(field.mesh)}) does not agree with that "
            f"specified in section header ({point_count})."
        )

    return GaussianEspData(charge, multiplicity, molecule_with_charges, dipole, quadrupole, field)


def _parse_prelude(lines: List[str]) -> Tuple[float, int, int]:
    assert len(lines) == 3

    # Line 1
    if lines[0] != " ESP FILE - ATOMIC UNITS":
        raise InputFormatError("Unexpected first line of .esp line.")

    # Line 2
    charge_and_multiplicity_re = re.compile(" CHARGE =\s+([-0-9.]+) - MULTIPLICITY =\s+([0-9.]+)")
    charge_and_multiplicity = charge_and_multiplicity_re.match(lines[1])
    if charge_and_multiplicity is None:
        raise InputFormatError("Failed parsing line 2 (charge and multiplicity expected).")

    charge = float(charge_and_multiplicity.group(1))
    multiplicity = int(charge_and_multiplicity.group(2))

    # Line 3
    atom_count_re = re.compile(" ATOMIC COORDINATES AND ESP CHARGES. #ATOMS =\s+([0-9.]+)")
    atom_count = atom_count_re.match(lines[2])
    if atom_count is None:
        raise InputFormatError("Failed parsing line 3 (molecule header and atom count).")

    return charge, multiplicity, int(atom_count.group(1))


def _parse_atom(line: str) -> Tuple[Atom, Charge]:
    line_split = line.split()
    return (
        Atom.from_symbol(
            line_split[0],
            make_coords(*(coord.replace('D', 'E') for coord in line_split[1:4]))
        ),
        make_charge(line_split[4].replace('D', 'E'))
    )


def _parse_dipole(line: str) -> Dipole:
    dipole_line_re = re.compile(" X=\s+([-0-9.D]+) Y=\s+([-0-9.D]+) Z=\s+([-0-9.D]+) Total=\s+([-0-9.D]+)")
    dipole_line_match = dipole_line_re.match(line)
    if dipole_line_match is None:
        raise InputFormatError("Failed parsing dipole specification.")
    # TODO: Are the units atomic or Debye?
    return Dipole(
        make_dipole_moment(dipole_line_match.group(1).replace('D', 'E')),
        make_dipole_moment(dipole_line_match.group(2).replace('D', 'E')),
        make_dipole_moment(dipole_line_match.group(3).replace('D', 'E'))
    )


def _parse_quadrupole(lines: List[str]) -> Quadrupole:
    assert len(lines) == 2

    line1_components = ("XX", "YY", "ZZ")
    line2_components = ("XY", "XZ", "YZ")
    get_line_re = lambda components: re.compile(
        "   {}=\s+([-0-9.D]+)   {}=\s+([-0-9.D]+)   {}=\s+([-0-9.D]+)".format(*components)
    )

    line1_match = get_line_re(line1_components).match(lines[0])
    line2_match = get_line_re(line2_components).match(lines[1])

    if line1_match is None or line2_match is None:
        raise InputFormatError("Failed parsing quadrupole specification.")

    # TODO: Are the units atomic or other?
    return Quadrupole(
        make_quadrupole_moment(line1_match.group(1).replace('D', 'E')),
        make_quadrupole_moment(line1_match.group(2).replace('D', 'E')),
        make_quadrupole_moment(line1_match.group(3).replace('D', 'E')),
        make_quadrupole_moment(line2_match.group(1).replace('D', 'E')),
        make_quadrupole_moment(line2_match.group(2).replace('D', 'E')),
        make_quadrupole_moment(line2_match.group(3).replace('D', 'E'))
    )


def _parse_esp_points(f: TextIO) -> NumericField[Esp]:
    points = []
    values = []
    for line in f:
        line_split = [val.replace('D', 'E') for val in line.split()]
        points.append(make_coords(*line_split[1:4]))
        values.append(make_esp(line_split[0]))

    return NumericField(
        NonGridMesh(points),
        values
    )


def parse_resp_esp(f: TextIO) -> EspData:

    get_line = lambda: f.readline().rstrip('\n')

    atom_and_point_count = get_line().split()

    if len(atom_and_point_count) != 2:
        raise InputFormatError(
            "Expected atom and point counts on the first line of .esp file in the `resp` format"
        )

    atom_count = int(atom_and_point_count[0])
    point_count = int(atom_and_point_count[1])

    atoms_coords = [make_coords(*get_line().split()) for _ in range(atom_count)]

    mesh_coords: List[Coords] = []
    esp_values: List[Esp] = []

    for _ in range(point_count):
        val, *coords = get_line().split()
        mesh_coords.append(make_coords(*coords))
        esp_values.append(make_esp(val))

    field = NumericField(
        NonGridMesh(
            mesh_coords
        ),
        esp_values
    )

    return EspData(
        atoms_coords,
        field
    )


def write_resp_esp(f: TextIO, esp_data: EspData):

    atoms_coords, field = esp_data.atoms_coords, esp_data.field

    # Numeric Fortran formats specified in resp input specification
    # http://upjv.q4md-forcefieldtools.org/RED/resp/#other3
    formats = {
        "header": "2I5",
        "atoms": "17X,3E16.7",
        "points": "1X,4E16.7",
    }

    f.write(
        FW(formats["header"]).write(
            [
                len(atoms_coords),
                len(field.mesh)
            ]
        ) + "\n"
    )

    for atom_coords in atoms_coords:
        f.write(FW(formats["atoms"]).write(atom_coords) + "\n")

    for point_coords, esp_val in zip(field.mesh.points(), field.values):
        f.write(
            FW(formats["points"]).write(
                [esp_val] + list(point_coords)
            ) + "\n"
        )
