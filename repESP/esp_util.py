"""Parsing and writing .esp file format describing molecular ESP field"""

from repESP.charges import AtomWithCoordsAndCharge, Charge, DipoleMoment, DipoleMomentValue
from repESP.charges import QuadrupoleMoment, QuadrupoleMomentValue
from repESP.fields import Esp, Field, Mesh
from repESP.exceptions import InputFormatError
from repESP.types import AtomWithCoords, Coords, Molecule
from repESP._util import get_line

from dataclasses import dataclass
from fortranformat import FortranRecordWriter as FW, FortranRecordReader as FR
from typing import Callable, cast, List, Pattern, TextIO, Tuple, Type, TypeVar
import re


@dataclass
class GaussianEspData:
    """Dataclass representing the .esp file in the Gaussian format

    Parameters
    ----------
    charge : int
        The total charge of the molecule.
    multiplicity : int
        The multiplicity of the molecule's electron configuration.
    molecule : Molecule[AtomWithCoordsAndCharge]
        The molecule with atom coordinates and partial charges specified.
    dipole_moment : DipoleMoment
        The dipole moment of the molecule.
    quadrupole_moment : QuadrupoleMoment
        The quadrupole moment of the molecule.
    field : Field[Esp]
        The ESP field described by the file.

    Attributes
    ----------
    charge
        See initialization parameter
    multiplicity
        See initialization parameter
    molecule
        See initialization parameter
    dipole_moment
        See initialization parameter
    quadrupole_moment
        See initialization parameter
    field
        See initialization parameter
    """
    charge: int
    multiplicity: int
    molecule: Molecule[AtomWithCoordsAndCharge]
    dipole_moment: DipoleMoment
    quadrupole_moment: QuadrupoleMoment
    field: Field[Esp]


EspDataT = TypeVar('EspDataT', bound='EspData')
@dataclass
class EspData:
    """Dataclass representing the .esp file in the ``resp`` format

    Parameters
    ----------
    atoms_coords : typing.List[Coords]
        The positions of atoms in space. Note that the identities of the atoms
        are not known and thus this member cannot be described with the `Atom` class.
    field : Field[Esp]
        The ESP field around the molecule.

    Attributes
    ----------
    atoms_coords
        See initialization parameter
    field
        See initialization parameter
    """
    atoms_coords: List[Coords]
    field: Field[Esp]

    @classmethod
    def from_gaussian(cls: Type[EspDataT], gaussian_esp_data: GaussianEspData) -> EspDataT:
        """Alternative initialization from .esp file in Gaussian format

        Note that the conversion in the opposite direction is not possible as
        the conversion is lossy.

        Parameters
        ----------
        gaussian_esp_data : GaussianEspData
            A dataclass representing an .esp file in the Gaussian format.
        """
        return cls(
            [atom.coords for atom in gaussian_esp_data.molecule.atoms],
            gaussian_esp_data.field
        )


def parse_gaussian_esp(f: TextIO) -> GaussianEspData:
    """Parse a file in the Gaussian .esp file format

    Parameters
    ----------
    f : TextIO
        File object opened in read mode containing the .esp file to be parsed.
        The file can be generated with Gaussian by specifying the ``IOp(6/50=1)``
        override.

    Raises
    ------
    InputFormatError
        Raised when the file does not follow the expected format. Note that
        this function has only been tested with the output of Gaussian 09.

    Returns
    -------
    GaussianEspData
        A dataclass representing the information in the given .esp file.
    """

    charge, multiplicity, atom_count = _parse_prelude([get_line(f) for i in range(3)])

    molecule = Molecule([_parse_atom(get_line(f)) for _ in range(atom_count)])

    if get_line(f) != " DIPOLE MOMENT:":
        raise InputFormatError("Expected dipole moment section header.")

    dipole_moment = _parse_dipole(get_line(f))

    if get_line(f) != " TRACELESS QUADRUPOLE MOMENT:":
        raise InputFormatError("Expected quadrupole moment section header.")

    quadrupole_moment = _parse_quadrupole([get_line(f), get_line(f)])

    points_header_re = re.compile(" ESP VALUES AND GRID POINT COORDINATES. #POINTS =\s+([0-9]+)")
    points_header_match = points_header_re.match(get_line(f))

    if points_header_match is None:
        raise InputFormatError("Expected ESP points section header.")

    point_count = int(points_header_match.group(1))
    field = _parse_esp_points(f)

    if len(field.mesh) != point_count:
        raise InputFormatError(
            f"The number of ESP points ({len(field.mesh)}) does not agree with that "
            f"specified in section header ({point_count})."
        )

    return GaussianEspData(charge, multiplicity, molecule, dipole_moment, quadrupole_moment, field)


def _parse_prelude(lines: List[str]) -> Tuple[int, int, int]:
    assert len(lines) == 3

    # Line 1
    if lines[0] != " ESP FILE - ATOMIC UNITS":
        raise InputFormatError("Unexpected first line of .esp line.")

    # Line 2
    charge_and_multiplicity_re = re.compile(" CHARGE =\s+([-0-9.]+) - MULTIPLICITY =\s+([0-9.]+)")
    charge_and_multiplicity = charge_and_multiplicity_re.match(lines[1])
    if charge_and_multiplicity is None:
        raise InputFormatError("Failed parsing line 2 (charge and multiplicity expected).")

    charge = int(charge_and_multiplicity.group(1))
    multiplicity = int(charge_and_multiplicity.group(2))

    # Line 3
    atom_count_re = re.compile(" ATOMIC COORDINATES AND ESP CHARGES. #ATOMS =\s+([0-9.]+)")
    atom_count = atom_count_re.match(lines[2])
    if atom_count is None:
        raise InputFormatError("Failed parsing line 3 (molecule header and atom count).")

    return charge, multiplicity, int(atom_count.group(1))


def _parse_atom(line: str) -> AtomWithCoordsAndCharge:
    line_split = line.split()
    return AtomWithCoordsAndCharge.from_symbol(
        line_split[0],
        Coords(tuple(coord.replace('D', 'E') for coord in line_split[1:4])),
        Charge(line_split[4].replace('D', 'E'))
    )


def _parse_dipole(line: str) -> DipoleMoment:
    dipole_line_re = re.compile(" X=\s+([-+0-9.D]+) Y=\s+([-+0-9.D]+) Z=\s+([-+0-9.D]+) Total=\s+([-+0-9.D]+)")
    dipole_line_match = dipole_line_re.match(line)
    if dipole_line_match is None:
        raise InputFormatError("Failed parsing dipole specification.")
    return DipoleMoment(
        DipoleMomentValue(dipole_line_match.group(1).replace('D', 'E')),
        DipoleMomentValue(dipole_line_match.group(2).replace('D', 'E')),
        DipoleMomentValue(dipole_line_match.group(3).replace('D', 'E'))
    )


def _parse_quadrupole(lines: List[str]) -> QuadrupoleMoment:
    assert len(lines) == 2

    line1_components = ("XX", "YY", "ZZ")
    line2_components = ("XY", "XZ", "YZ")
    get_line_re: Callable[[Tuple[str, str, str]], Pattern[str]] = lambda components: re.compile(
        "   {}=\s+([-+0-9.D]+)   {}=\s+([-+0-9.D]+)   {}=\s+([-+0-9.D]+)".format(*components)
    )

    line1_match = get_line_re(line1_components).match(lines[0])
    line2_match = get_line_re(line2_components).match(lines[1])

    if line1_match is None or line2_match is None:
        raise InputFormatError("Failed parsing quadrupole specification.")

    return QuadrupoleMoment(
        QuadrupoleMomentValue(line1_match.group(1).replace('D', 'E')),
        QuadrupoleMomentValue(line1_match.group(2).replace('D', 'E')),
        QuadrupoleMomentValue(line1_match.group(3).replace('D', 'E')),
        QuadrupoleMomentValue(line2_match.group(1).replace('D', 'E')),
        QuadrupoleMomentValue(line2_match.group(2).replace('D', 'E')),
        QuadrupoleMomentValue(line2_match.group(3).replace('D', 'E'))
    )


def _parse_esp_points(f: TextIO) -> Field[Esp]:
    points = []
    values = []
    for line in f:
        line_split = [val.replace('D', 'E') for val in line.split()]
        points.append(Coords(line_split[1:4]))
        values.append(Esp(line_split[0]))

    return Field(
        Mesh(points),
        values
    )


def parse_resp_esp(f: TextIO) -> EspData:
    """Parse a file in the .esp file format defined by ``resp``

    Parameters
    ----------
    f : TextIO
        File object opened in read mode containing the .esp file to be parsed.

    Raises
    ------
    InputFormatError
        Raised when the file does not follow the expected format.

    Returns
    -------
    EspData
        A dataclass representing the information in the given .esp file.
    """

    atom_and_point_count = get_line(f).split()

    if len(atom_and_point_count) != 2:
        raise InputFormatError(
            "Expected atom and point counts on the first line of .esp file in the `resp` format"
        )

    atom_count = int(atom_and_point_count[0])
    point_count = int(atom_and_point_count[1])

    atoms_coords = [Coords(get_line(f).split()) for _ in range(atom_count)]

    mesh_coords: List[Coords] = []
    esp_values: List[Esp] = []

    for _ in range(point_count):
        val, *coords = get_line(f).split()
        mesh_coords.append(Coords(coords))
        esp_values.append(Esp(val))

    field = Field(
        Mesh(
            mesh_coords
        ),
        esp_values
    )

    return EspData(
        atoms_coords,
        field
    )


def write_resp_esp(f: TextIO, esp_data: EspData) -> None:
    """Write a ``resp`` .esp file described by the given input data

    Parameters
    ----------
    f : TextIO
        File object to which the supplied data is to be saved. Must be opened
        in write mode.
    esp_data : EspData
        The dataclass containing the information needed to create a .esp file.
    """

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

    for point_coords, esp_val in zip(field.mesh.points, field.values):
        f.write(
            FW(formats["points"]).write(
                cast(List[float], [esp_val]) + list(point_coords)
            ) + "\n"
        )
