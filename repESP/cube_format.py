"""Parsing and writing Gaussian "cube" format describing molecular fields"""

from repESP.charges import AtomWithCoordsAndCharge, Charge
from repESP.fields import Ed, Esp, Field, FieldValue, GridMesh
from repESP.exceptions import InputFormatError
from repESP.types import Coords, Coords, Molecule

from dataclasses import dataclass
from typing import Callable, Generic, List, TextIO, Tuple


@dataclass
class Cube(Generic[FieldValue]):
    """Dataclass representing information in a Gaussian "cube" file

    This class is generic in the type of the values of the field described by
    the cube file. This type can be any type, as can `fields.FieldValue`.

    Parameters
    ----------
    info : Info
        Additional, less structured information about the cube file.
    molecule : Molecule[AtomWithCoordsAndCharge]
        The molecule which field is described by the cube file. Atoms have
        coordinates and partial charges specified.
    field : Field[FieldValue]
        The field described by the cube file.

    Attributes
    ----------
    info
        See initialization parameter
    molecule
        See initialization parameter
    field
        See initialization parameter
    """

    @dataclass
    class Info:
        """Additional, less structured information about the cube file.

        Parameters
        ----------
        input_line: str
            The top line of a cube file, typically containing information
            regarding the Gaussian input parameters used.
        title_line: str
            The second line of a cube file, containing a free-form title.

        Attributes
        ----------
        input_line
            See initialization parameter
        title_line
            See initialization parameter
        """

        input_line: str
        title_line: str

    info: Info
    molecule: Molecule[AtomWithCoordsAndCharge]
    field: Field[FieldValue]


@dataclass
class _GridPrelude:
    atom_count: int
    origin: Coords
    nval: float


def _parse_grid_prelude(line: str) -> _GridPrelude:
    line_split = line.split()
    if len(line_split) in (4, 5) :
        atom_count, *origin_coords = line_split[:4]
        nval = line_split[4] if len(line_split) == 5 else "1"
    else:
        raise InputFormatError(
            f"Cube file incorrectly formatted! Expected four or five fields "
            "(atom count, 3*origin coordinates, [NVal]) on line 3, found "
            "{len(line_split)} fields."
        )

    return _GridPrelude(
        int(atom_count),
        Coords(origin_coords),
        float(nval)
    )


def _parse_grid(origin: Coords, lines: List[str]) -> GridMesh:

    def parse_axis(line: str) -> GridMesh.Axis:
        # `cubegen` documentation, which contains the only definition of
        # the cube format, specifies that a negative `point_count` may be used
        # in the *input* to signify atomic units. However, this is not an
        # option in the output of `cubegen` and values are always in atomic units.
        point_count, *vector_components = line.split()
        return GridMesh.Axis(
            Coords(vector_components),
            int(point_count)
        )

    assert len(lines) == 3

    return GridMesh(
        origin,
        GridMesh.Axes(tuple(  # type: ignore # (asserted len is 3)
            parse_axis(line) for line in lines
        ))
    )


def _parse_atom(line: str) -> AtomWithCoordsAndCharge:
    atomic_number, cube_charge, *coords = line.split()
    return AtomWithCoordsAndCharge(int(atomic_number), Coords(coords), Charge(cube_charge))


def parse_cube(
    f: TextIO,
    make_value: Callable[[Cube.Info, str], FieldValue]
) -> Cube[FieldValue]:
    """Parse a file in the Gaussian "cube" file format

    You probably mean to use `parse_ed_cube` or `parse_esp_cube` unless
    your cube file is of neither of those types.

    Note that the values are expected to be space separated. If your cube file
    comes from elsewhere than Gaussian, you should ensure that the coordinates
    are given in bohr.

    Parameters
    ----------
    f : TextIO
        File object opened in read mode containing the cube file to be parsed.
    make_value : Callable[[Cube.Info, str], FieldValue]
        A function taking two parameters: the cube information and a string
        representing the field value. The function should parse the field value
        into the desired internal representation, for example an `Esp` object.
        The cube information is provided in case verification of the cube file
        type is required.

        Example
        -------

        In the simplest case this could be::

            lambda _, str_: float(str_)

        which ignores the cube information (thus performing no verification)
        and simply parses the string value as a float.

    Raises
    ------
    InputFormatError
        Raised when the file does not follow the expected format.

    Returns
    -------
    Cube[FieldValue]
        Data from the parsed cube file.
    """

    get_line = lambda: f.readline().rstrip('\n')

    # Lines 1-2
    info = Cube.Info(input_line=get_line(), title_line=get_line())

    # Line 3
    grid_prelude = _parse_grid_prelude(get_line())

    if float(grid_prelude.nval) != 1:
        # I don't know what NVal means, haven't seen it to be different than 1.
        raise InputFormatError("NVal is different than 1.")

    # Lines 4-6
    grid = _parse_grid(grid_prelude.origin, [get_line() for i in range(3)])

    # Molecule
    molecule = Molecule([
        _parse_atom(get_line()) for i in range(grid_prelude.atom_count)
    ])

    # Field values
    value_ctor = lambda x: make_value(info, x)
    values = [value_ctor(x) for x in f.read().split()]

    return Cube(
        info,
        molecule,
        # The implicit assumption here is that the order of points in `grid`
        # is the same as the order of `values`. This is correct, as the order
        # of points in a GridMesh is the same as that in a cube file.
        Field(grid, values)
    )


def _parse_cube_by_title_common(
    expected_title_start: str,
    value_ctor: Callable[[str], FieldValue],
    verify_title: bool
) -> Callable[[Cube.Info, str], FieldValue]:

    def make_value(info: Cube.Info, value: str) -> FieldValue:
        check_title = lambda title: title.startswith(expected_title_start)
        if verify_title and not check_title(info.title_line):
            raise InputFormatError(
                f'Title of cube file does not start with "{expected_title_start}".'
            )
        return value_ctor(value)

    return make_value


def parse_esp_cube(f: TextIO, verify_title=True) -> Cube[Esp]:
    """Parse a Gaussian "cube" file describing an ESP field

    If your cube file comes from elsewhere than Gaussian, you should ensure
    that the coordinates are given in bohr and ESP values in atomic units.

    Parameters
    ----------
    f : TextIO
        File object opened in read mode containing the cube file to be parsed.
    verify_title : bool, optional
        If this flag is set to True (default), an `InputFormatError` will
        be raised if the cube title does not start with the string
        ``" Electrostatic potential"``.

    Returns
    -------
    Cube[Esp]
        Data from the parsed cube file.
    """
    return parse_cube(
        f,
        _parse_cube_by_title_common(" Electrostatic potential", Esp, verify_title)
    )


def parse_ed_cube(f: TextIO, verify_title=True) -> Cube[Ed]:
    """Parse a Gaussian "cube" file describing electron density field

    If your cube file comes from elsewhere than Gaussian, you should ensure
    that the coordinates are given in bohr and electron density values in
    atomic units.

    Parameters
    ----------
    f : TextIO
        File object opened in read mode containing the cube file to be parsed.
    verify_title : bool, optional
        If this flag is set to True (default), an `InputFormatError` will
        be raised if the cube title does not start with the string
        ``" Electron density"``.

    Returns
    -------
    Cube[Ed]
        Data from the parsed cube file.
    """
    return parse_cube(
        f,
        _parse_cube_by_title_common(" Electron density", Ed, verify_title)
    )

def write_cube(f: TextIO, cube: Cube[Field.NumericValue]) -> None:
    """Write a Gaussian "cube" file described by the given input data

    Parameters
    ----------
    f : TextIO
        File object to which the supplied data is to be saved. Must be opened
        in write mode.
    cube : Cube[Field.NumericValue]
        The dataclass containing the information needed to create a cube
        file. Note that only cube files describing fields with values matching
        `Field.NumericValue` are supported.
    """

    f.write(f"{cube.info.input_line}\n{cube.info.title_line}\n")

    assert isinstance(cube.field.mesh, GridMesh)

    f.write(' {0:4}   {1: .6f}   {2: .6f}   {3: .6f}    1\n'.format(
        len(cube.molecule.atoms),
        *cube.field.mesh.origin
    ))

    for axis in cube.field.mesh.axes:
        f.write(' {0:4}   {1: .6f}   {2: .6f}   {3: .6f}\n'.format(
            axis.point_count,
            *axis.vector
        ))

    for atom in cube.molecule.atoms:
        f.write(' {0:4}   {1: .6f}   {2: .6f}   {3: .6f}   {4: .6f}\n'.format(
            atom.atomic_number,
            atom.charge,
            *atom.coords
        ))

    i = 1
    for value in cube.field.values:
        f.write(' {0: .5E}'.format(value))
        if not i % 6:
            f.write('\n')
        if not i % cube.field.mesh.axes[2].point_count:
            f.write('\n')
            i = 1
        else:
            i += 1
