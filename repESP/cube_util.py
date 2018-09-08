from .types import Atom, Ed, Esp, Coords, Field, GridMesh, Molecule
from .types import make_coords, make_ed, make_esp
from .exceptions import InputFormatError

from enum import Enum
from dataclasses import dataclass
from typing import Any, Callable, Generic, List, Mapping, NewType, Optional, TextIO, Tuple, TypeVar, Union


FieldValue = TypeVar('FieldValue')


@dataclass
class CubeInfo:
    input_line: str
    title_line: str


@dataclass
class Cube(Generic[FieldValue]):
    cube_info: CubeInfo
    molecule: Molecule
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
        make_coords(*origin_coords),
        float(nval)
    )


def _parse_grid(origin: Coords, lines: List[str]) -> GridMesh:

    def parse_axis(line: str) -> GridMesh.Axis:
        point_count, *vector_components = line.split()
        return GridMesh.Axis(
            make_coords(*vector_components),
            int(point_count)
        )

    assert len(lines) == 3

    return GridMesh(
        origin,
        GridMesh.Axes(tuple(  # type: ignore # (asserted len is 3)
            parse_axis(line) for line in lines
        ))
    )


def _parse_atom(line: str) -> Atom:
    identity, _cube_charge, *coords = line.split()
    return Atom(int(identity), make_coords(*coords))


def parse_cube(
    f: TextIO,
    make_value: Callable[[CubeInfo, str], FieldValue]
) -> Cube[FieldValue]:
    # Assumption: coordinates in bohr

    get_line = lambda: f.readline().rstrip('\n')

    # Lines 1-2
    cube_info = CubeInfo(input_line=get_line(), title_line=get_line())

    # Line 3
    grid_prelude = _parse_grid_prelude(get_line())

    if float(grid_prelude.nval) != 1:
        # I don't know what NVal means, haven't seen it to be different than 1.
        raise InputFormatError("NVal is different than 1.")

    # Lines 4-6
    grid = _parse_grid(grid_prelude.origin, [get_line() for i in range(3)])

    # Molecule
    molecule = Molecule(list(_parse_atom(get_line()) for i in range(grid_prelude.atom_count)))

    # Field values
    value_ctor = lambda x: make_value(cube_info, x)
    # TODO: this may be unfeasible for very large cubes
    values = [value_ctor(x) for x in f.read().split()]

    return Cube(
        cube_info,
        molecule,
        Field(grid, values)
    )


def _parse_cube_by_title_common(
    expected_title_start: str,
    value_ctor: Callable[[str], FieldValue],
    verify_title: bool
) -> Callable[[CubeInfo, str], FieldValue]:

    def make_value(cube_info: CubeInfo, value: str) -> FieldValue:
        check_title = lambda title: title.startswith(expected_title_start)
        if verify_title and not check_title(cube_info.title_line):
            raise InputFormatError(
                f'Title of cube file does not start with "{expected_title_start}".'
            )
        return value_ctor(value)

    return make_value


def parse_esp_cube(f: TextIO, verify_title=True) -> Cube[Esp]:
    return parse_cube(
        f,
        _parse_cube_by_title_common(" Electrostatic potential", make_esp, verify_title)
    )


def parse_ed_cube(f: TextIO, verify_title=True) -> Cube[Ed]:
    return parse_cube(
        f,
        _parse_cube_by_title_common(" Electron density", make_ed, verify_title)
    )
