from .types import AtomWithCoords, Ed, Esp, Coords, Field, GridMesh, Molecule
from .types import make_coords, make_ed, make_esp
from .exceptions import InputFormatError

from dataclasses import dataclass
from typing import Callable, Generic, List, TextIO, Tuple


@dataclass
class Cube(Generic[Field.NumericValue]):

    @dataclass
    class Info:
        input_line: str
        title_line: str

    info: Info
    molecule: Molecule[AtomWithCoords]
    electrons_on_atoms: List[float]
    field: Field[Field.NumericValue]


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


def _parse_atom(line: str) -> Tuple[AtomWithCoords, float]:
    identity, cube_charge, *coords = line.split()
    return (
        AtomWithCoords(int(identity), make_coords(*coords)),
        float(cube_charge)
    )


def parse_cube(
    f: TextIO,
    make_value: Callable[[Cube.Info, str], Field.NumericValue]
) -> Cube[Field.NumericValue]:
    # Assumption: coordinates in bohr

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
    atoms_with_electrons = [_parse_atom(get_line()) for i in range(grid_prelude.atom_count)]
    atoms, electrons_on_atoms = zip(*atoms_with_electrons)
    molecule = Molecule(list(atoms))

    # Field values
    value_ctor = lambda x: make_value(info, x)
    # TODO: this may be unfeasible for very large cubes
    values = [value_ctor(x) for x in f.read().split()]

    return Cube(
        info,
        molecule,
        list(electrons_on_atoms),
        Field(grid, values)
    )


def _parse_cube_by_title_common(
    expected_title_start: str,
    value_ctor: Callable[[str], Field.NumericValue],
    verify_title: bool
) -> Callable[[Cube.Info, str], Field.NumericValue]:

    def make_value(info: Cube.Info, value: str) -> Field.NumericValue:
        check_title = lambda title: title.startswith(expected_title_start)
        if verify_title and not check_title(info.title_line):
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

def write_cube(f: TextIO, cube: Cube):

    f.write(f"{cube.info.input_line}\n{cube.info.title_line}\n")

    assert isinstance(cube.field.mesh, GridMesh)

    f.write(' {0:4}   {1: .6f}   {2: .6f}   {3: .6f}    1\n'.format(
        len(cube.molecule.atoms),
        *cube.field.mesh._origin
    ))

    for axis in cube.field.mesh._axes:
        f.write(' {0:4}   {1: .6f}   {2: .6f}   {3: .6f}\n'.format(
            axis.point_count,
            *axis.vector
        ))

    for atom, electron_count in zip(cube.molecule.atoms, cube.electrons_on_atoms):
        f.write(' {0:4}   {1: .6f}   {2: .6f}   {3: .6f}   {4: .6f}\n'.format(
            atom.identity,
            electron_count,
            *atom.coords
        ))

    i = 1
    for value in cube.field.values:
        f.write(' {0: .5E}'.format(value))
        if not i % 6:
            f.write('\n')
        if not i % cube.field.mesh._axes[2].point_count:
            f.write('\n')
            i = 1
        else:
            i += 1
