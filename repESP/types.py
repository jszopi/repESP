from exceptions import InputFormatError

from abc import ABC, abstractmethod
from typing import Collection, Generic, Iterator, List, NamedTuple, Tuple, TypeVar

import functools
import math
import operator


class _Coord(float):

    # TODO: Should probably implement __gt__ because otherwise non-equality
    # comparisons will not honour approximate equality.
    def __eq__(self, other):
        print("Cmp: {} v {}".format(self, other))
        return math.isclose(self, other, abs_tol=1e-6)


class Coords(tuple):

    def __new__(cls, x: float, y: float, z: float):
        return super().__new__(cls, (_Coord(x), _Coord(y), _Coord(z)))


class Atom(NamedTuple):
    identity: int  # Atomic number
    coords: Coords


class Molecule(NamedTuple):
    atoms: List[Atom]


class Charges:

    def __init__(self, molecule: Molecule, charge_list: Collection[float]) -> None:
        if len(molecule.atoms) != len(charge_list):
            raise InputFormatError(
                "Construction of Charges failed due to mismatch between the "
                "number of atoms in molecule ({}) and the number of charges ({})".format(
                    len(molecule.atoms),
                    len(charge_list)
                )
            )

        self.values = list(charge_list)
        self.molecule = molecule


class Mesh(ABC):

    @abstractmethod
    def points(self) -> Iterator[Coords]:
        pass

    @abstractmethod
    def __eq__(self, other):
        pass

    def __ne__(self, other):
        return not self == other

    @abstractmethod
    def __len__(self) -> int:
        pass


class NonGridMesh(Mesh):

    def __init__(self, points_coords: Collection[Coords]) -> None:
        self._points_coords = list(points_coords)

    def points(self) -> Iterator[Coords]:
        return iter(self._points_coords)

    def __eq__(self, other) -> bool:
        return self._points_coords == other._points_coords

    def __len__(self) -> int:
        return len(self._points_coords)


class GridAxis(NamedTuple):
    vector: Coords  # Unit vector in xyz coordinates
    interval: float  # In multiples of unit vector
    point_count: int


GridAxes = Tuple[GridAxis, GridAxis, GridAxis]

class GridMesh(Mesh):

    def __init__(self, origin: Coords, axes: GridAxes) -> None:
        if axes != (Coords(1, 0, 0), Coords(0, 1, 0), Coords(0, 0, 1)):
            raise NotImplementedError(
                "GridMesh cannot currently be constructed with axes not aligned"
                " to coordinate axes. The provided axes are: {}".format(axes)
            )

        self.origin = origin
        self.axes = axes


    def points(self) -> Iterator[Coords]:
        for i in range(self.axes[0].point_count):
            for j in range(self.axes[1].point_count):
                for k in range(self.axes[2].point_count):
                    yield (
                        self.origin[0] + i*self.axes[0].interval,
                        self.origin[1] + j*self.axes[1].interval,
                        self.origin[2] + k*self.axes[2].interval
                    )

    def __eq__(self, other) -> bool:
        return self.origin == other.origin and self.axes == other.axes

    def __len__(self) -> int:
        return functools.reduce(
            operator.mul,
            (axis.point_count for axis in self.axes)
        )

FieldValue = TypeVar('FieldValue')

class Field(Generic[FieldValue]):

    def __init__(self, mesh: Mesh, values: Collection[FieldValue]) -> None:

        if len(values) != len(mesh):
            raise InputFormatError(
                "Construction of a Field failed due to mismatch between the "
                "number of points ({}) and the number of values ({})".format(
                    len(mesh),
                    len(values)
                )
            )

        self._mesh = mesh
        self._values = values

    def __eq__(self, other):

        self._mesh == other._mesh
        self._values == other._values
