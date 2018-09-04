from .exceptions import InputFormatError

from abc import ABC, abstractmethod
from typing import Callable, Collection, Generic, Iterator, List, NewType, NamedTuple, Tuple, TypeVar

import functools
import math
import operator


Coord = NewType("Coord", float)
Coords = NewType("Coords", Tuple[Coord, Coord, Coord])


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


FieldValue = TypeVar('FieldValue')


class Field(Generic[FieldValue]):

    def __init__(self, mesh: 'Mesh', values: Collection[FieldValue]) -> None:

        if len(values) != len(mesh):
            raise InputFormatError(
                "Construction of a Field failed due to mismatch between the "
                "number of points ({}) and the number of values ({})".format(
                    len(mesh),
                    len(values)
                )
            )

        self.mesh = mesh
        self.values = values


class Mesh(ABC):

    @abstractmethod
    def points(self) -> Iterator[Coords]:
        pass

    @abstractmethod
    def __len__(self) -> int:
        pass

    def calc_field(
        self,
        field_at_point: Callable[[Coords], FieldValue]
    ) -> Field[FieldValue]:
        """Calculate values at points to a function"""

        # This is a default, inefficient implementation
        values: List[FieldValue] = []

        for point in self.points():
            values.append(
                field_at_point(point)
            )

        return Field(self, values)


class NonGridMesh(Mesh):

    def __init__(self, points_coords: Collection[Coords]) -> None:
        self._points_coords = list(points_coords)

    def points(self) -> Iterator[Coords]:
        return iter(self._points_coords)

    def __len__(self) -> int:
        return len(self._points_coords)


class GridMesh(Mesh):

    class GridMeshAxis(NamedTuple):
        vector: Coords  # Unit vector in xyz coordinates
        point_count: int

    GridMeshAxes = Tuple[GridMeshAxis, GridMeshAxis, GridMeshAxis]

    def __init__(self, origin: Coords, axes: GridMeshAxes) -> None:
        # TODO: Remove this assumption (affects implementation of self.points)
        if (not self._axes_are_aligned_to_coordinate_axes(axes)):
            raise NotImplementedError(
                "GridMesh cannot currently be constructed with axes not aligned"
                " to coordinate axes. The provided axes are: {}".format(axes)
            )

        self.origin = origin
        self.axes = axes

    def _axes_are_aligned_to_coordinate_axes(self, axes: GridMeshAxes) -> bool:
        return functools.reduce(
            operator.and_,
            (math.isclose(vector_component, Coord(0)) for vector_component in [
                axes[0].vector[1],
                axes[0].vector[2],
                axes[1].vector[0],
                axes[1].vector[2],
                axes[2].vector[0],
                axes[2].vector[1],
            ])
        )

    def points(self) -> Iterator[Coords]:
        for i in range(self.axes[0].point_count):
            for j in range(self.axes[1].point_count):
                for k in range(self.axes[2].point_count):
                    yield Coords((
                        Coord(self.origin[0] + i*self.axes[0].vector[0]),
                        Coord(self.origin[1] + j*self.axes[1].vector[1]),
                        Coord(self.origin[2] + k*self.axes[2].vector[2])
                    ))

    def __len__(self) -> int:
        return functools.reduce(
            operator.mul,
            (axis.point_count for axis in self.axes)
        )
