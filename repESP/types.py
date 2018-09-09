from .exceptions import InputFormatError

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any, Callable, Collection, Generic, Iterator, List, NewType, Tuple, TypeVar

import functools
import math
import operator


# TODO: It would be neat if the __repr__ of these could be overriden, so that
# it's different from float. May not be possible without creating a class though.
Dist = NewType("Dist", float)  # Distance [bohr]
Coords = NewType("Coords", Tuple[Dist, Dist, Dist])

# NewType allows to avoid the overhead of creating a class, but that's at the
# cost of having custom constructors.
# https://github.com/python/typing/issues/415#issuecomment-297401553
def make_dist(x: Any) -> Dist:
    return Dist(float(x))


def make_coords(x: Any, y: Any, z: Any) -> Coords:
    return Coords((make_dist(x), make_dist(y), make_dist(z)))


Esp = NewType("Esp", float)  # Electrostatic potential [atomic units]
Ed = NewType("Ed", float)  # Electron density [atomic units]


def make_esp(x: Any) -> Esp:
    return Esp(float(x))


def make_ed(x: Any) -> Ed:
    return Ed(float(x))


@dataclass
class Atom:
    identity: int  # Atomic number
    coords: Coords


@dataclass
class Molecule:
    atoms: List[Atom]


Charge = NewType("Charge", float)  # Atomic charge [elementary charge]


def make_charge(x: Any) -> Charge:
    return Charge(float(x))


@dataclass(init=False)
class MoleculeWithCharges:

    molecule: Molecule
    charges: List[Charge]

    def __init__(self, molecule: Molecule, charges: Collection[Charge]) -> None:
        if len(molecule.atoms) != len(charges):
            raise InputFormatError(
                "Construction of MoleculeWithCharges failed due to mismatch between the "
                "number of atoms in molecule ({}) and the number of charges ({})".format(
                    len(molecule.atoms),
                    len(charges)
                )
            )

        self.charges = list(charges)
        self.molecule = molecule


FieldValue = TypeVar('FieldValue')


@dataclass(init=False)
class Field(Generic[FieldValue]):

    mesh: 'Mesh'
    values: Collection[FieldValue]

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
        # A default, possibly inefficient implementation
        return Field(self, [field_at_point(point) for point in self.points()])


class NonGridMesh(Mesh):

    def __init__(self, points_coords: Collection[Coords]) -> None:
        self._points_coords = list(points_coords)

    def points(self) -> Iterator[Coords]:
        return iter(self._points_coords)

    def __len__(self) -> int:
        return len(self._points_coords)


class GridMesh(Mesh):

    @dataclass
    class Axis:
        vector: Coords  # Unit vector in xyz coordinates
        point_count: int

    Axes = NewType("Axes", Tuple[Axis, Axis, Axis])

    def __init__(self, origin: Coords, axes: Axes) -> None:
        # TODO: Remove this assumption (affects implementation of self.points)
        if (not self._axes_are_aligned_to_coordinate_axes(axes)):
            raise NotImplementedError(
                "GridMesh cannot currently be constructed with axes not aligned"
                " to coordinate axes. The provided axes are: {}".format(axes)
            )

        self._origin = origin
        self._axes = axes

    @staticmethod
    def _axes_are_aligned_to_coordinate_axes(axes: Axes) -> bool:
        return functools.reduce(
            operator.and_,
            (math.isclose(vector_component, Dist(0)) for vector_component in [
                axes[0].vector[1],
                axes[0].vector[2],
                axes[1].vector[0],
                axes[1].vector[2],
                axes[2].vector[0],
                axes[2].vector[1],
            ])
        )

    def points(self) -> Iterator[Coords]:
        for i in range(self._axes[0].point_count):
            for j in range(self._axes[1].point_count):
                for k in range(self._axes[2].point_count):
                    yield Coords((
                        Dist(self._origin[0] + i*self._axes[0].vector[0]),
                        Dist(self._origin[1] + j*self._axes[1].vector[1]),
                        Dist(self._origin[2] + k*self._axes[2].vector[2])
                    ))

    def __len__(self) -> int:
        return functools.reduce(
            operator.mul,
            (axis.point_count for axis in self._axes)
        )
