from .exceptions import InputFormatError
from .util import _elements, _get_symbol, _get_atomic_number

from abc import ABC, abstractmethod
from dataclasses import dataclass, field, InitVar
from typing import Any, Callable, cast, Collection, Generic, Iterator, List, NewType, Tuple, TypeVar

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

    def __post_init__(self):
        if self.identity < 1 or self.identity >= len(_elements):
            raise ValueError("Atomic number is not within expected bounds.")

    @property
    def symbol(self) -> str:
        return _get_symbol(self.identity)

    @classmethod
    def from_symbol(cls, symbol: str, *args, **kwargs):
        # Generic type annotations as per:
        # https://github.com/python/typing/issues/58#issuecomment-326240794)
        # don't seem to be working for a dataclass, and that's even before
        # getting them to work with args and kwargs.
        return cls(_get_atomic_number(symbol), *args, **kwargs)  # type: ignore # (args, kwargs)


@dataclass
class AtomWithCoords(Atom):
    coords: Coords


GenericAtom = TypeVar('GenericAtom', bound=Atom)


@dataclass
class Molecule(Generic[GenericAtom]):
    atoms: List[GenericAtom]


FieldValue = TypeVar('FieldValue')

@dataclass
class Field(Generic[FieldValue]):

    mesh: 'Mesh'
    values_: InitVar[Collection[FieldValue]]
    values: List[FieldValue] = field(init=False)

    def __post_init__(self, values_) -> None:

        if len(values_) != len(self.mesh):
            raise InputFormatError(
                f"Construction of a Field failed due to mismatch between the "
                f"number of points ({len(self.mesh)}) and the number of values ({len(values_)})"
            )

        self.values = list(values_)

    # TODO: This would ideally be extended to numbers.Number but mypy throws errors.
    NumericValue = TypeVar('NumericValue', bound=float)

    def __add__(self: 'Field[NumericValue]', other: 'Field[NumericValue]') -> 'Field[NumericValue]':

        if not isinstance(other, Field):
            raise TypeError(
                "unsupported operand type(s) for +: 'Field' and 'type(other)"
            )

        if self.mesh != other.mesh:
            raise ValueError(
                "Cannot add or subtract Fields with different meshes."
            )

        return Field(
            self.mesh,
            [
                cast(Field.NumericValue, value_self + value_other)
                for value_self, value_other in zip(self.values, other.values)
            ]
        )

    def __neg__(self: 'Field[NumericValue]') -> 'Field[NumericValue]':
        return Field(
            self.mesh,
            [
                cast(Field.NumericValue, -value)
                for value in self.values
            ]
        )

    def __sub__(self: 'Field[NumericValue]', other: 'Field[NumericValue]') -> 'Field[NumericValue]':
        return self + (-other)

    # TODO: Could add div and sub but it's not needed at the moment.
    # TODO: __iadd__ can be added as an optimization (unless we decide to
    # freeze the dataclass.


class Mesh(ABC):

    @property
    @abstractmethod
    def points(self) -> Iterator[Coords]:
        pass

    @abstractmethod
    def __len__(self) -> int:
        pass


@dataclass
class NonGridMesh(Mesh):

    points_: InitVar[Collection[Coords]]
    _points: List[Coords] = field(init=False)

    def __post_init__(self, points_) -> None:
        self._points = list(points_)

    @property
    def points(self) -> Iterator[Coords]:
        return iter(self._points)

    def __len__(self) -> int:
        return len(self._points)


@dataclass
class GridMesh(Mesh):

    @dataclass
    class Axis:
        vector: Coords  # Unit vector in xyz coordinates
        point_count: int

    Axes = NewType("Axes", Tuple[Axis, Axis, Axis])

    origin: Coords
    axes: Axes

    def __post_init__(self) -> None:
        # TODO: Remove this assumption (affects implementation of self.points)
        if (not self._axes_are_aligned_to_coordinate_axes(self.axes)):
            raise NotImplementedError(
                f"GridMesh cannot currently be constructed with axes not aligned"
                f" to coordinate axes. The provided axes are: {self.axes}"
            )

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

    @property
    def points(self) -> Iterator[Coords]:
        for i in range(self.axes[0].point_count):
            for j in range(self.axes[1].point_count):
                for k in range(self.axes[2].point_count):
                    yield Coords((
                        Dist(self.origin[0] + i*self.axes[0].vector[0]),
                        Dist(self.origin[1] + j*self.axes[1].vector[1]),
                        Dist(self.origin[2] + k*self.axes[2].vector[2])
                    ))

    def __len__(self) -> int:
        return functools.reduce(
            operator.mul,
            (axis.point_count for axis in self.axes)
        )
