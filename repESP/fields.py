"""Types used to describe fields as values at selected points in space"""

from .exceptions import InputFormatError
from .types import Coords, Dist

from abc import ABC, abstractmethod
from dataclasses import dataclass, field, InitVar
from typing import Any, cast, Collection, Generic, Iterator, List, NewType, Tuple, TypeVar

import functools
import math
import operator


class Esp(float):

    """Electrostatic potential [atomic units]"""

    __slots__ = ()

    def __new__(cls, x: Any):
        return super().__new__(cls, float(x))  # type: ignore # (Too many arguments for "__new__" of "object")


class Ed(float):

    """Electron density [atomic units]"""

    __slots__ = ()

    def __new__(cls, x: Any):
        return super().__new__(cls, float(x))  # type: ignore # (Too many arguments for "__new__" of "object")


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
