"""Types used to describe fields as values at selected points in space

Attributes
----------
"""

from .types import Coords, Dist

from abc import ABC, abstractmethod
from dataclasses import dataclass, field, InitVar
from typing import Any, cast, Collection, Generic, Iterator, List, NewType, Tuple, TypeVar

import functools
import math
import operator


class Esp(float):
    """Electrostatic potential value in atomic units (:math:`E_h/e`)

    Parameters
    ----------
    value : Any
        Any value convertible to float representing the value in atomic units.
    """

    __slots__ = ()

    def __new__(cls, value: Any):
        return super().__new__(cls, float(value))  # type: ignore # (Too many arguments for "__new__" of "object")

    def __repr__(self) -> str:
        return f"Esp({super().__repr__()})"

    def __str__(self) -> str:
        return f"{super().__str__()} a.u."


class Ed(float):
    """Electron density value in atomic units (\ :math:`e / \mathrm{a}_0^3`\ )

    Parameters
    ----------
    value : Any
        Any value convertible to float representing the value in atomic units.
    """

    __slots__ = ()

    def __new__(cls, value: Any):
        return super().__new__(cls, float(value))  # type: ignore # (Too many arguments for "__new__" of "object")

    def __repr__(self) -> str:
        return f"Ed({super().__repr__()})"

    def __str__(self) -> str:
        return f"{super().__str__()} a.u."


class AbstractMesh(ABC):
    """Abstract base class for collections of points in space

    Calling ``len`` on instances of this class will return the number of points.
    """

    @property
    @abstractmethod
    def points(self) -> Iterator[Coords]:
        """Coordinates of points of which the mesh consists

        Yields
        ------
        Iterator[Coords]
            Iterator over the point coordinates
        """
        pass

    @abstractmethod
    def __len__(self) -> int:
        pass


@dataclass
class Mesh(AbstractMesh):
    """Collection of points in space without assumptions regarding structure

    This class stores all the points given on initialization and hence it's
    memory footprint is linear in the number of points.

    Parameters
    ----------
    points_ : Collection[Coords]
        The coordinates of points to be stored
    """
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
class GridMesh(AbstractMesh):
    """Collection of points in space organized in a grid

    This class only stores information regarding the grid which underlies the
    spatial organization of the points and thus it's memory footprint is constant
    with respect to the number of points it describes. This is at the cost of
    a small CPU cost whenever a point is retrieved from the ``points`` iterator.

    Parameters
    ----------
    origin : Coords
        The coordinates of the coordinates system origin.
    axes : Axes
        The coordinates of the coordinates axes' origin.

    Attributes
    ----------
    Axes : Tuple[GridMesh.Axis, GridMesh.Axis, GridMesh.Axis]
        Type alias for a tuple of three ``Axis`` objects. Note that currently
        only axes aligned with the coordinate system axes are supported, i.e.
        the axes' vectors are expected to be ((1, 0, 0), (0, 1, 0), (0, 0, 1)).
    """

    @dataclass
    class Axis:
        """Dataclass describing a coordinate system axis

        Parameters
        ----------
        vector : Coords
            Unit vector in xyz coordinates, e.g. (1, 0, 0) for a typical x-axis.
        point_count : int
            The number of points that the grid places along this axis.

        Attributes
        ----------
        vector
            See initialization parameter
        point_count
            See initialization parameter
        """
        vector: Coords
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


FieldValue = TypeVar('FieldValue')
"""typing.TypeVar : The generic type for the values of the `Field` class.

There are no restrictions on what this type can be.
"""

@dataclass
class Field(Generic[FieldValue]):
    """Dataclass representing values of a field at a "mesh" of points in space

    This class is generic in the type of the field value, which can be of any
    type. Classes where `FieldValue` matches `NumericValue`, additionally
    support arithmetic operations (currently only addition and subtraction).

    Parameters
    ----------
    mesh : AbstractMesh
        A "mesh" of points in space at which the field has values
    values\_ : Collection[FieldValue]
        A collection of values corresponding to the points in space given in
        the same order as the `AbstractMesh.points` iterator.

    Attributes
    ----------
    mesh
        See initialization parameter
    values : typing.List[FieldValue]
        Converted from the `values_` initialization parameter
    NumericValue : typing.TypeVar
        Generic type specifying a subset of FieldValue types for which arithmetic
        operations are defined. This can be any type matching "bound=float".
    """

    mesh: AbstractMesh
    values_: InitVar[Collection[FieldValue]]
    values: List[FieldValue] = field(init=False)

    def __post_init__(self, values_) -> None:

        if len(values_) != len(self.mesh):
            raise ValueError(
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
