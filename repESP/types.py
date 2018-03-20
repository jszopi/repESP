from abc import ABC, abstractmethod
from typing import Iterable, Iterator, List, NamedTuple, Tuple

import math


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


class Points(ABC):

    @abstractmethod
    def points(self) -> Iterator[Coords]:
        pass

    @abstractmethod
    def __eq__(self, other):
        pass

    def __ne__(self, other):
        return not self == other


class Mesh(Points):

    def __init__(self, points_coords: Iterable[Coords]) -> None:
        self._points = list(points_coords)

    def __eq__(self, other) -> bool:
        return self._points == other._points

    def points(self) -> Iterator[Coords]:
        return iter(self._points)


# class Grid(Points):
#     pass
#
