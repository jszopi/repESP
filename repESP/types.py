from .util import _elements, _get_symbol, _get_atomic_number

from dataclasses import dataclass
from typing import Any, Generic, List, NewType, Tuple, TypeVar

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


@dataclass
class Atom:
    atomic_number: int

    def __post_init__(self):
        if self.atomic_number < 1 or self.atomic_number >= len(_elements):
            raise ValueError("Atomic number is not within expected bounds.")

    @property
    def symbol(self) -> str:
        return _get_symbol(self.atomic_number)

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
