"""Types used to describe partial charges and higher moments"""

from .exceptions import InputFormatError
from .types import Atom, AtomWithCoords, Coords
# from ._util import NoValue

from dataclasses import dataclass
# from enum import auto
from typing import Any, List, Collection, Tuple, Type


class Charge(float):

    """Atomic charge [elementary charge]"""

    __slots__ = ()

    def __new__(cls, value: Any):
        return super().__new__(cls, float(value))  # type: ignore # (Too many arguments for "__new__" of "object")

    def __repr__(self) -> str:
        return f"Charge({super().__repr__()})"

    def __str__(self) -> str:
        return f"{super().__str__()} e"


# Could be useful as a common denominator between modules, e.g. to translate
# a given charge type to its corresponding ChargesSectionParser. However, this
# is currently not necessary.
#
# class ChargeType(NoValue):
#     MULLIKEN = auto()
#     MK = auto()
#     CHELP = auto()
#     CHELPG = auto()
#     HLY = auto()
#     NPA = auto()
#     AIM = auto()


@dataclass
class AtomWithCharge(Atom):
    charge: Charge


@dataclass
class AtomWithCoordsAndCharge(AtomWithCharge, AtomWithCoords):
    # NOTE: mypy incorrectly infers the argument order for __init__ to be:
    # (atomic_number, charge, coords), resulting in loads of spurious errors.
    pass


class DipoleMomentValue(float):

    """Dipole moment [bohr * fundamental charge]"""

    __slots__ = ()

    def __new__(cls, value: Any):
        return super().__new__(cls, float(value))  # type: ignore # (Too many arguments for "__new__" of "object")

    # TODO: Add conversion to/from Debye (and printing with it).


@dataclass
class DipoleMoment:
    x: DipoleMomentValue
    y: DipoleMomentValue
    z: DipoleMomentValue


class QuadrupoleMomentValue(float):

    """Quadrupole moment [bohr^2 * fundamental charge]"""

    __slots__ = ()

    def __new__(cls, value: Any):
        return super().__new__(cls, float(value))  # type: ignore # (Too many arguments for "__new__" of "object")

    # TODO: Add conversion to/from Buckingham (and printing with it).

@dataclass
class QuadrupoleMoment:
    xx: QuadrupoleMomentValue
    yy: QuadrupoleMomentValue
    zz: QuadrupoleMomentValue
    xy: QuadrupoleMomentValue
    xz: QuadrupoleMomentValue
    yz: QuadrupoleMomentValue

    def __post__init__(self):
        # TODO: Only five of those are independent as the tensor is required to
        # be traceless. Add this check.
        pass
