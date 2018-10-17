"""Types used to describe partial charges and higher moments"""

from .exceptions import InputFormatError
from .types import Atom, AtomWithCoords, Coords
# from ._util import NoValue

from dataclasses import dataclass
# from enum import auto
from typing import Any, List, Collection, Tuple, Type


class Charge(float):
    """Partial (atomic) charge in units of elementary charge (*e*)

    Parameters
    ----------
    value : Any
        Any value convertible to float representing the value in units of *e*.
    """

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
    """Dataclass representing an atom with an associated partial charge

    Parameters
    ----------
    atomic_number : int
        The atomic number of the element
    charge: Charge
        The partial charge on the atom
    """
    charge: Charge


@dataclass
class AtomWithCoordsAndCharge(AtomWithCharge, AtomWithCoords):
    """Dataclass representing an atom in space, with a partial charge

    Parameters
    ----------
    atomic_number : int
        The atomic number of the element
    charge : Charge
        The partial charge on the atom
    coords : Coords
        The coordinates of the atom in space
    """
    # NOTE: mypy incorrectly infers the argument order for __init__ to be:
    # (atomic_number, charge, coords), resulting in loads of spurious errors.
    pass


class DipoleMomentValue(float):
    """Dipole moment value in atomic units (:math:`e\mathrm{a}_0`)

    Parameters
    ----------
    value : Any
        Any value convertible to float representing the value in atomic units.
    """

    __slots__ = ()

    def __new__(cls, value: Any):
        return super().__new__(cls, float(value))  # type: ignore # (Too many arguments for "__new__" of "object")

    # TODO: Add conversion to/from Debye (and printing with it).


@dataclass
class DipoleMoment:
    """Dipole moment vector described through its three components

    Parameters
    ----------
    x : DipoleMomentValue
        The value in the *x* direction
    y : DipoleMomentValue
        The value in the *y* direction
    z : DipoleMomentValue
        The value in the *z* direction
    """

    x: DipoleMomentValue
    y: DipoleMomentValue
    z: DipoleMomentValue


class QuadrupoleMomentValue(float):
    """Quadrupole moment in atomic units (:math:`e\mathrm{a}_0^2`)

    Parameters
    ----------
    value : Any
        Any value convertible to float representing the value in atomic units.
    """

    __slots__ = ()

    def __new__(cls, value: Any):
        return super().__new__(cls, float(value))  # type: ignore # (Too many arguments for "__new__" of "object")

    # TODO: Add conversion to/from Buckingham (and printing with it).

@dataclass
class QuadrupoleMoment:
    """Quadrupole moment tensor described through its six components

    Parameters
    ----------
    xx : QuadrupoleMomentValue
        The value in the *xx* direction.
    yy : QuadrupoleMomentValue
        The value in the *yy* direction.
    zz : QuadrupoleMomentValue
        The value in the *zz* direction.
    xy : QuadrupoleMomentValue
        The value in the *xy* direction.
    xz : QuadrupoleMomentValue
        The value in the *xz* direction.
    yz : QuadrupoleMomentValue
        The value in the *yz* direction.
    """

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
