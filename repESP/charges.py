from .exceptions import InputFormatError
from .types import Atom, AtomWithCoords, Coords
from .util import _NoValue, _get_atomic_number

from dataclasses import dataclass, make_dataclass
from enum import auto
from typing import Any, List, Collection, NewType, Tuple, Type


Charge = NewType("Charge", float)  # Atomic charge [elementary charge]


def make_charge(x: Any) -> Charge:
    return Charge(float(x))


class ChargeType(_NoValue):
    MULLIKEN = auto()
    MK = auto()
    CHELP = auto()
    CHELPG = auto()
    HLY = auto()
    NPA = auto()
    AIM = auto()


@dataclass
class AtomWithCharge(Atom):
    charge: Charge

    @classmethod
    def from_symbol(  # type: ignore
        cls: 'AtomWithCharge',
        symbol: str,
        charge: Charge
    ) -> 'AtomWithCharge':
        return cls(_get_atomic_number(symbol), charge)


@dataclass
class AtomWithCoordsAndCharge(AtomWithCharge, AtomWithCoords):

    @classmethod
    def from_symbol(  # type: ignore
        cls: 'AtomWithCoordsAndCharge',
        symbol: str,
        coords: Coords,
        charge: Charge
    ) -> 'AtomWithCoordsAndCharge':
        # TODO: mypy incorrectly infers the argument order for __init__ to be:
        # (identity, charge, coords). Not ignoring here as it's also failing
        # the tests, so leaving this one as reference.
        return cls(_get_atomic_number(symbol), coords, charge)


DipoleMoment = NewType("DipoleMoment", float)  # Dipole moment [bohr * fundamental charge]


def make_dipole_moment(x: Any) -> DipoleMoment:
    return DipoleMoment(float(x))


@dataclass
class Dipole:
    x: DipoleMoment
    y: DipoleMoment
    z: DipoleMoment


QuadrupoleMoment = NewType("QuadrupoleMoment", float)  # Quadrupole moment [bohr^2 * fundamental charge]


def make_quadrupole_moment(x: Any) -> QuadrupoleMoment:
    return QuadrupoleMoment(float(x))


@dataclass
class Quadrupole:
    xx: QuadrupoleMoment
    yy: QuadrupoleMoment
    zz: QuadrupoleMoment
    xy: QuadrupoleMoment
    xz: QuadrupoleMoment
    yz: QuadrupoleMoment
