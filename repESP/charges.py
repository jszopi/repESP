from .exceptions import InputFormatError
from .types import Molecule
from .util import _NoValue

from dataclasses import dataclass
from enum import auto
from typing import Any, List, Collection, NewType


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
