from .exceptions import InputFormatError
from .types import Molecule

from dataclasses import dataclass
from enum import Enum, auto
from typing import Any, List, Collection, NewType


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
