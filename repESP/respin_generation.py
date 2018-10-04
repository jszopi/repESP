from .charges import Charge
from .respin_util import Respin, Equivalence
from .types import Atom, Molecule
from .util import _zip_exact

from abc import ABC, abstractmethod
from typing import List, Optional


class TypeSpecificRespinGenerator(ABC):

    @abstractmethod
    def resp_type(self) -> str:
        pass

    @abstractmethod
    def get_ivary(self) -> Respin.Ivary:
        pass

    @abstractmethod
    def get_cntrl(self) -> Respin.Cntrl:
        pass


class RespRespinGenerator(TypeSpecificRespinGenerator):

    def __init__(self, ivary: Respin.Ivary) -> None:
        self.ivary = ivary

    def get_ivary(self) -> Respin.Ivary:
        return self.ivary

    @classmethod
    def from_respin(cls, respin: Respin):
        return cls(respin.ivary)

    @classmethod
    @abstractmethod
    def from_methyl_and_methylene(
        cls,
        equivalence: Equivalence,
        methyl_methylene_mask: List[bool]
    ):
        pass


class RespStage1RespinGenerator(RespRespinGenerator):

    def resp_type(self) -> str:
        return "RESP stage 1"

    def get_cntrl(self) -> Respin.Cntrl:
        return Respin.Cntrl(
            qwt=0.0005
        )

    @classmethod
    def from_methyl_and_methylene(
        cls,
        equivalence: Equivalence,
        methyl_methylene_mask: List[bool]
    ):
        return cls(Respin.Ivary([
            0 if is_methyl_or_methylene else ivary_from_equivalence_value
            for ivary_from_equivalence_value, is_methyl_or_methylene
            in _zip_exact(
                Respin.Ivary.from_equivalence(equivalence).values,
                methyl_methylene_mask
            )
        ]))


class RespStage2RespinGenerator(RespRespinGenerator):

    def resp_type(self) -> str:
        return "RESP stage 2"

    def get_cntrl(self) -> Respin.Cntrl:
        return Respin.Cntrl(
            iqopt=2,
            qwt=0.001
        )

    @classmethod
    def from_methyl_and_methylene(
        cls,
        equivalence: Equivalence,
        methyl_methylene_mask: List[bool]
    ):
        return cls(Respin.Ivary([
            ivary_from_equivalence_value if is_methyl_or_methylene else -1
            for ivary_from_equivalence_value, is_methyl_or_methylene
            in _zip_exact(
                Respin.Ivary.from_equivalence(equivalence).values,
                methyl_methylene_mask
            )
        ]))


class FitHydrogensOnlyRespinGenerator(TypeSpecificRespinGenerator):

    def __init__(self, equivalence: Equivalence, molecule: Molecule[Atom]) -> None:
        self.equivalence = equivalence
        self.molecule = molecule

    def resp_type(self) -> str:
        return "fitting of hydrogen atoms"

    def get_ivary(self) -> Respin.Ivary:
        return Respin.Ivary([
            -1 if atom.identity != 1 else ivary
            for atom, ivary in zip(
                self.molecule.atoms,
                Respin.Ivary.from_equivalence(self.equivalence).values
            )
        ])

    def get_cntrl(self) -> Respin.Cntrl:
        return Respin.Cntrl(
            iqopt=2,
            ihfree=0,  # Not just hydrogens...
            qwt=0.0  # ...because there are no restraints.
        )


class EquivalenceOnlyRespinGenerator(TypeSpecificRespinGenerator):

    def __init__(self, equivalence: Equivalence) -> None:
        self.equivalence = equivalence

    def resp_type(self) -> str:
        return "atom equivalencing"

    def get_ivary(self) -> Respin.Ivary:
        return Respin.Ivary.from_equivalence(self.equivalence)

    def get_cntrl(self) -> Respin.Cntrl:
        return Respin.Cntrl(
            ihfree=0,  # Not just hydrogens...
            qwt=0.0  # ...because there are no restraints.
        )


class FrozenAtomsRespinGenerator(TypeSpecificRespinGenerator):

    def __init__(self, equivalence: Equivalence, frozen_atoms: List[int]) -> None:

        for atom_label in frozen_atoms:
            if atom_label < 0 or atom_label >= len(equivalence.values):
                raise ValueError(
                    f"Label of atom to be frozen ({atom_label}) is outside of "
                    f"expected range i.e. [0, {len(equivalence.values)})."
                )

        self.equivalence = equivalence
        self.frozen_atoms = frozen_atoms

    def resp_type(self) -> str:
        return "fitting with selected atom charges frozen"

    def get_ivary(self) -> Respin.Ivary:
        return Respin.Ivary([
            -1 if i in self.frozen_atoms else ivary
            for i, ivary in enumerate(
                Respin.Ivary.from_equivalence(self.equivalence).values
            )
        ])

    def get_cntrl(self) -> Respin.Cntrl:
        return Respin.Cntrl(
            iqopt=2,
            ihfree=0,  # Not just hydrogens...
            qwt=0.0  # ...because there are no restraints.
        )


def prepare_respin(
    respin_generator: TypeSpecificRespinGenerator,
    total_charge: int,
    molecule: Molecule[Atom],
    title: Optional[str]=None,
    subtitle: Optional[str]=None,
    read_charges: bool=False
) -> Respin:

    cntrl = respin_generator.get_cntrl()
    if read_charges:
        cntrl.iqopt = 2

    default_title = f"Respin file prepared by `repESP` to perform {respin_generator.resp_type()}."
    default_subtitle = "Resp charges for organic molecule"

    return Respin(
        title=title if title is not None else default_title,
        cntrl=cntrl,
        wtmol=1.0,
        subtitle=subtitle if subtitle is not None else default_subtitle,
        charge=total_charge,
        iuniq=len(molecule.atoms),
        molecule=molecule,
        ivary=respin_generator.get_ivary()
    )
