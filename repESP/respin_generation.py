from .charges import Charge
from .respin_util import Respin, Equivalence

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


class RespStage1RespinGenerator(TypeSpecificRespinGenerator):

    def __init__(self, respin1: Respin) -> None:
        self.respin1 = respin1

    def resp_type(self) -> str:
        return "RESP stage 1"

    def get_ivary(self) -> Respin.Ivary:
        return self.respin1.ivary

    def get_cntrl(self) -> Respin.Cntrl:
        return Respin.Cntrl(
            qwt=0.0005
        )


class RespStage2RespinGenerator(TypeSpecificRespinGenerator):

    def __init__(self, respin1: Respin) -> None:
        self.respin1 = respin1

    def resp_type(self) -> str:
        return "RESP stage 2"

    def get_ivary(self) -> Respin.Ivary:
        return self.respin1.ivary

    def get_cntrl(self) -> Respin.Cntrl:
        return Respin.Cntrl(
            iqopt=2,
            qwt=0.001
        )

class FitHydrogensOnlyRespinGenerator(TypeSpecificRespinGenerator):

    def __init__(self, equivalence: Equivalence, atomic_numbers: List[int]) -> None:
        self.equivalence = equivalence
        self.atomic_numbers = atomic_numbers

    def resp_type(self) -> str:
        return "fitting of hydrogen atoms"

    def get_ivary(self) -> Respin.Ivary:
        return Respin.Ivary([
            -1 if atomic_number != 1 else ivary
            for atomic_number, ivary in zip(
                self.atomic_numbers,
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
        # This should check compatibility between equivalence and frozen_atoms.
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
    atomic_numbers: List[int],
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
        iuniq=len(atomic_numbers),
        atomic_numbers=atomic_numbers,
        ivary=respin_generator.get_ivary()
    )
