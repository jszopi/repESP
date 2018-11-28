"""Creating ``resp`` program input for predefined ESP fitting types

This module may be useful if you require custom instructions for the ``resp``
program. Otherwise, the functions in the `resp_wrapper` module should be
sufficient and more convenient.
"""

from repESP.charges import Charge
from repESP.equivalence import Equivalence
from repESP.respin_format import Respin
from repESP.types import Atom, Molecule
from repESP._util import zip_exact

from abc import ABC, abstractmethod
from typing import List, Optional


class RespinGenerator(ABC):
    """Interface for creators of ``resp`` fitting instructions

    This interface will be implemented differently for different types of
    fitting with the ``resp`` program. This module provides implementations
    for the most common fitting types.

    An implementation of this interface is required by the `prepare_respin`
    function.
    """

    @abstractmethod
    def resp_type(self) -> str:
        """Return the human-readable description of the type of fitting"""
        pass

    @abstractmethod
    def get_ivary(self) -> Respin.Ivary:
        """Create the "ivary" section of the "respin" file"""
        pass

    @abstractmethod
    def get_cntrl(self) -> Respin.Cntrl:
        """Create the "cntrl" section of the "respin" file"""
        pass


class RespRespinGenerator(RespinGenerator):
    """Base class providing "respin" creation methods for RESP charges"""

    def __init__(self, ivary: Respin.Ivary) -> None:
        self.ivary = ivary

    def get_ivary(self) -> Respin.Ivary:
        return self.ivary

    @classmethod
    def from_respin(cls, respin: Respin):
        """Alternative initialization from a "respin" file

        This alternative initializer expresses the fact that this object can
        be created directly from an existing "respin" file, when the contained
        "ivary" section expresses the desired per-atom fitting instructions.

        Parameters
        ----------
        respin : Respin
            An object representing pre-existing instructions to the ``resp`` program.
        """
        return cls(respin.ivary)

    @classmethod
    @abstractmethod
    def from_methyl_and_methylene(
        cls,
        equivalence: Equivalence,
        methyl_methylene_mask: List[bool]
    ):
        """Alternative initialization from chemical information

        The per-atom fitting instructions in the RESP charge optimization depend
        on the chemical equivalence and whether a given carbon or hydrogen atom
        belongs to a methyl or methylene group.

        Parameters
        ----------
        equivalence : Equivalence
            The chemical equivalence relations between atoms in the molecule.
        methyl_methylene_mask : typing.List[bool]
            A list, in which each value describes an atom of the molecule. The
            value should be True when the atom is a carbon or hydrogen atom in
            a methyl or methylene group and False otherwise.

        Raises
        ------
        ValueError
            Raised when the length of the `values` attribute of the
            `equivalence` argument is not equal to the length of the
            `methyl_methylene_mask` argument.
        """
        pass


class RespStage1RespinGenerator(RespRespinGenerator):
    """Creator of "respin" files for 1st stage RESP fitting

    The original RESP procedure is conducted in two stages and this object
    prepares the input for the first stage. In this stage, carbon and hydrogen
    atoms in methyl and methylene groups vary freely, while other atoms need to
    obey chemical equivalence relations. The magnitudes of charges on all atoms
    are constrained with hyperbolic constraints with a weight factor of 0.0005.

    Parameters
    ----------
    ivary : Respin.Ivary
        The "ivary" section of a "respin" file for 1st stage RESP fitting.
    """

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
        """See documentation in base class (`RespRespinGenerator`)"""
        return cls(Respin.Ivary([
            0 if is_methyl_or_methylene else ivary_from_equivalence_value
            for ivary_from_equivalence_value, is_methyl_or_methylene
            in zip_exact(
                Respin.Ivary.from_equivalence(equivalence).values,
                methyl_methylene_mask
            )
        ]))


class RespStage2RespinGenerator(RespRespinGenerator):
    """Creator of "respin" files for 2nd stage RESP fitting

    The original RESP procedure is conducted in two stages and this object
    prepares the input for the second stage. In this stage, carbon and hydrogen
    atoms in methyl and methylene groups are varied subject to equivalence
    relations and the magnitudes of the charges are constrained with hyperbolic
    constraints with a weight factor of 0.001. Other atoms are frozen at
    initial charge values, which should come from 1st stage fitting.

    Parameters
    ----------
    ivary : Respin.Ivary
        The "ivary" section of a "respin" file for 2nd stage RESP fitting.
    """

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
        """See documentation in base class (`RespRespinGenerator`)"""
        return cls(Respin.Ivary([
            ivary_from_equivalence_value if is_methyl_or_methylene else -1
            for ivary_from_equivalence_value, is_methyl_or_methylene
            in zip_exact(
                Respin.Ivary.from_equivalence(equivalence).values,
                methyl_methylene_mask
            )
        ]))


class FitHydrogensOnlyRespinGenerator(RespinGenerator):
    """Supports creation of "respin" files for fitting only hydrogen atoms

    The generated "respin" file instructs ``resp`` to perform fitting in which
    only hydrogen atoms are free to vary, subject to equivalence relations.

    Parameters
    ----------
    equivalence : Equivalence
        The equivalence relations between atoms of the molecule.
    molecule : Molecule[Atom]
        The molecule which charges are being fitted. Only atom identities are required.
    """

    def __init__(self, equivalence: Equivalence, molecule: Molecule[Atom]) -> None:
        self.equivalence = equivalence
        self.molecule = molecule

    def resp_type(self) -> str:
        return "fitting of hydrogen atoms"

    def get_ivary(self) -> Respin.Ivary:
        return Respin.Ivary([
            -1 if atom.atomic_number != 1 else ivary
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


class EquivalenceOnlyRespinGenerator(RespinGenerator):
    """Supports creation of "respin" files with no constraints beyond equivalence

    The generated "respin" file instructs ``resp`` to perform fitting in which
    atoms are free to vary, subject only to equivalence relations.

    Parameters
    ----------
    equivalence : Equivalence
        The equivalence relations between atoms of the molecule.
    """

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


class FrozenAtomsRespinGenerator(RespinGenerator):
    """Supports creation of "respin" files excluding selected atoms from fitting

    The generated "respin" file instructs ``resp`` to perform fitting in which
    atoms are free to vary (subject to equivalence relations) unless they are
    explicitly excluded from the fitting by being listed in the `frozen_atoms`
    argument.

    Parameters
    ----------
    equivalence : Equivalence
        The equivalence relations between atoms of the molecule.
    frozen_atoms : typing.List[int]
        List of zero-based indices of atoms in the molecule which are to be
        excluded from the fitting.
    """

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
    respin_generator: RespinGenerator,
    total_charge: int,
    molecule: Molecule[Atom],
    title: Optional[str]=None,
    subtitle: Optional[str]=None,
    read_charges: bool=False
) -> Respin:
    """Create ``resp`` program instructions based on input options

    This function allows to create ``resp`` program instructions from higher-level
    abstractions.

    Parameters
    ----------
    respin_generator : RespinGenerator
        An object implementing the `RespinGenerator` interface.
        This object expresses a given type of charge fitting through
        instructions understood by ``resp``. Implementations of this interface
        for common fitting types are provided in this module.
    total_charge : int
        The total charge of the molecule.
    molecule : Molecule[Atom]
        The molecule which charges are being fitted. Only atom identities are required.
    title : Optional[str], optional
        The title for the optimization. If set to None (default), a default
        title will be used, referencing the `repESP` library and the type of
        the fitting as given by the `respin_generator` used.
    subtitle : Optional[str], optional
        The subttitle for the fitted structure. If set to None (default), a
        default generic subtitle will be used.
    read_charges : bool, optional
        If this option is set to True, ``resp`` will require initial charges
        for the fitting. Defaults to False.

        .. note::
            It is the responsibility of the `respin_generator` to instruct
            ``resp`` to read initial charges when they affect the result of the
            optimization. This parameter should only be used when the user
            wishes to aid the optimization by supplying an initial guess.

        .. It could also be used in case the ``resp`` algorithm converges to
            different charges depending on the initial charges but I am not
            aware of such an issue with the ``resp`` algorithm).

    Returns
    -------
    Respin
        An object representing the fitting instructions for the ``resp`` program.
    """

    cntrl = respin_generator.get_cntrl()
    if read_charges:
        cntrl.iqopt = 2

    default_title = f"Respin file prepared by `repESP` to perform {respin_generator.resp_type()}."
    default_subtitle = "Resp charges for organic molecule"

    return Respin(
        title=title if title is not None else default_title,
        cntrl=cntrl,
        subtitle=subtitle if subtitle is not None else default_subtitle,
        charge=total_charge,
        molecule=molecule,
        ivary=respin_generator.get_ivary()
    )
