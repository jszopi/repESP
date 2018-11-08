"""Parsing and writing ``resp`` program input format ("respin" files)"""

from dataclasses import dataclass, asdict
from fortranformat import FortranRecordWriter as FW
from itertools import zip_longest
import io
import math
import re
import sys
from typing import Dict, List, Optional, TextIO, Tuple, TypeVar, Union

from .exceptions import InputFormatError
from .equivalence import Equivalence
from .types import Atom, Molecule
from ._util import zip_exact


@dataclass
class Respin:
    """Dataclass describing the ``resp`` program input

    Note that the functionality is currently limited to a single molecule and
    a single structure.

    Parameters
    ----------
    title : str
        The title of the calculation to be performed.
    cntrl : Cntrl
        Dataclass representing the "cntrl" section of the input.
    subtitle : str
        Subtitle describing the considered molecular structure.
    charge : int
        The total charge of the molecule.
    molecule : Molecule[Atom]
        The molecule as described by the positions of constituent atoms.
    ivary : Ivary
        The "ivary" values for fitting the considered structure. These determine
        how the charge on each atom is allowed to vary during the fitting.

    Attributes
    ----------
    title
        See initialization parameter
    cntrl
        See initialization parameter
    subtitle
        See initialization parameter
    charge
        See initialization parameter
    molecule
        See initialization parameter
    ivary
        See initialization parameter
    """

    _ValueType = TypeVar("_ValueType", int, float, str)

    @staticmethod
    def _check_value(
        attr_name: str,
        attr_value: _ValueType,
        allowed_values: List[_ValueType]
    ) -> None:
        if attr_value not in allowed_values:
            raise ValueError(
                f"Invalid value for `{attr_name}`: {attr_value}."
            )

    @dataclass
    class Cntrl:
        """Dataclass describing the "cntrl" section of a "respin" file

        See ``resp`` program documentation for more details.

        Parameters
        ----------
        inopt : int, optional
            If equal to 1, ``resp`` will cycle through different "qwt" values
            from the file specified with the ``-w`` option. Defaults to 0.
        ioutopt : int, optional
            If equal to 1, ``resp`` will write restart info of new ESP field
            to the file specified with the ``-s`` option. Defaults to 0.
        iqopt : int, optional
            Controls setting initial charges. If equal to 1 (default), all
            initial charges will be set to zero. If equal to 2, initial charges
            are read from the file specified with the ``-q`` option. If equal
            to 3, charges are read as with the previous option and will
            additionally be averaged according to "ivary" values (normally not
            used).
        ihfree : int, optional
            If equal to 0, the charge magnitude restraint is applied to all
            charges. If equal to 1 (default), the restraint does not apply to
            hydrogen atoms.
        irstrnt : int, optional
            Controls the type of charge magnitude restraint. If equal to 0,
            harmonic restraints are used (old-style). If equal to 1 (default),
            hyperbolic restraints are used. If equal to 2, no charge fitting
            is carried out and only analysis of input charges is performed.
        qwt : float, optional
            The weight of the charge magnitude restraint to be used during
            the fitting. Defaults to 0.0 (no charge magnitude restraint).

            .. warning::
                The default used here is different from the default used by ``resp``.

                That the ``resp`` documentation specifies that it uses the
                Amber force field values by default. However, it is not clear how
                it can determine the fitting stage. Thus, to remove the ambiguity,
                this dataclass assumes a weight of zero by default.

            .. note::
                Amber force fields use values of 0.0005 and 0.001 for
                stages 1 and 2, respectively. The Glycam force field is derived with
                one stage fitting with a value of 0.01.


        Attributes
        ----------
        inopt
            See initialization parameter
        ioutopt
            See initialization parameter
        iqopt
            See initialization parameter
        ihfree
            See initialization parameter
        irstrnt
            See initialization parameter
        qwt
            See initialization parameter
        """
        inopt: int = 0
        ioutopt: int = 0
        iqopt: int = 1
        ihfree: int = 1
        irstrnt: int = 1
        qwt: float = 0

        @property
        def nmol(self) -> int:
            """Number of structures in a multiple structure fit.

            With the current implementation this will always be equal to 1.
            """
            return 1

        def __post_init__(self) -> None:
            Respin._check_value("inopt", self.inopt, [0, 1])
            Respin._check_value("ioutopt", self.ioutopt, [0, 1])
            Respin._check_value("iqopt", self.iqopt, [1, 2, 3])
            Respin._check_value("ihfree", self.ihfree, [0, 1])
            Respin._check_value("irstrnt", self.irstrnt, [0, 1, 2])
            if self.qwt < 0:
                raise ValueError(f"Invalid value for `qwt`: {self.qwt}.")

    @dataclass
    class Ivary:
        """Dataclass representing per atom fitting instructions for ``resp``

        The fitting instructions are represented as a list of values stored in
        the `values` attribute. The length of this list must be the same as the
        number of atoms in the molecule it describes. Consecutive values refer
        to consecutive atoms of the molecule.

        The values determine how the charge on the atom can vary during the
        fitting and the allowed values are:

        * -1, meaning that the atom's charge is "frozen" at the initial value
        * 0, meaning that this atom will be varied freely
        * Larger than zero, representing the 1-based index of the atom in the
          molecule to which this atom is to be equivalenced.

        Example
        -------

        Consider fitting RESP charges in a molecule of methylamine:

        >>> methylamine = Molecule([Atom(atomic_number) for atomic_number in [6, 1, 1, 1, 7, 1, 1]])

        Fitting RESP charges consists of two stages. The ivary section for the
        second stage of the fitting for the methylamine molecule should be as
        follows:

        >>> ivary = Respin.Ivary([0, 0, 2, 2, -1, -1, -1])

        The carbon atom is free to vary during the fitting. The first of the methyl
        hydrogens is equivalenced to the remaining two but they haven't been
        specified yet, so it also has a value of 0. These two hydrogen atoms
        are equivalenced to the first one, and thus are assigned its one-based
        index in the molecule, i.e. 2 (meaning "equivalenced to the second atom
        of the molecule"). The nitrogen atom and the two hydrogen atoms attached
        to it are frozen during the second stage of the fitting and are thus
        assigned values of -1.

        Parameters
        ----------
        values : List[int]
            The per-atom instructions for the ``resp`` program.

        Attributes
        ----------
        values
            See initialization parameter
        """

        values: List[int]

        def __post_init__(self):
            for i, elem in enumerate(self.values):
                if elem < -1 or elem > len(self.values):
                    raise ValueError(
                        f"Value number {i} passed as `ivary` with value {elem}, "
                        f"which is either lower than 0 or outside the list length."
                    )

        def describe(self, molecule: Optional[Molecule[Atom]]=None) -> str:
            """Verbosely report the "ivary" actions

            Example
            -------

            >>> print(ivary.describe(methylamine))
            Atom (C) number 1
            Atom (H) number 2
            Atom (H) number 3, equivalenced to atom 2
            Atom (H) number 4, equivalenced to atom 2
            Atom (N) number 5, frozen
            Atom (H) number 6, frozen
            Atom (H) number 7, frozen

            Parameters
            ----------
            molecule : Optional[Molecule[Atom]], optional
                The molecule to which the ivary information refers. This
                argument is optional and defaults to None. If it is provided,
                atom identities will be included in the output.

            Raises
            ------
            ValueError
                Raised when the number of atoms in the molecule does not match
                the length of the list of values in this object.

            Returns
            -------
            str
                A verbose description of the "ivary" instructions.
            """

            if molecule is not None and len(molecule.atoms) != len(self.values):
                raise ValueError(
                    f"The number of atoms ({len(molecule.atoms)} is not the same "
                    f"as the number of ivary values ({len(self.values)}."
                )

            zipped = zip_longest(self.values, molecule.atoms if molecule is not None else [])

            f = io.StringIO()
            for i, (ivary, atom) in enumerate(zipped):
                atomic_number = atom.symbol if molecule is not None else None
                id_str = f" ({atomic_number})" if atomic_number is not None else ""

                if ivary < 0:
                    # TODO: This could also report at what value if charges are provided
                    ivary_str = ", frozen"
                elif ivary > 0:
                    ivary_str = f", equivalenced to atom {ivary}"
                else:
                    ivary_str = ""

                print(f"Atom{id_str} number {i+1}{ivary_str}", file=f)

            return f.getvalue()

        @classmethod
        def from_equivalence(cls, equivalence: Equivalence):
            """Alternative initialization from equivalence information

            .. note:
                The resulting ivary instructions will correspond to fitting
                with equivalent atoms assigned identical charges. This may not
                be the type of fitting that you want to perform.

            Parameters
            ----------
            equivalence : Equivalence
                Information about chemical equivalence of atoms in a molecule.
            """
            return cls([0 if val is None else val + 1 for val in equivalence.values])


    title: str
    cntrl: Cntrl
    subtitle: str
    charge: int
    molecule: Molecule[Atom]
    ivary: Ivary

    @property
    def wtmol(self) -> float:
        """Relative weight of the structure in a multistructure fitting.

        A value of 1.0 is always returned in the current implementation.
        """
        return 1.0

    @property
    def iuniq(self) -> int:
        """The number of atoms in the fitted structure"""
        return len(self.molecule.atoms)

    def __post_init__(self) -> None:
        if len(self.molecule.atoms) != len(self.ivary.values):
            raise ValueError(
                f"Number of atoms ({len(self.molecule.atoms)}) does not match number "
                f"of ivary values ({len(self.ivary.values)})."
            )


def _get_equivalence_from_ivary(ivary: Respin.Ivary) -> Equivalence:
    """Get atom equivalence information from an `Respin.Ivary` object

    This function is private as users probably mean to use the
    `get_equivalence` function instead.

    `Ivary` objects are specific to ``resp`` program input and thus may not
    provide information about atom equivalence. The "respin" file may have been
    generated to perform any custom fitting with ``resp``. Only use this
    function when you're certain that the "respin" file contains all the
    equivalence information that you need.
    """
    return Equivalence([
        None if ivary_value == 0 else ivary_value - 1
        for ivary_value in ivary.values
    ])


def get_equivalence_from_two_stage_resp_ivary(ivary1: Respin.Ivary, ivary2: Respin.Ivary) -> Equivalence:
    """Get atom equivalence from two input for two 2-stage ``resp``

    Derive atom equivalence based on the data in two "respin" files
    (represented by the `Respin` objects) created for the purpose of two-stage
    fitting with the ``resp`` program. The files can be generated with the
    ``respgen`` program with the following commands::

        respgen -i methane.ac -o methane.respin1 -f resp1
        respgen -i methane.ac -o methane.respin2 -f resp2

    .. warning::
        The correctness of this function relies on:

        1. Antechamber and ``respgen`` correctly recognizing the symmetry
           relations between atoms. Fast-exchanging atoms may not be identified.
        2. The author's understanding of how ``respgen`` generates the "respin"
           files for two-stage RESP fitting.

        Thus it is advised to always check that the result of this function
        agrees with the domain knowledge about the studied molecule.
    """
    # The equivalence logic is explained somewhat inconsistently in the RESP
    # papers but I've additionally re-engineered the ``resp`` program's logic
    # to be sure that reading both the ``respin`` files will give the desired
    # behaviour. In fact, it's pretty simple. In the first stage atoms of the
    # methyl and methylene groups are free, while all the others are
    # equivalenced. In the second stage the former are equivalenced, while all
    # the others are frozen.

    return _get_equivalence_from_ivary(Respin.Ivary([
        max(ivary1_value, ivary2_value)
        for ivary1_value, ivary2_value in zip_exact(ivary1.values, ivary2.values)
    ]))


def _parse_cntrl(f: TextIO) -> Respin.Cntrl:
    line_re = re.compile(" (\w+) = ([0-9.]+)")
    kwargs: Dict[str, Union[int, float]] = {}
    for line in f:
        if line.rstrip('\n') == " &end":
            break
        if line.rstrip('\n') == "":
            continue

        line_match = line_re.match(line)
        if line_match is None:
            raise InputFormatError(
                f"Failed parsing cntrl section of respin file:\n{line}"
            )
        key = line_match.group(1)
        value = line_match.group(2)

        kwargs[key] = float(value) if key == "qwt" else int(value)

    # nmol is not a parameter of Cntrl.__init__ and must be equal to 1.
    nmol = kwargs.pop("nmol", None)
    if nmol is not None and nmol != 1:
        raise InputFormatError("Parsing multiple structures is not supported")

    return Respin.Cntrl(**kwargs)  # type: ignore # (not sure why not recognized)


def parse_respin(f) -> Respin:
    """Parse a file in the "respin" format (input format of ``resp``)

    Note that only files describing a single structure fit are currently supported.

    Parameters
    ----------
    f : TextIO
        File object opened in read mode containing the "respin" file.

    Raises
    ------
    InputFormatError
        Raised when the file does not follow the expected format.

    Returns
    -------
    Respin
        Object representing the fitting instructions for the ``resp`` program.
    """

    get_line = lambda: f.readline().rstrip('\n')
    title = get_line()
    while get_line() != " &cntrl":
        pass

    cntrl = _parse_cntrl(f)

    wtmol = get_line().strip()
    if not math.isclose(float(wtmol), 1.0, rel_tol=0, abs_tol=1e-6):
        raise InputFormatError(
            f"Encountered value of `wtmol` different from 1.0 ({wtmol}) but "
            f"parsing is supported only for single-structure respin files."
        )

    subtitle = get_line()

    charge_and_iuniq = get_line()
    if len(charge_and_iuniq.split()) != 2:
        raise InputFormatError(
            f"Expected two ints for the line specifying charge and iuniq, found:\n{charge_and_iuniq}"
        )

    charge = int(charge_and_iuniq.split()[0])
    iuniq = int(charge_and_iuniq.split()[1])

    atoms: List[Atom] = []
    ivary = Respin.Ivary([])

    for line in f:
        if line.rstrip('\n') == "":
            break
        if len(line.split()) != 2:
            raise InputFormatError(
                f"Expected two ints for the line specifying atom and ivary, found:\n{line}"
            )

        atoms.append(Atom(int(line.split()[0])))
        ivary_value = int(line.split()[1])
        # `respgen` uses a value of -99 but internally we use -1 as per resp spec.
        ivary.values.append(ivary_value if ivary_value != -99 else -1)

    if len(atoms) != iuniq:
        raise InputFormatError(
            f"The value of `iuniq` ({iuniq}) is different from the number of"
            f"atoms in the described molecule ({len(atoms)})."
        )

    return Respin(
        title,
        cntrl,
        subtitle,
        charge,
        Molecule(atoms),
        ivary
    )


def _write_cntrl(f: TextIO, cntrl: Respin.Cntrl, skip_defaults: bool) -> None:

    default_cntrl: Dict[str, Union[int, float]] = asdict(Respin.Cntrl())
    default_cntrl["nmol"] = 1

    dict_: Dict[str, Union[int, float]] = asdict(cntrl)
    dict_["nmol"] = cntrl.nmol

    print(" &cntrl\n", file=f)
    for key, value in dict_.items():
        if key == "qwt":
            print(" {} = {:.5f},".format(key, value), file=f)
        else:
            if not skip_defaults or value != default_cntrl[key]:
                print(" {} = {},".format(key, value), file=f)
    print("\n &end", file=f)


def write_respin(f: TextIO, respin: Respin, skip_cntrl_defaults: bool=True) -> None:
    """Write a "respin" file described by the given input data

    Parameters
    ----------
    f : TextIO
        The file object to which the instructions are to be saved. The file
        must be opened for writing.
    respin : Respin
        The dataclass representing all the instructions needed by the ``resp``
        program. Note that only single-structure fitting is currently supported.
    skip_cntrl_defaults : bool, optional
        When this option is set to True (default), fitting options with default
        values will not be written to the file.
    """

    print(respin.title, file=f)
    print(file=f)
    _write_cntrl(f, respin.cntrl, skip_cntrl_defaults)
    print(FW("F7.1").write([respin.wtmol]), file=f)
    print(respin.subtitle, file=f)
    print(FW("2I5").write([respin.charge, respin.iuniq]), file=f)
    for atom, ivary in zip(respin.molecule.atoms, respin.ivary.values):
        print(FW("2I5").write([atom.atomic_number, ivary]), file=f)
    # According to the spec a blank line is only for multi-structures but `resp`
    # fails without it.
    print(file=f)
