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
from .types import Atom, Molecule
from ._util import zip_exact


@dataclass
class Equivalence:
    """Dataclass representing chemical equivalence relations between atoms

    Atoms in a molecule are considered equivalent for the purpose of fitting
    partial charges when they are symmetry-related or fast-exchanging.

    The equivalence information is represented as a list of values stored in
    the `values` attribute. The length of this list must be the same as the
    number of atoms in the molecule it describes. Consecutive values refer
    to consecutive atoms of the molecule. Each value is either None, if the
    described atom is not equivalenced to any other atom, or a zero-based index
    of the atom to which the described atom is equivalenced.

    Example
    -------

    Consider a methane molecule defined as follows:

    >>> methane = Molecule([Atom(atomic_number) for atomic_number in [6, 1, 1, 1, 1]])

    The corresponding equivalence information would be:

    >>> equivalence = Equivalence([None, None, 1, 1, 1])

    There are no atoms equivalent to the carbon atom, and thus its value is
    None. The first hydrogen atom is equivalent to the other three but those
    have not been specified yet, so its value is also None. The remaining
    three hydrogen atoms are equivalent to the first hydrogen atom, which
    zero-based index in the `methane.atoms` list is 1.

    Parameters
    ----------
    values : typing.List[Optional[int]]
        The list of equivalence values.

    Raises
    ------
    ValueError
        Raised when any of the values are outside of the expected bounds. In
        the future the initialization may further verify that there are no cyclic
        relations.

    Attributes
    ----------
    values
        See initialization parameter
    """
    values: List[Optional[int]]

    def __post_init__(self):
        for i, elem in enumerate(self.values):
            if elem is not None and (elem < 0 or elem >= len(self.values)):
                raise ValueError(
                    f"Value number {i} is not valid as equivalence information "
                    f"in a molecule of {len(self.values)}."
                )

    def describe(self, molecule: Optional[Molecule[Atom]]=None) -> str:
        """Verbosely report the equivalence information

        Example
        -------

        >>> print(equivalence.describe(methane))
        Atom (C) number 1
        Atom (H) number 2
        Atom (H) number 3, equivalenced to atom 1
        Atom (H) number 4, equivalenced to atom 2
        Atom (H) number 5, equivalenced to atom 2

        Parameters
        ----------
        molecule : Optional[Molecule[Atom]], optional
            The molecule to which the equivalence information refers. This
            argument is optional and defaults to None. If it is provided, atom
            identities will be included in the output.

        Raises
        ------
        ValueError
            Raised when the number of atoms in the molecule does not match
            the length of the list of values in this object.

        Returns
        -------
        str
            A verbose description of the equivalence information.
        """
        if molecule is not None and len(molecule.atoms) != len(self.values):
            raise ValueError(
                f"The number of atoms ({len(molecule.atoms)} is not the same "
                f"as the number of equivalence values ({len(self.values)}."
            )

        zipped = zip_longest(self.values, molecule.atoms if molecule is not None else [])

        f = io.StringIO()
        for i, (equivalence, atom) in enumerate(zipped):
            atomic_number = atom.symbol if molecule is not None else None
            id_str = f" ({atomic_number})" if atomic_number is not None else ""
            equivalence_str = f", equivalenced to atom {equivalence+1}" if equivalence is not None else ""
            print(f"Atom{id_str} number {i+1}{equivalence_str}", file=f)

        return f.getvalue()

    @classmethod
    def from_ivary(cls, ivary: "Respin.Ivary"):
        """Get atom equivalence information from an `Respin.Ivary` object

        .. note::
            You probably don't mean to use this function. Use the
            `get_equivalence` function instead.

            `Ivary` objects are specific to ``resp`` program input and thus may
            not provide information about atom equivalence. The "respin" file
            may have been generated to perform any custom fitting with
            ``resp``. Only use this function when you're certain that the
            "respin" file contains all the equivalence information that you need.
        """
        return cls([
            # Whether `ivary - 1` fulfills other preconditions will be checked in __post_init__
            None if ivary_value == 0 else ivary_value - 1
            for ivary_value in ivary.values
        ])


def get_equivalence(ivary1: "Respin.Ivary", ivary2: "Respin.Ivary") -> Equivalence:
    """Get atom equivalence from two input for two 2-stage ``resp``

    Derive atom equivalence based on the data in two "respin" files
    (represented by the `Respin` objects) created for the purpose of two-stage
    fitting with the ``resp`` program. The files can be generated with the
    ``respgen`` program with the following commands::

        respgen -i methane.ac -o methane.respin1 -f resp1
        respgen -i methane.ac -o methane.respin2 -f resp2
    """
    # The equivalence logic is explained somewhat inconsistently in the RESP
    # papers but I've additionally re-engineered the ``resp`` program's logic
    # to be sure that reading both the ``respin`` files will give the desired
    # behaviour. In fact, it's pretty simple. In the first stage atoms of the
    # methyl and methylene groups are free, while all the others are
    # equivalenced. In the second stage the former are equivalenced, while all
    # the others are frozen.

    return Equivalence.from_ivary(Respin.Ivary([
        max(ivary1_value, ivary2_value)
        for ivary1_value, ivary2_value in zip_exact(ivary1.values, ivary2.values)
    ]))


def get_equivalence_from_ac(ac) -> Equivalence:
    # TODO: This should wrap the two respgen calls:
    #     respgen -i methane.ac -o methane.respin1 -f resp1
    #     respgen -i methane.ac -o methane.respin2 -f resp2
    # and call get_equivalence on them.
    raise NotImplementedError()


@dataclass
class Respin:

    # Limited to a single molecule and structure at the moment

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
        inopt: int = 0
        ioutopt: int = 0
        iqopt: int = 1
        ihfree: int = 1
        irstrnt: int = 1
        qwt: float = 0

        @property
        def nmol(self) -> int:
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
        # Zero refers to no equivalencing, -1 for frozen and one-indexed value
        # for equivalencing to a particular atom in the molecule.
        values: List[int]

        def __post_init__(self):
            for i, elem in enumerate(self.values):
                if elem < -1 or elem > len(self.values):
                    raise ValueError(
                        f"Value number {i} passed as `ivary` with value {elem}, "
                        f"which is either lower than 0 or outside the list length."
                    )

        def describe(self, molecule: Optional[Molecule[Atom]]=None) -> str:
            """Verbosely report the ``ivary`` actions"""
            # I'm undecided about printing functions. Can `file` point to managed log object?
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
            """Created ivary instructions based on equivalence information

            Note: the resulting ivary instructions will correspond to fitting
            with equivalent atoms assigned identical charges.
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
        return 1.0

    @property
    def iuniq(self) -> int:
        return len(self.molecule.atoms)

    def __post_init__(self) -> None:
        if len(self.molecule.atoms) != len(self.ivary.values):
            raise ValueError(
                f"Number of atoms ({len(self.molecule.atoms)}) does not match number "
                f"of ivary values ({len(self.ivary.values)})."
            )


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
