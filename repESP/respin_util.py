from dataclasses import dataclass, asdict
from fortranformat import FortranRecordWriter as FW
from itertools import zip_longest
import re
import sys
from typing import Dict, List, Optional, TextIO, Tuple, TypeVar, Union

from .exceptions import InputFormatError
from .util import _get_symbol, _zip_exact


@dataclass
class Equivalence:
    # Zero-indexed, None if not equivalenced to any atom. Note similarities in
    # implementation with Respin.Ivary - perhaps should be refactored.
    values: List[Optional[int]]

    def __post_init__(self):
        for i, elem in enumerate(self.values):
            if elem is not None and (elem < 0 or elem >= len(self.values)):
                raise ValueError(
                    f"Value number {i} is not valid as equivalence information "
                    f"in a molecule of {len(self.values)}."
                )

    def describe(self, atomic_numbers: Optional[List[int]]=None, file=sys.stdout):
        """Verbosely report the equivalence information."""
        if atomic_numbers is not None and len(atomic_numbers) != len(self.values):
            raise ValueError(
                f"The number of atoms ({len(atomic_numbers)} is not the same "
                f"as the number of equivalence values ({len(self.values)}."
            )

        zipped = zip_longest(self.values, atomic_numbers if atomic_numbers is not None else [])

        for i, (equivalence, atomic_number) in enumerate(zipped):
            identity = _get_symbol(atomic_number) if atomic_numbers is not None else None
            id_str = f" ({identity})" if identity is not None else ""
            equivalence_str = f", equivalenced to atom {equivalence+1}" if equivalence is not None else ""
            print(f"Atom{id_str} number {i+1}{equivalence_str}", file=file)

    @classmethod
    def from_ivary(cls, ivary: "Respin.Ivary"):
        """Get atom equivalence information from an Respin.Ivary object.

        Warning: You probably don't mean to use this function. Use the
        get_equivalence function instead.

        Ivary objects are specific to `resp` program input and thus may
        not provide information about atom equivalence. The `respin` file may
        have been generated to perform any custom fitting by RESP. Only use
        this function when you're certain that the `respin` file contains all
        the equivalence information that you need.
        """
        return cls([
            # Whether `ivary - 1` fulfills other preconditions will be checked in __post_init__
            None if ivary_value == 0 else ivary_value - 1
            for ivary_value in ivary.values
        ])


def get_equivalence(ivary1: "Respin.Ivary", ivary2: "Respin.Ivary") -> Equivalence:
    """Get atom equivalence from two input for two 2-stage `resp`

    Derive atom equivalence based on the data in two `.respin` files (represented
    by the `Respin` objects) created for the purpose two-stage fitting with the
    `resp` program. The files can be generated with the `respgen` program with
    the following commands:

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
        for ivary1_value, ivary2_value in _zip_exact(ivary1.values, ivary2.values)
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
        nmol: int = 1
        ihfree: int = 1
        irstrnt: int = 1
        qwt: float = 0

        def __post_init__(self) -> None:
            Respin._check_value("inopt", self.inopt, [0, 1])
            Respin._check_value("ioutopt", self.ioutopt, [0, 1])
            Respin._check_value("iqopt", self.iqopt, [1, 2, 3])
            Respin._check_value("nmol", self.nmol, [1])
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

        def describe(self, atomic_numbers: Optional[List[int]]=None, file: TextIO=sys.stdout):
            """Verbosely report the ``ivary`` actions"""
            # I'm undecided about printing functions. Can `file` point to managed log object?
            if atomic_numbers is not None and len(atomic_numbers) != len(self.values):
                raise ValueError(
                    f"The number of atoms ({len(atomic_numbers)} is not the same "
                    f"as the number of ivary values ({len(self.values)}."
                )

            zipped = zip_longest(self.values, atomic_numbers if atomic_numbers is not None else [])

            for i, (ivary, atomic_number) in enumerate(zipped):
                identity = _get_symbol(atomic_number) if atomic_number is not None else None
                id_str = f" ({identity})" if identity is not None else ""

                if ivary < 0:
                    # TODO: This could also report at what value if charges are provided
                    ivary_str = ", frozen"
                elif ivary > 0:
                    ivary_str = f", equivalenced to atom {ivary}"
                else:
                    ivary_str = ""

                print(f"Atom{id_str} number {i+1}{ivary_str}", file=file)

        @classmethod
        def from_equivalence(cls, equivalence: Equivalence):
            """Created ivary instructions based on equivalence information

            Note: the resulting ivary instructions will correspond to fitting
            with equivalent atoms assigned identical charges.
            """
            return cls([0 if val is None else val + 1 for val in equivalence.values])


    title: str
    cntrl: Cntrl
    wtmol: float
    subtitle: str
    charge: int
    iuniq: int
    atomic_numbers: List[int]
    ivary: Ivary

    def __post_init__(self) -> None:
        Respin._check_value("wtmol", "{:.1f}".format(self.wtmol), ["1.0"])

        if len(self.atomic_numbers) != len(self.ivary.values):
            raise ValueError(
                f"Number of atoms ({self.atomic_numbers}) does not match number "
                f"of ivary values ({self.ivary.values})."
            )

        if self.iuniq != len(self.ivary.values):
            raise ValueError(
                f"`iuniq` value ({self.iuniq}) doesn't match the number of atoms "
                f"({self.atomic_numbers})."
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

    return Respin.Cntrl(**kwargs)  # type: ignore # (not sure why not recognized)


def _parse_respin(f) -> Respin:

    get_line = lambda: f.readline().rstrip('\n')
    title = get_line()
    while get_line() != " &cntrl":
        pass

    cntrl = _parse_cntrl(f)
    wtmol = float(get_line().strip())
    subtitle = get_line()

    charge_and_iuniq = get_line()
    if len(charge_and_iuniq.split()) != 2:
        raise InputFormatError(
            f"Expected two ints for the line specifying charge and iuniq, found:\n{charge_and_iuniq}"
        )

    charge = int(charge_and_iuniq.split()[0])
    iuniq = int(charge_and_iuniq.split()[1])

    atomic_numbers: List[int] = []
    ivary = Respin.Ivary([])

    for line in f:
        if line.rstrip('\n') == "":
            break
        if len(line.split()) != 2:
            raise InputFormatError(
                f"Expected two ints for the line specifying atom and ivary, found:\n{line}"
            )

        atomic_numbers.append(int(line.split()[0]))
        ivary_value = int(line.split()[1])
        # `respgen` uses a value of -99 but internally we use -1 as per resp spec.
        ivary.values.append(ivary_value if ivary_value != -99 else -1)

    return Respin(
        title,
        cntrl,
        wtmol,
        subtitle,
        charge,
        iuniq,
        atomic_numbers,
        ivary
    )


def _write_cntrl(f: TextIO, cntrl: Respin.Cntrl, skip_defaults: bool) -> None:
    default_cntrl: Dict[str, Union[int, float]] = asdict(Respin.Cntrl())
    print(" &cntrl\n", file=f)
    for key, value in asdict(cntrl).items():  # type: ignore
        if key == "qwt":
            print(" {} = {:.5f},".format(key, value), file=f)
        else:
            if not skip_defaults or value != default_cntrl[key]:
                print(" {} = {},".format(key, value), file=f)
    print("\n &end", file=f)


def _write_respin(f: TextIO, respin: Respin, skip_cntrl_defaults: bool=True) -> None:

    print(respin.title, file=f)
    print(file=f)
    _write_cntrl(f, respin.cntrl, skip_cntrl_defaults)
    print(FW("F7.1").write([respin.wtmol]), file=f)
    print(respin.subtitle, file=f)
    print(FW("2I5").write([respin.charge, respin.iuniq]), file=f)
    for atomic_number, ivary in zip(respin.atomic_numbers, respin.ivary.values):
        print(FW("2I5").write([atomic_number, ivary]), file=f)
    # According to the spec a blank line is only for multi-structures but `resp`
    # fails without it.
    print(file=f)
