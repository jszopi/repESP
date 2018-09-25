from dataclasses import dataclass, asdict
from fortranformat import FortranRecordWriter as FW
import re
from typing import Dict, List, TextIO, Tuple, TypeVar, Union

from .exceptions import InputFormatError


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
