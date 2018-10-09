"""Parsing and writing charges to format handled by the `resp` program"""

from .charges import Charge

from fortranformat import FortranRecordWriter as FW, FortranRecordReader as FR
from functools import reduce
from operator import add
from typing import List, TextIO


def parse_resp_charges(f: TextIO) -> List[Charge]:
    formatter = FR("8F10.6")
    return list(map(
        Charge,
        filter(
            lambda elem: elem is not None,
            reduce(add, [formatter.read(line) for line in f], [])
        )
    ))


def write_resp_charges(f: TextIO, charges: List[Charge]) -> None:
    print(FW("8F10.6").write(charges), file=f)
