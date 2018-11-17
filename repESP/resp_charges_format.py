"""Parsing and writing charges to format handled by the ``resp`` program

This is the file expected and created by the ``resp`` program (``qin`` and
``qout`` arguments). It is specified with the following Fortran formatting:
"8F10.6" i.e. eight floating point numbers per line, with 10 characters per
number and 6 digits after the decimal point.

Example
-------

The following file describes the charges of a molecule with 14 atoms::

    -0.383354  0.192503  0.192503  0.192503 -0.383354  0.192503  0.192503  0.192503
    -0.383354  0.192503  0.192503  0.192503  0.129430  0.288107

"""

from repESP.charges import Charge
from repESP.exceptions import InputFormatError

from fortranformat import FortranRecordWriter as FW, FortranRecordReader as FR
from functools import reduce
from operator import add
from typing import List, TextIO


def parse_resp_charges(f: TextIO) -> List[Charge]:
    """Parse a file in the ``resp`` charges format

    Parameters
    ----------
    f : TextIO
        File object opened in read mode containing charges in the ``resp`` format.

    Raises
    ------
    InputFormatError
        Raised when the file does not follow the expected format.

    Returns
    -------
    typing.List[Charge]
        List of charges described in the given input file.
    """
    formatter = FR("8F10.6")
    try:
        return list(map(
            Charge,
            filter(
                lambda elem: elem is not None,
                reduce(add, [formatter.read(line) for line in f], [])
            )
        ))
    except ValueError as e:
        raise InputFormatError(e)


def write_resp_charges(f: TextIO, charges: List[Charge]) -> None:
    """Write given charges to a file in the ``resp`` charges format.

    Parameters
    ----------
    f : TextIO
        File object to which the supplied charges is to be saved. Must be opened
        in write mode.
    charges : typing.List[Charge]
        List of charges to be written to the file.
    """
    print(FW("8F10.6").write(charges), file=f)
