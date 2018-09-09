from .charges import Charge, ChargeType, make_charge
from .exceptions import InputFormatError
from .types import Atom, Coords, Molecule

from typing import List, Optional, TextIO, Tuple


handled_charge_types = [
    ChargeType.MULLIKEN,
    ChargeType.MK,
    ChargeType.CHELP,
    ChargeType.CHELPG,
    ChargeType.HLY,
    ChargeType.NPA,
]

esp_charge_types = [
    ChargeType.MK,
    ChargeType.CHELP,
    ChargeType.CHELPG,
    ChargeType.HLY,
]

esp_charge_identifier_line = {
    ' Merz-Kollman atomic radii used.': ChargeType.MK,
    ' Francl (CHELP) atomic radii used.': ChargeType.CHELP,
    ' Breneman (CHELPG) radii used.': ChargeType.CHELPG,
    ' Generate Potential Derived Charges using the Hu-Lu-Yang model:': ChargeType.HLY,
}

esp_charges_section_header = ' Charges from ESP fit,'

section_headers = {
    ChargeType.MULLIKEN: ' Mulliken charges:',
    ChargeType.MK: esp_charges_section_header,
    ChargeType.CHELP: esp_charges_section_header,
    ChargeType.CHELPG: esp_charges_section_header,
    ChargeType.HLY: esp_charges_section_header,
    ChargeType.NPA: ' Summary of Natural Population Analysis:'
}


def get_charges_from_log(
    f: TextIO,
    charge_type: ChargeType,
    verify_against: Optional[Molecule]=None,
    occurrence: int=-1
) -> List[Charge]:
    """Extract charges from the charges section in Gaussian output

    Parameters
    ----------
    f : TextIO
        The file object corresponding to the Gaussian `.log`/`.out` output file
        from which the charges are to be extracted.
    charge_type: ChargeType
        Charge type to be extracted. Supported options are: {handled_charge_types}.
    verify_against: Molecule, optional
        Molecule against which the output is to be verified. Defaults to None.
        Note that currently the verification only involves comparing the number
        of extracted charges against the number of atoms. In the future this may
        be extended to verifying the atom identities or Gaussian labels (TODO).
    occurrence: int, optional
        Determines which charges section to use for extracting the charges.
        Defaults to -1 i.e. the last section.

        On occasion, the output can contain multiple charge sections regarding
        the same charge type. AFAIR, this happens when multiple Gaussian jobs
        were processed. Thus, the sensible default is to select the final one,
        as this is the final optimization. Other values can be specified,
        starting with 0 for the first occurrence.

    Returns
    -------
    List[float]
        List of charges in order of occurrence in output file.

    Raises
    ------
    InputFormatError
    """

    if charge_type not in handled_charge_types:
        raise InputFormatError(
            f"Extracting charge type {charge_type} from Gaussian output file is not supported."
        )

    charge_section = _get_charge_section(f, charge_type, occurrence)
    charges = _parse_charges_section(charge_type, charge_section)

    if verify_against is not None:
        if len(verify_against.atoms) != len(charges):
            raise InputFormatError(
                "Charges from log file failed verification against given molecule."
            )

    return charges


def _get_charge_section(
    f: TextIO,
    charge_type: ChargeType,
    occurrence: int=-1
) -> List[str]:
    """Go to the selected occurrence of input about charges in a Gaussian output file.

    The function includes extra logic to ensure that, if an ESP charge type was
    requested, all ESP charge sections were recognized correctly. This is meant
    as a safeguard in case the parsing fails without other symptoms.
    """
    current_offset = 0
    charge_section_offsets = []
    esp_types = []

    for line in f:
        current_offset += len(line)
        line = line.rstrip('\n')

        # All ESP charges are added here, as they cannot be distinguished just
        # by the header
        if section_headers[charge_type] in line.rstrip():
            charge_section_offsets.append(current_offset)
        # The information about the type of ESP charges is gathered separately
        if charge_type in esp_charge_types and line in esp_charge_identifier_line:
            esp_types.append(esp_charge_identifier_line[line])

    if charge_type in esp_charge_types:
        # Verify if all ESP charge output has been recognized correctly
        if len(esp_types) != len(charge_section_offsets):
            raise InputFormatError('Information about the type of some '
                                   'ESP charges was not recognized.')
        # Filter only the requested ESP charge type
        charge_section_offsets = [elem for i, elem in enumerate(charge_section_offsets) if
                  esp_types[i] == charge_type]

    if not charge_section_offsets:
        raise InputFormatError("Output about charge type '{0}' not found."
                               .format(charge_type))

    try:
        f.seek(charge_section_offsets[occurrence])
    except IndexError:
        raise IndexError(
            "Cannot find occurrence '{0}' in a list of recognized pieces of "
            "output about charges, which length is {1}.".format(
                occurrence,
                len(charge_section_offsets)
            )
        )

    # Skip unnecessary lines
    lines_count = 1
    if charge_type == ChargeType.NPA:
        lines_count = 5
    if charge_type in esp_charge_types:
        lines_count = 2

    for counter in range(lines_count):
        f.readline()

    # Read charge section until termination line
    charge_lines: List[str] = []
    for line in f:
        if _is_section_termination(line, charge_type):
            break
        charge_lines.append(line)

    return charge_lines


def _is_section_termination(line: str, charge_type: ChargeType) -> bool:

    if charge_type == ChargeType.NPA:
        return line.startswith(' =======')
    else:
        # The latter happens when IOp(6/50=1) is requested
        return line.startswith(' Sum of ') or line.startswith(' Charges')


def _parse_charge_line(line: str, charge_type: ChargeType) -> Tuple[int, str, float]:
    if charge_type == ChargeType.NPA:
        symbol, label, charge, *other = line.split()
    else:
        label, symbol, charge = line.split()
    # TODO: The labels and element symbols could be used to further verify
    # against the passed molecule
    return int(label), symbol, float(charge)


def _parse_charges_section(
    charge_type: ChargeType,
    charge_section: List[str]
) -> List[Charge]:
    return [
        make_charge(_parse_charge_line(line, charge_type)[2])
        for line in charge_section
    ]
