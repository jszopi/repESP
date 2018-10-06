from .charges import Charge, ChargeType
from .fields import Esp
from .exceptions import InputFormatError
from .types import Coords, Molecule

from dataclasses import dataclass
import re
from typing import List, Mapping, Optional, TextIO, Tuple


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

charge_type_to_identifier_line: Mapping[ChargeType, str] = {
    ChargeType.MULLIKEN: " Mulliken charges:",
    ChargeType.MK: " Merz-Kollman atomic radii used.",
    ChargeType.CHELP: " Francl (CHELP) atomic radii used.",
    ChargeType.CHELPG: " Breneman (CHELPG) radii used.",
    ChargeType.HLY: " Generate Potential Derived Charges using the Hu-Lu-Yang model",
    ChargeType.NPA: " Summary of Natural Population Analysis:                  "
}

identifier_line_to_charge_type = {
    v: k
    for k, v in charge_type_to_identifier_line.items()
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
    List[Charge]
        List of charges in order of occurrence in output file.

    Raises
    ------
    InputFormatError
        Raised when parsing of the input file fails.
    IndexError
        Raised whenn the requested occurence of the charges section cannot
        be found could not be found in the output file.
    NotImplementedError
        Raised when finding and parsing the charges section of the requested
        charge_type is not implemented.
    """
    return _get_charges_section_from_log(f, charge_type, verify_against, occurrence).charges


def get_esp_fit_stats_from_log(
    f: TextIO,
    charge_type: ChargeType,
    verify_against: Optional[Molecule]=None,
    occurrence: int=-1
) -> Tuple[Esp, float]:
    """Extract ESP fit statistics from charges in Gaussian output

    See documentation of `get_charges_from_log` for parameters and raised exceptions.

    Returns
    -------
    Tuple[Esp, float]
        RMS and RRMS.
    """
    if charge_type not in esp_charge_types:
        raise ValueError(
            f"ESP fit statistics can only be extracted for ESP charges, got "
            f"{charge_type} instead."
        )

    charges_section = _get_charges_section_from_log(f, charge_type, verify_against, occurrence)

    assert isinstance(charges_section, _EspChargesSectionData)
    return (charges_section.rms, charges_section.rrms)


@dataclass
class _ChargesSectionData:
    charge_type: ChargeType
    charges: List[Charge]


@dataclass
class _EspChargesSectionData(_ChargesSectionData):
    rms: Esp
    rrms: float
    # TODO: Could be extended with number of points and molecule


def _get_charges_section_from_log(
    f: TextIO,
    charge_type: ChargeType,
    verify_against: Optional[Molecule],
    occurrence: int
) -> _ChargesSectionData:

    if charge_type not in handled_charge_types:
        raise InputFormatError(
            f"Extracting charge type {charge_type} from Gaussian output file is not supported."
        )

    charges_sections = _get_charges_sections(f, charge_type)

    parsed_charges_sections = [
        _parse_charges_section(charges_section, charge_type)
        for charges_section in charges_sections
    ]

    charges_sections_of_correct_type = [
        charges_section
        for charges_section in parsed_charges_sections
        if charges_section.charge_type == charge_type
    ]

    try:
        selected_charges_section = charges_sections_of_correct_type[occurrence]
    except IndexError:
        raise IndexError(
            f"Cannot find occurrence {occurrence} in a list of recognized charges "
            "sections in the output. Check if the program that produced the output "
            "finished without errors. If you're sure that the charge section "
            "should appear at least the specified number of times, please "
            "submit a bug report attaching the file that failed parsing."
        )

    if not _verify_charges_section(selected_charges_section, verify_against):
        raise InputFormatError(
            "Charges from log file failed verification against given molecule."
        )

    return selected_charges_section


def _verify_charges_section(
    charges_section: _ChargesSectionData,
    verify_against: Optional[Molecule]
) -> bool:
    # TODO: This could be extended to check atom identities if those get parsed
    if verify_against is None:
        return True
    else:
        return len(verify_against.atoms) == len(charges_section.charges)


def _get_charges_sections(
    f: TextIO,
    charge_type: ChargeType
) -> List[List[str]]:
    """Extract all charges sections which *may* be of the given type

    Further verification of charge type is necessary based on parsing the section.
    """
    charges_sections: List[List[str]] = []
    current_section: Optional[List[str]] = None
    for line in f:
        line = line.rstrip('\n')
        if _is_section_start(line, charge_type):
            if current_section is not None:
                raise InputFormatError(
                    "Encountered start of new charge section start while "
                    "parsing a charge section. Please submit a bug report "
                    "attaching the input file that failed parsing."
                )
            current_section = []

        if current_section is not None:
            current_section.append(line)
            # Section end lines are less generic, hence we're only checking for
            # them when inside a section.
            if _is_section_end(line, charge_type):
                charges_sections.append(current_section)
                current_section = None

    return charges_sections


def _parse_charges_section(section: List[str], charge_type: ChargeType) -> _ChargesSectionData:

    if charge_type in esp_charge_types:
        return _parse_esp_charge_section(section)
    elif charge_type == ChargeType.NPA:
        charges = []
        for line in section[6:-1]:
            _symbol, _label, charge, *_other = line.split()
            charges.append(Charge(charge))
        return _ChargesSectionData(
            ChargeType.NPA,
            charges
        )
    elif charge_type == ChargeType.MULLIKEN:
        charges = []
        for line in section[2:-1]:
            _label, _symbol, charge = line.split()
            charges.append(Charge(charge))
        return _ChargesSectionData(
            ChargeType.MULLIKEN,
            charges
        )
    else:
        raise NotImplementedError(
            f"Parsing of charge section for charge type {charge_type} is not implemented."
        )


def _parse_esp_charge_section(section: List[str]) -> _EspChargesSectionData:

    try:
        charge_type = identifier_line_to_charge_type[section[0]]
    except KeyError:
        # This shouldn't happen as the start line already matched in _is_section_start
        raise InputFormatError(
            "Failed parsing first line of ESP charges section. Please submit "
            " a bug report attaching the input file that failed parsing."
        )

    charges_and_stats_re = re.compile(" Charges from ESP fit, RMS=\s+(\d+\.\d+) RRMS=\s+(\d+\.\d+):$")

    for i, line in enumerate(section):
        matched_charges_and_stats = charges_and_stats_re.match(line)
        if matched_charges_and_stats is not None:
            rms = Esp(matched_charges_and_stats.group(1))
            rrms = float(matched_charges_and_stats.group(2))
            break

    charges = []
    for line in section[i+3:]:
        if _is_section_end(line, charge_type):
            break
        _label, _symbol, charge = line.split()
        charges.append(Charge(charge))

    return _EspChargesSectionData(
        charge_type,
        charges,
        rms,
        rrms
    )


def _is_section_start(line: str, charge_type: ChargeType) -> bool:

    try:
        return line == charge_type_to_identifier_line[charge_type]
    except KeyError:
        raise NotImplementedError(
            f"Detection of start line for charge type {charge_type} is not implemented."
        )


def _is_section_end(line: str, charge_type: ChargeType) -> bool:

    if charge_type in esp_charge_types or charge_type == ChargeType.MULLIKEN:
        # TODO: This will not work with IOp(6/50=1), where all I know is that
        # the termination line starts with ' Charges' but that yields false
        # positives, so I have to find out a more specific regex.
        return line.startswith(' Sum of ')
    elif charge_type == ChargeType.NPA:
        return line.startswith(' =======')
    else:
        raise NotImplementedError(
            f"Detection of termination line for charge type {charge_type} is not implemented."
        )
