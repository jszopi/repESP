"""Parsing the Gaussian output format (.log or .out files)"""

from .charges import Charge
from .fields import Esp
from .exceptions import InputFormatError
from .types import Coords, Molecule

from abc import ABC, abstractmethod
from dataclasses import dataclass
import re
from typing import List, Mapping, Optional, TextIO, Tuple


@dataclass
class ChargesSectionData:
    charges: List[Charge]
    # TODO: Could be extended with Molecule[Atom]


class ChargesSectionParser(ABC):

    @abstractmethod
    def is_section_start(self, line: str) -> bool:
        pass

    @abstractmethod
    def is_section_end(self, line: str) -> bool:
        pass

    @abstractmethod
    def parse_section(self, section: List[str]) -> ChargesSectionData:
        pass


@dataclass
class EspChargesSectionData(ChargesSectionData):
    rms: Esp
    rrms: float
    # TODO: Could be extended with number of points


class EspChargesSectionParser(ChargesSectionParser):

    def is_section_end(self, line: str) -> bool:
        # TODO: This will not work with IOp(6/50=1), where all I know is that
        # the termination line starts with ' Charges' but that yields false
        # positives, so I have to find out a more specific regex.
        return line.startswith(' Sum of ')

    def parse_section(self, section: List[str]) -> EspChargesSectionData:

        charges_and_stats_re = re.compile(" Charges from ESP fit, RMS=\s+(\d+\.\d+) RRMS=\s+(\d+\.\d+):$")

        for i, line in enumerate(section):
            matched_charges_and_stats = charges_and_stats_re.match(line)
            if matched_charges_and_stats is not None:
                rms = Esp(matched_charges_and_stats.group(1))
                rrms = float(matched_charges_and_stats.group(2))
                break

        charges = []
        for line in section[i+3:]:
            if self.is_section_end(line):
                break
            _label, _symbol, charge = line.split()
            charges.append(Charge(charge))

        return EspChargesSectionData(
            charges,
            rms,
            rrms
        )


class MullikenChargeSectionParser(ChargesSectionParser):

    def is_section_start(self, line: str) -> bool:
        return line == " Mulliken charges:"

    def is_section_end(self, line: str) -> bool:
        # TODO: Same as this function in EspChargesSectionParser
        return line.startswith(' Sum of ')

    def parse_section(self, section: List[str]) -> ChargesSectionData:
        charges = []
        for line in section[2:-1]:
            _label, _symbol, charge = line.split()
            charges.append(Charge(charge))
        return ChargesSectionData(charges)


class MkChargeSectionParser(EspChargesSectionParser):

    def is_section_start(self, line: str) -> bool:
        return line == " Merz-Kollman atomic radii used."


class ChelpChargeSectionParser(EspChargesSectionParser):

    def is_section_start(self, line: str) -> bool:
        return line == " Francl (CHELP) atomic radii used."


class ChelpgChargeSectionParser(EspChargesSectionParser):

    def is_section_start(self, line: str) -> bool:
        return line == " Breneman (CHELPG) radii used."


class HlyChargeSectionParser(EspChargesSectionParser):

    def is_section_start(self, line: str) -> bool:
        return line == " Generate Potential Derived Charges using the Hu-Lu-Yang model"


class NpaChargeSectionParser(ChargesSectionParser):

    def is_section_start(self, line: str) -> bool:
        return line == " Summary of Natural Population Analysis:                  "

    def is_section_end(self, line: str) -> bool:
        return line.startswith(' =======')

    def parse_section(self, section: List[str]) -> ChargesSectionData:
        charges = []
        for line in section[6:-1]:
            _symbol, _label, charge, *_other = line.split()
            charges.append(Charge(charge))
        return ChargesSectionData(charges)


def get_charges_from_log(
    f: TextIO,
    charges_section_parser: ChargesSectionParser,
    verify_against: Optional[Molecule]=None,
    occurrence: int=-1
) -> List[Charge]:
    """Extract charges from the charges section in Gaussian output

    Parameters
    ----------
    f : TextIO
        The file object corresponding to the Gaussian `.log`/`.out` output file
        from which the charges are to be extracted.
    charges_section_parser: ChargesSectionParser
        Class implementing the `ChargesSectionParser` interface for the desired
        charge type, e.g. `MullikenChargeSectionParser`.
    verify_against: Molecule, optional
        Molecule against which the output is to be verified. Defaults to None.
        Note that currently the verification only involves comparing the number
        of extracted charges against the number of atoms. In the future this may
        be extended to verifying the atom identities (TODO).
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
        Raised when the requested occurence of the charges section cannot
        be found could not be found in the output file.
    """
    charges_section = _get_charges_section_from_log(f, charges_section_parser, occurrence)
    parsed_charges_section = charges_section_parser.parse_section(charges_section)
    _verify_charges_section(parsed_charges_section, verify_against)

    return parsed_charges_section.charges


def get_esp_fit_stats_from_log(
    f: TextIO,
    charges_section_parser: EspChargesSectionParser,
    verify_against: Optional[Molecule]=None,
    occurrence: int=-1
) -> Tuple[Esp, float]:
    """Extract ESP fit statistics from charges section in Gaussian output

    See documentation of `get_charges_from_log` for parameters and raised exceptions.

    Returns
    -------
    Tuple[Esp, float]
        RMS and RRMS.
    """
    charges_section = _get_charges_section_from_log(f, charges_section_parser, occurrence)
    parsed_charges_section = charges_section_parser.parse_section(charges_section)
    _verify_charges_section(parsed_charges_section, verify_against)

    return (parsed_charges_section.rms, parsed_charges_section.rrms)


def _get_charges_section_from_log(
    f: TextIO,
    charges_section_parser: ChargesSectionParser,
    occurrence: int
) -> List[str]:

    charges_sections = _get_charges_sections(f, charges_section_parser)

    try:
        selected_charges_section = charges_sections[occurrence]
    except IndexError:
        raise IndexError(
            f"Cannot find occurrence {occurrence} in a list of recognized charges "
            f"sections in the output. Check if the program that produced the output "
            f"finished without errors. If you're sure that the charge section "
            f"should appear at least the specified number of times, please "
            f"submit a bug report attaching the file that failed parsing."
        )

    return selected_charges_section


def _verify_charges_section(
    charges_section: ChargesSectionData,
    verify_against: Optional[Molecule]
) -> None:
    # TODO: This could be extended to check atom identities if those get parsed
    if verify_against is None:
        return
    elif len(verify_against.atoms) != len(charges_section.charges):
        raise InputFormatError(
            "Charges from log file failed verification against given molecule."
        )
    else:
        return


def _get_charges_sections(
    f: TextIO,
    charges_section_parser: ChargesSectionParser
) -> List[List[str]]:
    """Extract all charges sections which *may* be of the given type

    Further verification of charge type is necessary based on parsing the section.
    """
    charges_sections: List[List[str]] = []
    current_section: Optional[List[str]] = None
    for line in f:
        line = line.rstrip('\n')
        if charges_section_parser.is_section_start(line):
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
            if charges_section_parser.is_section_end(line):
                charges_sections.append(current_section)
                current_section = None

    return charges_sections
