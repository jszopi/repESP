from repESP.charges import *
from repESP.fields import *
from repESP.types import *
from repESP.cube_util import parse_esp_cube
from repESP.gaussian_util import get_charges_from_log, get_esp_fit_stats_from_log
from repESP.gaussian_util import ChargesSectionParser, EspChargesSectionParser
from repESP.gaussian_util import MullikenChargeSectionParser, MkChargeSectionParser
from repESP.gaussian_util import ChelpgChargeSectionParser, NpaChargeSectionParser

from my_unittest import TestCase

from io import StringIO
from typing import TextIO


class TestFromLog(TestCase):

    def setUp(self) -> None:
        with open("data/methane/methane_esp.cub") as f:
            cube = parse_esp_cube(f)
            self.molecule = cube.molecule

    @staticmethod
    def concatenate(filenames: List[str]) -> TextIO:
        files_contents = []
        for filename in filenames:
            with open(filename) as f:
                files_contents.append(f.read())

        return StringIO("".join(files_contents))

    @staticmethod
    def replace_first_occurrence(
        f: TextIO,
        line_to_replace: str,
        line_to_replace_with: str
    ) -> StringIO:
        result_lines = []
        is_first_occurrence = True
        for line in f:
            if line == line_to_replace and is_first_occurrence:
                result_lines.append(line_to_replace_with)
                is_first_occurrence = False
            else:
                result_lines.append(line)
        return StringIO("".join(result_lines))


class TestGetChargesFromLog(TestFromLog):

    def common(
        self,
        f: TextIO,
        charges_section_parser: ChargesSectionParser,
        expected: List[float],
        occurrence: int=-1
    ) -> None:
        charges = get_charges_from_log(
            f,
            charges_section_parser,
            verify_against=self.molecule,
            occurrence=occurrence
        )

        expected_charges = [Charge(x) for x in expected]
        self.assertListsAlmostEqual(charges, expected_charges)

    def test_mulliken(self) -> None:
        with open("data/methane/methane_mk.log") as f:
            self.common(
                f,
                MullikenChargeSectionParser(),
                [-0.437226, 0.109307, 0.109307, 0.109307, 0.109307]
            )

    def test_mk(self) -> None:
        with open("data/methane/methane_mk.log") as f:
            self.common(
                f,
                MkChargeSectionParser(),
                [-0.500314, 0.125323, 0.124834, 0.124834, 0.125323]
            )

    def test_chelpg(self) -> None:
        with open("data/methane/methane_chelpg.log") as f:
            self.common(
                f,
                ChelpgChargeSectionParser(),
                [-0.344877, 0.086219, 0.086219, 0.086219, 0.086219]
            )

    def test_npa(self) -> None:
        with open("data/methane/methane_nbo.log") as f:
            self.common(
                f,
                NpaChargeSectionParser(),
                [-0.79151, 0.19788, 0.19788, 0.19788, 0.19788]
            )

    def test_mk_from_combined_mk_and_npa(self) -> None:
        self.common(
            self.concatenate(["data/methane/methane_mk.log", "data/methane/methane_nbo.log"]),
            MkChargeSectionParser(),
            [-0.500314, 0.125323, 0.124834, 0.124834, 0.125323]
        )

    def test_npa_from_combined_mk_and_npa(self) -> None:
        self.common(
            self.concatenate(["data/methane/methane_mk.log", "data/methane/methane_nbo.log"]),
            NpaChargeSectionParser(),
            [-0.79151, 0.19788, 0.19788, 0.19788, 0.19788]
        )

    def test_mk_from_combined_mk_and_chelpg(self) -> None:
        self.common(
            self.concatenate(["data/methane/methane_mk.log", "data/methane/methane_chelpg.log"]),
            MkChargeSectionParser(),
            [-0.500314, 0.125323, 0.124834, 0.124834, 0.125323]
        )

    def test_chelpg_from_combined_mk_and_chelpg(self) -> None:
        self.common(
            self.concatenate(["data/methane/methane_mk.log", "data/methane/methane_chelpg.log"]),
            ChelpgChargeSectionParser(),
            [-0.344877, 0.086219, 0.086219, 0.086219, 0.086219]
        )

    def test_two_mk_occurrences(self) -> None:
        f = self.replace_first_occurrence(
            self.concatenate(["data/methane/methane_mk.log", "data/methane/methane_mk.log"]),
            "     1  C   -0.500314\n",
            "     1  C    0.123456\n"
        )
        self.common(
            f,
            MkChargeSectionParser(),
            [0.123456, 0.125323, 0.124834, 0.124834, 0.125323],
            occurrence=0,
        )
        f.seek(0)
        self.common(
            f,
            MkChargeSectionParser(),
            [-0.500314, 0.125323, 0.124834, 0.124834, 0.125323],
            occurrence=1
        )


class TestGetEspFitStatsFromLog(TestFromLog):

    def common(
        self,
        f: TextIO,
        charges_section_parser: EspChargesSectionParser,
        expected: Tuple[float, float],
        occurrence: int=-1
    ) -> None:
        rms, rrms = get_esp_fit_stats_from_log(
            f,
            charges_section_parser,
            verify_against=self.molecule,
            occurrence=occurrence
        )

        expected_rms = Esp(expected[0])
        expected_rrms = expected[1]

        self.assertAlmostEqual(rms, expected_rms)
        self.assertAlmostEqual(rrms, expected_rrms)

    def test_mk(self) -> None:
        with open("data/methane/methane_mk.log") as f:
            self.common(
                f,
                MkChargeSectionParser(),
                (0.00069, 0.35027)
            )

    def test_chelpg(self) -> None:
        with open("data/methane/methane_chelpg.log") as f:
            self.common(
                f,
                ChelpgChargeSectionParser(),
                (0.00121, 0.62228)
            )

    def test_mk_from_combined_mk_and_npa(self) -> None:
        self.common(
            self.concatenate(["data/methane/methane_mk.log", "data/methane/methane_nbo.log"]),
            MkChargeSectionParser(),
            (0.00069, 0.35027)
        )

    def test_mk_from_combined_mk_and_chelpg(self) -> None:
        self.common(
            self.concatenate(["data/methane/methane_mk.log", "data/methane/methane_chelpg.log"]),
            MkChargeSectionParser(),
            (0.00069, 0.35027)
        )

    def test_chelpg_from_combined_mk_and_chelpg(self) -> None:
        self.common(
            self.concatenate(["data/methane/methane_mk.log", "data/methane/methane_chelpg.log"]),
            ChelpgChargeSectionParser(),
            (0.00121, 0.62228)
        )

    def test_two_mk_occurrences(self) -> None:
        f = self.replace_first_occurrence(
            self.concatenate(["data/methane/methane_mk.log", "data/methane/methane_mk.log"]),
            " Charges from ESP fit, RMS=   0.00069 RRMS=   0.35027:\n",
            " Charges from ESP fit, RMS=   1.23456 RRMS=   7.65432:\n"
        )
        self.common(
            f,
            MkChargeSectionParser(),
            (1.23456, 7.65432),
            occurrence=0,
        )
        f.seek(0)
        self.common(
            f,
            MkChargeSectionParser(),
            (0.00069, 0.35027),
            occurrence=1
        )
