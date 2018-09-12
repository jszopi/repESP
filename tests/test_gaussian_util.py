from repESP.types import *
from repESP.charges import *
from repESP.cube_util import parse_esp_cube
from repESP.gaussian_util import get_charges_from_log, get_esp_fit_stats_from_log

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


class TestGetChargesFromLog(TestFromLog):

    def common(
        self,
        filenames: List[str],
        charge_type: ChargeType,
        expected: List[float]
    ) -> None:
        charges = get_charges_from_log(
            self.concatenate(filenames),
            charge_type,
            verify_against=self.molecule
        )

        expected_charges = [make_charge(x) for x in expected]
        self.assertListsAlmostEqual(charges, expected_charges)

    def test_mulliken(self) -> None:
        self.common(
            ["data/methane/methane_mk.log"],
            ChargeType.MULLIKEN,
            [-0.437226, 0.109307, 0.109307, 0.109307, 0.109307]
        )

    def test_mk(self) -> None:
        self.common(
            ["data/methane/methane_mk.log"],
            ChargeType.MK,
            [-0.500314, 0.125323, 0.124834, 0.124834, 0.125323]
        )

    def test_chelpg(self) -> None:
        self.common(
            ["data/methane/methane_chelpg.log"],
            ChargeType.CHELPG,
            [-0.344877, 0.086219, 0.086219, 0.086219, 0.086219]
        )

    def test_npa(self) -> None:
        self.common(
            ["data/methane/methane_nbo.log"],
            ChargeType.NPA,
            [-0.79151, 0.19788, 0.19788, 0.19788, 0.19788]
        )

    def test_mk_from_combined_mk_and_npa(self) -> None:
        self.common(
            ["data/methane/methane_mk.log", "data/methane/methane_nbo.log"],
            ChargeType.MK,
            [-0.500314, 0.125323, 0.124834, 0.124834, 0.125323]
        )

    def test_npa_from_combined_mk_and_npa(self) -> None:
        self.common(
            ["data/methane/methane_mk.log", "data/methane/methane_nbo.log"],
            ChargeType.NPA,
            [-0.79151, 0.19788, 0.19788, 0.19788, 0.19788]
        )

    def test_mk_from_combined_mk_and_chelpg(self) -> None:
        self.common(
            ["data/methane/methane_mk.log", "data/methane/methane_chelpg.log"],
            ChargeType.MK,
            [-0.500314, 0.125323, 0.124834, 0.124834, 0.125323]
        )

    def test_chelpg_from_combined_mk_and_chelpg(self) -> None:
        self.common(
            ["data/methane/methane_mk.log", "data/methane/methane_chelpg.log"],
            ChargeType.CHELPG,
            [-0.344877, 0.086219, 0.086219, 0.086219, 0.086219]
        )


class TestGetEspFitStatsFromLog(TestFromLog):

    def common(
        self,
        filenames: List[str],
        charge_type: ChargeType,
        expected: Tuple[float, float]
    ) -> None:
        rms, rrms = get_esp_fit_stats_from_log(
            self.concatenate(filenames),
            charge_type,
            verify_against=self.molecule
        )

        expected_rms = make_esp(expected[0])
        expected_rrms = expected[1]

        self.assertAlmostEqual(rms, expected_rms)
        self.assertAlmostEqual(rrms, expected_rrms)

    def test_mk(self) -> None:
        self.common(
            ["data/methane/methane_mk.log"],
            ChargeType.MK,
            (0.00069, 0.35027)
        )

    def test_chelpg(self) -> None:
        self.common(
            ["data/methane/methane_chelpg.log"],
            ChargeType.CHELPG,
            (0.00121, 0.62228)
        )

    def test_mk_from_combined_mk_and_npa(self) -> None:
        self.common(
            ["data/methane/methane_mk.log", "data/methane/methane_nbo.log"],
            ChargeType.MK,
            (0.00069, 0.35027)
        )

    def test_mk_from_combined_mk_and_chelpg(self) -> None:
        self.common(
            ["data/methane/methane_mk.log", "data/methane/methane_chelpg.log"],
            ChargeType.MK,
            (0.00069, 0.35027)
        )

    def test_chelpg_from_combined_mk_and_chelpg(self) -> None:
        self.common(
            ["data/methane/methane_mk.log", "data/methane/methane_chelpg.log"],
            ChargeType.CHELPG,
            (0.00121, 0.62228)
        )
