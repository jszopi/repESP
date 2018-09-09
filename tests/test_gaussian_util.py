from repESP.types import *
from repESP.charges import *
from repESP.cube_util import parse_esp_cube
from repESP.gaussian_util import get_charges_from_log

from my_unittest import TestCase

class TestGetChargesFromLog(TestCase):

    def setUp(self) -> None:
        with open("data/methane/methane_esp.cub") as f:
            cube = parse_esp_cube(f)
            self.molecule = cube.molecule

    def common(self, filename: str, charge_type: ChargeType, expected: List[float]) -> None:
        with open(f"data/methane/{filename}", 'r') as f:
            charges = get_charges_from_log(
                f,
                charge_type,
                verify_against=self.molecule
            )

        expected_charges = [make_charge(x) for x in expected]
        self.assertListsAlmostEqual(charges, expected_charges)

    def test_mulliken_charges(self) -> None:
        self.common(
            "methane_mk.log",
            ChargeType.MULLIKEN,
            [-0.437226, 0.109307, 0.109307, 0.109307, 0.109307]
        )

    def test_mk_charges(self) -> None:
        self.common(
            "methane_mk.log",
            ChargeType.MK,
           [-0.500314, 0.125323, 0.124834, 0.124834, 0.125323] 
        )

    def test_chelpg_charges(self) -> None:
        self.common(
            "methane_chelpg.log",
            ChargeType.CHELPG,
           [-0.344877, 0.086219, 0.086219, 0.086219, 0.086219]
        )

    def test_npa_charges(self) -> None:
        self.common(
            "methane_nbo.log",
            ChargeType.NPA,
           [-0.79151, 0.19788, 0.19788, 0.19788, 0.19788]
        )
