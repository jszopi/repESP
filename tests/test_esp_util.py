from repESP.types import *
from repESP.charges import *
from repESP.esp_util import GaussianEspData, parse_gaussian_esp
from repESP.esp_util import EspData, parse_resp_esp, write_resp_esp

from my_unittest import TestCase

from io import StringIO


gaussian_esp_data = GaussianEspData(
    0,
    1,
    MoleculeWithCharges(
        Molecule([
            Atom(6, make_coords(0, 0, 0)),
            Atom(1, make_coords(1.23, 0.456, 0.0789))
        ]),
        [
            make_charge(-0.50031415),
            make_charge( 0.12532268)
        ]
    ),
    Dipole(
        make_dipole_moment( 0.38811727e-15),
        make_dipole_moment( 0.42690461e-16),
        make_dipole_moment(-0.29029513e-15)
    ),
    Quadrupole(
        make_quadrupole_moment(-0.26645353e-14),
        make_quadrupole_moment( 0.35527137e-14),
        make_quadrupole_moment(-0.88817842e-15),
        make_quadrupole_moment(-0.13868301e-15),
        make_quadrupole_moment(-0.97158067e-16),
        make_quadrupole_moment( 0.72144168e-15),
    ),
    Field(
        NonGridMesh(
            [
                make_coords( 0.00000000,  0.0000000, 3.9684249),
                make_coords(-0.99210622,  1.7183784, 3.4367568),
                make_coords( 0.99210622, -1.7183784, 3.4367568)
            ],
        ),
        [
            make_esp(-0.26293556e-2),
            make_esp(-0.28665426e-2),
            make_esp(-0.28665426e-2)
        ]
    )
)


class TestGaussianEsp(TestCase):

    def test_parsing(self) -> None:

        with open("tests/test_gaussian.esp") as f:
            parsed_gaussian_esp_data = parse_gaussian_esp(f)

        self.assertAlmostEqualRecursive(
            gaussian_esp_data,
            parsed_gaussian_esp_data
        )


class TestRespEsp(TestCase):

    def test_writing(self) -> None:

        esp_data = EspData.from_gaussian(gaussian_esp_data)

        written = StringIO()
        write_resp_esp(written, esp_data)
        written.seek(0)

        with open("tests/test_resp.esp") as f:
            self.assertListEqual(f.readlines(), written.readlines())

    def test_parsing(self) -> None:

        with open("tests/test_resp.esp") as f:
            esp_data = parse_resp_esp(f)

        expected_esp_data = EspData.from_gaussian(gaussian_esp_data)

        self.assertAlmostEqualRecursive(
            esp_data,
            expected_esp_data,
            places=6
        )
