from repESP.types import *
from repESP.charges import *
from repESP.esp_util import GaussianEspData, parse_gaussian_esp

from my_unittest import TestCase

from typing import TextIO


class TestGaussianEspData(TestCase):

    def setUp(self) -> None:

        self.expected_gaussian_esp_data = GaussianEspData(
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
            NumericField(
                NonGridMesh(
                    [
                        make_coords(3.9684249, 0, 0),
                        make_coords(3.4367568, -0.99210622,  1.7183784),
                        make_coords(3.4367568,  0.99210622, -1.7183784)
                    ],
                ),
                [
                    make_esp(-0.26293556e-2),
                    make_esp(-0.28665426e-2),
                    make_esp(-0.28665426e-2)
                ]
            )
        )


    def test_parsing(self) -> None:

        with open("tests/test_gaussian.esp") as f:
            gaussian_esp = parse_gaussian_esp(f)
