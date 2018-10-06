from repESP.charges import Charge
from repESP.resp_charges_util import parse_resp_charges, write_resp_charges

from my_unittest import TestCase

from io import StringIO

class TestRespCharges(TestCase):

    def setUp(self) -> None:
        self.charges = list(map(
            Charge,
            [
                -0.162076,  0.040051,  0.042711,  0.033893, -0.001304, 0.012177,
                 0.011888,  0.200011, -0.001915, -0.223809,  0.044483, 0.048844,
                 0.058568, -0.315388,  0.081265,  0.062636,  0.067967
            ]
        ))

    def test_parsing(self) -> None:
        with open("tests/test_resp_charges.qout") as f:
            parsed_charges = parse_resp_charges(f)
        self.assertListsAlmostEqual(self.charges, parsed_charges)

    def test_writing(self) -> None:

        written = StringIO()
        write_resp_charges(written, self.charges)
        written.seek(0)

        with open("tests/test_resp_charges.qout") as f:
            self.assertListEqual(f.readlines(), written.readlines())
