from repESP.esp_util import EspData, parse_gaussian_esp
from repESP.types import *
from repESP.resp_util import run_resp, Equivalence
from repESP.respin_util import Respin

from my_unittest import TestCase

from io import StringIO


class TestEquivalence(TestCase):

    def test_init_validates_values(self) -> None:

        Equivalence([None, None, 0, 1, 1])

        with self.assertRaises(ValueError):
            Equivalence([0, 0, 2, -1, 2])

        with self.assertRaises(ValueError):
            Equivalence([5, 0, 2, 2, 2])

    def test_init_from_ivary(self) -> None:

        equivalence = Equivalence.from_ivary(Respin.Ivary([0, 0, 2, 2, 2]))
        self.assertListEqual(equivalence.values, [None, None, 1, 1, 1])

        with self.assertRaises(ValueError):
            Equivalence.from_ivary(Respin.Ivary([0, 0, 2, -1, 2]))

        with self.assertRaises(ValueError):
            Equivalence.from_ivary(Respin.Ivary([6, 0, 2, 2, 2]))

    def test_description(self) -> None:

        equivalence = Equivalence([None, None, 0, 1, 1])

        output = StringIO()
        expected_lines = [
            "Atom number 1\n",
            "Atom number 2\n",
            "Atom number 3, equivalenced to atom 1\n",
            "Atom number 4, equivalenced to atom 2\n",
            "Atom number 5, equivalenced to atom 2\n",
        ]

        equivalence.describe(file=output)
        output.seek(0)
        self.assertListEqual(expected_lines, output.readlines())

    def test_description_with_molecule(self) -> None:

        equivalence = Equivalence([None, None, 0, 1, 1])

        output = StringIO()
        expected_lines = [
            "Atom (C) number 1\n",
            "Atom (H) number 2\n",
            "Atom (H) number 3, equivalenced to atom 1\n",
            "Atom (H) number 4, equivalenced to atom 2\n",
            "Atom (H) number 5, equivalenced to atom 2\n",
        ]

        equivalence.describe(atomic_numbers=[6, 1, 1, 1, 1], file=output)
        output.seek(0)
        self.assertListEqual(expected_lines, output.readlines())


class TestResp(TestCase):

    def setUp(self) -> None:
        with open("data/methane/methane_mk.esp", 'r') as f:
            self.esp_data = EspData.from_gaussian(parse_gaussian_esp(f))

        # First stage RESP
        self.respin = Respin(
            title="File generated for unit tests only (can been removed).",
            cntrl=Respin.Cntrl(
                nmol=1,
                ihfree=1,
                ioutopt=1,
                qwt=0.0005,
            ),
            wtmol=1.0,
            subtitle="Resp charges for organic molecule",
            charge=0,
            iuniq=5,
            atomic_numbers=[6, 1, 1, 1, 1],
            ivary=Respin.Ivary([0, 0, 0, 0, 0])
        )

        self.result_charges = run_resp(
            self.esp_data,
            self.respin
        )

    def test_resp_charges(self) -> None:
        self.assertListsAlmostEqual(
            self.result_charges,
            # Expected values from resp calculations done a while back
            [-0.407205, 0.101907, 0.101695, 0.101695, 0.101907]
        )
