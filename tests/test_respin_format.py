from repESP.fields import *
from repESP.types import *
from repESP.respin_format import Respin, Equivalence, get_equivalence, parse_respin, write_respin

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

        equivalence.describe(
            molecule=Molecule([Atom(atomic_number) for atomic_number in [6, 1, 1, 1, 1]]),
            file=output
        )
        output.seek(0)
        self.assertListEqual(expected_lines, output.readlines())


class TestRespin(TestCase):

    def setUp(self) -> None:

        self.respin = Respin(
            title="RESP input of type '2' generated by the repESP program",
            cntrl=Respin.Cntrl(
                nmol=1,
                ihfree=1,
                ioutopt=1,
                qwt=0.001,
                iqopt=2
            ),
            wtmol=1.0,
            subtitle="Resp charges for organic molecule",
            charge=0,
            iuniq=5,
            molecule=Molecule([Atom(atomic_number) for atomic_number in [6, 1, 1, 1, 1]]),
            ivary=Respin.Ivary([0, 0, 2, 2, 2])
        )


class TestParsingAndWriting(TestRespin):

    def setUp(self) -> None:
        super().setUp()

    def test_parsing(self) -> None:
        with open("tests/test_respin_format.respin") as f:
            parsed_respin = parse_respin(f)
        self.assertAlmostEqualRecursive(self.respin, parsed_respin)

    def test_writing(self) -> None:
        written = StringIO()
        write_respin(written, self.respin)
        written.seek(0)

        with open("tests/test_respin_format.respin") as f:
            self.assertListEqual(f.readlines(), written.readlines())


class TestIvary(TestRespin):

    def setUp(self) -> None:
        super().setUp()

    def test_ivary_init_validates_values(self) -> None:

        Respin.Ivary([0, 0, 1, -1, 2])

        with self.assertRaises(ValueError):
            Respin.Ivary([0, 0, 2, -99, 2])

        with self.assertRaises(ValueError):
            Respin.Ivary([6, 0, 2, 2, 2])

    def test_ivary_description(self) -> None:

        output = StringIO()
        expected_lines = [
            "Atom number 1\n",
            "Atom number 2\n",
            "Atom number 3, equivalenced to atom 2\n",
            "Atom number 4, equivalenced to atom 2\n",
            "Atom number 5, equivalenced to atom 2\n",
        ]

        self.respin.ivary.describe(file=output)
        output.seek(0)
        self.assertListEqual(expected_lines, output.readlines())

    def test_ivary_description_with_molecule(self) -> None:

        output = StringIO()
        expected_lines = [
            "Atom (C) number 1\n",
            "Atom (H) number 2\n",
            "Atom (H) number 3, equivalenced to atom 2\n",
            "Atom (H) number 4, equivalenced to atom 2\n",
            "Atom (H) number 5, equivalenced to atom 2\n",
        ]

        self.respin.ivary.describe(
            molecule=Molecule([Atom(atomic_number) for atomic_number in [6, 1, 1, 1, 1]]),
            file=output
        )
        output.seek(0)
        self.assertListEqual(expected_lines, output.readlines())

    def test_ivary_from_equivalence(self) -> None:
        self.assertListEqual(
            Respin.Ivary([0, 0, 1, 2, 2]).values,
            Respin.Ivary.from_equivalence(
                Equivalence([None, None, 0, 1, 1])
            ).values
        )


class TestRespinAndEquivalence(TestCase):

    def test_getting_equivalence(self) -> None:
        # MeSO4 [S O O (bridging) O C H H H]
        ivary1 = Respin.Ivary([0, 0, 0, 2, 2, 0, 0, 0, 0])
        ivary2 = Respin.Ivary([-1, -1, -1, -1, -1, 0, 0, 7, 7])
        equivalence = get_equivalence(ivary1, ivary2)
        expected = Equivalence([None, None, None, 1, 1, None, None, 6, 6])
        self.assertListEqual(equivalence.values, expected.values)