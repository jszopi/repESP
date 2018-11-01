from repESP.equivalence import Equivalence
from repESP.types import *

from my_unittest import TestCase


class TestEquivalence(TestCase):

    def test_init_validates_values(self) -> None:

        Equivalence([None, None, 0, 1, 1])

        with self.assertRaises(ValueError):
            Equivalence([0, 0, 2, -1, 2])

        with self.assertRaises(ValueError):
            Equivalence([5, 0, 2, 2, 2])

    def test_description(self) -> None:

        equivalence = Equivalence([None, None, 0, 1, 1])

        expected_lines = [
            "Atom number 1",
            "Atom number 2",
            "Atom number 3, equivalenced to atom 1",
            "Atom number 4, equivalenced to atom 2",
            "Atom number 5, equivalenced to atom 2",
        ]

        result = equivalence.describe()
        self.assertListEqual(expected_lines, result.splitlines())

    def test_description_with_molecule(self) -> None:

        equivalence = Equivalence([None, None, 0, 1, 1])

        expected_lines = [
            "Atom (C) number 1",
            "Atom (H) number 2",
            "Atom (H) number 3, equivalenced to atom 1",
            "Atom (H) number 4, equivalenced to atom 2",
            "Atom (H) number 5, equivalenced to atom 2",
        ]

        result = equivalence.describe(
            molecule=Molecule([Atom(atomic_number) for atomic_number in [6, 1, 1, 1, 1]]),
        )
        self.assertListEqual(expected_lines, result.splitlines())
