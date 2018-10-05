from repESP.charges import *
from repESP.types import *

from my_unittest import TestCase

from copy import copy


class TestAtom(TestCase):

    def setUp(self) -> None:
        self.atomic_number = 1
        self.coords = make_coords(1, 1, 2)
        self.charge = make_charge(0.6)

    def test_ctor(self) -> None:
        atom = Atom(self.atomic_number)
        self.assertEqual(atom.atomic_number, self.atomic_number)

    def test_ctor_fails_with_invalid_atomic_number(self) -> None:
        with self.assertRaises(ValueError):
            Atom(-1)
        with self.assertRaises(ValueError):
            Atom(300)

        # Checking that the check in __post_init__ is inherited correctly:
        with self.assertRaises(ValueError):
            AtomWithCoords(-1, self.coords)
        with self.assertRaises(ValueError):
            AtomWithCharge(-1, self.charge)
        with self.assertRaises(ValueError):
            AtomWithCoordsAndCharge(-1, self.coords, self.charge)

    def test_with_coords(self) -> None:
        atom = AtomWithCoords(self.atomic_number, self.coords)
        self.assertEqual(atom.atomic_number, self.atomic_number)
        self.assertEqual(atom.coords, self.coords)

    def test_with_charges(self) -> None:
        atom = AtomWithCharge(self.atomic_number, self.charge)
        self.assertEqual(atom.atomic_number, self.atomic_number)
        self.assertEqual(atom.charge, self.charge)

    def test_with_coords_and_charges(self) -> None:
        atom = AtomWithCoordsAndCharge(self.atomic_number, self.coords, self.charge)
        self.assertEqual(atom.atomic_number, self.atomic_number)
        self.assertEqual(atom.coords, self.coords)
        self.assertEqual(atom.charge, self.charge)


class TestMolecule(TestCase):

    def test_ctors(self) -> None:

        Molecule([Atom(1), Atom(3)])

        Molecule(
            [
                AtomWithCoords(1, make_coords(1, 1, 2)),
                AtomWithCoords(2, make_coords(2, 0, 2)),
            ]
        )

        Molecule(
            [
                AtomWithCharge(1, make_charge(0.6)),
                AtomWithCharge(2, make_charge(-0.2)),
            ]
        )

        Molecule(
            [
                AtomWithCoordsAndCharge(1, make_charge(0.6), make_coords(1, 1, 2)),
                AtomWithCoordsAndCharge(2, make_charge(-0.2), make_coords(2, 0, 2)),
            ]
        )
