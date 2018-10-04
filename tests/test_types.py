from repESP.charges import *
from repESP.types import *

from my_unittest import TestCase

from copy import copy


class TestAtom(TestCase):

    def setUp(self) -> None:
        self.identity = 1
        self.coords = make_coords(1, 1, 2)
        self.charge = make_charge(0.6)

    def test_ctor(self) -> None:
        atom = Atom(self.identity)
        self.assertEqual(atom.identity, self.identity)

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
        atom = AtomWithCoords(self.identity, self.coords)
        self.assertEqual(atom.identity, self.identity)
        self.assertEqual(atom.coords, self.coords)

    def test_with_charges(self) -> None:
        atom = AtomWithCharge(self.identity, self.charge)
        self.assertEqual(atom.identity, self.identity)
        self.assertEqual(atom.charge, self.charge)

    def test_with_coords_and_charges(self) -> None:
        atom = AtomWithCoordsAndCharge(self.identity, self.coords, self.charge)
        self.assertEqual(atom.identity, self.identity)
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


class TestNonGridMesh(TestCase):

    def setUp(self) -> None:

        self.mesh = NonGridMesh([
            make_coords(1, 1, 1),
            make_coords(-1, 0, -0.9)
        ])

    def test_points(self) -> None:

        points = self.mesh.points
        self.assertListsAlmostEqual(next(points), make_coords(1, 1, 1))
        self.assertListsAlmostEqual(next(points), make_coords(-1, 0, -0.9))


class TestGridMesh(TestCase):

    def setUp(self) -> None:
        self.origin = make_coords(0.1, 0.2, 0.3)
        self.axes = GridMesh.Axes((
            GridMesh.Axis(
                vector=make_coords(0.2, 0, 0),
                point_count=3
            ),
            GridMesh.Axis(
                vector=make_coords(0, 0.3, 0),
                point_count=3
            ),
            GridMesh.Axis(
                vector=make_coords(0, 0, 0.4),
                point_count=3
            ),
        ))

        self.mesh = GridMesh(origin=self.origin, axes=self.axes)

    def test_construction_fails_with_misaligned_axes(self) -> None:
        axes = GridMesh.Axes((
            GridMesh.Axis(
                vector=make_coords(0.2, 0, 0.0001),
                point_count=3
            ),
            self.axes[1],
            self.axes[2]
        ))

        with self.assertRaises(NotImplementedError):
            GridMesh(origin=self.origin, axes=axes)

    def test_points(self) -> None:
        points = [
            make_coords(*coords) for coords in [
                (0.1, 0.2, 0.3),
                (0.1, 0.2, 0.7),
                (0.1, 0.2, 1.1),
                (0.1, 0.5, 0.3),
                (0.1, 0.5, 0.7),
                (0.1, 0.5, 1.1),
                (0.1, 0.8, 0.3),
                (0.1, 0.8, 0.7),
                (0.1, 0.8, 1.1),
                (0.3, 0.2, 0.3),
                (0.3, 0.2, 0.7),
                (0.3, 0.2, 1.1),
                (0.3, 0.5, 0.3),
                (0.3, 0.5, 0.7),
                (0.3, 0.5, 1.1),
                (0.3, 0.8, 0.3),
                (0.3, 0.8, 0.7),
                (0.3, 0.8, 1.1),
                (0.5, 0.2, 0.3),
                (0.5, 0.2, 0.7),
                (0.5, 0.2, 1.1),
                (0.5, 0.5, 0.3),
                (0.5, 0.5, 0.7),
                (0.5, 0.5, 1.1),
                (0.5, 0.8, 0.3),
                (0.5, 0.8, 0.7),
                (0.5, 0.8, 1.1),
            ]
        ]

        self.assertAlmostEqualRecursive(list(self.mesh.points), points)


class TestField(TestCase):

    def setUp(self) -> None:

        self.mesh = NonGridMesh([
            make_coords(1, 1, 1),
            make_coords(-1, 0, -0.9)
        ])

        self.values = [make_esp(0.5), make_esp(-0.7)]

    def test_construction(self) -> None:
        Field(self.mesh, self.values)

    def test_construction_fails_when_lengths_mismatched(self) -> None:
        with self.assertRaises(InputFormatError):
            Field(self.mesh, [make_esp(0.5)])

    def test_addition_fails_for_different_meshes(self) -> None:

        field1 = Field(self.mesh, self.values)
        field2 = Field(
            NonGridMesh([
                make_coords(1, 1, 1),
                make_coords(-1, 0, 0.9)
            ]),
            self.values
        )

        with self.assertRaises(ValueError):
            field3 = field1 + field2

    def test_addition(self) -> None:

        field1 = Field(self.mesh, self.values)
        field2 = Field(self.mesh, [make_esp(0.1), make_esp(1) ])
        field3 = field1 + field2

        self.assertEqual(field1.mesh, field3.mesh)
        self.assertListsAlmostEqual(field3.values, [make_esp(0.6), make_esp(0.3)])

    def test_subtraction(self) -> None:

        field1 = Field(self.mesh, self.values)
        field2 = Field(self.mesh, [make_esp(0.1), make_esp(1) ])
        field3 = field1 - field2

        self.assertEqual(field1.mesh, field3.mesh)
        self.assertListsAlmostEqual(field3.values, [make_esp(0.4), make_esp(-1.7)])
