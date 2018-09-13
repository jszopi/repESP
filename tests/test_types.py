from repESP.charges import *
from repESP.types import *

from my_unittest import TestCase


class TestMoleculeWithCharges(TestCase):

    def setUp(self) -> None:
        self.molecule = Molecule(
            atoms=[
                Atom(1, make_coords(1, 1, 2)),
                Atom(2, make_coords(2, 0, 2)),
            ]
        )

    def test_construction(self) -> None:
        values = [make_charge(x) for x in [0, 1]]
        MoleculeWithCharges(molecule=self.molecule, charges=values)

    def test_construction_fails_for_mismatched_data(self) -> None:
        values = [make_charge(x) for x in [0, 1, 2, 4]]
        with self.assertRaises(InputFormatError):
            MoleculeWithCharges(molecule=self.molecule, charges=values)


class TestNonGridMesh(TestCase):

    def setUp(self) -> None:

        self.mesh = NonGridMesh([
            make_coords(1, 1, 1),
            make_coords(-1, 0, -0.9)
        ])

    def test_points(self) -> None:

        points = self.mesh.points()
        self.assertListsAlmostEqual(next(points), make_coords(1, 1, 1))
        self.assertListsAlmostEqual(next(points), make_coords(-1, 0, -0.9))

    def test_calc_field(self) -> None:
        # For testing convenience, the field function is a sum of coordinates.
        func = lambda coords: sum(coords)
        self.assertListsAlmostEqual(
            self.mesh.calc_field(func).values,
            [3, -1.9]
        )


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
            make_coords(0.1, 0.2, 0.3),
            make_coords(0.1, 0.2, 0.7),
            make_coords(0.1, 0.2, 1.1),
            make_coords(0.1, 0.5, 0.3),
            make_coords(0.1, 0.5, 0.7),
            make_coords(0.1, 0.5, 1.1),
            make_coords(0.1, 0.8, 0.3),
            make_coords(0.1, 0.8, 0.7),
            make_coords(0.1, 0.8, 1.1),
            make_coords(0.3, 0.2, 0.3),
            make_coords(0.3, 0.2, 0.7),
            make_coords(0.3, 0.2, 1.1),
            make_coords(0.3, 0.5, 0.3),
            make_coords(0.3, 0.5, 0.7),
            make_coords(0.3, 0.5, 1.1),
            make_coords(0.3, 0.8, 0.3),
            make_coords(0.3, 0.8, 0.7),
            make_coords(0.3, 0.8, 1.1),
            make_coords(0.5, 0.2, 0.3),
            make_coords(0.5, 0.2, 0.7),
            make_coords(0.5, 0.2, 1.1),
            make_coords(0.5, 0.5, 0.3),
            make_coords(0.5, 0.5, 0.7),
            make_coords(0.5, 0.5, 1.1),
            make_coords(0.5, 0.8, 0.3),
            make_coords(0.5, 0.8, 0.7),
            make_coords(0.5, 0.8, 1.1),
        ]

        self.assertListsAlmostEqualRecursive(list(self.mesh.points()), points)

    def test_calc_field(self) -> None:
        # For testing convenience, the field function is a sum of coordinates.
        func = lambda coords: sum(coords)
        field = self.mesh.calc_field(func)

        for coords, value in zip(self.mesh.points(), field.values):
            self.assertAlmostEqual(sum(coords), value)


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
