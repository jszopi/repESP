from repESP.types import *

from my_unittest import TestCase


class TestCoords(TestCase):

    def test_equality(self):
        a = Coords(1, 2, 3)
        b = Coords(1, 2, 3)
        self.assertEqual(a, b)

    def test_inequality(self):
        a = Coords(1, 2, 3)
        b = Coords(10, 2, 3)
        self.assertNotEqual(a, b)

    def test_approximate_equality(self):
        a = Coords(1.0000001, 2.0000001, 3.0000001)
        b = Coords(1.0000002, 2.0000002, 3.0000002)
        self.assertEqual(a, b)


class TestNonGridMesh(TestCase):

    def setUp(self):

        self.mesh = NonGridMesh([
            Coords(1, 1, 1),
            Coords(-1, 0, -0.9)
        ])

    def test_points(self):

        points = self.mesh.points()
        self.assertEqual(next(points), Coords(1, 1, 1))
        self.assertEqual(next(points), Coords(-1, 0, -0.9))

    def test_calc_field(self):
        func = lambda coords: sum(coords)
        self.assertEqual(
            self.mesh.calc_field(func).values,
            [3, -1.9]
        )


class TestGridMesh(TestCase):

    def setUp(self):
        origin = Coords(0.1, 0.2, 0.3)
        axes = (
            GridMeshAxis(
                vector=Coords(0.2, 0, 0),
                point_count=3
            ),
            GridMeshAxis(
                vector=Coords(0, 0.3, 0),
                point_count=3
            ),
            GridMeshAxis(
                vector=Coords(0, 0, 0.4),
                point_count=3
            ),
        )

        self.mesh = GridMesh(origin=origin, axes=axes)

    def test_points(self):
        points = [
            Coords(0.1, 0.2, 0.3),
            Coords(0.1, 0.2, 0.7),
            Coords(0.1, 0.2, 1.1),
            Coords(0.1, 0.5, 0.3),
            Coords(0.1, 0.5, 0.7),
            Coords(0.1, 0.5, 1.1),
            Coords(0.1, 0.8, 0.3),
            Coords(0.1, 0.8, 0.7),
            Coords(0.1, 0.8, 1.1),
            Coords(0.3, 0.2, 0.3),
            Coords(0.3, 0.2, 0.7),
            Coords(0.3, 0.2, 1.1),
            Coords(0.3, 0.5, 0.3),
            Coords(0.3, 0.5, 0.7),
            Coords(0.3, 0.5, 1.1),
            Coords(0.3, 0.8, 0.3),
            Coords(0.3, 0.8, 0.7),
            Coords(0.3, 0.8, 1.1),
            Coords(0.5, 0.2, 0.3),
            Coords(0.5, 0.2, 0.7),
            Coords(0.5, 0.2, 1.1),
            Coords(0.5, 0.5, 0.3),
            Coords(0.5, 0.5, 0.7),
            Coords(0.5, 0.5, 1.1),
            Coords(0.5, 0.8, 0.3),
            Coords(0.5, 0.8, 0.7),
            Coords(0.5, 0.8, 1.1),
        ]

        self.assertListsEqualWhenSorted(list(self.mesh.points()), points)

    def test_calc_field(self):
        func = lambda coords: sum(coords)
        field = self.mesh.calc_field(func)

        for coords, value in zip(self.mesh.points(), field.values):
            self.assertEqual(sum(coords), value)
