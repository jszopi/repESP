from repESP.types import *

from my_unittest import TestCase, makeCoords


class TestNonGridMesh(TestCase):

    def setUp(self):

        self.mesh = NonGridMesh([
            makeCoords(1, 1, 1),
            makeCoords(-1, 0, -0.9)
        ])

    def test_points(self):

        points = self.mesh.points()
        self.assertListsAlmostEqual(next(points), makeCoords(1, 1, 1))
        self.assertListsAlmostEqual(next(points), makeCoords(-1, 0, -0.9))

    def test_calc_field(self):
        func = lambda coords: sum(coords)
        self.assertAlmostEqual(
            self.mesh.calc_field(func).values,
            [3, -1.9]
        )


class TestGridMesh(TestCase):

    def setUp(self):
        origin = makeCoords(0.1, 0.2, 0.3)
        axes = (
            GridMesh.GridMeshAxis(
                vector=makeCoords(0.2, 0, 0),
                point_count=3
            ),
            GridMesh.GridMeshAxis(
                vector=makeCoords(0, 0.3, 0),
                point_count=3
            ),
            GridMesh.GridMeshAxis(
                vector=makeCoords(0, 0, 0.4),
                point_count=3
            ),
        )

        self.mesh = GridMesh(origin=origin, axes=axes)

    def test_points(self):
        points = [
            makeCoords(0.1, 0.2, 0.3),
            makeCoords(0.1, 0.2, 0.7),
            makeCoords(0.1, 0.2, 1.1),
            makeCoords(0.1, 0.5, 0.3),
            makeCoords(0.1, 0.5, 0.7),
            makeCoords(0.1, 0.5, 1.1),
            makeCoords(0.1, 0.8, 0.3),
            makeCoords(0.1, 0.8, 0.7),
            makeCoords(0.1, 0.8, 1.1),
            makeCoords(0.3, 0.2, 0.3),
            makeCoords(0.3, 0.2, 0.7),
            makeCoords(0.3, 0.2, 1.1),
            makeCoords(0.3, 0.5, 0.3),
            makeCoords(0.3, 0.5, 0.7),
            makeCoords(0.3, 0.5, 1.1),
            makeCoords(0.3, 0.8, 0.3),
            makeCoords(0.3, 0.8, 0.7),
            makeCoords(0.3, 0.8, 1.1),
            makeCoords(0.5, 0.2, 0.3),
            makeCoords(0.5, 0.2, 0.7),
            makeCoords(0.5, 0.2, 1.1),
            makeCoords(0.5, 0.5, 0.3),
            makeCoords(0.5, 0.5, 0.7),
            makeCoords(0.5, 0.5, 1.1),
            makeCoords(0.5, 0.8, 0.3),
            makeCoords(0.5, 0.8, 0.7),
            makeCoords(0.5, 0.8, 1.1),
        ]

        self.assertListsAlmostEqualRecursive(list(self.mesh.points()), points)

    def test_calc_field(self):
        # For testing convenience, the field function is a sum of coordinates.
        func = lambda coords: sum(coords)
        field = self.mesh.calc_field(func)

        for coords, value in zip(self.mesh.points(), field.values):
            self.assertAlmostEqual(sum(coords), value)
