from repESP.calc_fields import *
from repESP.types import *

from my_unittest import TestCase, makeCoords

nonGridMesh = NonGridMesh([
    makeCoords(1, 1, 1),
    makeCoords(-1, 0, -0.9)
])

gridMesh = GridMesh(
    origin=makeCoords(0.1, 0.2, 0.3),
    axes=(
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
)

molecule = Molecule(
    atoms=[
        Atom(identity=1, coords=makeCoords(0, 1, 0.5)),
        Atom(identity=1, coords=makeCoords(-0.4, 0.2, 0.5))
    ]
)

class TestEspFromCharges(TestCase):

    def setUp(self) -> None:
        self.charges = Charges(molecule, [0.5, -0.9])

    def test_non_grid_esp(self) -> None:

        result = esp_from_charges(nonGridMesh, self.charges)
        expected = Field(
            nonGridMesh,
            [
                Esp(-0.08590039),
                Esp(-0.33459064)
            ]
        )

        self.assertEqual(result.mesh, expected.mesh)
        self.assertListsAlmostEqual(result.values, expected.values)

    def test_grid_esp(self) -> None:

        expected_floats = [
            -1.06932877, -1.06932877, -0.65481332, -0.54712186, -0.54712186,
            -0.44070511,  0.55035405,  0.55035405, -0.13294273, -0.66644219,
            -0.66644219, -0.49727391, -0.33189403, -0.33189403, -0.33066481,
             0.25868003,  0.25868003, -0.10389610, -0.45771121, -0.45771121,
            -0.38483669, -0.24786530, -0.24786530, -0.26261985,  0.05220646,
             0.05220646, -0.10743320
        ]

        expected = Field(gridMesh, [Esp(x) for x in expected_floats])
        result = esp_from_charges(gridMesh, self.charges)

        self.assertEqual(result.mesh, expected.mesh)
        self.assertListsAlmostEqual(result.values, expected.values)


class TestVoronoi(TestCase):

    def test_non_grid_esp(self) -> None:

        result = voronoi(nonGridMesh, molecule)
        expected = Field(
            nonGridMesh,
            [
                (0, Dist(1.11803398)),
                (1, Dist(1.53622914))
            ]
        )

        self.assertEqual(result.mesh, expected.mesh)
        self.assertListsAlmostEqualRecursive(result.values, expected.values)

    def test_grid_esp(self) -> None:

        expected_floats = [
            0.53851648, 0.53851648, 0.78102496, 0.54772255, 0.54772255,
            0.78740078, 0.30000000, 0.30000000, 0.64031242, 0.72801098,
            0.72801098, 0.92195444, 0.61644140, 0.61644140, 0.83666002,
            0.41231056, 0.41231056, 0.70000000, 0.92195444, 0.92195444,
            1.08166538, 0.73484692, 0.73484692, 0.92736184, 0.57445626,
            0.57445626, 0.80622577
        ]

        expected_closest_atom = [
            1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1,
            1, 0, 0, 0, 0, 0, 0
        ]

        expected = Field(
            gridMesh,
            list(zip(
                expected_closest_atom,
                [Dist(x) for x in expected_floats]
            ))
        )

        result = voronoi(gridMesh, molecule)

        self.assertEqual(result.mesh, expected.mesh)
        self.assertListsAlmostEqualRecursive(result.values, expected.values)
