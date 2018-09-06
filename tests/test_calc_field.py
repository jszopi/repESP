from repESP.calc_fields import *
from repESP.types import *

from my_unittest import TestCase, makeCoords


class TestEspFromCharges(TestCase):

    def setUp(self) -> None:

        molecule = Molecule(
            atoms=[
                Atom(identity=1, coords=makeCoords(0, 1, -2)),
                Atom(identity=1, coords=makeCoords(-1, 0.2, 0.5))
            ]
        )

        self.charges = Charges(molecule, [0.5, -0.9])

        self.nonGridMesh = NonGridMesh([
            makeCoords(1, 1, 1),
            makeCoords(-1, 0, -0.9)
        ])

    def test_non_grid_esp(self) -> None:

        result = esp_from_charges(self.nonGridMesh, self.charges)
        expected = Field(
            self.nonGridMesh,
            [
                Esp(-0.248880185),
                Esp(-0.357323317)
            ]
        )

        self.assertEqual(result.mesh, expected.mesh)
        self.assertListsAlmostEqual(result.values, expected.values)
