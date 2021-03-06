from repESP.calc_fields import esp_from_charges, voronoi, calc_rms_error, calc_relative_rms_error
from repESP.charges import *
from repESP.esp_util import parse_gaussian_esp
from repESP.fields import *
from repESP.types import *
from repESP.gaussian_format import get_charges_from_log, MkChargeSectionParser

from my_unittest import TestCase


class SmallTestCase(TestCase):

    def setUp(self) -> None:

        self.mesh = Mesh([
            Coords((1, 1, 1)),
            Coords((-1, 0, -0.9))
        ])

        self.gridMesh = GridMesh(
            origin=Coords((0.1, 0.2, 0.3)),
            axes=GridMesh.Axes((
                GridMesh.Axis(
                    vector=Coords((0.2, 0, 0)),
                    point_count=3
                ),
                GridMesh.Axis(
                    vector=Coords((0, 0.3, 0)),
                    point_count=3
                ),
                GridMesh.Axis(
                    vector=Coords((0, 0, 0.4)),
                    point_count=3
                ),
            ))
        )

        self.molecule = Molecule(
            atoms=[
                AtomWithCoords(atomic_number=1, coords=Coords((0, 1, 0.5))),
                AtomWithCoords(atomic_number=1, coords=Coords((-0.4, 0.2, 0.5)))
            ]
        )


class TestEspFromCharges(SmallTestCase):

    def setUp(self) -> None:
        super().setUp()
        self.molecule_with_charges = Molecule([
            AtomWithCoordsAndCharge(
                atom.atomic_number,
                atom.coords,
                charge
            )
            for atom, charge in zip(
                self.molecule.atoms,
                [Charge(x) for x in [0.5, -0.9]]
            )
        ])

    def test_non_grid_esp(self) -> None:

        result = esp_from_charges(self.mesh, self.molecule_with_charges)
        expected = Field(
            self.mesh,
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

        expected = Field(self.gridMesh, [Esp(x) for x in expected_floats])
        result = esp_from_charges(self.gridMesh, self.molecule_with_charges)

        self.assertEqual(result.mesh, expected.mesh)
        self.assertListsAlmostEqual(result.values, expected.values)


class TestVoronoi(SmallTestCase):

    def setUp(self) -> None:
        super().setUp()

    def test_non_grid_esp(self) -> None:

        result = voronoi(self.mesh, self.molecule)
        expected = Field(
            self.mesh,
            [
                (0, Dist(1.11803398)),
                (1, Dist(1.53622914))
            ]
        )

        self.assertAlmostEqualRecursive(expected, result)

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
            self.gridMesh,
            list(zip(
                expected_closest_atom,
                [Dist(x) for x in expected_floats]
            ))
        )

        result = voronoi(self.gridMesh, self.molecule)

        self.assertAlmostEqualRecursive(expected, result)


class TestCalcStats(TestCase):

    def setUp(self) -> None:

        with open("data/methane/methane_mk.esp") as f:
            gaussian_esp = parse_gaussian_esp(f)

        self.esp_values = gaussian_esp.field.values
        molecule = gaussian_esp.molecule

        with open("data/methane/methane_mk.log") as f:
            charges = get_charges_from_log(f, MkChargeSectionParser(), verify_against=molecule)

        molecule_with_charges = Molecule([
            AtomWithCoordsAndCharge(
                atom.atomic_number,
                atom.coords,
                charge
            )
            for atom, charge in zip(molecule.atoms, charges)
        ])
        self.rep_esp_values = esp_from_charges(gaussian_esp.field.mesh, molecule_with_charges).values

    def test_rms(self) -> None:
        self.assertAlmostEqual(
            calc_rms_error(self.esp_values, self.rep_esp_values),
            0.00069,  # from Gaussian log
            places=5

        )

    def test_rrms(self) -> None:
        self.assertAlmostEqual(
            calc_relative_rms_error(self.esp_values, self.rep_esp_values),
            0.35027,  # from Gaussian log
            places=5
        )
