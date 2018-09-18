from repESP.calc_fields import esp_from_charges, voronoi, calc_rms_error, calc_relative_rms_error
from repESP.charges import *
from repESP.types import *
from repESP.gaussian_util import get_charges_from_log, get_esp_fit_stats_from_log
from repESP.esp_util import parse_gaussian_esp

from my_unittest import TestCase


class SmallTestCase(TestCase):

    def setUp(self) -> None:

        self.nonGridMesh = NonGridMesh([
            make_coords(1, 1, 1),
            make_coords(-1, 0, -0.9)
        ])

        self.gridMesh = GridMesh(
            origin=make_coords(0.1, 0.2, 0.3),
            axes=GridMesh.Axes((
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
        )

        self.molecule = Molecule(
            atoms=[
                Atom(identity=1, coords=make_coords(0, 1, 0.5)),
                Atom(identity=1, coords=make_coords(-0.4, 0.2, 0.5))
            ]
        )


class TestEspFromCharges(SmallTestCase):

    def setUp(self) -> None:
        super().setUp()
        self.molecule_with_charges = MoleculeWithCharges(
            self.molecule,
            [make_charge(x) for x in [0.5, -0.9]]
        )

    def test_non_grid_esp(self) -> None:

        result = esp_from_charges(self.nonGridMesh, self.molecule_with_charges)
        expected = Field(
            self.nonGridMesh,
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

        result = voronoi(self.nonGridMesh, self.molecule)
        expected = Field(
            self.nonGridMesh,
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
            self.gridMesh,
            list(zip(
                expected_closest_atom,
                [Dist(x) for x in expected_floats]
            ))
        )

        result = voronoi(self.gridMesh, self.molecule)

        self.assertEqual(result.mesh, expected.mesh)
        self.assertListsAlmostEqualRecursive(result.values, expected.values)


class TestCalcStats(TestCase):

    def setUp(self) -> None:

        with open("data/methane/methane_mk.esp") as f:
            gaussian_esp = parse_gaussian_esp(f)

        self.esp_field = gaussian_esp.field
        molecule = gaussian_esp.molecule_with_charges.molecule

        with open("data/methane/methane_mk.log") as f:
            self.charges = get_charges_from_log(f, ChargeType.MK, verify_against=molecule)
            f.seek(0)
            self.esp_fit_stats = get_esp_fit_stats_from_log(f, ChargeType.MK, verify_against=molecule)

        self.molecule_with_charges = MoleculeWithCharges(molecule, self.charges)
        self.rep_esp_field = esp_from_charges(self.esp_field.mesh, self.molecule_with_charges)

    def test_rms(self) -> None:
        self.assertAlmostEqual(
            calc_rms_error(self.esp_field, self.rep_esp_field),
            self.esp_fit_stats[0],
            places=5  # Gaussian log output precision

        )

    def test_rrms(self) -> None:
        self.assertAlmostEqual(
            calc_relative_rms_error(self.esp_field, self.rep_esp_field),
            self.esp_fit_stats[1],
            places=5  # Gaussian log output precision
        )
