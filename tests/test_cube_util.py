from repESP.types import *
from repESP.cube_util import parse_ed_cube, write_cube

from io import StringIO
from my_unittest import TestCase, make_coords

class TestCubeParser(TestCase):

    def setUp(self) -> None:
        with open("tests/test_mol_den.cub", 'r') as f:
            self.cube = parse_ed_cube(f)

    def test_cube_data(self) -> None:
        self.assertEqual(
            self.cube.cube_info.input_line,
            " Test molecule density=mp2"
        )
        self.assertEqual(
            self.cube.cube_info.title_line,
            " Electron density from Total MP2 Density"
        )

    def test_grid(self) -> None:

        self.assertListsAlmostEqual(
            [0.1, 0.2, 0.3],
            self.cube.field.mesh._origin  # type: ignore # (accessing subclass attribute)
        )

        self.assertListsAlmostEqual(
            [0.1, 0.2, 0.3],
            self.cube.field.mesh._origin  # type: ignore # (accessing subclass attribute)
        )

        expected_axis_vecs = [
            make_coords(0.2, 0, 0),
            make_coords(0, 0.3, 0),
            make_coords(0, 0, 0.4)
        ]

        for axis, expected_axis_vec in zip(
            self.cube.field.mesh._axes,  # type: ignore # (accessing subclass attribute)
            expected_axis_vecs
        ):
            self.assertEqual(axis.point_count, 3)
            self.assertListsAlmostEqual(axis.vector, expected_axis_vec)

    def test_molecule(self) -> None:

        self.assertEqual(len(self.cube.molecule.atoms), 1)
        self.assertEqual(self.cube.molecule.atoms[0].identity, 1)
        self.assertEqual(self.cube.molecule.atoms[0].coords, make_coords(0.1, 0.2, 0.4))

    def test_electron_counts(self) -> None:
        self.assertListEqual(self.cube.electrons_on_atoms.values, [0.9])

    def test_field_values(self) -> None:

        values = [
            1.36229e-08, 2.28426e-08, 3.69340e-08, 5.75781e-08, 8.65339e-08,
            1.25366e-07, 1.75080e-07, 2.35719e-07, 3.06022e-07, 3.83225e-07,
            4.63120e-07, 5.40378e-07, 6.09116e-07, 6.63619e-07, 6.99085e-07,
            7.12274e-07, 7.01958e-07, 6.69101e-07, 6.16729e-07, 5.49515e-07,
            4.73127e-07, 3.93467e-07, 3.15935e-07, 2.44838e-07, 1.83067e-07,
            1.32029e-07, 9.18219e-08
        ]

        ed_values = [Ed(value) for value in values]

        self.assertEqual(
            ed_values,
            self.cube.field.values
        )

class TestCubeWriter(TestCase):

    def setUp(self) -> None:
        with open("tests/test_mol_den.cub", 'r') as f:
            self.cube = parse_ed_cube(f)
            f.seek(0)
            self.input = f.readlines()

    def test_compare_input_to_written_cube(self) -> None:
        stringIO = StringIO()
        write_cube(stringIO, self.cube)
        stringIO.seek(0)
        output = stringIO.readlines()
        self.assertListEqual(self.input, output)
