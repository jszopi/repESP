import my_unittest
from repESP import cube_helpers

# I'm not sure whether this should be here but it seems fine, it's only
# evaluated once and will be available to all the TestCases in this module
cube = cube_helpers.Cube("tests/test_mol_den.cub")
atom = cube.molecule[0]
grid = cube.field.grid

# This is getting complicated so these initialization lines should be a
# separate test. That likely requires a custom test suite so that they're only
# executed once and before the others.
dist_field = cube.molecule.calc_field(grid, 'dist')
repESP_field = cube.molecule.calc_field(grid, 'rep_esp', ['cube'])[0]
# Values calculated with an independent ad-hoc script
dist_result = [0.1, 0.3, 0.7, 0.316227766, 0.424264068, 0.761577310,
               0.608276253, 0.670820393, 0.921954445, 0.223606797, 0.360555127,
               0.728010988, 0.374165738, 0.469041575, 0.787400787, 0.640312423,
               0.7, 0.943398113, 0.412310562, 0.5, 0.806225774, 0.509901951,
               0.583095189, 0.860232526, 0.728010988, 0.781024967, 1.004987562]

# Isovalue selected so that only two points fall above it
ed_isoval = 7.01e-07
edt_field = cube.field.distance_transform(ed_isoval)
# This was calculated with the aforementioned script as the smaller elements of
# the distance matrices from those two points.
edt_result = [0.63245553, 0.63245553, 0.74833147, 0.36055512, 0.36055512,
              0.53851648, 0.2, 0.2, 0.44721359, 0.6, 0.6, 0.72111025, 0.3, 0.3,
              0.5, 0.0, 0.0, 0.4, 0.63245553, 0.63245553, 0.74833147,
              0.36055512, 0.36055512, 0.53851648, 0.2, 0.2, 0.44721359]


class TestCube(my_unittest.TestCase):

    def test_title(self):
        self.assertEqual(cube.title,
                         " Electron density from Total MP2 Density")

    def test_atom_count(self):
        self.assertEqual(cube.atom_count, 1)

    def test_type(self):
        self.assertEqual(cube.cube_type, 'ed')


class TestGrid(my_unittest.TestCase):

    def test_origin_coords(self):
        self.assertListAlmostEqual(grid.origin_coords, [0.1, 0.2, 0.3])

    def test_dir_intervals(self):
        self.assertListAlmostEqual(grid.dir_intervals, [0.2, 0.3, 0.4])

    def test_aligned(self):
        self.assertTrue(grid.aligned_to_coord)

    def test_points_on_axes(self):
        self.assertListAlmostEqual(grid.points_on_axes, [3, 3, 3])


class TestMolecule(my_unittest.TestCase):

    def test_label(self):
        self.assertEqual(atom.label, 1)

    def test_atomic_no(self):
        self.assertEqual(atom.atomic_no, 1)

    def test_coords(self):
        self.assertListAlmostEqual(atom.coords, [0.1, 0.2, 0.4])

    def test_cube_charges(self):
        self.assertAlmostEqual(atom.charges['cube'], 0.9)

    def test_atom_count(self):
        self.assertEqual(len(cube.molecule), 1)


class TestDistField(my_unittest.TestCase):

    def test_min_atom(self):
        self.assertListEqual(list(dist_field[0].values.flatten()), [1]*27)

    def test_min_dist(self):
        self.assertListAlmostEqual(list(dist_field[1].values.flatten()),
                                   dist_result)

    def test_field_type(self):
        self.assertEqual(dist_field[0].field_type, 'parent_atom')
        self.assertEqual(dist_field[1].field_type, 'dist')

    def test_atom_field_info(self):
        self.assertListEqual(dist_field[0].field_info, ['Voronoi', 'own'])


class TestRepESPField(my_unittest.TestCase):

    def test_values(self):
        self.assertListAlmostEqual(list(repESP_field.values.flatten()),
                                   [0.9/dist for dist in dist_result])

    def test_field_type(self):
        self.assertEqual(repESP_field.field_type, 'rep_esp')

    def test_atom_field_info(self):
        self.assertEqual(repESP_field.field_info, 'cube')


class TestField(my_unittest.TestCase):

    def test_edt_type(self):
        self.assertEqual(edt_field.field_type, 'dist')

    def test_edt_info(self):
        self.assertEqual(edt_field.field_info[0], 'ed')
        self.assertAlmostEqual(edt_field.field_info[1], ed_isoval)

    def test_edt_values(self):
        self.assertListAlmostEqual(list(edt_field.values.flatten()),
                                   edt_result)
