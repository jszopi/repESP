import my_unittest
from repESP import cube_helpers

# I'm not sure whether this should be here but it seems fine, it's only
# evaluated once and will be available to all the TestCases in this module
cube = cube_helpers.Cube("tests/test_mol_den.cub")


class TestCube(my_unittest.TestCase):

    def test_title(self):
        self.assertEqual(cube.title,
                         " Electron density from Total MP2 Density")

    def test_atom_count(self):
        self.assertEqual(cube.atom_count, 1)

    def test_type(self):
        self.assertEqual(cube.cube_type, 'ed')
