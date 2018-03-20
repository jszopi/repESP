from repESP.types import *

import unittest


class TestCoords(unittest.TestCase):

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
