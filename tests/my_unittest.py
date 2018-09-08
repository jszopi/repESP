from repESP.types import *

from typing import Any, List

import unittest


def make_coords(x: float, y: float, z: float) -> Coords:
    return Coords((Dist(x), Dist(y), Dist(z)))


class TestCase(unittest.TestCase):
    '''Extend TestCase to include own assertion methods'''

    def assertListsEqualWhenSorted(
            self,
            first: List[Any],
            second: List[Any],
            msg=None
    ):

        self.assertEqual(len(first), len(second), msg)
        for a, b in zip(sorted(first), sorted(second)):
            self.assertEqual(a, b, msg)

    def assertListsAlmostEqual(
            self,
            first: Collection[float],
            second: Collection[float],
            places=None,
            msg=None,
            delta=None
    ):
        '''Element-wise float collection comparison

        http://stackoverflow.com/a/8312110
        '''
        self.assertEqual(len(first), len(second))
        for a, b in zip(first, second):
            self.assertAlmostEqual(a, b, places, msg, delta)

    def assertListsAlmostEqualRecursive(
            self,
            first: Collection[Any],
            second: Collection[Any],
            places=None,
            msg=None,
            delta=None
    ):
        self.assertEqual(len(first), len(second))
        for a, b in zip(first, second):
            try:
                self.assertAlmostEqual(a, b, places, msg, delta)
            except TypeError:
                self.assertListsAlmostEqual(a, b, places, msg, delta)
