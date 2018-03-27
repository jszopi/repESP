import unittest


class TestCase(unittest.TestCase):
    '''Extend TestCase to include own assertion methods'''

    def assertListsEqualWhenSorted(self, first, second):

        self.assertEqual(len(first), len(second))
        for a, b in zip(sorted(first), sorted(second)):
            self.assertEqual(a, b)
