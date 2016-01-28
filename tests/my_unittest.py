import unittest


class TestCase(unittest.TestCase):
    '''Extend TestCase to include own assertion methods'''

    def assertListAlmostEqual(self, first, second, places=None, msg=None,
                              delta=None):
        '''Element-wise float collection comparison

        http://stackoverflow.com/a/8312110
        '''
        self.assertEqual(len(first), len(second))
        for a, b in zip(first, second):
            self.assertAlmostEqual(a, b, places, msg, delta)
