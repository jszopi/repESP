from repESP.charge_util import average
from repESP.charges import Charge
from repESP.respin_util import Equivalence

from my_unittest import TestCase


class TestAveraging(TestCase):

    def setUp(self) -> None:
        self.charges = [Charge(x) for x in [1, 2, 3, 4, 5, 6]]

    def test_no_equivalence(self) -> None:
        equivalence = Equivalence([None]*6)
        result = average(self.charges, equivalence)
        self.assertListsAlmostEqual(
            self.charges,
            result
        )

    def test_basic_equivalence(self) -> None:
        equivalence = Equivalence([None, 0] + [None]*4)
        result = average(self.charges, equivalence)
        self.assertListsAlmostEqual(
            [1.5, 1.5, 3, 4, 5, 6],
            result
        )

    def test_typical_euqivalence(self) -> None:
        equivalence = Equivalence([None, None, 0, 1, 0, 1])
        result = average(self.charges, equivalence)
        self.assertListsAlmostEqual(
            [3, 4, 3, 4, 3, 4],
            result
        )

    def test_cyclic_euqivalence(self) -> None:
        # It may be later decided that this should fail Equivalence validation
        # but until it does, functions in the library should cope with it.
        equivalence = Equivalence([5, 0, 1, 2, 3, 4])
        result = average(self.charges, equivalence)
        self.assertListsAlmostEqual(
            [3.5]*6,
            result
        )
