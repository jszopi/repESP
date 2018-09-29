from repESP.util import _list_from_dict, _mask_from_list

from typing import Dict, List, Mapping

from my_unittest import TestCase


class TestListFromDictionary(TestCase):

    def test_example(self) -> None:
        dictionary: Dict[int, bool] = {1: True, 3: False}
        expected = [None, True, None, False]
        result = _list_from_dict(
            dictionary,
            4,
            default=None,
            one_indexed=False
        )

        self.assertListEqual(expected, result)

    def test_one_indexed(self) -> None:
        dictionary: Dict[int, bool] = {1: True, 3: False}
        expected = [True, None, False, None]
        result = _list_from_dict(
            dictionary,
            4,
            default=None,
            one_indexed=True
        )

        self.assertListEqual(expected, result)


class TestMaskFromList(TestCase):

    def test_with_default_values(self) -> None:
        expected = [False, True, False, True]
        result: List[bool] = _mask_from_list([1, 3], 4, one_indexed=False)
        self.assertListEqual(expected, result)

    def test_with_custom_values(self) -> None:
        expected = [0, 1, 0, 1]
        result: List[int]  = _mask_from_list(
            [1, 3],
            4,
            value_if_present=1,
            value_if_absent=0,
            one_indexed=False
        )
        self.assertListEqual(expected, result)

    def test_one_indexed(self) -> None:
        expected = [True, False, True, False]
        result: List[bool] = _mask_from_list([1, 3], 4, one_indexed=True)
        self.assertListEqual(expected, result)
