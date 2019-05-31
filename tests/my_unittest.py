from typing import Any, Callable, Collection, List, Optional, overload

import unittest
import dataclasses


class TestCase(unittest.TestCase):
    '''Extend TestCase to include own assertion methods'''

    def assertListsEqualWhenSorted(
            self,
            first: List[Any],
            second: List[Any],
            msg: Optional[str]=None
    ) -> None:

        self.assertEqual(len(first), len(second), msg)
        for a, b in zip(sorted(first), sorted(second)):
            self.assertEqual(a, b, msg)

    @overload
    def assertListsAlmostEqual(
            self,
            first: Collection[float],
            second: Collection[float],
            places: Optional[int]=None,
            msg: Optional[str]=None
    ) -> None: ...

    @overload
    def assertListsAlmostEqual(
            self,
            first: Collection[float],
            second: Collection[float],
            *,
            msg: Optional[str]=None,
            delta: Optional[float]=None
    ) -> None: ...

    def assertListsAlmostEqual(
            self,
            first: Collection[float],
            second: Collection[float],
            places: Optional[int]=None,
            msg: Optional[str]=None,
            delta: Optional[float]=None
    ) -> None:
        '''Element-wise float collection comparison

        http://stackoverflow.com/a/8312110
        '''
        self.assertEqual(len(first), len(second))
        for a, b in zip(first, second):
            self.assertAlmostEqual(a, b, places, msg, delta)  # type: ignore # (offending arguments forwarded from identically overloaded definition)

    def assertEqualRecursive(self) -> None:
        # TODO: Would be useful for symmetry.
        raise NotImplementedError()

    @overload
    def _assertDataclassesAlmostEqual(
            self,
            first: Any,
            second: Any,
            places: Optional[int]=None,
            msg: Optional[str]=None
    ) -> None: ...

    @overload
    def _assertDataclassesAlmostEqual(
            self,
            first: Any,
            second: Any,
            *,
            msg: Optional[str]=None,
            delta: Optional[float]=None
    ) -> None: ...

    def _assertDataclassesAlmostEqual(
            self,
            first: Any,
            second: Any,
            places: Optional[int]=None,
            msg: Optional[str]=None,
            delta: Optional[float]=None
    ) -> None:
        is_dataclass: Callable[[Any], bool] = lambda obj: dataclasses.is_dataclass(obj) and not isinstance(obj, type)

        if not is_dataclass(first) or not is_dataclass(second):
            raise TypeError("At least one of the supplied objects is not a dataclass.")

        if not isinstance(first, type(second)) or not isinstance(second, type(first)):
            self.fail(f"Dataclass types differ: {type(first)} v. {type(second)}")

        self.assertAlmostEqualRecursive(  # type: ignore # (offending arguments forwarded from identically overloaded definition)
            dataclasses.astuple(first),
            dataclasses.astuple(second),
            places,
            msg,
            delta
        )

    @overload
    def assertAlmostEqualRecursive(
            self,
            first: Any,
            second: Any,
            places: Optional[int]=None,
            msg: Optional[str]=None,
    ) -> None: ...

    @overload
    def assertAlmostEqualRecursive(
            self,
            first: Any,
            second: Any,
            *,
            msg: Optional[str]=None,
            delta: Optional[float]=None
    ) -> None: ...

    def assertAlmostEqualRecursive(
            self,
            first: Any,
            second: Any,
            places: Optional[int]=None,
            msg: Optional[str]=None,
            delta: Optional[float]=None
    ) -> None:

        try:
            self.assertAlmostEqual(first, second, places, msg, delta)  # type: ignore # (offending arguments forwarded from identically overloaded definition)
        except TypeError:
            try:
                self._assertDataclassesAlmostEqual(first, second, places, msg, delta)  # type: ignore # (offending arguments forwarded from identically overloaded definition)
            except TypeError:
                try:
                    if len(first) != len(second):
                        self.fail(
                            f"Comparison failed due to sequence lengths ({len(first)} v. {len(second)}):"
                            f"\n{first}\nv.\n{second}"
                        )
                    self.assertEqual(len(first), len(second))
                    zipped = zip(first, second)
                except TypeError:
                    self.assertEqual(first, second, msg)

                for a, b in zipped:
                    self.assertAlmostEqualRecursive(a, b, places, msg, delta)  # type: ignore # (recursive call with offending arguments forwarded)
