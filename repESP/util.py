"""Constants and convenience functions for interacting with the library"""

from typing import Dict, List, TypeVar, Union

"""
Attributes
----------
"""

angstrom_per_bohr = 0.5291772086
"""float : The conversion factor from Bohr radii (a₀) to angstrom (Å)

Value used by Gaussian 09, based on:
P. J. Mohr, B. N. Taylor, and D. B. Newell, “CODATA Recommended Values of the
Fundamental Physical Constants: 2006,” *Rev. Mod. Phys.*, **80** (2008)
633-730. DOI: `10.1103/RevModPhys.80.633 <http://dx.doi.org/10.1103/RevModPhys.80.633>`_.
"""


T1 = TypeVar('T1')
T2 = TypeVar('T2')


def list_from_dict(
    dict_: Dict[int, T1],
    length: int,
    default: T2=None,
    one_indexed: bool=False
) -> List[Union[T1, T2]]:
    """Convenience function for creating lists from dictionaries

    For sparse data or certain user input it may be easier to represent
    properties of a list (e.g. equivalence information for each atom in a
    molecule) as a dictionary. However, many interfaces of this library expect
    lists and thus this convenience function is provided for the conversion.

    Example
    -------
    >>> list_from_dict({1: "foo"}, 3, "bar")
    ["bar", "foo", "bar"]

    Parameters
    ----------
    dict\_ : typing.Dict[int, T1]
        The dictionary mapping the index of the desired list to its property,
        which is of the generic type `T1`. If this argument contains keys
        outside the range of valid indices in the desired list, the input
        is not valid; no errors will be raised but these keys will be ignored.

        Valid indices are integer values of `i` that fulfil ``0 <= i <
        length`` if `one_indexed` is False and ``0 < i <= length`` otherwise.
    length : int
        The length of the desired list.
    default : T2, optional
        The value to be used in the output list when the index is not included
        in the `dict_` argument. The type of this value is not required to
        be of the same type `T1` as the values of `dict_`. Defaults to None.
    one_indexed : bool, optional
        Whether the list indices in the keys of the input ``dict_`` are counting
        from one. Defaults to False, meaning they are counting from zero.

    Returns
    -------
    typing.List[Union[T1, T2]]
        A list of items of the given `length` with values specified in `dict_`.
    """
    offset = 1 if one_indexed else 0
    return [
        # mypy infers the result to be a Union including None, whereas I see it
        # as falling under T2.
        dict_[i+offset] if i+offset in dict_ else default  # type: ignore
        for i in range(length)
    ]


def mask_from_list(  # type: ignore # (defaults with generics: https://github.com/python/mypy/issues/3737)
    list_: List[int],
    length: int,
    value_if_present: T1=True,
    value_if_absent: T2=False,
    one_indexed: bool=False
) -> List[Union[T1, T2]]:
    """Convenience function for creating masks from lists

    For sparse data or certain user input it may be easier to represent
    properties of a list (e.g. whether each atom belongs to a methyl or
    methylene group) as a list of elements which have a certain value of the
    property. However, some interfaces of this library expect complete lists
    ('masks') and thus this convenience function is provided for the conversion.

    Example
    -------
    >>> mask_from_list([1], 3)
    [False, True, False]

    Note that the same output can always be achieved with the `list_from_dict`
    function.

    >>> list_from_dict({1: True}, 3, False)
    [False, True, False]

    The difference is that `list_from_dict` allows more than two values in the
    output list. When only two values are needed, `mask_from_list` has a more
    convenient input format.

    Parameters
    ----------
    list\_ : typing.List[int]
        The list of indices, for which the corresponding values in the resulting
        list are to have `value_if_present` as the value (and `value_if_absent`
        otherwise).

        If this argument contains values outside the range of valid indices in
        the desired list, the input is not valid; no errors will be raised
        but these keys will be ignored.

        Valid indices are integer values of `i` that fulfil ``0 <= i <
        length`` if `one_indexed` is False and ``0 < i <= length`` otherwise.
    length : int
        The length of the desired list.
    value_if_present : T1, optional
        The value to be assigned to the element of the output list if the
        corresponding argument is present in the `list_` argument. Defaults
        to True.
    value_if_absent : T2, optional
        The value to be assigned to the element of the output list if the
        corresponding argument is absent from the `list_` argument. Defaults
        to False.
    one_indexed : bool, optional
        Whether the list indices in the keys of the input `dict_` are counting
        from one. Defaults to False, meaning they are counting from zero.

    Returns
    -------
    typing.List[Union[T1, T2]]
        A list of length `length` with values equal to `value_if_present`
        if the corresponding item is in the input `list_` and `value_if_absent`
        otherwise.
    """
    offset = 1 if one_indexed else 0
    return [value_if_present if i+offset in list_ else value_if_absent for i in range(length)]
