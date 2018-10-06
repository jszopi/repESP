from typing import Dict, List, TypeVar, Union


T = TypeVar('T')
U = TypeVar('U')


def list_from_dict(
    dictionary: Dict[int, T],
    length: int,
    default: U,
    one_indexed: bool=False
) -> List[Union[T, U]]:
    """Convenience function for creating lists from dictionaries

    For sparse data or certain user input it may be easier to represent properties of a
    list (e.g. equivalence information for each atom in a molecule) as a dictionary.
    However, many interfaces of this library expect lists and thus this convenience
    function is provided for the conversion.
    """
    offset = 1 if one_indexed else 0
    return [dictionary[i+offset] if i+offset in dictionary else default for i in range(length)]


def mask_from_list(  # type: ignore # (defaults with generics: https://github.com/python/mypy/issues/3737)
    l: List[int],
    length: int,
    value_if_present: T=True,
    value_if_absent: U=False,
    one_indexed: bool=False
) -> List[Union[T, U]]:
    """Convenience function for creating masks from lists

    For sparse data or certain user input it may be easier to represent
    properties of a list (e.g. whether each atom belongs to a methyl or
    methylene group) as a list of elements which have a certain value of the
    property. However, some interfaces of this library expect complete lists
    ('masks') and thus this convenience function is provided for the conversion.
    """
    offset = 1 if one_indexed else 0
    return [value_if_present if i+offset in l else value_if_absent for i in range(length)]
