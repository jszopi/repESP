from enum import Enum
from typing import Collection, Dict, Iterable, List, Tuple, TypeVar, Union


# As per Python docs
class _NoValue(Enum):
    def __repr__(self):
        return '<%s.%s>' % (self.__class__.__name__, self.name)


# TODO: this should be handled by a library
_elements = [('H', 'Hydrogen'),
            ('He', 'Helium'),
            ('Li', 'Lithium'),
            ('Be', 'Beryllium'),
            ('B', 'Boron'),
            ('C', 'Carbon'),
            ('N', 'Nitrogen'),
            ('O', 'Oxygen'),
            ('F', 'Fluorine'),
            ('Ne', 'Neon'),
            ('Na', 'Sodium'),
            ('Mg', 'Magnesium'),
            ('Al', 'Aluminum'),
            ('Si', 'Silicon'),
            ('P', 'Phosphorus'),
            ('S', 'Sulfur'),
            ('Cl', 'Chlorine'),
            ('Ar', 'Argon')]

_symbol_to_atomic_nummber = {v[0]: i+1 for i, v in enumerate(_elements)}


def _get_symbol(atomic_number: int) -> str:
    return _elements[atomic_number-1][0]


def _get_atomic_number(symbol: str) -> int:
    return _symbol_to_atomic_nummber[symbol]


T = TypeVar('T')
U = TypeVar('U')

def _zip_exact(first: Collection[T], second: Collection[U]) -> Iterable[Tuple[T, U]]:

    if len(first) != len(second):
        raise ValueError(
            f"Failed to zip arguments due to unequal lengths: {len(first)} v. {len(second)}."
        )

    return zip(first, second)


def _list_from_dict(
    dictionary: Dict[int, T],
    length: int,
    default: U,
    one_indexed: bool=False
) -> List[Union[T, U]]:
    """Convenience function for creating lists from dictionaries

    For sparse data or certain user input it may be easier to represent properties of a
    list (e.g. equivalence information for each atom in a molecule) as a dictionary.
    However, many interfaces of this library expect lists and thus this convenience
    function is provided for the converstion.
    """
    offset = 1 if one_indexed else 0
    return [dictionary[i+offset] if i+offset in dictionary else default for i in range(length)]


def _mask_from_list(  # type: ignore # (defaults with generics: https://github.com/python/mypy/issues/3737)
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
    ('masks') and thus this convenience function is provided for the converstion.
    """
    offset = 1 if one_indexed else 0
    return [value_if_present if i+offset in l else value_if_absent for i in range(length)]
