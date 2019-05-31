from enum import Enum
from typing import Collection, Iterable, TextIO, Tuple, TypeVar


# As per Python docs
class NoValue(Enum):
    def __repr__(self) -> str:
        return '<%s.%s>' % (self.__class__.__name__, self.name)


# TODO: this should be handled by a library
elements = [
    ('H', 'Hydrogen'),
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
    ('Ar', 'Argon')
]

symbol_to_atomic_nummber = {v[0]: i+1 for i, v in enumerate(elements)}


def get_symbol(atomic_number: int) -> str:
    return elements[atomic_number-1][0]


def get_atomic_number(symbol: str) -> int:
    return symbol_to_atomic_nummber[symbol]


T = TypeVar('T')
U = TypeVar('U')

def zip_exact(first: Collection[T], second: Collection[U]) -> Iterable[Tuple[T, U]]:

    if len(first) != len(second):
        raise ValueError(
            f"Failed to zip arguments due to unequal lengths: {len(first)} v. {len(second)}."
        )

    return zip(first, second)

def get_line(f: TextIO) -> str:
    return f.readline().rstrip('\n')
