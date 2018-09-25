from enum import Enum


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


def get_symbol(atomic_number: int) -> str:
    return _elements[atomic_number-1][0]


def get_atomic_number(symbol: str) -> int:
    return _symbol_to_atomic_nummber[symbol]
