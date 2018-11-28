"""Fundamental types used to describe molecules in space"""

from repESP._util import elements, get_symbol, get_atomic_number
from repESP.util import angstrom_per_bohr

from dataclasses import dataclass
from typing import Any, Generic, List, Tuple, TypeVar


# NewType had many limitations: not supported in sphinx, not possible to
# override __repr__, and requirement for standalone helper constructors, like:
# https://github.com/python/typing/issues/415#issuecomment-297401553
class Dist(float):
    """Distance in atomic units i.e. Bohr radii (a\ :sub:`0`\ )

    Bohr radius is the measure of length in atomic units and angstrom is a unit
    commonly used in computational chemistry. The string representations are
    thus implemented as follows:

    >>> d = Dist(2.3)
    >>> d
    Dist(2.3)
    >>> print(d)
    2.3 a₀ (1.2 Å)

    Parameters
    ----------
    value : Any
        Any value convertible to float representing the value in units of Bohr
        radii, e.g. 1.23 or "2.34".
    """

    __slots__ = ()

    def __new__(cls, x: Any):
        return super().__new__(cls, float(x))  # type: ignore # (Too many arguments for "__new__" of "object")

    @classmethod
    def from_angstrom(cls, value: Any):
        """Alternative initialization from value in angstrom (Å)

        Parameters
        ----------
        value : Any
            Any value convertible to float representing the value in units of
            angstrom.
        """
        return cls(float(value)/angstrom_per_bohr)

    def angstrom(self) -> float:
        """Value in angstrom (Å)"""
        return angstrom_per_bohr*self

    def __repr__(self) -> str:
        return f"Dist({super().__repr__()})"

    def __str__(self) -> str:

        def get_decimals_in_str(str_: str) -> int:
            return len(str_) - str_.find(".") - 1

        # Angstroms are displayed with the same decimal precision as a.u. value
        angstrom_format = f"{{:.{get_decimals_in_str(super().__str__())}f}}"
        str_ = f"{super().__str__()} a₀ ({angstrom_format} Å)"
        return str_.format(self.angstrom())


class Coords(tuple):
    """Three-dimensional coordinates of a point in space, in atomic units

    Parameters
    ----------
    values : Tuple[Any, Any, Any]
        Tuple of three values of any type convertible to float representing the
        coordinates in units of Bohr radii. The values will be converted to `Dist`.
    """

    __slots__ = ()

    # A constructor from individual distances would be more convenient but that
    # was causing issues due to libraries and built-ins assuming the same
    # interface as `tuple`.
    def __new__(cls, values: Tuple[Any, Any, Any]):
        # To allow generators (required by `dataclasses.astuple`)
        values = tuple(values)  # type: ignore # (`values` is now a variadic-size tuple)
        return super().__new__(cls, tuple((Dist(values[0]), Dist(values[1]), Dist(values[2]))))

    # TODO: Implement __str__ to avoid falling back on Dist.__repr__


@dataclass
class Atom:
    """Dataclass representing an atom of a certain element

    If more information should be associated with an atom, it's easy to inherit
    from this class.

    Parameters
    ----------
    atomic_number : int
        The atomic number of the element

    Raises
    ------
    ValueError
        Raised in case the given atomic number is not recognized

    Attributes
    ----------
    atomic_number
        See initialization parameter
    """
    atomic_number: int

    def __post_init__(self):
        if self.atomic_number < 1 or self.atomic_number >= len(elements):
            raise ValueError(
                f"Atomic number is not within expected bounds: {self.atomic_number}."
            )

    @property
    def symbol(self) -> str:
        """The chemical symbol of this atom"""
        return get_symbol(self.atomic_number)

    @classmethod
    def from_symbol(cls, symbol: str, *args, **kwargs):
        """Alternative initialization from chemical symbol

        Despite not being documented, this constructor is inherited by the
        subclasses of this class. `*args` and `**kwargs` will be forwarded to
        the initializer of the subclass.

        Parameters
        ----------
        symbol : str
            The chemical symbol of the element, e.g. "Be".

        Raises
        ------
        ValueError
            Raised in case the given chemical symbol is not recognized
        """
        # Generic type annotations as per:
        # https://github.com/python/typing/issues/58#issuecomment-326240794)
        # don't seem to be working for a dataclass, and that's even before
        # getting them to work with args and kwargs.
        return cls(get_atomic_number(symbol), *args, **kwargs)  # type: ignore # (args, kwargs)


@dataclass
class AtomWithCoords(Atom):
    """Dataclass representing an atom of a certain element in space

    Parameters
    ----------
    atomic_number : int
        The atomic number of the element
    coords : Coords
        The coordinates of the atom in space

    Raises
    ------
    ValueError
        Raised in case the given atomic number is not recognized

    Attributes
    ----------
    atomic_number
        See initialization parameter
    coords
        See initialization parameter
    """
    coords: Coords


GenericAtom_co = TypeVar('GenericAtom_co', bound=Atom, covariant=True)

@dataclass
class Molecule(Generic[GenericAtom_co]):
    """A basic dataclass representing a molecule

    "Basic" because no characteristics of the molecule, like bonding, are
    stored in this class in addition to the list of constituent atoms.

    This class is generic in the type of `GenericAtom_co`, which must be a subtype
    of `Atom`. This allows the molecule to store `Atom` objects as well as
    other objects, e.g. `AtomWithCoords`.

    Parameters
    ----------
    atoms : typing.List[Atom]
        The constituent atoms of the molecule. Any `Atom` subtype is valid.

    Attributes
    ----------
    atoms
        See initialization parameter
    """
    atoms: List[GenericAtom_co]
