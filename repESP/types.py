"""Fundamental types used to describe molecules in space"""

from repESP._util import elements, get_symbol, get_atomic_number
from repESP.util import angstrom_per_bohr

from dataclasses import dataclass
from typing import Any, Generic, List, Iterable, Tuple, Type, TypeVar


# NewType had many limitations: not supported in sphinx, not possible to
# override __repr__, and requirement for standalone helper constructors, like:
# https://github.com/python/typing/issues/415#issuecomment-297401553
DistT = TypeVar('DistT', bound='Dist')
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

    def __new__(cls: Type[DistT], x: Any) -> DistT:
        return super().__new__(cls, float(x))  # type: ignore # (Too many arguments for "__new__" of "object")

    @classmethod
    def from_angstrom(cls: Type[DistT], value: Any) -> DistT:
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


CoordsT = TypeVar('CoordsT', bound='Coords')
class Coords(Tuple[Dist, Dist, Dist]):
    """Three-dimensional coordinates of a point in space, in atomic units

    Parameters
    ----------
    values : Iterable[Any]
        An iterable yielding three values of any type convertible to float
        representing the coordinates in units of Bohr radii. The values will be
        converted to `Dist`.
    """

    __slots__ = ()

    # A constructor requiring the element type to be Dist would be more strict
    # but less convenient, to be discussed in a future library revision.
    def __new__(cls: Type[CoordsT], values: Iterable[Any]) -> CoordsT:
        self_to_be = super().__new__(cls, (Dist(value) for value in values))  # type: ignore # https://github.com/python/mypy/issues/8541
        if len(self_to_be) != 3:
            raise ValueError("Coords constructor expected an iterable yielding three elements.")
        return self_to_be

    # TODO: Consider implementing __str__ to avoid falling back on Dist.__repr__.
    # Needs to be tackled more broadly than just this type, but I'm not sure
    # what to do about dataclasses, which would require overriding the whole
    # __str__ even if only one member has an interesting __str__ implementation.


AtomT = TypeVar('AtomT', bound='Atom')
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

    def __post_init__(self) -> None:
        if self.atomic_number < 1 or self.atomic_number >= len(elements):
            raise ValueError(
                f"Atomic number is not within expected bounds: {self.atomic_number}."
            )

    @property
    def symbol(self) -> str:
        """The chemical symbol of this atom"""
        return get_symbol(self.atomic_number)

    @classmethod
    def from_symbol(cls: Type[AtomT], symbol: str, *args: Any, **kwargs: Any) -> AtomT:
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
