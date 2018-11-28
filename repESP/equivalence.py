"""Chemical equivalence relations between atoms of a molecule"""

from dataclasses import dataclass, asdict
from itertools import zip_longest
import io
from typing import List, Optional

from repESP.types import Atom, Molecule


@dataclass
class Equivalence:
    """Dataclass representing chemical equivalence relations between atoms

    Atoms in a molecule are considered equivalent for the purpose of fitting
    partial charges when they are symmetry-related or fast-exchanging.

    The equivalence information is represented as a list stored in the `values`
    attribute. The length of this list must be the same as the number of atoms
    in the molecule it describes. Consecutive values refer to consecutive atoms
    of the molecule. Each value is either None, if the described atom is not
    equivalenced to any other atom, or a zero-based index of the atom to which
    the described atom is equivalenced.

    Example
    -------

    Consider a methane molecule defined as follows:

    >>> methane = Molecule([Atom(atomic_number) for atomic_number in [6, 1, 1, 1, 1]])

    The corresponding equivalence information would be:

    >>> equivalence = Equivalence([None, None, 1, 1, 1])

    There are no atoms equivalent to the carbon atom, and thus its value is
    None. The first hydrogen atom is equivalent to the other three but those
    have not been specified yet, so its value is also None. The remaining
    three hydrogen atoms are equivalent to the first hydrogen atom, which
    zero-based index in the `methane.atoms` list is 1.

    Parameters
    ----------
    values : typing.List[Optional[int]]
        The list of equivalence values.

    Raises
    ------
    ValueError
        Raised when any of the values are outside of the expected bounds. In
        the future the initialization may further verify that there are no cyclic
        relations.

    Attributes
    ----------
    values
        See initialization parameter
    """
    values: List[Optional[int]]

    def __post_init__(self) -> None:
        for i, elem in enumerate(self.values):
            if elem is not None and (elem < 0 or elem >= len(self.values)):
                raise ValueError(
                    f"Value number {i} is not valid as equivalence information "
                    f"in a molecule of {len(self.values)}."
                )

    def describe(self, molecule: Optional[Molecule[Atom]]=None) -> str:
        """Verbosely report the equivalence information

        Example
        -------

        >>> print(equivalence.describe(methane))
        Atom (C) number 1
        Atom (H) number 2
        Atom (H) number 3, equivalenced to atom 1
        Atom (H) number 4, equivalenced to atom 2
        Atom (H) number 5, equivalenced to atom 2

        Parameters
        ----------
        molecule : Optional[Molecule[Atom]], optional
            The molecule to which the equivalence information refers. This
            argument is optional and defaults to None. If it is provided, atom
            identities will be included in the output.

        Raises
        ------
        ValueError
            Raised when the number of atoms in the molecule does not match
            the length of the list of values in this object.

        Returns
        -------
        str
            A verbose description of the equivalence information.
        """
        if molecule is not None and len(molecule.atoms) != len(self.values):
            raise ValueError(
                f"The number of atoms ({len(molecule.atoms)} is not the same "
                f"as the number of equivalence values ({len(self.values)}."
            )

        zipped = zip_longest(self.values, molecule.atoms if molecule is not None else [])

        f = io.StringIO()
        for i, (equivalence, atom) in enumerate(zipped):
            atomic_number = atom.symbol if molecule is not None else None
            id_str = f" ({atomic_number})" if atomic_number is not None else ""
            equivalence_str = f", equivalenced to atom {equivalence+1}" if equivalence is not None else ""
            print(f"Atom{id_str} number {i+1}{equivalence_str}", file=f)

        return f.getvalue()


def _get_equivalence_from_ac(ac) -> Equivalence:
    # TODO: This should wrap the two respgen calls:
    #     respgen -i methane.ac -o methane.respin1 -f resp1
    #     respgen -i methane.ac -o methane.respin2 -f resp2
    # and call get_equivalence on them.
    raise NotImplementedError()
