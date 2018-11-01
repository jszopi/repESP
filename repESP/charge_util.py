"""Utility functions for working with partial charges"""

from .charges import Charge
from .equivalence import Equivalence

import numpy as np
from typing import Dict, List


def average(charges: List[Charge], equivalence: Equivalence) -> List[Charge]:
    """Average values of charges on chemically equivalent atoms

    Parameters
    ----------
    charges : List[Charge]
        The charges to be averaged.
    equivalence : Equivalence
        Information about the chemical equivalence of atoms (symmetry-related
        and fast-exchanging atoms).

    Returns
    -------
    List[Charge]
        The list of averaged charges.
    """

    if len(charges) != len(equivalence.values):
        raise ValueError(
            f"Mismatched lengths of charges ({len(charges)}) and equivalence "
            f"({len(equivalence.values)}) arguments."
        )

    equivalent_groups: List[List[int]] = []
    equivalent_group_of_atom: Dict[int, int] = {}

    for referrer, referred in enumerate(equivalence.values):
        if referred is not None and referred in equivalent_group_of_atom:
            equivalent_group = equivalent_group_of_atom[referred]
            equivalent_groups[equivalent_group].append(referrer)
            equivalent_group_of_atom[referrer] = equivalent_group
        else:
            equivalent_groups.append([referrer])
            equivalent_group_of_atom[referrer] = len(equivalent_groups) - 1

    return [
        Charge(
            np.mean([charges[j] for j in equivalent_groups[equivalent_group_of_atom[i]]])
        ) for i in range(len(charges))
    ]
