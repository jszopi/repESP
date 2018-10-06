from .respin_util import Equivalence
from .charges import Charge

import numpy as np
from typing import Dict, List


def average(charges: List[Charge], equivalence: Equivalence) -> List[Charge]:

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
