"""Functions producing or acting on molecular fields"""

from .fields import Esp, Field, AbstractMesh
from .exceptions import InputFormatError
from .charges import AtomWithCoordsAndCharge
from .types import AtomWithCoords, Coords, Dist, Molecule

from scipy.spatial.distance import euclidean
import numpy as np
from typing import cast, Collection, List, Optional, Tuple, TypeVar


def _esp_from_charges_at_point(coords: Coords, molecule: Molecule[AtomWithCoordsAndCharge]) -> Esp:

    return Esp(sum(
        atom.charge/(euclidean(coords, atom.coords))
        for atom in molecule.atoms
    ))


def esp_from_charges(mesh: AbstractMesh, molecule: Molecule[AtomWithCoordsAndCharge]) -> Field[Esp]:
    """Calculate ESP value at given point due to charges on atoms"""
    return Field(
        mesh,
        [_esp_from_charges_at_point(coords, molecule) for coords in mesh.points]
    )


def _voronoi_at_point(coords: Coords, molecule: Molecule[AtomWithCoords]) -> Tuple[Optional[int], Dist]:
    """For a given point, find the closest atom and its distance"""
    min_dist = Dist(float('inf'))
    min_atom = None
    for atom_index, atom in enumerate(molecule.atoms):
        dist = euclidean(coords, atom.coords)
        if dist < min_dist:
            min_dist = dist
            min_atom = atom_index
    return (min_atom, min_dist)


def voronoi(mesh: AbstractMesh, molecule: Molecule[AtomWithCoords]) -> Field[Tuple[Optional[int], Dist]]:
    # Voronoi means closest-atom in molecular partitioning lingo
    return Field(
        mesh,
        [_voronoi_at_point(coords, molecule) for coords in mesh.points]
    )


# Meant to mirror fields.Field.NumericValue, and similarly the bound should be
# numbers.Number but mypy throws errors.
NumericValue = TypeVar('NumericValue', bound=float)


def calc_rms_value(values: Collection[NumericValue]) -> NumericValue:

    return np.sqrt(np.mean(np.square(values)))


def calc_rms_error(
    values1: Collection[NumericValue],
    values2: Collection[NumericValue]
) -> NumericValue:
    return calc_rms_value([cast(NumericValue, value1 - value2) for value1, value2 in zip(values1, values2)])


def calc_relative_rms_error(
    values1: Collection[NumericValue],
    values2: Collection[NumericValue]
) -> float:
    return calc_rms_error(values1, values2)/calc_rms_value(values1)
