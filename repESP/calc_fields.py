from .types import Charges, Coords, Dist, Esp, Field, Mesh, Molecule

from scipy.spatial.distance import euclidean
from typing import List, NewType, Optional, Tuple


def _esp_from_charges_at_point(coords: Coords, charges: Charges) -> Esp:

    return Esp(sum(
        charge/(euclidean(coords, atom.coords))
        for atom, charge in zip(charges.molecule.atoms, charges.values)
    ))


def esp_from_charges(mesh: Mesh, charges: Charges) -> Field[Esp]:
    """Calculate ESP value at given point due to charges on atoms"""
    return mesh.calc_field(
        lambda coords: _esp_from_charges_at_point(coords, charges)
    )


def _voronoi_at_point(coords: Coords, molecule: Molecule) -> Tuple[Optional[int], Dist]:
    """For a given point, find the closest atom and its distance"""
    min_dist = Dist(float('inf'))
    min_atom = None
    for atom_index, atom in enumerate(molecule.atoms):
        dist = euclidean(coords, atom.coords)
        if dist < min_dist:
            min_dist = dist
            min_atom = atom_index
    return (min_atom, min_dist)


def voronoi(mesh: Mesh, molecule: Molecule) -> Field[Tuple[Optional[int], Dist]]:
    # Voronoi means closest-atom in molecular partitioning lingo
    return mesh.calc_field(
        lambda coords: _voronoi_at_point(coords, molecule)
    )
