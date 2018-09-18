from .charges import MoleculeWithCharges
from .types import Coords, Dist, Esp, Field, Mesh, NumericField, NumericFieldValue, Molecule
from .types import make_esp

from scipy.spatial.distance import euclidean
from typing import List, NewType, Optional, Tuple


def _esp_from_charges_at_point(coords: Coords, molecule_with_charges: MoleculeWithCharges) -> Esp:

    return make_esp(sum(
        charge/(euclidean(coords, atom.coords))
        for atom, charge in zip(molecule_with_charges.molecule.atoms, molecule_with_charges.charges)
    ))


def esp_from_charges(mesh: Mesh, molecule_with_charges: MoleculeWithCharges) -> NumericField[Esp]:
    """Calculate ESP value at given point due to charges on atoms"""
    return NumericField(
        mesh,
        [_esp_from_charges_at_point(coords, molecule_with_charges) for coords in mesh.points()]
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
    return Field(
        mesh,
        [_voronoi_at_point(coords, molecule) for coords in mesh.points()]
    )
