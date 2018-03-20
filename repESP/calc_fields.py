from .types import Charges, Coords, Field, FieldValue, Mesh, Molecule

from scipy.spatial.distance import euclidean
from typing import Callable, List, Tuple


def calc_field(
    mesh: Mesh,
    field_at_point: Callable[[Coords], FieldValue]
) -> Field[FieldValue]:

    values: List[FieldValue] = []

    for point in mesh.points():
        values.append(
            field_at_point(point)
        )

    return Field(mesh, values)


def _esp_from_charges_at_point(coords: Coords, charges: Charges) -> float:

    return sum(
        # NOTE: I've removed the conversion from angstrom to bohr
        charge/(euclidean(coords, atom.coords))
        for atom, charge in zip(charges.molecule.atoms, charges.values)
    )


def esp_from_charges(mesh: Mesh, charges: Charges) -> Field[float]:

    return calc_field(
        mesh,
        lambda coords: _esp_from_charges_at_point(coords, charges)
    )


def _voronoi_at_point(coords: Coords, molecule: Molecule) -> Tuple[int, float]:
    """For a given point, find the closest atom and its distance"""
    min_dist = float('inf')
    min_atom = None
    for atom_index, atom in enumerate(molecule.atoms):
        dist = euclidean(coords, atom.coords)
        if dist < min_dist:
            min_dist = dist
            min_atom = atom_index
    return (min_atom, min_dist)


def voronoi(mesh: Mesh, molecule: Molecule):
    # Voronoi means closest-atom in molecular partitioning lingo
    return calc_field(
        mesh,
        lambda coords: _voronoi_at_point(coords, molecule)
    )
