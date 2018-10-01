from .exceptions import InputFormatError
from .charges import MoleculeWithCharges
from .types import Coords, Dist, Esp, Field, Mesh, Molecule
from .types import make_esp

from scipy.spatial.distance import euclidean
import numpy as np
from typing import cast, List, NewType, Optional, Tuple


def _esp_from_charges_at_point(coords: Coords, molecule_with_charges: MoleculeWithCharges) -> Esp:

    return make_esp(sum(
        charge/(euclidean(coords, atom.coords))
        for atom, charge in zip(molecule_with_charges.molecule.atoms, molecule_with_charges.charges)
    ))


def esp_from_charges(mesh: Mesh, molecule_with_charges: MoleculeWithCharges) -> Field[Esp]:
    """Calculate ESP value at given point due to charges on atoms"""
    return Field(
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


def calc_rms_value(field: Field[Field.NumericValue]) -> Field.NumericValue:
    return np.sqrt(np.mean(np.square(field.values)))


def calc_rms_error(
    field1: Field[Field.NumericValue],
    field2: Field[Field.NumericValue]
) -> Field.NumericValue:

    if field1.mesh != field2.mesh:
        raise InputFormatError(
            "Calculating RMS requires the underlying grids to be the same for both fields."
        )

    return calc_rms_value(field1-field2)


def calc_relative_rms_error(
    field1: Field[Field.NumericValue],
    field2: Field[Field.NumericValue]
) -> float:
    return calc_rms_error(field1, field2)/calc_rms_value(field1)
