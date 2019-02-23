"""Functions producing or acting on molecular fields

Attributes
----------
"""

from repESP.fields import Esp, Field, AbstractMesh
from repESP.charges import AtomWithCoordsAndCharge
from repESP.types import AtomWithCoords, Coords, Dist, Molecule

from scipy.spatial.distance import euclidean
import numpy as np
from typing import cast, Collection, List, Optional, Tuple, TypeVar


def esp_from_charges(mesh: AbstractMesh, molecule: Molecule[AtomWithCoordsAndCharge]) -> Field[Esp]:
    """Calculate ESP value at specified points due to charges on atoms

    Parameters
    ----------
    mesh : AbstractMesh
        The points at which the ESP values are to be calculated.
    molecule : Molecule[AtomWithCoordsAndCharge]
        A molecule with atom coordinates and partial charges specified.

    Returns
    -------
    Field[Esp]
        The ESP field at the specified points reproduced from the partial
        charges on atoms of the given molecule.
    """
    value_at_point = lambda coords: Esp(sum(
        atom.charge/(euclidean(coords, atom.coords))
        for atom in molecule.atoms
    ))

    return Field(
        mesh,
        [value_at_point(coords) for coords in mesh.points]
    )


def voronoi(mesh: AbstractMesh, molecule: Molecule[AtomWithCoords]) -> Field[Tuple[Optional[int], Dist]]:
    """Find the atom closest to each point and its distance

    Example
    -------
    Imagine a molecule of carbon monoxide placed along the x-axis.

    >>> mesh = Mesh([Coords([0, 0, 0])])
    >>> molecule = Molecule([
        AtomWithCoords(6, Coords([-2, 0, 0])),
        AtomWithCoords(8, Coords([0.13, 0, 0]))
    ])
    >>> print(voronoi(mesh, molecule).values)
    [(1, 0.13)]

    The result means that for the only point of the mesh (0, 0, 0), atom of
    index 1 (i.e. the second atom, oxygen) is closer than any other atom
    (carbon) and the distance to it is 0.13 a.u.

    Parameters
    ----------
    mesh : AbstractMesh
        The points at which the ESP values are to be calculated.
    molecule : Molecule[AtomWithCoords]
        A molecule consisting of atoms with the coordinates specified.

    Returns
    -------
    Field[Tuple[Optional[int], Dist]]
        A `Field` object specifying for each point the atom to which the
        point is nearest (represented as ordinal, zero-based index into the
        molecule) and the distance from that atom.
    """
    def value_at_point(coords: Coords) -> Tuple[Optional[int], Dist]:
        min_dist = Dist(float('inf'))
        min_atom = None
        for atom_index, atom in enumerate(molecule.atoms):
            dist = euclidean(coords, atom.coords)
            if dist < min_dist:
                min_dist = dist
                min_atom = atom_index
        return (min_atom, min_dist)


    return Field(
        mesh,
        [value_at_point(coords) for coords in mesh.points]
    )


# Meant to mirror fields.Field.NumericValue, and similarly the bound should be
# numbers.Number but mypy throws errors.
NumericValue = TypeVar('NumericValue', bound=float)
"""TypeVar : Generic type for numeric values

This type matches that of `fields.Field.NumericValue` and similarly it can be
any type matching "bound=float".
"""


def calc_rms_value(values: Collection[NumericValue]) -> NumericValue:
    """Calculate root-mean-square (RMS) value of a collection of values

    Parameters
    ----------
    values : Collection[NumericValue]
        The values which RMS is to be calculated.

    Returns
    -------
    NumericValue
        RMS value fo the given values.
    """
    return np.sqrt(np.mean(np.square(values)))


def calc_rms_error(
    values1: Collection[NumericValue],
    values2: Collection[NumericValue]
) -> NumericValue:
    """Calculate RMS error between two collections of values

    Parameters
    ----------
    values1 : Collection[NumericValue]
        One of the collections of values to be used for the calculation.
    values2 : Collection[NumericValue]
        The other collection of values to be used for the calculation.

    Returns
    -------
    NumericValue
        RMS error between the two given collections of values
    """
    # TODO: The cast wouldn't be required if subtraction was implemented to
    # yield the correct type.
    return calc_rms_value([cast(NumericValue, value1 - value2) for value1, value2 in zip(values1, values2)])


def calc_relative_rms_error(
    values1: Collection[NumericValue],
    values2: Collection[NumericValue]
) -> float:
    """Calculate relative RMS error of two collections of values

    This is calculated as the RMS error between the two collections of values
    (as given by `calc_rms_error`) divided by the RMS error of the first
    collection of values.

    Parameters
    ----------
    values1 : Collection[NumericValue]
        One of the collections of values which are to be used for the
        calculation. The RMS error will be reported relative to the RMS value
        of this collection.
    values2 : Collection[NumericValue]
        The other collection of values to be used for the calculation.

    Returns
    -------
    float
        Relative RMS error of the field values.
    """
    return calc_rms_error(values1, values2)/calc_rms_value(values1)
