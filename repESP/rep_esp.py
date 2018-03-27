from scipy.spatial.distance import euclidean
import numpy as np

from .resp_helpers import Points, NonGridField
from .cube_helpers import GridField, Grid, angstrom_per_bohr

au_per_debye = 0.393430307

    # PERFORMANCE CONSIDERATIONS
    # This will likely be a bottleneck and, while a mathematical trick to
    # reduce the complexity from O(n*g^3) may exist, I don't know one.
    # Still, this method should be optimized but **only once its proven a
    # bottleneck** by profiling.
    # (1) Some straightforward optimization has been performed. The reason
    # why this method expects a field_func returning more than one value is
    # that calling field_func only once prevents the reevaluation of some
    # common parts of its body. For _rep_esp_func and _dist_func, it is the
    # `euclidean' distance, and for the latter also the logic behind
    # choosing the closest atom.
    # (2) However, `euclidean' still gets reevaluated for different calls
    # of this method. It could be memoized as a distance field resulting
    # from a new Molecule method. TODO
    # (3) np.array should be preferred over the intermediate `results` list
    # (4) Then it may be worth mapping field_func onto array elements
    # instead of iterating over it but that's a disputable topic:
    # https://wiki.python.org/moin/PythonSpeed/PerformanceTips#Loops
    # https://www.python.org/doc/essays/list2str/
    # (5) Finally, iterating the points has a good potential for
    # parallelization.

def calc_dipole(molecule, charge_type):
    """3D components of dipole moment in Debye"""
    # Wouldn't this function fit better to another module?
    dipole = []
    for coord in range(3):
        component = 0
        for atom in molecule:
            component += (atom.charges[charge_type]*atom.coords[coord] /
                          angstrom_per_bohr / au_per_debye)
        dipole.append(component)
    return dipole
