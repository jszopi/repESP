from repESP import cube_helpers, resp
from repESP.rep_esp import calc_grid_field
from repESP.field_comparison import difference, filter_by_dist, filter_by_atom
from repESP.charges import update_with_charges, _update_molecule_with_charges

import copy
import os
import pickle
import sys

# NOTE: This ad-hoc script has been replaced with the more general field_diff.py

molecule_name = 'MMIm_plus'
path = '../../../Results/cations/{0}/'.format(molecule_name)
output_path = path + 'compare_diffs/'
os.mkdir(output_path)

# log and esp files
fn = lambda ext, charge_type: path + molecule_name + "_{0}.{1}".format(
    charge_type, ext)

esp_cube = cube_helpers.Cube(path + molecule_name + '_esp.cub')
ed_cube = cube_helpers.Cube(path + molecule_name + '_den.cub')
molecule = esp_cube.molecule

isoval = 0.01
dist = ed_cube.field.distance_transform(isoval)
exclusion_dist = 0

check_ivary = True
for charge_type in ['mk', 'chelpg']:

    esp_equiv_molecule = resp.run_resp(
        path, output_path + 'unrest_' + charge_type, resp_type='unrest',
        esp_fn=molecule_name + '_' + charge_type + '.esp',
        check_ivary=check_ivary)

    check_ivary = False
    charges = [atom.charges['resp'] for atom in esp_equiv_molecule]
    _update_molecule_with_charges(molecule, charges, charge_type + '_equiv')

    pickled_fn = path + "compromise_nbo_and_{0}/molecule.p".format(charge_type)
    with open(pickled_fn, "rb") as f:
        pickled_mol = pickle.load(f)

    charges = [atom.charges[charge_type] for atom in pickled_mol]
    _update_molecule_with_charges(molecule, charges, 'compr_' + charge_type)

charges = [atom.charges['nbo_equiv'] for atom in pickled_mol]
_update_molecule_with_charges(molecule, charges, 'nbo_equiv')

all_charges = ['chelpg_equiv',
               'mk_equiv',
               'nbo_equiv',
               'compr_chelpg',
               'compr_mk'
               ]

for charge_type in all_charges:

    print('\n', charge_type.upper())
    for atom in molecule:
        atom.print_with_charge(charge_type)
    sys.stdout.flush()

    rep = calc_grid_field(molecule, esp_cube.field.grid, 'rep_esp',
                          [charge_type])[0]
    # Change details of calculating difference here (absolute, relative)
    diff = difference(esp_cube.field, rep)

    # Write cube
    diff.write_cube(output_path + charge_type + '_diff.cub', molecule,
                    charge_type)
    # and with elements within the ED isosurface excluded
    diff_filtered = copy.deepcopy(diff)
    diff_filtered.check_nans = False
    _dist, diff_filtered.values = filter_by_dist(exclusion_dist, dist, diff)
    diff_filtered.write_cube(output_path + charge_type + '_diff_filtered-iso_'
                             + str(isoval) + '.cub', molecule, charge_type)
