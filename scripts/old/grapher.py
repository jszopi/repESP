from repESP import cube_helpers
from repESP.charges import update_with_charges
from repESP.field_comparison import difference, filter_by_dist, filter_by_atom
from repESP.field_comparison import rms_and_rrms
from repESP.graphs import plot
from repESP.rep_esp import calc_grid_field

import os
import copy
import numpy as np

# NOTE: This ad-hoc script has been replaced with the more general fit_points.py

molecule_name = 'methane'
path = '../data/' + molecule_name + '/'
output_path = path + 'grapher_output/'
os.mkdir(output_path)

if not os.path.isdir(path):
    raise FileNotFoundError("The input directory does not exist! Currently "
                            "this script assumes a certain directory tree and "
                            "should be run from its own `scripts` directory.")

esp_cube = cube_helpers.Cube(path + molecule_name + '_esp.cub')
ed_cube = cube_helpers.Cube(path + molecule_name + '_den.cub')

molecule = esp_cube.molecule

# SET OPTIONS HERE
isoval = 0.01
dist = ed_cube.field.distance_transform(isoval)
exclusion_dist = 0
rand_skim = 1
axes_limits = [[0, 4]]

charge_types = {'aim': '.sumviz',
                # 'aim': 'ACF.dat',
                'mulliken': '.log',
                'chelpg': '_chelpg.log',
                'mk': '_mk.log',
                'nbo': '_nbo.log',
                }

# Use the same color span in all plots (true ESP doesn't depend on charges so
# common to all charge_types and can be outside of loop):
_dist, esp_values = filter_by_dist(exclusion_dist, dist, esp_cube.field)
color_span = [np.nanmin(esp_values), np.nanmax(esp_values)]
# A symmetric span can also be used:
# color_limit = max(abs(np.nanmin(esp_values)), abs(np.nanmax(esp_values)))
# color_span = [-color_limit, color_limit]

# Note that the y-axis is not set for different charge_types! If this is
# desired behaviour, it could be set in a similar manner here.

for charge_type in charge_types.keys():
    charge_dir = output_path + charge_type
    # raises OSError if directory exists
    os.mkdir(charge_dir)
    if charge_types[charge_type][0] in ['.', '_']:
        filename = path + molecule_name + charge_types[charge_type]
    else:
        filename = path + charge_types[charge_type]

    update_with_charges(charge_type, filename, molecule)

    print("\n{0} charges:".format(charge_type.upper()))
    for atom in molecule:
        atom.print_with_charge(charge_type)
    # The same but to file (would be better with my Tee class from featsel)
    with open(charge_dir + '/charges.txt', 'w') as f:
        for atom in esp_cube.molecule:
            atom.print_with_charge(charge_type, f=f)

    # This is costly but was designed to be easy for many charge types, so
    # should be moved outside of the loop
    rep = calc_grid_field(esp_cube.molecule, esp_cube.field.grid, 'rep_esp',
                          [charge_type])[0]
    # Change details of calculating difference here (absolute, relative)
    diff = difference(esp_cube.field, rep)

    # Write cube
    diff.write_cube(charge_dir + '/diff.cub', molecule, charge_type)
    # and with elements within the ED isosurface excluded
    diff_filtered = copy.deepcopy(diff)
    diff_filtered.check_nans = False
    _dist, diff_filtered.values = filter_by_dist(exclusion_dist, dist, diff)
    diff_filtered.write_cube(charge_dir + '/diff_filtered-iso_' + str(isoval) +
                             '.cub', molecule, charge_type)

    rep_filtered = copy.deepcopy(rep)
    rep_filtered.check_nans = False
    _dist, rep_filtered.values = filter_by_dist(exclusion_dist, dist, rep)
    rms_val, rrms_val = rms_and_rrms(esp_cube.field, rep_filtered,
                                     ignore_nans=True)
    text = "RMS on filtered cube points:"
    # Report both RMS and RRMS
    text = "\n {0} {1:.7f}\nR{0} {2:.7f}".format(text, rms_val, rrms_val)
    print(text)
    # TODO: my Tee object from my previous project would help DRY
    with open(charge_dir + '/charges.txt', 'a') as f:
        print(text, file=f)

    # Whole molecule plot
    title = "Not filtered (whole molecule)"
    # Use given limits and save the limits that were actually used in this
    # graph to reuse in all the plots.
    # Actually this is redundant because the y-axis is scaled to the whole data
    # set in each case.
    limits = []
    plot(dist, diff, color=esp_cube.field, color_span=color_span,
         dist_field_filter=dist, exclusion_dist=exclusion_dist,
         rand_skim=rand_skim, axes_limits=axes_limits,
         save_to=charge_dir+'/whole_mol.png', title=title, get_limits=limits)

    for basin, basin_name in zip(('dist', 'qtaim'), ('Voronoi', 'QTAIM')):

        basin_dir = charge_dir + '/' + basin_name
        os.mkdir(basin_dir)

        for atom in molecule:
            atom_id = str(atom.label) + atom.atomic_number
            atom_fn = basin_dir + '/' + atom_id + '.png'
            title = (charge_type.upper() + ", filtered by " + basin_name +
                     " basin of " + atom_id)

            atom_filter = lambda *fields: filter_by_atom(molecule, atom.label,
                                                         basin, *fields)[1:]

            diff_atom = copy.deepcopy(diff)
            diff_atom.check_nans = False
            _atom, diff_atom.values = filter_by_atom(molecule, atom.label,
                                                     basin, diff_atom)
            diff_atom.write_cube(basin_dir + '/diff_filtered-atom_' +
                                 atom_id + '.cub', molecule, charge_type)

            plot(dist, diff, color=esp_cube.field, color_span=color_span,
                 dist_field_filter=dist, exclusion_dist=exclusion_dist,
                 rand_skim=rand_skim, extra_filter=atom_filter,
                 # Uses limits from the whole-molecule plot
                 axes_limits=limits, save_to=atom_fn, title=title)
