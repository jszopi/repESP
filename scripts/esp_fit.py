from repESP import resp, resp_helpers, rep_esp, charges, graphs
from repESP.field_comparison import _check_grids, difference, rms_and_rep
from repESP.resp import get_atom_signature

from numpy import mean, sqrt, square, linspace, meshgrid
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import os
import shutil
import sys

import numpy as np
import pickle

esp_charge_type = 'mk'
alt_esp_charge_type = 'chelpg'
# esp_charge_type = 'chelpg'
# alt_esp_charge_type = 'mk'

molecule_name = 'NMe3H_plus'
# molecule_name = 'methane'

# Parameters to change when changing the molecule
# Methane/water:
vary_label1 = 1
# charge_dict = lambda a1: {1: a1}
xlim1 = (-1, 0)
xlim2 = (-0.5, 1)
# NMe3H_plus, NMe3:
vary_label2 = 13
charge_dict = lambda c, n: {1: c, 5: c, 9: c, 13: n}
# NMe4_plus:
# vary_label2 = 17
# charge_dict = lambda c, n: {1: c, 5: c, 9: c, 13: c, 17: n}
sampling_num = 21
# SLICING plots
# Methane:
# cut_through = 1, 2, 3
# NMe3_plus:
cut_through = 1, 13, 14
# NMe4_plus:
# cut_through = 1, 5, 17
# NMe3:
# cut_through = 1, 5, 9

path = '../data/' + molecule_name + '/'
output_path = path + "esp_fit" + '_' + esp_charge_type + '/'
resp_output_path = output_path + 'resp_calcs/'

charge_type = 'nbo'
charge_log_fn = path + molecule_name + "_" + charge_type + ".log"
temp_output_path = output_path + 'ratio/'

common_fn = path + molecule_name + "_" + esp_charge_type
log_fn = common_fn + ".log"
input_esp = common_fn + ".esp"
esp_fn = molecule_name + "_" + esp_charge_type + ".esp"
alt_esp_fn = molecule_name + "_" + alt_esp_charge_type + ".esp"
output_esp = common_fn + ".esp"

print("esp_fit.py script")
print("\nMolecule:             ", molecule_name)
print("Charge type:          ", esp_charge_type.upper())
print("Non-ESP charge type:  ", charge_type.upper())
print("Alt. ESP charge type: ", alt_esp_charge_type.upper(), '\n')

g = resp_helpers.G09_esp(input_esp)

# Write the .esp file in the correct format expected by the `resp` program
if False:
    g.field.write_to_file(output_esp, g.molecule)

# Division into Voronoi basins:
# parent_atom, dist = rep_esp.calc_non_grid_field(g.molecule, g.field.points,
#                                                 'dist')

title = molecule_name + " " + esp_charge_type.upper()
# Plot the grid in 3 and 2D:
if False:
    color_span = [min(g.field.values), max(g.field.values)]
    for dimension in (3, 2):
        graphs.plot_points(
            g.field, dimension, title=title, molecule=g.molecule,
            plane_eqn=graphs.plane_through_atoms(g.molecule, *cut_through),
            dist_thresh=0.5, axes_limits=[(-5, 5)]*dimension,
            color_span=color_span)
        graphs.plot_points(
            g.field, dimension, title=title, molecule=g.molecule,
            plane_eqn=[1, 0, 0, 0], dist_thresh=0.5,
            axes_limits=[(-5, 5)]*dimension, color_span=color_span)

# Calculate RMS for various charges
import copy
from itertools import chain
molecule = copy.deepcopy(g.molecule)


def polygon_area(vertices):
    x, y = vertices[:, 0], vertices[:, 1]
    # http://stackoverflow.com/a/30408825
    # Tested against TADlib.polygon.shoelace (python2 only)
    return 0.5*np.abs(np.dot(x, np.roll(y, 1))-np.dot(y, np.roll(x, 1)))


def contour_vertices(contour_isovalue, contour_plot, levels):
    contour_index = levels.index(contour_isovalue)
    # http://stackoverflow.com/a/5666461
    return contour_plot.collections[contour_index].get_paths()[0].vertices


def interpret(molecule, charge_dict, vary_label1, vary_label2=None):
    dictio = charge_dict(1, 2)  # Example numbers to get the dict
    print("\nCharges on these atoms will be varied:")
    for vary_label in vary_label1, vary_label2:
        if vary_label is None:
            break
        print('*', molecule[vary_label-1])
        equiv = [label for label in dictio if
                 dictio[label] == dictio[vary_label] and label != vary_label]
        if equiv:
            print("  with the following atoms equivalenced to it:")
            for equiv_label in sorted(equiv):
                print("  -", molecule[equiv_label-1])
    print("\nSee below for equivalence information of other atoms.")


class Result(object):

    def __init__(self, sampling_num, xlim1, xlim2):
        self.sampling_num = sampling_num
        self.xlim1 = xlim1
        self.xlim2 = xlim2
        self.inp1, self.inp2 = self.get_meshgrid(xlim1, xlim2, sampling_num)
        self.rrms = []

    @staticmethod
    def get_meshgrid(xlim1, xlim2, sampling_num):
        charges1 = linspace(xlim1[0], xlim1[1], num=sampling_num)
        charges2 = linspace(xlim2[0], xlim2[1], num=sampling_num)
        return meshgrid(charges1, charges2)


# Calculation --- set to True if running the 'one charge variation' version.
# Also set to True when performing the calculations for the 'two charges'
# version. Then you must also switch on the other calculation section.
if True:
    os.mkdir(output_path)
    os.mkdir(resp_output_path)

    print("\nRunning unrestrained RESP to fit ESP with equivalence:")
    esp_equiv_molecule = resp.run_resp(
        path, resp_output_path + 'unrest', resp_type='unrest',
        esp_fn=esp_fn)
    # Equivalence alternative charge as well (i.e. unrest RESP on its own grid)
    alt_esp_equiv_molecule = resp.run_resp(
        path, resp_output_path + 'alt_unrest', resp_type='unrest',
        esp_fn=alt_esp_fn, check_ivary=False)

    charges.update_with_charges(esp_charge_type, log_fn, g.molecule)
    # This should actually be called esp_charge_rms
    charge_rms, charge_rrms = rms_and_rep(g.field, g.molecule,
                                          esp_charge_type)[:2]
    resp_rms, resp_rrms = rms_and_rep(g.field, esp_equiv_molecule, 'resp')[:2]
    # Note that, crucially, the equivalenced alternative charges are evaluated
    # on the same grid as the original charges, i.e. `g.field`
    alt_resp_rms, alt_resp_rrms = rms_and_rep(g.field, alt_esp_equiv_molecule,
                                              'resp')[:2]

    print("\nThe molecule with {0} charges:".format(esp_charge_type.upper()))
    print(" RMS: {0:.5f}".format(charge_rms))
    print("RRMS: {0:.5f}".format(charge_rrms))
    print("RMSV: {0:.5f}".format(charge_rms/charge_rrms))
    for atom in g.molecule:
        atom.print_with_charge(esp_charge_type)

    print("\nThe molecule with equivalenced {0} charges (unrestrained RESP):"
          .format(esp_charge_type.upper()))
    print(" RMS: {0:.5f}".format(resp_rms))
    print("RRMS: {0:.5f}".format(resp_rrms))
    print("RMSV: {0:.5f}".format(resp_rms/resp_rrms))
    for atom in esp_equiv_molecule:
        atom.print_with_charge('resp')

    print("\nChecking differences between raw and equivalenced charges ...")
    print(charges.compare_charges(esp_charge_type, 'resp', g.molecule,
          esp_equiv_molecule))

    print("\nThe molecule with equivalenced {0} charges (unrestrained RESP) "
          "evaluated on the {1} grid:".format(alt_esp_charge_type.upper(),
                                              esp_charge_type.upper()))
    print(" RMS: {0:.5f}".format(alt_resp_rms))
    print("RRMS: {0:.5f}".format(alt_resp_rrms))
    print("RMSV: {0:.5f}".format(alt_resp_rms/alt_resp_rrms))
    print("Percentage increase over equivalenced {0} charges: {1:.2f} %"
          .format(esp_charge_type.upper(),
                  100*(alt_resp_rrms-resp_rrms)/resp_rrms))
    for atom in alt_esp_equiv_molecule:
        atom.print_with_charge('resp')

# One charge (e.g. methane, water): 2D plot
if False:
    charges = linspace(xlim1[0], xlim1[1], num=sampling_num)

    result = []
    print("\nOne-dimensional scan:")
    check_ivary = True
    for i, charge in enumerate(charges):
        inp_charges = resp.charges_from_dict(charge_dict(charge),
                                             len(molecule))
        updated_molecule = resp.run_resp(
            path, resp_output_path + "{0}{1:+.3f}".format(
                get_atom_signature(molecule, vary_label1), charge),
            resp_type='dict', inp_charges=inp_charges, esp_fn=esp_fn,
            check_ivary=check_ivary)
        # check_ivary is supposed to be True only on the first run
        if check_ivary:
            check_ivary = False
            print()
        rrms_val = rms_and_rep(g.field, updated_molecule, 'resp')[1]
        result.append(rrms_val)
        sys.stdout.write("\rSampling progress: {0:.2f} %".format(
            100*(i+1)/sampling_num))

    min_charge = esp_equiv_molecule[vary_label1-1].charges['resp']

    plt.title(title)
    plt.xlabel("Charge on " + get_atom_signature(molecule, vary_label1))
    plt.ylabel("RRMS at fitting points")
    plt.plot(charges, result)
    # Guiding line at zero y
    plt.plot((-1.2, 1.2), (0, 0), 'r--')
    # Guiding line at zero x
    plt.plot((0, 0), (0, 1.2*max(result)), 'r--')
    # Location of minimum
    plt.plot((min_charge, min_charge), (0, resp_rrms), 'g--')
    plt.plot((-1.2, min_charge), (resp_rrms, resp_rrms), 'g--')

    axes = plt.gca()
    axes.set_xlim(xlim1)
    axes.set_ylim([0, max(result)])

    save_to = output_path + molecule_name + "_rrms_" + esp_charge_type
    plt.savefig(save_to + ".pdf", format='pdf')
    plt.show()
    plt.close()


def plot_common(x_atom_label, y_atom_label, molecule, title):
    if title is not None:
        plt.title(title)
    plt.ylabel("Charge on " + get_atom_signature(molecule, x_atom_label))
    plt.xlabel("Charge on " + get_atom_signature(molecule, y_atom_label))

# 2 charges (NMe3H_plus, NMe4_plus etc.): Calculation. Note that the previous
# calculation section must also be switched on.
if True:
    print("\nTwo-dimensional scan:")
    new_result = Result(sampling_num, xlim1, xlim2)
    i = 0
    check_ivary = True
    print("\nEvaluating the meshgrid of H-only RESPs.")
    interpret(g.molecule, charge_dict, vary_label1, vary_label2)
    for c, n in zip(new_result.inp1.flat, new_result.inp2.flat):
        inp_charges = resp.charges_from_dict(charge_dict(c, n), len(molecule))
        updated_molecule = resp.run_resp(
            path, resp_output_path + "{0}{1:+.3f}-{2}{3:+.3f}".format(
                get_atom_signature(molecule, vary_label1), c,
                get_atom_signature(molecule, vary_label2), n),
            resp_type='dict', inp_charges=inp_charges, esp_fn=esp_fn,
            check_ivary=check_ivary)

        # check_ivary is supposed to be True only on the first run
        if check_ivary:
            check_ivary = False
            print()

        rrms_val = rms_and_rep(g.field, updated_molecule, 'resp')[1]
        new_result.rrms.append(rrms_val)

        i += 1
        sys.stdout.write("\rMeshgrid progress: {0:.2f} %".format(
            100*i/sampling_num**2))
        sys.stdout.flush()

    new_result.rrms = np.array(new_result.rrms)
    new_result.rrms.resize([sampling_num, sampling_num])
    new_result.resp_rms = resp_rms
    new_result.resp_rrms = resp_rrms
    new_result.esp_equiv_molecule = esp_equiv_molecule
    new_result.alt_esp_equiv_molecule = alt_esp_equiv_molecule

    with open(output_path + "result.p", "wb") as f:
        pickle.dump(new_result, f)

# 2 charges: Presentation
if True:
    with open(output_path + "result.p", "rb") as f:
        read_result = pickle.load(f)

    rel_rrms = [100*(elem-read_result.resp_rrms)/read_result.resp_rrms for
                elem in read_result.rrms]
    rel_rrms = np.array(rel_rrms)
    rel_rrms.resize([read_result.sampling_num, read_result.sampling_num])

    # Non-ESP charge and its minimized ratio
    charges.update_with_charges(charge_type, charge_log_fn, molecule)
    equiv_start_charges = resp.equivalence(molecule, charge_type, path)
    charge_type += '_equiv'
    charges._update_molecule_with_charges(molecule, equiv_start_charges,
                                          charge_type)
    os.mkdir(temp_output_path)
    # Scan roughly various ratios to find bracket for minimization
    print("\nScanning roughly various ratios. This shouldn't take long.")
    heavy_args = (g.field, path, temp_output_path, esp_fn, False)
    heavy_result, indicator_charge, ratio_values = resp.eval_ratios(
        'heavy', (0, 2), equiv_start_charges, 10, vary_label2, heavy_args,
        first_verbose=True)
    # Minimization
    print("\nStarting minimization of charge ratio.")
    heavy_args = (equiv_start_charges, g.field, path, temp_output_path,
                  esp_fn, True)  # True for optimization
    heavy_min_ratio, heavy_min_ratio_rrms = resp.minimize_ratio(
        'heavy', ratio_values, heavy_result, heavy_args)
    shutil.rmtree(temp_output_path)

    # Presentation: 3D
    if False:
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        surf = ax.plot_surface(read_result.inp2, read_result.inp1,
                               read_result.rrms, cmap=cm.plasma, rstride=1,
                               cstride=1)
        fig.colorbar(surf, label="RRMS at fitting points")
        plot_common(vary_label2, vary_label1, molecule, title)
        plt.show()
        plt.close()

    # Presentation: 2D contour
    if True:
        levels = [1, 5, 10, 20, 30, 50, 100]
        CS = plt.contour(read_result.inp2, read_result.inp1, rel_rrms,
                         levels, rstride=1, ctride=1, inline=1, colors='k',
                         zorder=1)

        print("Number of charges sampled along each axis:", sampling_num)
        print("Reporting contour surface areas:")
        for contour_isovalue in levels:
            vertices = contour_vertices(contour_isovalue, CS, levels)
            print("{0:>3}% contour area: {1:.5f}".format(
                contour_isovalue, polygon_area(vertices)))

        plt.clabel(CS, fmt="%1.0f", inline=1, fontsize=15, colors='b')
        axes = plt.gca()

        axes.set_xlim(xlim2)
        axes.set_ylim(xlim1)
        plt.axes().set_aspect('equal')

        # Add non-esp point
        new_point = (molecule[vary_label2-1].charges[charge_type],
                     molecule[vary_label1-1].charges[charge_type])
        plt.scatter(*new_point, zorder=2, s=50, lw=0.3)
        # Non-esp ratio charge point
        plt.scatter(heavy_min_ratio*equiv_start_charges[vary_label2-1],
                    heavy_min_ratio*equiv_start_charges[vary_label1-1],
                    marker='D', zorder=2, s=50, lw=0.3)
        # Add ratio line
        y_coord = axes.get_xlim()[0]*new_point[1]/new_point[0]
        plt.plot((axes.get_xlim()[0], 0), (y_coord, 0), 'r--', zorder=1)

        plt.scatter(
            read_result.esp_equiv_molecule[vary_label2-1].charges['resp'],
            read_result.esp_equiv_molecule[vary_label1-1].charges['resp'],
            marker='x', zorder=2, s=50, lw=0.8)

        plt.scatter(
            read_result.alt_esp_equiv_molecule[vary_label2-1].charges['resp'],
            read_result.alt_esp_equiv_molecule[vary_label1-1].charges['resp'],
            marker='+', zorder=2, s=50, lw=0.8)

        plot_common(vary_label2, vary_label1, molecule, title)
        save_to = output_path + molecule_name + "_" + esp_charge_type
        plt.savefig(save_to + ".pdf", format='pdf')
        plt.close()
