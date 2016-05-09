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
# The first label is the more important atom, for which a 1D scan will be done
# NMe3H_plus, NMe3:
vary_label1 = 13
charge_dict_1D = lambda a1: {13: a1}
vary_label2 = 1
charge_dict_2D = lambda a1, a2: {1: a2, 5: a2, 9: a2, 13: a1}
# This shouldn't be more than 5 because that's the max number of colors in plot
labels_to_monitor = [1, 2, 14]
# NMe4_plus:
# vary_label1 = 17
# charge_dict_1D = lambda a1: {17: a1}
# charge_dict_2D = lambda a2, a1: {1: a2, 5: a2, 9: a2, 13: a2, 17: a1}
xlim1 = (-0.5, 1)
xlim2 = (-1, 0)
sampling_num = 21
swap_axes = False

# SLICING plots
# Note that slicing plots in planes different than those of the coordinate
# system do not preserve distances. TODO
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
opt_output_path = output_path + 'opt/'

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
    if vary_label2 is None:
        dictio = charge_dict(1)  # Example number to get the dict
    else:
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


# CALCULATIONS COMMON TO 1 AND 2D VARIATIONS
if True:
    os.mkdir(output_path)
    os.mkdir(resp_output_path)
    levels = [1, 5, 10, 20, 30, 50, 100]

    print("\nRunning unrestrained RESP to fit ESP with equivalence:")
    esp_equiv_molecule = resp.run_resp(
        path, resp_output_path + 'unrest', resp_type='unrest', esp_fn=esp_fn)
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


# ONE CHARGE VARIATION
if True:
    charge_vals = linspace(xlim1[0], xlim1[1], num=sampling_num)
    result = []
    all_charges_result = []
    print("\nOne-dimensional scan:")
    check_ivary = True
    interpret(g.molecule, charge_dict_1D, vary_label1)
    resp_args = [g.field, path, resp_output_path, esp_fn, molecule,
                 vary_label1, charge_dict_1D]
    for i, charge in enumerate(charge_vals):
        rrms_val, all_charges = resp.eval_one_charge_resp(charge, *resp_args,
                                                          check_ivary)
        result.append(rrms_val)
        all_charges_result.append(all_charges)
        # check_ivary is supposed to be True only on the first run
        if check_ivary:
            check_ivary = False
            print()
        sys.stdout.write("\rSampling progress: {0:.2f} %".format(
            100*(i+1)/sampling_num))

    min_charge = esp_equiv_molecule[vary_label1-1].charges['resp']

    os.mkdir(opt_output_path)
    resp_args[2] = opt_output_path
    flex_limits = []
    print("\n\nFlexibility limits on", molecule[vary_label1-1])
    for level in levels:
        if (1+level/100)*resp_rrms > min(result[0], result[1]):
            # Otherwise find_flex would throw an error
            print("Further limits were beyond the sampled range.")
            break

        sol1, sol2, charges1, charges2 = resp.find_flex(
            (1+level/100)*resp_rrms, charge_vals, result, resp_args)
        flex_limits.append([sol1, sol2])
        if level == 10:
            charges1_at_10 = charges1
            charges2_at_10 = charges2
        # Difference also shown as percentage of charge on that atom
        print("{0:>3}% limits: {1: .5f}, {2: .5f}, diff: {3:.5f} ({4:.1f}%)"
              .format(level, sol1, sol2, sol2-sol1,
                      100*abs((sol2-sol1)/min_charge)))

    print("\nThe monitored charges at 10% flexibility limits:")
    i = 1
    for charge1, charge2 in zip(charges1_at_10, charges2_at_10):
        if i in labels_to_monitor:
            print("{0:>3}: {1: .5f}, {2: .5f}, diff: {3:.5f}".format(
                get_atom_signature(molecule, i), charge1, charge2,
                abs(charge1 - charge2)))
        i += 1

    shutil.rmtree(opt_output_path)

    fig, ax1 = plt.subplots()
    ax1.set_xlabel("Charge on " + get_atom_signature(molecule, vary_label1))
    ax1.plot(charge_vals, result)
    # Make the y-axis label and tick labels match the line color
    ax1.set_ylabel("RRMS", color='b')
    for tl in ax1.get_yticklabels():
        tl.set_color('b')

    ax2 = ax1.twinx()
    colors = iter(['g', 'r', 'c', 'm', 'y'])
    # Plot charges on other atoms
    ax2.set_ylabel("Charges on other atoms (see legend)")
    for i, atom_charge in enumerate(zip(*all_charges_result)):
        if i+1 in labels_to_monitor:
            ax2.plot(charge_vals, atom_charge, color=next(colors),
                     label=get_atom_signature(molecule, i+1))
    ax2.legend()
    # Guiding line at zero y2
    ax2.plot((-1.2, 1.2), (0, 0), 'k:')
    # Guiding line at zero x
    ax1.plot((0, 0), (0, 1.2*max(result)), 'k:')
    # Location of minimum
    ax1.plot((min_charge, min_charge), (0, resp_rrms), 'k--')
    ax1.plot((-1.2, min_charge), (resp_rrms, resp_rrms), 'k--')
    # 10% flexibility limits
    ax1.axvspan(flex_limits[2][0], flex_limits[2][1], alpha=0.2, color='grey')

    ax1.set_xlim(xlim1)
    ax1.set_ylim([0, max(result)])

    save_to = output_path + molecule_name + "_" + esp_charge_type + "_1D"
    plt.savefig(save_to + ".pdf", format='pdf')
    plt.show()
    plt.close()


def plot_common(x_atom_label, y_atom_label, molecule, title):
    if title is not None:
        plt.title(title)
    plt.xlabel("Charge on " + get_atom_signature(molecule, x_atom_label))
    plt.ylabel("Charge on " + get_atom_signature(molecule, y_atom_label))


# TWO CHARGE VARIATION (switch off for molecules with less than 2 heavy atoms
# or when this is not of interest)
if True:
    print("\nTwo-dimensional scan:")
    new_result = Result(sampling_num, xlim1, xlim2)
    i = 0
    check_ivary = True
    interpret(g.molecule, charge_dict_2D, vary_label1, vary_label2)
    for a1, a2 in zip(new_result.inp1.flat, new_result.inp2.flat):
        inp_charges = resp.charges_from_dict(charge_dict_2D(a1, a2),
                                             len(molecule))
        updated_molecule = resp.run_resp(
            path, resp_output_path + "{0}{1:+.3f}-{2}{3:+.3f}".format(
                get_atom_signature(molecule, vary_label1), a1,
                get_atom_signature(molecule, vary_label2), a2),
            resp_type='dict', inp_charges=inp_charges, esp_fn=esp_fn,
            check_ivary=check_ivary)
        rrms_val = rms_and_rep(g.field, updated_molecule, 'resp')[1]
        new_result.rrms.append(rrms_val)
        # check_ivary is supposed to be True only on the first run
        if check_ivary:
            check_ivary = False
            print()
        # Show progress
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
    os.mkdir(opt_output_path)
    # Scan roughly various ratios to find bracket for minimization
    print("\nScanning roughly various ratios. This shouldn't take long.")
    heavy_args = (g.field, path, opt_output_path, esp_fn, False)
    heavy_result, indicator_charge, ratio_values = resp.eval_ratios(
        'heavy', (0, 2), equiv_start_charges, 10, vary_label1, heavy_args,
        first_verbose=True)
    # Minimization
    print("\nStarting minimization of charge ratio.")
    heavy_args = (equiv_start_charges, g.field, path, opt_output_path,
                  esp_fn, True)  # True for optimization
    heavy_min_ratio, heavy_min_ratio_rrms = resp.minimize_ratio(
        'heavy', ratio_values, heavy_result, heavy_args)
    shutil.rmtree(opt_output_path)

    # Presentation: 3D
    if False:
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        surf = ax.plot_surface(read_result.inp1, read_result.inp2,
                               read_result.rrms, cmap=cm.plasma, rstride=1,
                               cstride=1)
        fig.colorbar(surf, label="RRMS at fitting points")
        plot_common(vary_label1, vary_label2, molecule, title)
        plt.show()
        plt.close()

    # Presentation: 2D contour
    if True:
        inp1 = read_result.inp1
        inp2 = read_result.inp2
        if swap_axes:
            xlim1, xlim2 = xlim2, xlim1
            vary_label1, vary_label2 = vary_label2, vary_label1
            inp1, inp2 = inp2, inp1

        # NOTE: the x and y axes may need to be swapped for these plots to look
        # better if desired. These have been designed with NMe3H in mind.
        CS = plt.contour(inp1, inp2, rel_rrms, levels, rstride=1, ctride=1,
                         inline=1, colors='k', zorder=1)

        print("Number of charges sampled along each axis:", sampling_num)
        print("Reporting contour surface areas:")
        for contour_isovalue in levels:
            vertices = contour_vertices(contour_isovalue, CS, levels)
            print("{0:>3}% contour area: {1:.5f}".format(
                contour_isovalue, polygon_area(vertices)))

        plt.clabel(CS, fmt="%1.0f", inline=1, fontsize=15, colors='b')
        axes = plt.gca()

        axes.set_xlim(xlim1)
        axes.set_ylim(xlim2)
        plt.axes().set_aspect('equal')

        # Add non-esp point
        new_point = (molecule[vary_label1-1].charges[charge_type],
                     molecule[vary_label2-1].charges[charge_type])
        plt.scatter(*new_point, zorder=2, s=50, lw=0.3)
        # Non-esp ratio charge point
        plt.scatter(heavy_min_ratio*equiv_start_charges[vary_label1-1],
                    heavy_min_ratio*equiv_start_charges[vary_label2-1],
                    marker='D', zorder=2, s=50, lw=0.3)
        # Add ratio line
        y_coord = axes.get_xlim()[0]*new_point[1]/new_point[0]
        plt.plot((axes.get_xlim()[0], 0), (y_coord, 0), 'r--', zorder=1)

        plt.scatter(
            read_result.esp_equiv_molecule[vary_label1-1].charges['resp'],
            read_result.esp_equiv_molecule[vary_label2-1].charges['resp'],
            marker='x', zorder=2, s=50, lw=0.8)

        plt.scatter(
            read_result.alt_esp_equiv_molecule[vary_label1-1].charges['resp'],
            read_result.alt_esp_equiv_molecule[vary_label2-1].charges['resp'],
            marker='+', zorder=2, s=50, lw=0.8)

        plot_common(vary_label1, vary_label2, molecule, None)
        save_to = output_path + molecule_name + "_" + esp_charge_type + '_2D'
        plt.savefig(save_to + ".pdf", format='pdf')
        plt.close()
