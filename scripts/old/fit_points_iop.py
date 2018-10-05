from repESP import resp_helpers, charges, graphs
from repESP.field_comparison import rms_and_rep, difference
from repESP.esp_fit_calc import FitCalc, IOpCalcSet

import shutil
import os
import copy
import pickle
import numpy as np
import matplotlib.pyplot as plt

# The purpose of the script is investigating the change in generated fitting
# points upon varying IOp 41-43 settings. The first part of this script creates
# input files and pickles filenames. The user should then run all .com files
# with Gaussian and run the second part of this script, which produces graphs.


def calc_min_max(alist):
    return min(alist), max(alist)


def smart_range(min_max, low=0.9, high=1.1):
    # Something like this:
    # color_span = [0.8*min_max[0], 1.2*min_max[1]]
    # but allowing for different signs
    color_span = []
    for elem, sign in zip(min_max, [-1, 1]):
        if np.sign(elem) == sign:
            color_span.append(high*elem)
        else:
            color_span.append(low*elem)
    return color_span


def plot_range(alist, margin=0.05):
    min_max = calc_min_max(alist)
    diff = abs(min_max[1] - min_max[0])
    return min_max[0] - margin*diff, min_max[1] + margin*diff


def calc_plot(calcs, to_plot, title, set_lim=False, save_to=None):
    plt.scatter(list(range(len(calcs))), to_plot)
    plt.xticks(list(range(len(calcs))), [elem[-6:] for elem in calcs],
               rotation='vertical')
    plt.title(title)
    axes = plt.gca()
    if set_lim:
        axes.set_ylim(plot_range(to_plot))
    graphs._save_or_display(save_to)


def check_color_span(values, color_span, default=None):
    min_max = calc_min_max(values)
    if color_span == []:
        if default is None:
            color_span = smart_range(min_max)
        else:
            color_span = default
    if min_max[0] < color_span[0] or min_max[1] > color_span[1]:
        # The extent of this is unlikely to be large, since both MK and
        # CHelpG use a fixed exclusion radius (or do they? That's up to
        # Gaussian's implementation and is to be investigated).
        print("WARNING: Values on graph (min_max = {0:.4f}, {1:.4f}) will "
              "be outside of color scale ({2:.4f}, {3:.4f})".format(
                  *min_max, *color_span))
    return color_span


charge_type = 'mk'
path = '../data/methane/fit_points-1/'

# PART 1
if True:
    os.mkdir(path)
    shutil.copy('../data/methane/input/methane.chk', path)

    # A simple investigation --- varying only IOp42
    calc_set = IOpCalcSet(iop42=list(range(1, 6)))
    calcs = []
    print(calc_set.create_param_list())
    for iop41, iop42, iop43 in zip(*calc_set.create_param_list()):
        calc = FitCalc(path, 'methane', 'MP2/6-31+G(d,p)', charge_type, 0, 1,
                       iop41, iop42, iop43)
        calc.create_input()
        calcs.append(calc)
        print("Created file: ", calc.filename)

    # Pickle filenames for part 2
    with open(path + "fit_points.p", 'wb') as f:
        pickle.dump([elem.filename for elem in calcs], f)

# PART 2 --- run when the Gaussian calculations have been completed
if False:
    with open(path + "fit_points.p", 'rb') as f:
        calcs = pickle.load(f)

    rms_list = []
    charges_dict = {}

    color_span = []
    error_color_span = []
    for calc in calcs:
        g = resp_helpers.G09_esp(path + calc + '.esp')
        charges.update_with_charges(charge_type, path + calc + '.log',
                                    g.molecule)
        with open(path + calc + "-charges.txt", "a") as fc:
            for atom in g.molecule:
                atom.print_with_charge(charge_type, fc)
                if atom.label in charges_dict:
                    charges_dict[atom.label].append(atom.charges[charge_type])
                else:
                    charges_dict[atom.label] = [atom.charges[charge_type]]

            min_rms, min_rrms, rep_esp_field = rms_and_rep(g.field, g.molecule,
                                                           charge_type)
            rms_list.append(min_rms)
            print("\n", min_rms, file=fc)

        # Default given as extremal values of methane CHelpG
        color_span = check_color_span(g.field.values, color_span,
                                      default=[-0.0045, 0.011])
        diff_field = difference(g.field, rep_esp_field)
        error_color_span = check_color_span(diff_field.values,
                                            error_color_span,
                                            default=[-0.0012, 0.019])

        graphs.plot_points(
            g.field, 2, title=calc, molecule=g.molecule,
            plane_eqn=graphs.plane_through_atoms(g.molecule, 1, 2, 3),
            dist_thresh=0.5, axes_limits=[(-5, 5)]*2, color_span=color_span,
            save_to=path + calc[-6:] + '_V.pdf')

        graphs.plot_points(
            diff_field, 2, title=calc + " Errors", molecule=g.molecule,
            plane_eqn=graphs.plane_through_atoms(g.molecule, 1, 2, 3),
            dist_thresh=0.5, axes_limits=[(-5, 5)]*2,
            color_span=error_color_span, save_to=path + calc[-6:] + '_E.pdf')

    save_to = path + "RMS.pdf"
    calc_plot(calcs, rms_list, charge_type.upper() + " RMS value",
              set_lim=True, save_to=save_to)

    for atom in g.molecule:
        save_to = path + atom.atomic_number + str(atom.label) + "_charge.pdf"
        title = "Charge on " + atom.atomic_number + str(atom.label)
        calc_plot(calcs, charges_dict[atom.label], title, save_to=save_to)
