from repESP import resp, rep_esp, charges
from repESP.field_comparison import _check_grids, difference

from numpy import mean, sqrt, square, linspace, meshgrid
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

import numpy as np
import pickle

charge_type = 'mk'
# charge_type = 'chelpg'
# charge_type = 'hly'

# molecule_name = 'methane'
molecule_name = 'tma'
path = '../data/' + molecule_name + '/'
input_path = path + 'input/'

common_fn = input_path + molecule_name + "_" + charge_type
log_fn = common_fn + ".log"
input_esp = common_fn + "_resp.esp"
output_esp = common_fn + "_resp_reformatted.esp"

print("To see a demonstration of all the capabilities of the script, change "
      "the hard-coded conditional values to True. You can also change the "
      "charges type between MK and CHelp(G).")
print("\nMolecule:    ", molecule_name.capitalize())
print("Charge type: ", charge_type.upper(), '\n')

g = resp.G09_esp(input_esp)

# Write the .esp file in the correct format expected by the `resp` program
if False:
    g.field.write_to_file(output_esp, g.molecule)

charges.update_with_charges(charge_type, log_fn, g.molecule)
for atom in g.molecule:
    atom.print_with_charge(charge_type)


# Reproduce ESP values at those points
rep_esp_field = rep_esp.calc_non_grid_field(g.molecule, g.field.points,
                                            'rep_esp', [charge_type])[0]
# Division into Voronoi basins:
# parent_atom, dist = rep_esp.calc_non_grid_field(g.molecule, g.field.points,
#                                                 'dist')


diff = difference(g.field, rep_esp_field).values
rel_diff = difference(g.field, rep_esp_field, relative=True).values

# RMS value -- can compare with than in log file
min_rms = sqrt(mean(square(diff)))
print("\nRMS: {0:.6f}".format(min_rms))
# Trying to reverese-engineer RRMS
if False:
    rrms_by_mean = sqrt(mean(square(diff)))/mean(g.field.values)
    print("\nRRMS by mean:            {0:6f}".format(rrms_by_mean))
    mean_val = mean([abs(elem) for elem in g.field.values])
    rrms_by_mean_val = sqrt(mean(square(diff)))/mean_val
    print("RRMS by mean value:      {0:6f}".format(rrms_by_mean_val))
    rrms_by_rel_err = sqrt(mean(square(rel_diff)))
    print("RRMS by relative error:  {0:6f}".format(rrms_by_rel_err))


# Plot the grid in 3D
if False:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    color = g.field.values
    cmap = plt.get_cmap('plasma')

    image = ax.scatter(*list(zip(*g.field.points)), c=color, cmap=cmap)
    cbar = fig.colorbar(image, label="ESP value")

    plt.show()
    plt.close()

# Calculate RMS for various charges
import copy
from itertools import chain
molecule = copy.deepcopy(g.molecule)
# Could also do a 2D plot for tetramethylammonium.


def change_charges(charge, molecule, vary_atom_label):
    # Assumes neutral molecule
    charge_other = -charge/(len(molecule)-1)
    for i, atom in enumerate(molecule):
        if i+1 == vary_atom_label:
            atom.charges['temp'] = charge
        else:
            atom.charges['temp'] = charge_other


def change_charges2(net_charge, charge_dict, molecule):
    # e.g. {17: n, 1: c, 5: c, 9: c, 13: c}
    described_labels = list(charge_dict.keys())
    inferred_labels = list(set(range(1, len(molecule)+1))
                           - set(described_labels))

    left_charge = net_charge
    for label in charge_dict:
        left_charge -= charge_dict[label]
    inferred_charge = left_charge/len(inferred_labels)

    for i, atom in enumerate(molecule):
        if i+1 in described_labels:
            atom.charges['temp'] = charge_dict[i+1]
        else:
            atom.charges['temp'] = inferred_charge


def calc_rms(charge, molecule, vary_atom_label=1):
    # Reproduce ESP values at those points
    change_charges(charge, molecule, vary_atom_label)
    rep_esp_field = rep_esp.calc_non_grid_field(molecule, g.field.points,
                                                'rep_esp', ['temp'])[0]
    diff = difference(g.field, rep_esp_field).values
    rms = sqrt(mean(square(diff)))
    return rms


def calc_rms2(net_charge, charge_dict, molecule):
    # Reproduce ESP values at those points
    change_charges2(net_charge, charge_dict, molecule)
    rep_esp_field = rep_esp.calc_non_grid_field(molecule, g.field.points,
                                                'rep_esp', ['temp'])[0]
    diff = difference(g.field, rep_esp_field).values
    rms = sqrt(mean(square(diff)))
    return rms

# One charge: 2D plot
if False:
    vary_atom_label = 1
    xlim = (-1, 0.5)
    num = 150
    charges = linspace(xlim[0], xlim[1], num=num)
    result = []
    for i, charge in enumerate(charges):
        if not i % 10:
            print("{0:.2f}%".format(100*i/num))
        result.append(calc_rms(charge, molecule, vary_atom_label))

    min_charge = molecule[0].charges[charge_type]

    plt.title(molecule_name.capitalize() + " " + charge_type.upper())
    plt.xlabel("Charge on " + molecule[vary_atom_label-1].identity)
    plt.ylabel("RMSE at fitting points")
    plt.plot(charges, result)
    plt.plot((-1.2, 1.2), (0, 0), 'r--')
    plt.plot((0, 0), (0, 1.2*max(result)), 'r--')
    plt.plot((min_charge, min_charge), (0, min_rms), 'g--')
    plt.plot((-1.2, min_charge), (min_rms, min_rms), 'g--')

    axes = plt.gca()
    axes.set_xlim(xlim)
    axes.set_ylim([0, max(result)])

    save_to = path + molecule_name + "_rms_" + charge_type
    plt.savefig(save_to + ".pdf", format='pdf')
    plt.show()
    plt.close()


def plot_common():
    plt.title(molecule_name.capitalize() + " " + charge_type.upper())
    plt.xlabel("Charge on N")
    plt.ylabel("Charge on C")


class Result(object):

    def __init__(self, num, c_xlim, n_xlim):
        self.num = num
        self.c_xlim = c_xlim
        self.n_xlim = n_xlim
        self.c_inp, self.n_inp = self.get_meshgrid(c_xlim, n_xlim, num)
        self.rms = []

    @staticmethod
    def get_meshgrid(c_xlim, n_xlim, num):
        c_charges = linspace(c_xlim[0], c_xlim[1], num=num)
        n_charges = linspace(n_xlim[0], n_xlim[1], num=num)
        return meshgrid(c_charges, n_charges)

# 2 charges: Calculation
if True:
    # Change input values here
    net_charge = 1
    charges = lambda n, c: {17: n, 1: c, 5: c, 9: c, 13: c}
    c_xlim = (-1, 0.5)
    n_xlim = (-0.5, 1)
    num = 51
    new_result = Result(num, c_xlim, n_xlim)

    i = 0
    for n, c in zip(new_result.n_inp.flat, new_result.c_inp.flat):
        new_result.rms.append(calc_rms2(1, charges(n, c), molecule))
        for atom in molecule:
            atom.print_with_charge('temp')
        i += 1
        print("{0:.2f} %".format(100*i/num**2))

    new_result.rms = np.array(new_result.rms)
    new_result.rms.resize([num, num])

    with open(molecule_name + "_" + charge_type + "_result.p", "wb") as f:
        pickle.dump(new_result, f)

# 2 charges: Presentation
if True:
    with open(molecule_name + "_" + charge_type + "_result.p", "rb") as f:
        read_result = pickle.load(f)

    rel_rms = [100*(elem-min_rms)/min_rms for elem in read_result.rms]
    rel_rms = np.array(rel_rms)
    rel_rms.resize([read_result.num, read_result.num])

    # Presentation: 3D
    if False:
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        surf = ax.plot_surface(read_result.n_inp, read_result.c_inp,
                               read_result.rms, cmap=cm.plasma, rstride=1,
                               cstride=1)
        fig.colorbar(surf, label="RMSE at fitting points")
        plot_common()
        plt.show()
        plt.close()

    # Presentation: 2D contour
    if True:
        levels = [1, 5, 10, 20, 30, 50, 100]
        CS = plt.contour(read_result.n_inp, read_result.c_inp, rel_rms, levels,
                         rstride=1, ctride=1, inline=1, colors='k')
        plt.clabel(CS, fmt="%1.0f", inline=1, fontsize=10, colors='b')
        axes = plt.gca()

        # Change axes here:
        axes.set_ylim([-1, 0])
        plt.axes().set_aspect('equal')

        plot_common()
        save_to = molecule_name + "_" + charge_type
        plt.savefig(save_to + ".pdf", format='pdf')
        plt.show()
        plt.close()
