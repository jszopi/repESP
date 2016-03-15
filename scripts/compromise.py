from repESP import resp, charges
from repESP.field_comparison import rms_and_rep

import numpy as np
import os
import matplotlib.pyplot as plt

# This was necessary to prevent title from being cut-off when it's shifted up
# due to the second x-axis label.
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

esp_charge_type = 'mk'
# esp_charge_type = 'chelpg'

charge_type = 'nbo'
# charge_type = 'aim'

molecule_name = 'methane'
indicator_label = 1

path = '../data/methane/input/'
output_path = path + "compromise" + '_' + charge_type + '/'
os.mkdir(output_path)

zero_net_charge = True
if "plus" in molecule_name or "minus" in molecule_name:
    zero_net_charge = False

print("\nThe molecule was found {0}to be neutral based on its name. You should"
      " check if this is correct.".format("" if zero_net_charge else "NOT "))

resp_output_path = output_path + 'resp_calcs/'
os.mkdir(resp_output_path)

log_fn = path + molecule_name + "_" + charge_type + ".log"
esp_log_fn = path + molecule_name + "_" + esp_charge_type + ".log"
g = resp.G09_esp(path + molecule_name + '_' + esp_charge_type + '.esp')

charges.update_with_charges(esp_charge_type, esp_log_fn, g.molecule)
charges.update_with_charges(charge_type, log_fn, g.molecule)

print("\nThe molecule with {0} charges:".format(charge_type.upper()))
for atom in g.molecule:
    atom.print_with_charge(charge_type)

print("\nThe molecule with {0} charges:".format(esp_charge_type.upper()))
for atom in g.molecule:
    atom.print_with_charge(esp_charge_type)

start_charges = [atom.charges[charge_type] for atom in g.molecule]

min_rms, min_rrms, rep_esp_field = rms_and_rep(g.field, g.molecule,
                                               esp_charge_type)

num = 50
heavy_result = []
result = []
ratio_limits = (0, 1.5)
ratio_values = np.linspace(*ratio_limits, num=num)
indicator_charge = []

for ratio in ratio_values:
    inp_charges = [charge*ratio for charge in start_charges]
    indicator_charge.append(inp_charges[indicator_label-1])

    updated_molecule = resp.run_resp(
        path, resp_output_path + "ratio{0:+.3f}".format(ratio),
        resp_type='h_only', inp_charges=inp_charges, esp_fn=molecule_name + "_"
        + esp_charge_type + '.esp')
    rrms_val = rms_and_rep(g.field, updated_molecule, 'resp')[1]
    heavy_result.append(rrms_val)

    print("\nHEAVY: RATIO: {0:.3f}, RRMS: {1:.3f}".format(ratio, rrms_val))
    for atom in updated_molecule:
        atom.print_with_charge('resp')

    if zero_net_charge:
        # Scaling all charges is only possible with neutral molecules as
        # otherwise in this case there's no free hydrogens to compensate as in
        # the 'heavy_only' version
        charges._update_molecule_with_charges(g.molecule, inp_charges, 'temp')
        rrms_val = rms_and_rep(g.field, g.molecule, 'temp')[1]
        result.append(rrms_val)

        print("\nREGULAR: RATIO: {0:.3f}, RRMS: {1:.3f}".format(ratio,
                                                                rrms_val))
        for atom in g.molecule:
            atom.print_with_charge('temp')


def plot(result_list, title):
    fig, ax1 = plt.subplots()
    ax1.plot(ratio_values, result_list)
    ax2 = ax1.twiny()
    ax2.plot(indicator_charge, result_list)
    if start_charges[indicator_label-1] < 0:
        ax2.invert_xaxis()

    ax1.set_xlabel(charge_type.upper() + " ratio")
    ax2.set_xlabel("Charge on " + g.molecule[indicator_label-1].identity +
                   str(indicator_label))
    ax1.set_ylabel("RRMS")

    ax1.set_ylim(0, ax1.get_ylim()[1])

    ax1.set_xlim(*ratio_limits)
    ax2.set_xlim(indicator_charge[0], indicator_charge[-1])
    # The lines should overlap so one of the lines can be removed:
    ax2.lines[0].remove()
    # NOTE: if the plots don't look right, try disabling the above option and
    # see whether you get two different lines. However, hard-coding the x-axis
    # limits should ensure that the lines do overlap.
    ax1.plot((1, 1), (0, ax1.get_ylim()[1]), 'g--')
    ax1.plot(ratio_limits, (min_rrms, min_rrms), 'r--')

    plt.title(title, y=1.15)
    plt.show()
    plt.close()

title = "RRMS values for various {0} charge ratios".format(charge_type.upper())
if zero_net_charge:
    plot(result, title)

title += "\nset on heavy atoms ONLY (H free to improve fit)"
plot(heavy_result, title)
