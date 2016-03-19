from repESP import resp, resp_helpers, graphs
from repESP.field_comparison import rms_and_rep
from repESP.charges import update_with_charges, _update_molecule_with_charges

import numpy as np
import os
import matplotlib.pyplot as plt
import math

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
esp_fn = molecule_name + "_" + esp_charge_type + '.esp'

zero_net_charge = True
if "plus" in molecule_name or "minus" in molecule_name:
    zero_net_charge = False

print("\nThe molecule was found {0}to be neutral based on its name. You should"
      " check if this is correct.".format("" if zero_net_charge else "NOT "))

resp_output_path = output_path + 'resp_calcs/'
os.mkdir(resp_output_path)

log_fn = path + molecule_name + "_" + charge_type + ".log"
esp_log_fn = path + molecule_name + "_" + esp_charge_type + ".log"
g = resp_helpers.G09_esp(path + esp_fn)

# Both the Gaussian ESP fitting methods and other charge assignment methods may
# not yield equivalent charges. As equivalent charges make more sense for force
# field development, they will be used. The ESP charges are equivalenced by
# performing unrestrained RESP, which will be used as a reference for the fit
# minimum. Charges from the other method will be equivalenced manually by my
# averaging function `resp.equivalence`. They will be scaled to obtain
# different ratio charges. All the charges are calculated and printed at the
# start for reference.
update_with_charges(esp_charge_type, esp_log_fn, g.molecule)
update_with_charges(charge_type, log_fn, g.molecule)
equiv_charges = resp.equivalence(g.molecule, charge_type, path)
_update_molecule_with_charges(g.molecule, equiv_charges, charge_type+'_equiv')
print("\nNOTE: Running unrestrained RESP to fit ESP with equivalence:")
esp_equiv_molecule = resp.run_resp(
    path, resp_output_path + 'unrest', resp_type='unrest', esp_fn=esp_fn)

charge_rrms = rms_and_rep(g.field, g.molecule, charge_type)[1]
equiv_charge_rrms = rms_and_rep(g.field, g.molecule, charge_type + '_equiv')[1]
esp_charge_rrms = rms_and_rep(g.field, g.molecule, esp_charge_type)[1]
resp_charge_rrms = rms_and_rep(g.field, esp_equiv_molecule, 'resp')[1]

print("\nThe molecule with {0} charges:".format(charge_type.upper()))
print("RRMS: {0:.5f}".format(charge_rrms))
for atom in g.molecule:
    atom.print_with_charge(charge_type)

print("\nThe molecule with equivalenced {0} charges:".format(
      charge_type.upper()))
print("RRMS: {0:.5f}".format(equiv_charge_rrms))
for atom in g.molecule:
    atom.print_with_charge(charge_type + '_equiv')

print("\nThe molecule with {0} charges:".format(esp_charge_type.upper()))
print("RRMS: {0:.5f}".format(esp_charge_rrms))
for atom in g.molecule:
    atom.print_with_charge(esp_charge_type)

print("\nThe molecule with equivalenced {0} charges (unrestrained RESP):"
      .format(esp_charge_type.upper()))
print("RRMS: {0:.5f}".format(resp_charge_rrms))
for atom in esp_equiv_molecule:
    atom.print_with_charge('resp')

start_charges = [atom.charges[charge_type + '_equiv'] for atom in g.molecule]

num = 50
heavy_result = []
result = []
ratio_limits = (0, 1.5)
ratio_values = np.linspace(*ratio_limits, num=num)
indicator_charge = []

for ratio in ratio_values:
    indicator_charge.append(ratio*start_charges[indicator_label-1])

    rrms_val = resp.eval_heavy_ratio(ratio, start_charges, g.field, path,
                                     resp_output_path, esp_fn, verbose=True)
    heavy_result.append(rrms_val)

    if zero_net_charge:
        # Scaling all charges is only possible with neutral molecules as
        # otherwise in this case there's no free hydrogens to compensate as in
        # the 'heavy_only' version
        rrms_val = resp.eval_ratio(ratio, start_charges, g.molecule, g.field,
                                   verbose=True)
        result.append(rrms_val)


def plot(result_list, title):
    fig, ax1 = plt.subplots()
    ax1.plot(ratio_values, result_list)
    ax2 = ax1.twiny()
    ax2.plot(indicator_charge, result_list)
    if start_charges[indicator_label-1] < 0:
        ax2.invert_xaxis()

    ax1.set_xlabel(charge_type.upper() + " charge ratio")
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
    mark_charge = g.molecule[indicator_label-1].charges[charge_type]
    ax2.scatter(mark_charge, charge_rrms, marker='D')
    ax2.annotate(charge_type.upper(), xy=(mark_charge, charge_rrms),
                 textcoords='offset points', xytext=(5, -10))
    ax1.scatter(1, equiv_charge_rrms)

    mark_equiv_charge = g.molecule[indicator_label-1].charges[charge_type +
                                                              '_equiv']
    if (math.isclose(mark_charge, mark_equiv_charge, rel_tol=0.04) and
            math.isclose(charge_rrms, equiv_charge_rrms, rel_tol=0.04)):
        print("WARNING: The NBO and NBO (equiv) points overlap or are close to"
              " overlapping. Only one label is plotted.")
    else:
        ax1.annotate(charge_type.upper() + ' (equiv)',
                     xy=(1, equiv_charge_rrms), textcoords='offset points',
                     xytext=(5, -10))

    ax1.plot((1, 1), (0, ax1.get_ylim()[1]), 'g--')
    ax1.plot(ratio_limits, (resp_charge_rrms, resp_charge_rrms), 'r--')

    plt.title(title, y=1.15)
    plt.show()
    plt.close()

molecule_name = graphs.pretty_molecule_name(molecule_name)

title = "{0}: RRMS on {1} fitting points v. {2} ratio".format(
        molecule_name, esp_charge_type.upper(), charge_type.upper())

if zero_net_charge:
    plot(result,
         title + "\nset on ALL atoms (only possible for neutral molecules)")

plot(heavy_result, title + "\nset on heavy atoms ONLY (H free to improve fit)")
