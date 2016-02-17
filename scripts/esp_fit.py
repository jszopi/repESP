from repESP import resp, rep_esp, charges
from repESP.field_comparison import _check_grids, difference

from numpy import mean, sqrt, square, arange

charge_type = 'mk'
# charge_type = 'chelpg'

molecule_name = 'methane'
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
rms = sqrt(mean(square(diff)))
print("\nRMS: {0:.6f}".format(rms))
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
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    color = g.field.values
    cmap = plt.get_cmap('plasma')

    image = ax.scatter(*list(zip(*g.field.points)), c=color, cmap=cmap)
    cbar = fig.colorbar(image, label="ESP value")

    plt.show()
    plt.close()
