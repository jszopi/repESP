from repESP import resp_helpers, graphs

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm


esp_charge_type = 'mk'

molecule_name = 'NMe3H_plus'
# molecule_name = 'methane'

# SLICING plots
# Note that 2D slicing plots in planes different than those of the coordinate
# system distort distances (3D plots are fine). TODO
# Methane:
# cut_through = 1, 2, 3
# NMe3_plus:
cut_through = 1, 13, 14
# NMe4_plus:
# cut_through = 1, 5, 17
# NMe3:
# cut_through = 1, 5, 9

path = '../data/' + molecule_name + '/'

common_fn = path + molecule_name + "_" + esp_charge_type
input_esp = common_fn + ".esp"

g = resp_helpers.G09_esp(input_esp)

title = molecule_name + " " + esp_charge_type.upper()

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
