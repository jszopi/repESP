import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_pdf import PdfPages
import os

import field_comparison

DIR_LABELS = ['x', 'y', 'z']


def _plot_common(dimension):
    """Set up plot of correct dimensionality and return related objects"""
    fig = plt.figure()
    if dimension == 3:
        ax = fig.add_subplot(111, projection='3d')
    elif dimension == 2:
        ax = fig.add_subplot(111)
    else:
        raise NotImplementedError("Plotting of dimension {0} not implemented"
                                  .format(dimension))

    # Add a horizontal line at 0 for 2D plots
    if dimension == 2:
        ax.axhline(color='k', linestyle='--')

    return fig, ax


def plot(*fields, color=None, dist_field_filter=None, exclusion_dist=0,
         rand_skim=0.01, save_to=None):

    # Still plenty TODO:
    # * Color bar legend!
    # For easy comparisons between plots:
    # * Set ranges of axes. The axes should be created by the caller, perhaps
    #   using a function or class from this module, so it can be used as a
    #   parameter in calling this function for separate graphs.
    # * How to change color scale

    assert 2 <= len(fields) <= 3

    # Pack fields and color together
    if color is not None:
        fields_and_color = list(fields) + [color]
    else:
        fields_and_color = fields

    field_comparison._check_grids(*fields_and_color)
    field_comparison._check_fields_for_nans(*fields_and_color)

    fig, ax = _plot_common(len(fields))
    _set_axis_labels(ax, *fields)

    if dist_field_filter is not None:
        if dist_field_filter.field_type != 'dist':
            print("WARNING: The field selected for filtering is not of type "
                  "'dist' but ", dist_field_filter.field_type)
        dist_field_filter, *fields_and_color = field_comparison.filter_by_dist(
            exclusion_dist, *([dist_field_filter] + fields_and_color))
    elif exclusion_dist:
        print("WARNING: exclusion distance specified but no Field passed to "
              "filter by.")
    fields_and_color = field_comparison.skim(rand_skim, *fields_and_color)
    fields_and_color = list(map(field_comparison._flatten_no_nans,
                            fields_and_color))

    if color is not None:
        cmap = _get_cmap(len(fields), color.field_type)
        *fields, color = fields_and_color
        # This has to be inside of the conditional because an error occurs when
        # the kwarg ``c`` is set to None, even though it's the default value.
        ax.scatter(*fields, c=color, cmap=cmap)
    else:
        fields = fields_and_color
        ax.scatter(*fields)

    _save_or_display(save_to)


def _set_axis_labels(ax, *fields):
    """Set axis labels based on free-form names of Fields being plotted"""
    for field, dir_label in zip(fields, DIR_LABELS):
        getattr(ax, "set_" + dir_label + "label")(field.lookup_name())


def _get_cmap(dimension, field_type):
    """Return a color map based on plot dimensionality and field type"""
    if field_type == 'dist':
        if dimension != 3:
            # Shading by distance is more intuitive
            return plt.get_cmap('Blues_r')
        else:
            print("WARNING: Shading by distance doesn't look good on a 3D "
                  "plot. Colouring instead.")
    # This colour map is fairly intuitive, looks good in greysale and
    # to coloublind people (according to a simulation on:
    # http://www.color-blindness.com/coblis-color-blindness-simulator/
    return plt.get_cmap('gnuplot2')


def _save_or_display(save_to=None):
    """Save the plot or display it if save_to is None"""
    if save_to is None:
        plt.show()
    else:
        if type(save_to) is PdfPages:
            # Need to check the type first, because it may be a file object if
            # a pdf is to be created, see:
            # http://matplotlib.org/faq/howto_faq.html#save-multiple-plots-to-one-pdf-file
            plt.savefig(save_to, format="pdf")
        elif os.path.isfile(save_to):
            raise FileExistsError("File exists: " + save_to)
        else:
            # DPI may need to be increased
            plt.savefig(save_to)
    plt.close()
