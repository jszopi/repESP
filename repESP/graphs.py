import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_pdf import PdfPages
import os

# This was necesssary to prevent y-axis label from being cut off when plotting
# http://stackoverflow.com/a/17390833
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

from . import field_comparison

DIR_LABELS = ['x', 'y', 'z']


def _plot_common(dimension, title):
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

    if title is not None:
        plt.title(title)

    return fig, ax


def plot(*fields, color=None, color_span=None, dist_field_filter=None,
         exclusion_dist=0, rand_skim=0.01, extra_filter=None, save_to=None,
         axes_limits=None, title=None, get_limits=None):

    assert 2 <= len(fields) <= 3

    # Pack fields and color together
    if color is not None:
        fields_and_color = list(fields) + [color]
    else:
        fields_and_color = fields

    field_comparison._check_grids(*fields_and_color)
    field_comparison._check_fields_for_nans(*fields_and_color)
    # Necessary, as the original Field will be overwritten when filtering
    dist_field_filter_type = dist_field_filter.field_type

    fig, ax = _plot_common(len(fields), title)
    _set_axis_labels(ax, *fields)

    # This function got really fat due to all that filtering and it can still
    # handle only one additional filter. Some refactoring is due. TODO
    if extra_filter is not None:
        fields_and_color = extra_filter(*fields_and_color)
        # This filtering step changes all Fields to np.arrays. As a result, in
        # the next filtering step, by dist_field_filter, a mixture of np.arrays
        # and Fields is passed, which is not handled by the filters. While that
        # deficiency was intentional, I don't think there's a reason it should
        # not be handled (TODO). But for now, a kludge:
        if dist_field_filter is not None:
            dist_field_filter = dist_field_filter.values
    if dist_field_filter is not None:
        if dist_field_filter_type != 'dist':
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
        cmap_name = color.lookup_name()
        *fields, color = fields_and_color
        # ax.scatter has to be inside of the 'color is not None' conditional
        # because an error occurs when the kwarg ``c`` is explicitly set to
        # None, even though it's the default value.
        if color_span is None:
            image = ax.scatter(*fields, c=color, cmap=cmap)
        else:
            image = ax.scatter(*fields, c=color, cmap=cmap, vmin=color_span[0],
                               vmax=color_span[1])
        cbar = fig.colorbar(image, label=cmap_name)
    else:
        fields = fields_and_color
        ax.scatter(*fields)

    _set_axes_limits(len(fields), ax, axes_limits)
    _save_or_display(save_to)

    # Save limits to get_limits. This is useful when they are to be reused in
    # other plots. Saving the limits to an argument was more intuitive than
    # returning them.
    if get_limits is not None:
        # Code copied from _set_axes_limits (TODO: DRY)
        limits = []
        for dir_label in DIR_LABELS[:len(fields)]:
            # Get current limits
            limits.append(getattr(ax, "get_" + dir_label + "lim")())
        get_limits[:] = limits


def _set_axes_limits(dimension, ax, axes_limits):
    """Set axes limits"""
    if axes_limits is None:
        return
    # Smaller lengths are allowed, will be interpreted as the first few axes.
    # This should be an Exception not assertion though.
    assert len(axes_limits) <= dimension

    for axis_limits, dir_label in zip(axes_limits, DIR_LABELS):
        # Get current limits
        limits = list(getattr(ax, "get_" + dir_label + "lim")())
        for i, axis_limit in enumerate(axis_limits):
            if axis_limit is not None:
                limits[i] = axis_limit
        getattr(ax, "set_" + dir_label + "lim")(limits)

    # Although **not for my purposes at the moment** (I only want to set limits
    # so that different plots can be easily compared, so both axes will be
    # getting set), it would be nice to rescale the axes which were not
    # modified. However, when autoscaling, matplotlib always uses all the data.
    # ax.relim() with ax.autoscale_view() seemed to be relevant but they do not
    # easily operate on datapoints I think.


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
    # to coloublind people (except for those with the very rare tritanomaly),
    # according to a simulation on:
    # http://www.color-blindness.com/coblis-color-blindness-simulator/
    return plt.get_cmap('plasma')
    # 'gnuplot2' was also considered but in greyscale it yields black and white
    # at the ends and IMO contains unnecessarily many hues. 'plasma' is simpler
    # and uses similar colours.


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
