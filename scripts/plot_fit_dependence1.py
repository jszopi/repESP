#!/usr/bin/env python3

from repESP.cube_helpers import InputFormatError

import matplotlib.pyplot as plt

import argparse
import itertools
import pandas
import re

help_description = """
    Plot the dependence of the ESP fit and/or values on other charges as a
    function of charges on one atom, based on the output of the `fit_dependence`
    script.
    """

parser = argparse.ArgumentParser(
    description=help_description,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("scan_output",
                    help="output of the fit_dependence script (csv).",
                    metavar="FILENAME")

parser.add_argument("--varied_label",
                    help="label of the varied atom to be displayed as x-axis label",
                    type=str
)

parser.add_argument("--output",
                    help="""file to save the plot to (without extension, plot
                    will be saved as pdf). If not given, an interactive plot is shown.""",
                    type=str, metavar="FILENAME")

parser.add_argument("--plot_fit",
                    help="""plot the selected fit statistic (rms, rrms or
                    none). Note that "none" only makes sense if the `--charges`
                    options is specified.""",
                    type=str,
                    metavar="STATISTIC",
                    choices=["rms", "rrms", "none"],
                    default="rms")

charges_group = parser.add_argument_group(
    title="options regarding plotting charges on other atoms",
    description="""Optionally the charges on other atoms can be plotted
                   (alongside or instead the ESP fit)."""
)

charges_group.add_argument(
    "--charges",
    help="""labels of atoms which charges are to be plotted""",
    type=int,
    nargs="*",
    metavar="LABELS",
    default=[]
)

charges_group.add_argument(
    "--legend_labels",
    help="""custom legend labels of atoms given in `--charges`. If not given, the
    numeric labels will be used""",
    type=str,
    nargs="*",
    metavar="LABELS",
    default=[]
)

marker_group = parser.add_argument_group(
    title="options regarding marking additional features on the graph""",
    description="""These options were added for creating plots consistent
                with the original publication."""
)

marker_group.add_argument(
    "--mark_fit",
    help="""the charge and corresponding value of the selected
    fit statistic (see --plot_fit option) for a point to be
    marked on the fit dependence graph (intended to be the minimum)""",
    type=float,
    nargs=2,
    metavar=("CHARGE", "VALUE")
)

marker_group.add_argument(
    "--shade_region",
    help="region of charge values to be shaded (intended for flexibility range)",
    type=float,
    nargs=2,
    metavar=("LOWER", "UPPER")
)

marker_group.add_argument(
    "--title",
    help="graph title",
    type=str,
    metavar="TITLE"
)

args = parser.parse_args()


def plot_flexibility_func(axis, df, varied, statistic, mark_fit, flex_limits):

    x_values = df[get_col_header(varied)]
    y_values = df[statistic.upper()]
    axis.plot(x_values, y_values)

    # Make the y-axis label and tick labels match the line color
    axis.set_ylabel(statistic.upper(), color='b')
    for tl in axis.get_yticklabels():
        tl.set_color('b')

    if mark_fit:
        # Location of minimum
        axis.plot((mark_fit[0], mark_fit[0]), (0, mark_fit[1]), 'k--')
        axis.plot((x_values.min(), mark_fit[0]), (mark_fit[1], mark_fit[1]), 'k--')

    if flex_limits:
        # 10% flexibility limits
        axis.axvspan(flex_limits[0], flex_limits[1], alpha=0.2, color='grey')

    axis.set_ylim([0, y_values.max()])


def plot_response_func(axis, df, monitored, legend_labels, ylabel_see_legend=False):
    x_values = df[get_col_header(varied)]

    monitored_dict = {
        label: label if legend_label is None else legend_label
        for label, legend_label in itertools.zip_longest(monitored, legend_labels)
    }

    ylabel = "Charges on other atoms"
    if ylabel_see_legend:
        ylabel += " (see legend)"
    axis.set_ylabel(ylabel)

    y_min, y_max = float("+inf"), float("-inf")

    for label, color in zip(monitored_dict, iter(['g', 'r', 'c', 'm', 'y'])):

        try:
            y_values = df[get_col_header(label)]
        except KeyError:
            raise KeyError(
                "Requested charge on atom {} which is not in the scan output file".format(
                    label
                )
            )

        axis.plot(
            x_values,
            y_values,
            color=color,
            label=monitored_dict[label]
        )

        if y_values.min() < y_min:
            y_min = y_values.min()
        if y_values.max() > y_max:
            y_max = y_values.max()

    axis.legend()
    axis.set_ylim([y_min, y_max])

    # Guiding line at zero y2
    axis.plot(axis.get_xlim(), (0, 0), 'k:')


def get_label(col):
    regex = re.compile("Charge on ([0-9]*)")
    match = regex.match(col)
    if match is None:
        raise InputFormatError("Expected charge column but found {}".format(col))
    return match.group(1)


def get_col_header(label):
    return "Charge on {}".format(label)


def interpret_header(df):

    varied = get_label(df.columns.values[0])
    monitored = list(map(get_label, df.columns.values[3:]))
    _, rms, rrms, *_ = df.columns.values

    if rms != "RMS":
        raise InputFormatError("Expected RMS column but found {}".format(rms))
    if rrms != "RRMS":
        raise InputFormatError("Expected RRMS column but found {}".format(rrms))

    return varied, monitored


def line_at_zero_x(axis):
    ylim = axis.get_ylim()
    axis.plot((0, 0), ylim, 'k:')
    # Restore ylimit affected by the above
    axis.set_ylim(ylim)


def plot_common(df, varied, varied_label, title):

    fig, ax1 = plt.subplots()

    if title is not None:
        plt.title(title)

    x_values = df[get_col_header(varied)]
    ax1.set_xlabel("Charge on " + (varied if varied_label is None else varied_label))
    ax1.set_xlim(x_values.min(), x_values.max())
    line_at_zero_x(ax1)

    return ax1

if __name__ == "__main__":

    args.plot_fit = args.plot_fit if args.plot_fit != "none" else None

    if not args.plot_fit and not args.charges:
        raise KeyError(
            "Requested plotting neither fit quality (--plot_fit) nor other charges "
            "(--charges) to be plotted"
        )

    df = pandas.read_csv(args.scan_output)
    varied, _ = interpret_header(df)

    ax1 = plot_common(df, varied, args.varied_label, args.title)

    is_first_plot = True
    if args.plot_fit :

        plot_flexibility_func(ax1, df, varied, args.plot_fit, args.mark_fit, args.shade_region)
        is_first_plot = False

    if args.charges:

        if args.legend_labels and len(args.legend_labels) != len(args.charges):
            raise KeyError(
                "If specified, the length of the argument list passsed as "
                "`--legend_labels` must equal that of `--charges`"
            )

        plot_response_func(
            ax1 if is_first_plot else ax1.twinx(),
            df,
            args.charges,
            args.legend_labels,
            not is_first_plot
        )

    if args.output is None :
        plt.show()
        plt.close()
    else:
        plt.savefig(args.output + ".pdf", format='pdf')
