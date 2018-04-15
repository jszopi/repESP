#!/usr/bin/env python3

from repESP.cube_helpers import InputFormatError

import matplotlib.pyplot as plt

import argparse
import itertools
import pandas
import re

help_description = """
    Plot the dependence of the ESP fit and/or values on monitored charges as a
    function of charges on one atom, based on the output of the `fit_dependence`
    script.
    """

parser = argparse.ArgumentParser(
    description=help_description,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("scan_output",
                    help="output of the fit_dependence script (csv).",
                    metavar="FILENAME")

parser.add_argument("--output",
                    help="""file to save the plot to (without extension, plot
                    will be saved as pdf). If not given, an interactive plot is shown.""",
                    type=str, metavar="FILENAME")

parser.add_argument("--plot_fit",
                    help="""plot the selected fit statistic (y-axis). Allowed
                    values are (rms, rrms or none).""",
                    type=str,
                    metavar="STATISTIC",
                    choices=["rms", "rrms", "none"],
                    default="rms")

parser.add_argument("--monitored_charges",
                    help="""labels of atoms which charges were monitored and
                    are to be plotted (y-axis)""",
                    type=int,
                    nargs="*",
                    metavar="LABELS",
                    default=[])

plot_appearance_group = parser.add_argument_group(
    title="options regarding the appearance of the graph""",
    description="""These options were added to enable creating plots consistent
                with the original publication."""
)

plot_appearance_group.add_argument(
    "--title",
    help="graph title",
    type=str,
    metavar="TITLE"
)

plot_appearance_group.add_argument(
    "--varied_label",
    help="""label of the varied atom to be displayed as x-axis label. If not
    given, the numeric label will be used.""",
    type=str
)

plot_appearance_group.add_argument(
    "--legend_labels",
    help="""custom legend labels of atoms given in `--monitored_charges`.
    If not given, the numeric labels will be used""",
    type=str,
    nargs="*",
    metavar="LABELS",
    default=[]
)

plot_appearance_group.add_argument(
    "--plot_fit_limits",
    help="""override for graph limits for the fit statistic y-axis. If not
    specified the limits are from 0 to maximum value of the fit statistic""",
    type=float,
    nargs=2,
    metavar=("LOWER", "UPPER"),
)

plot_appearance_group.add_argument(
    "--monitored_charges_limits",
    help="""override for graph limits for the monitored charges y-axis. If not
    specified the values are from minimum to maximum value of charges on the
    monitored atoms""",
    type=float,
    nargs=2,
    metavar=("LOWER", "UPPER"),
)

plot_appearance_group.add_argument(
    "--mark_fit",
    help="""the charge and corresponding value of the selected
    fit statistic (see `--plot_fit` option) for a point to be
    marked on the fit dependence graph (intended to be the fit minimum)""",
    type=float,
    nargs=2,
    metavar=("CHARGE", "VALUE")
)

plot_appearance_group.add_argument(
    "--shade_region",
    help="region of charge values to be shaded (intended for flexibility range)",
    type=float,
    nargs=2,
    metavar=("LOWER", "UPPER")
)

args = parser.parse_args()


def plot_flexibility_func(axis, df, varied, statistic, mark_fit, ylim_override):

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

    if ylim_override is None:
        axis.set_ylim([0, y_values.max()])
    else:
        axis.set_ylim(ylim_override)


def plot_response_func(axis, df, monitored, legend_labels, ylim_override, ylabel_see_legend=False):
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

    if ylim_override is None:
        axis.set_ylim([y_min, y_max])
    else:
        axis.set_ylim(ylim_override)

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


def plot_common(df, varied, varied_label, title, shaded_region):

    fig, ax1 = plt.subplots()

    if title is not None:
        plt.title(title)

    x_values = df[get_col_header(varied)]
    ax1.set_xlabel("Charge on " + (varied if varied_label is None else varied_label))
    ax1.set_xlim(x_values.min(), x_values.max())
    line_at_zero_x(ax1)

    if shaded_region:
        # 10% flexibility limits
        ax1.axvspan(shaded_region[0], shaded_region[1], alpha=0.2, color='grey')

    return ax1

if __name__ == "__main__":

    args.plot_fit = args.plot_fit if args.plot_fit != "none" else None

    if not args.plot_fit and not args.monitored_charges:
        raise KeyError(
            "Requested plotting neither fit quality (`--plot_fit`) nor"
            "monitored charges (`--monitored_charges`) to be plotted"
        )

    df = pandas.read_csv(args.scan_output)
    varied, _ = interpret_header(df)

    ax1 = plot_common(df, varied, args.varied_label, args.title, args.shade_region)
    is_first_plot = True

    if args.plot_fit:

        plot_flexibility_func(
            ax1,
            df,
            varied,
            args.plot_fit,
            args.mark_fit,
            args.plot_fit_limits
        )

        is_first_plot = False

    if args.monitored_charges:

        if args.legend_labels and len(args.legend_labels) != len(args.monitored_charges):
            raise KeyError(
                "If specified, the length of the argument list passsed as "
                "`--legend_labels` must equal that of `--monitored_charges`"
            )

        plot_response_func(
            ax1 if is_first_plot else ax1.twinx(),
            df,
            args.monitored_charges,
            args.legend_labels,
            args.monitored_charges_limits,
            not is_first_plot
        )

    if args.output is None :
        plt.show()
        plt.close()
    else:
        plt.savefig(args.output + ".pdf", format='pdf')
