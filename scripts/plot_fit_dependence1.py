#!/usr/bin/env python3

from plot_fit_dependence_common import get_parser, preprocess_args

from repESP.cube_helpers import InputFormatError

import matplotlib.pyplot as plt

import pandas
import re


def add_specific_cli_args(parser, plot_appearance_group):

    plot_appearance_group.add_argument(
        "--monitored_atoms_display",
        help="""display labels for atoms given in `--monitored_atoms` to be used in
        the legend. If not given, the numeric labels will be used""",
        type=str,
        nargs="*",
        metavar="DISPLAY_LABELS",
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
        "--monitored_atoms_limits",
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

def plot_flexibility_func(axis, df, varied_atom, statistic, mark_fit, ylim_override):

    x_values = df[get_col_header(varied_atom)]
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


def plot_response_func(
        axis,
        df,
        monitored_atoms,
        monitored_atoms_display,
        ylim_override,
        ylabel_see_legend=False
):
    x_values = df[get_col_header(varied_atom)]

    monitored_dict = {
        label: label if display is None else display
        for label, display in zip(monitored_atoms, monitored_atoms_display)
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

    varied_atom = get_label(df.columns.values[0])
    monitored_atoms = list(map(get_label, df.columns.values[3:]))
    _, rms, rrms, *_ = df.columns.values

    if rms != "RMS":
        raise InputFormatError("Expected RMS column but found {}".format(rms))
    if rrms != "RRMS":
        raise InputFormatError("Expected RRMS column but found {}".format(rrms))

    return varied_atom, monitored_atoms


def line_at_zero_x(axis):
    ylim = axis.get_ylim()
    axis.plot((0, 0), ylim, 'k:')
    # Restore ylimit affected by the above
    axis.set_ylim(ylim)


def plot_common(df, varied_atom, varied_atom_display, title, shaded_region):

    fig, ax1 = plt.subplots()

    if title is not None:
        plt.title(title)

    x_values = df[get_col_header(varied_atom)]
    ax1.set_xlabel("Charge on " + (varied_atom if varied_atom_display is None else varied_atom_display))
    ax1.set_xlim(x_values.min(), x_values.max())
    line_at_zero_x(ax1)

    if shaded_region:
        # 10% flexibility limits
        ax1.axvspan(shaded_region[0], shaded_region[1], alpha=0.2, color='grey')

    return ax1

if __name__ == "__main__":

    parser, plot_appearance_group = get_parser(isTwoAtoms=False)
    add_specific_cli_args(parser, plot_appearance_group)

    args = parser.parse_args()
    preprocess_args(args)

    df = pandas.read_csv(args.scan_output)
    varied_atom, _ = interpret_header(df)

    ax1 = plot_common(df, varied_atom, args.varied_atom_display, args.title, args.shade_region)
    is_first_plot = True

    if args.plot_fit:

        plot_flexibility_func(
            ax1,
            df,
            varied_atom,
            args.plot_fit,
            args.mark_fit,
            args.plot_fit_limits
        )

        is_first_plot = False

    if args.monitored_atoms:

        if args.monitored_atoms_display and len(args.monitored_atoms_display) != len(args.monitored_atoms):
            raise KeyError(
                "If specified, the length of the argument list passsed as "
                "`--monitored_atoms_display` must equal that of `--monitored_atoms`"
            )

        plot_response_func(
            ax1 if is_first_plot else ax1.twinx(),
            df,
            args.monitored_atoms,
            args.monitored_atoms_display,
            args.monitored_atoms_limits,
            not is_first_plot
        )

    if args.output is None :
        plt.show()
        plt.close()
    else:
        plt.savefig(args.output + ".pdf", format='pdf')
