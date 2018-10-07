#!/usr/bin/env python3

from plot_fit_dependence_common import *

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.tri as mtri

import json
import pandas
import textwrap


def interpret_point_marker(string):

    try:
        coord1, coord2, *marker_split_json = string.split()
        coord1, coord2 = float(coord1), float(coord2)
    except ValueError as e:
        print(e)
        raise argparse.ArgumentTypeError(
            "Incorrect point marker specification. (Is the argument enclosed in "
            "single quotes and are the first two arguments floating point "
            "numbers?: {}".format(e)
        )

    if marker_split_json == []:
        marker_spec = {
            "marker": "x",
            "s": 50,
            "lw": 0.8
        }
    else:
        marker_json = " ".join(marker_split_json)
        try:
            marker_spec = json.loads(marker_json)
        except ValueError:
            raise argparse.ArgumentTypeError(
                "Failed parsing the JSON marker specification"
            )

    return coord1, coord2, marker_spec


def add_specific_cli_args(parser, plot_appearance_group):

    parser.add_argument(
        "--plot_relative",
        help="""value of the z-axis variable relative to which other values are
        to be scaled. As a result, the plotted z-axis values will be plotted as
        percentage deviations from the given value. If the ESP fit is plotted,
        a sensible value for this argument is the fit minimum.""",
        type=float,
        metavar="VALUE"
    )

    parser.add_argument(
        "--contours",
        help="""Specify this option if you want a 2D contour plot instead of a
        3D plot. The values given to this option will be used as the contour values.""",
        type=float,
        nargs="*",
        metavar="VALUES",
    )

    plot_appearance_group.add_argument(
        "--swap_axes",
        help="""Swap x and y axes. By default the first column in the input file is the
        x-axis and the second column is the y-axis.""",
        action="store_true"
    )

    plot_appearance_group.add_argument(
        "--draw_ratio_line",
        help="""Draw the line of constant charge ratio. Required are coordinates
        of one point on the line. Another point defining the line will be (0, 0)""",
        type=float,
        nargs=2,
        metavar=("COORD1", "COORD2")
    )

    plot_appearance_group.add_argument(
        "--mark_point",
        help="R|" + textwrap.dedent("""\
        If there are any extra points to be marked on the
        graph, their coordinates and optionally JSON
        specifications of the plot marker should be given here.
        Note that this argument group must be wrapped in single
        quotes (then double quotes can be freely used in the
        JSON specification without escaping). Multiple points
        can be marked by specifying this option multiple times.

        Coordinates should be given in the same order as the
        columns in the input file, i.e. first the charge on the
        atom in the first column. A simple example of a valid
        argument is: '0.5 -0.7'.

        The marker specification is optional and should be a
        JSON-encoded list of keyword arguments passed to the
        `matplotlib.pyplot.scatter` function. If the JSON
        specification is not given, the following default will
        be used:

        '{"marker": x", s": 50, "lw": 0.8}'

        where the "marker" key specifies the pictorial type of
        the marker, "s" specifies the size and "lw" the
        linewidth. For a full list of possible keyword
        arguments please refer to matplotlib documentation at:
        https://matplotlib.org/api/_as_gen/matplotlib.pyplot.scatter.html
        """),
        type=interpret_point_marker,
        action='append',
        metavar="'COORD1 COORD2 JSON_MARKER_SPEC'"
    )


def plot_common(varied_atoms, varied_atoms_display, swap_axes, title):

    if title is not None:
        plt.title(title)

    label = lambda firstOrSecond: (
        varied_atoms_display[firstOrSecond]
        if varied_atoms_display
        else varied_atoms[firstOrSecond]
    )

    labels = (label(0), label(1)) if not swap_axes else (label(1), label(0))

    plt.xlabel("Charge on " + labels[0])
    plt.ylabel("Charge on " + labels[1])


def plot_3d(x_axis, y_axis, z_axis, colorbar_label):

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    tri = mtri.Triangulation(x_axis, y_axis)
    surf = ax.plot_trisurf(x_axis, y_axis, z_axis, triangles=tri.triangles, cmap=cm.plasma)
    fig.colorbar(surf, label=colorbar_label)


def get_data(df, varied_atoms, plot_fit, plot_relative, monitored_atom, swap_axes):

    x_axis = df[get_col_header(varied_atoms[0])]
    y_axis = df[get_col_header(varied_atoms[1])]

    if swap_axes:
        x_axis, y_axis = y_axis, x_axis

    z_axis = df[get_col_header(monitored_atom)] if plot_fit is None else df[plot_fit.upper()]
    z_axis = z_axis if plot_relative is None else [100*(z/plot_relative-1) for z in z_axis]

    z_label = get_col_header(monitored_atom) if plot_fit is None else plot_fit.upper()
    z_label = "% increase in {} from {}".format(z_label, plot_relative) if plot_relative else z_label

    return x_axis, y_axis, z_axis, z_label


def plot_contours(x_axis, y_axis, z_axis, contours):

    fig, ax1 = plt.subplots()
    contours = plt.tricontour(x_axis, y_axis, z_axis, sorted(contours), colors='k', zorder=1)
    plt.clabel(contours, fmt="%1.0f", inline=1, fontsize=15, colors='b')
    ax1.set_aspect('equal')
    return ax1


def plot_ratio_line(ax, coords, swap_axes):
    tangent = coords[1]/coords[0] if not swap_axes else coords[0]/coords[1]
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    plt.plot(xlim, list(tangent*x for x in xlim), 'r--', zorder=1)
    ax.set_ylim(*ylim)


def mark_point(ax, coords, marker_spec, swap_axes):
    coords = coords if not swap_axes else (coords[1], coords[0])
    ax.scatter(coords[0], coords[1], **marker_spec)


if __name__ == "__main__":

    parser, plot_appearance_group = get_parser(isTwoAtoms=True)
    add_specific_cli_args(parser, plot_appearance_group)

    args = parser.parse_args()
    preprocess_args(args, isTwoAtoms=True)

    df = pandas.read_csv(args.scan_output)
    varied_atoms, _ = interpret_header(df, isTwoAtoms=True)

    *data, colorbar = get_data(
        df,
        varied_atoms,
        args.plot_fit,
        args.plot_relative,
        args.monitored_atom,
        args.swap_axes
    )

    if args.contours is None:
        plot_3d(*data, colorbar)
    else:
        ax = plot_contours(*data, args.contours)

        if args.draw_ratio_line:
            plot_ratio_line(ax, args.draw_ratio_line, args.swap_axes)

        if args.mark_point:
            for mark_points in args.mark_point:
                *coords, marker_spec = mark_points
                mark_point(ax, coords, marker_spec, args.swap_axes)

    plot_common(varied_atoms, args.varied_atoms_display, args.swap_axes, args.title)

    if args.output is None :
        plt.show()
        plt.close()
    else:
        plt.savefig(args.output + ".pdf", format='pdf')
