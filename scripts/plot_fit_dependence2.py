#!/usr/bin/env python3

from plot_fit_dependence_common import *

import pandas


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
        "--contour",
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


if __name__ == "__main__":

    parser, plot_appearance_group = get_parser(isTwoAtoms=True)
    add_specific_cli_args(parser, plot_appearance_group)

    args = parser.parse_args()
    preprocess_args(args)

    df = pandas.read_csv(args.scan_output)
    varied_atoms, _ = interpret_header(df, isTwoAtoms=True)
