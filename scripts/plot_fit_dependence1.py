#!/usr/bin/env python3

import argparse

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

parser.add_argument("--output",
                    help="file to save the plot to. If not given, an interactive plot is shown.",
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
    "--legend",
    help="""custom legend labels of atoms given in `--charges`. If not given, the
    numeric labels will be used""",
    type=int,
    nargs="*",
    metavar="LABELS",
    default=[]
)

marker_group = parser.add_argument_group(
    title="options regarding marking additional features on the graph""",
    description="""These options were added for creating plots consistent
                with the original publication."""
)

marker_group.add_argument("--mark_fit",
                    help="""the charge and corresponding value of the selected
                    fit statistic (see --plot_fit option) for a point to be
                    marked on the fit dependence graph (intended to be the minimum)""",
                    type=float,
                    nargs=2,
                    metavar=("CHARGE", "VALUE"))

marker_group.add_argument("--shade_region",
                    help="region of charge values to be shaded (intended for flexibility range)",
                    type=float,
                    nargs=2,
                    metavar=("LOWER", "UPPER"))

args = parser.parse_args()
