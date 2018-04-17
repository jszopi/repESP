import argparse

def get_parser(isTwoAtoms=False):

    help_description = """
        Plot the dependence of the ESP fit and/or values on monitored charges as a
        function of charges on one atom, based on the output of the `fit_dependence`
        script.
        """

    parser = argparse.ArgumentParser(
        description=help_description,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

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

    parser.add_argument("--monitored_atoms",
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
        "--varied_atom_display",
        help="""display label of the varied atom to be used in the x-axis label.
        If not given, the numeric label will be used.""",
        type=str,
        metavar="DISPLAY_LABEL",
    )

    return parser, plot_appearance_group
