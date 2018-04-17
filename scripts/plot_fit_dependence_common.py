import argparse


def preprocess_args(args):

    args.plot_fit = args.plot_fit if args.plot_fit != "none" else None

    if not args.plot_fit and not args.monitored_atoms:
        raise KeyError(
            "Requested plotting neither fit quality (`--plot_fit`) nor"
            "monitored charges (`--monitored_atoms`) to be plotted"
        )


def get_parser(isTwoAtoms=False):

    help_description = """
        Plot the dependence of the ESP fit and/or values on monitored charges as a
        function of charges on {}, based on the output of the `fit_dependence`
        script.
        """.format("two atoms" if isTwoAtoms else "one atom")

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
                        help="""plot the selected fit statistic ({}). Allowed
                        values are (rms, rrms or none).""".format(
                            "z-axis" if isTwoAtoms else "y-axis"
                        ),
                        type=str,
                        metavar="STATISTIC",
                        choices=["rms", "rrms", "none"],
                        default="rms")

    if not isTwoAtoms:
        parser.add_argument("--monitored_atoms",
                            help="""labels of atoms which charges were monitored and
                            are to be plotted (y-axis)""",
                            type=int,
                            nargs="*",
                            metavar="LABELS",
                            default=[])
    else:
        parser.add_argument("--monitored_atom",
                            help="""label of atom which charge was monitored and
                            is to be plotted (z-axis)""",
                            type=int,
                            metavar="LABEL")

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

    if not isTwoAtoms:
        plot_appearance_group.add_argument(
            "--varied_atom_display",
            help="""display label of the varied atom to be used in the x-axis label.
            If not given, the numeric label will be used.""",
            type=str,
            metavar="DISPLAY_LABEL",
        )
    else:
        plot_appearance_group.add_argument(
            "--varied_atoms_display",
            help="""display labels of the varied atoms to be used in the x-axis x- and y-axis labels.
            If not given, the numeric labels will be used.""",
            type=str,
            nargs=2,
            metavar=("DISPLAY_LABEL1", "DISPLAY_LABEL2"),
            default=[]
        )

    return parser, plot_appearance_group
