from repESP.cube_helpers import InputFormatError

import argparse
import re


def preprocess_args(args, isTwoAtoms=False):

    args.plot_fit = args.plot_fit if args.plot_fit != "none" else None

    monitored_plot = args.monitored_atom if isTwoAtoms else args.monitored_atoms

    if not args.plot_fit and not monitored_plot:
        raise ValueError(
            "Requested plotting neither fit quality (`--plot_fit`) nor"
            "monitored charges (`--monitored_atoms`) to be plotted"
        )

    if isTwoAtoms:
        if args.plot_fit is not None and args.monitored_atom is not None:
            raise ValueError('`--monitored_atom` can only be specified if `--plot_fit` is set to "none".')

        if args.monitored_atom is not None and args.plot_relative:
            raise ValueError('`--plot_relative` cannot be specified with `--monitored_atom`.')


# https://stackoverflow.com/a/22157136
class SmartFormatter(argparse.HelpFormatter):
    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)


def get_parser(isTwoAtoms=False):

    help_description = """
        Plot the dependence of the ESP fit and/or values on monitored charges as a
        function of charges on {}, based on the output of the `fit_dependence`
        script.
        """.format("two atoms" if isTwoAtoms else "one atom")

    parser = argparse.ArgumentParser(
        description=help_description,
        formatter_class=SmartFormatter
    )

    parser.add_argument("scan_output",
                        help="output of the fit_dependence script (csv).",
                        metavar="FILENAME")

    parser.add_argument("-o", "--output",
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


def get_label(col):
    regex = re.compile("Charge on ([0-9]*)")
    match = regex.match(col)
    if match is None:
        raise InputFormatError("Expected charge column but found {}".format(col))
    return match.group(1)


def get_col_header(label):
    return "Charge on {}".format(label)


def interpret_header(df, isTwoAtoms=False):

    rms_index = 2 if isTwoAtoms else 1
    rrms_index = rms_index + 1

    rms = df.columns.values[rms_index]
    rrms = df.columns.values[rrms_index]

    if rms != "RMS":
        raise InputFormatError("Expected RMS column but found {}".format(rms))
    if rrms != "RRMS":
        raise InputFormatError("Expected RRMS column but found {}".format(rrms))

    varied_atoms = list(map(get_label, df.columns.values[0:rms_index]))
    monitored_atoms = list(map(get_label, df.columns.values[rrms_index+1:]))

    return varied_atoms, monitored_atoms
