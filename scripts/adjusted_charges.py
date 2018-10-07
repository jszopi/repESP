import sys, os
sys.path.append(f"{os.path.dirname(__file__)}")
sys.path.append(f"{os.path.dirname(__file__)}/../repESP_old")

from resp import eval_ratios, minimize_ratio
import charges
import resp_helpers

import argparse
import charges_parser
import resp_parser

import shutil


def main():
    help_description = """Calculate adjusted rational charges for the given ESP
        fitting points. Note that the input charges should be averaged beforehand
        and charges passed as a 'list' type. Nevertheless, this script also works
        if raw charges are requested to be extracted from Gaussian or AIMAll output
        files."""

    parser = argparse.ArgumentParser(
        parents=[charges_parser.parser, resp_parser.parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=help_description)

    parser.add_argument("esp_filename",
                        help=resp_parser.esp_file_help,
                        metavar="ESP_FILE_NAME")

    parser.add_argument("--limits",
                        help="The interval for ratios to be considered",
                        nargs=2, type=float, default=(0, 1.5),
                        metavar=("LOWER", "UPPER"))

    parser.add_argument("--sampling",
                        help="The number of points sampled in the interval",
                        default=16, type=int, metavar="SAMPLING")

    # This will be used when plotting is available
    # parser.add_argument("--indicator",
    #                     help="The label of the atom whose charge is to be shown "
    #                     "on the plot.",
    #                     default=1, type=int, metavar="INDICATOR")

    parser.add_argument("--scale_all",
                        help="""Scale charges on all atoms. Note that when used
                        with non-neutral molecules this leads to net charge
                        scaling, which massively deteriorates fit quality. For
                        neutral molecules this approach has not been tested
                        thoroughly and presumably leads to a significantly larger
                        fit deterioration than the default adjusted rational charges.
                        """,
                        action="store_true")

    parser.add_argument("--save_resp_to",
                        help="save the files of the RESP calculation to the given "
                        "directory. Note that RESP calculations are not performed "
                        "if all charges are to be scaled.",
                        metavar="SAVE_RESP_TO")

    parser.add_argument("-o", "--output",
                        help="output file name",
                        default='adjusted_rational_charges.txt')

    args = parser.parse_args()

    if args.save_resp_to is None:
        save_resp_to = "adjusted_charges-temp_dir-dont_remove/"
    else:
        save_resp_to = args.save_resp_to + "/"
        if args.scale_all:
            raise ValueError("The option '--save_resp_to' should only be specified"
                             " when the '--scale_all' option is not set.")

    if os.path.exists(args.output):
        raise FileExistsError("Output file exists: " + args.output)
    if not args.scale_all:
        min_resp_calc_dir = save_resp_to + "/min_resp_calcs/"
        os.mkdir(save_resp_to)
        os.mkdir(min_resp_calc_dir)

    esp_file = resp_helpers.G09_esp(args.respin_location + '/' + args.esp_filename)
    molecule = esp_file.molecule

    # This will be used when plotting is available:
    # if args.indicator > len(molecule):
    #     raise ValueError("Indicator label {0} requested but the molecule has "
    #                      "only {1} atoms.".format(args.indicator, len(molecule)))
    # For now a patch:
    args.indicator = 1

    input_type = charges_parser.input_type(args.input_charge_type)
    charges._get_charges(args.input_charge_type, args.input_charge_file,
                         input_type, molecule)

    start_charges = [atom.charges[args.input_charge_type] for atom in molecule]
    print("\nEvaluating ratios...")

    if args.scale_all:
        regular_args = (molecule, esp_file.field)
        result, indicator_charge, ratio_values = eval_ratios(
            'regular', args.limits, start_charges, args.sampling, args.indicator,
            regular_args, first_verbose=True)
    else:
        func_args = (esp_file.field, args.respin_location, save_resp_to,
                     args.esp_filename, False)
        try:
            result, indicator_charge, ratio_values = eval_ratios(
                'heavy', args.limits, start_charges, args.sampling, args.indicator,
                func_args, first_verbose=True)
        except Exception as e:
            shutil.rmtree(save_resp_to)
            raise e

    print("\nOptimizing fit quality...")
    if args.scale_all:
        regular_args = (start_charges, molecule, esp_file.field)
        reg_min_ratio, reg_min_ratio_rrms, compr_charges = minimize_ratio(
            'regular', ratio_values, result, regular_args)
    else:
        func_args = (start_charges, esp_file.field, args.respin_location,
                     min_resp_calc_dir, args.esp_filename, True)
        try:
            min_ratio, min_ratio_rrms, compr_charges = minimize_ratio(
                'heavy', ratio_values, result, func_args)
        except Exception as e:
            shutil.rmtree(save_resp_to)
            raise e
        shutil.rmtree(min_resp_calc_dir)

    if not args.scale_all and args.save_resp_to is None:
        shutil.rmtree(save_resp_to)

    # TODO: plotting requires information about other charges. This would make this
    # script much bigger or require input from other scripts.

    charges._update_molecule_with_charges(molecule, compr_charges, "compr")
    for atom in molecule:
        atom.print_with_charge("compr")
    charges.dump_charges_to_qout(molecule, "compr", args.output)
