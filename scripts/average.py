import sys, os
sys.path.append(f"{os.path.dirname(__file__)}")
sys.path.append(f"{os.path.dirname(__file__)}/../repESP_old")

from resp import equivalence, _get_input_files, _read_respin
from resp import _check_ivary, run_resp
import charges
import resp_helpers

import argparse
import charges_parser
import resp_parser

import shutil
import copy


def main():
    help_description = """Average charges according to atom equivalence information
        from 'respin' files. Note that ESP-based charges should be equivalenced
        instead! This can be achieved by passing the fitting points to the
        --esp_file option. Equivalencing is then performed by running an
        unconstrained RESP calculation."""

    parser = argparse.ArgumentParser(
        parents=[charges_parser.parser, resp_parser.parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=help_description)

    parser.add_argument("--esp_file",
                        help=resp_parser.esp_file_help,
                        metavar="ESP_FILE")

    parser.add_argument("--thresh",
                        help="Inform about charges which change upon averaging "
                        "by more than THRESH",
                        default=0.05, type=float, metavar="THRESH")

    parser.add_argument("--save_resp_to",
                        help="save the files of the RESP calculation to the given "
                        "directory",
                        metavar="SAVE_RESP_TO")

    # TODO: This is a a useful hack to extract charges. This should probably be a
    # separate script but it would still require the information about the molecule
    # from a cube or respin files. That would be undesriable, so I don't want to
    # implement it before I have a better idea.
    parser.add_argument("--dump_raw",
                        help="dump extracted raw charges and exit",
                        metavar="DUMP_RAW")

    parser.add_argument("-o", "--output",
                        help="output file name",
                        default='averaged_charges.txt')

    args = parser.parse_args()

    # Check settings related to equivalencing
    if args.save_resp_to is None:
        save_resp_to = "averaging_temp_dir-dont_remove"
    else:
        save_resp_to = args.save_resp_to
        if args.esp_file is None:
            raise ValueError("The option '--save_resp_to' should only be specified"
                             " when the '--esp_file' option is set.")

    if os.path.exists(args.output):
        raise FileExistsError("Output file exists: " + args.output)
    if os.path.exists(save_resp_to):
        raise FileExistsError("Output directory exists: " + save_resp_to)

    input_type = charges_parser.input_type(args.input_charge_type)

    respin1, respin2 = _get_input_files(args.respin_location, respin1_fn="",
                                        respin2_fn="")
    # Note: the molecule is taken from one of the respin files...
    molecule = _read_respin(respin1)[-1]
    # ... but this will throw an error if the molecule in log is not the same one
    charges._get_charges(args.input_charge_type, args.input_charge_file,
                         input_type, molecule)

    if args.dump_raw:
        charges.dump_charges_to_qout(molecule, args.input_charge_type, args.dump_raw)
        import sys
        sys.exit(0)

    if args.esp_file is not None:
        try:
            averaged_molecule = run_resp(
                args.respin_location, save_resp_to, resp_type='unrest',
                check_ivary=True, esp_fn=args.esp_file)
        except Exception as e:
            shutil.rmtree(save_resp_to)
            raise e
        if args.save_resp_to is None:
            shutil.rmtree(save_resp_to)
    else:
        averaged_charges, ivary_list = equivalence(
            molecule, args.input_charge_type, args.respin_location, respin1_fn="",
            respin2_fn="")
        # Check equivalence information from respin:
        _check_ivary(True, molecule, ivary_list)
        # The new molecule and name 'resp' is for consistency with the
        # equivalencing code above
        averaged_molecule = copy.deepcopy(molecule)
        charges._update_molecule_with_charges(averaged_molecule, averaged_charges,
                                              "resp")

    message = charges.compare_charges(args.input_charge_type, "resp", molecule,
                                      averaged_molecule, thresh=args.thresh)
    thresh_message = " vary by more than the threshold of " + str(args.thresh)

    if message:
        print("\nThe following charges" + thresh_message)
        print(message)
    else:
        print("\nNo charges" + thresh_message)

    charges.dump_charges_to_qout(averaged_molecule, "resp", args.output)
