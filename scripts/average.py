#!/usr/bin/env python3

from repESP.resp import equivalence, _get_input_files, _read_respin
from repESP import charges

import argparse
import parent_parser
import os

help_description = """Average charges according to atom equivalence information
    from 'respin' files. Note that ESP-based charges should be equivalenced
    instead (use the 'equivalence.py' script)."""

parser = argparse.ArgumentParser(
    parents=[parent_parser.parent_parser],
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description=help_description)

parser.add_argument("--respin_location",
                    help="""The directory with the location of the *two* respin
                    files. Their names are expected to end with '.respin1' and
                    '.respin2'. They can be generated using the 'respgen'
                    program (see 'instructions.md').""",
                    default='.', metavar="RESPIN_LOCATION")

parser.add_argument("--thresh",
                    help="Inform about charges which change upon averaging "
                    "by more than THRESH",
                    default=0.05, type=float, metavar="THRESH")

parser.add_argument("-o", "--output",
                    help="output file name",
                    default='charges.txt')


args = parser.parse_args()

if os.path.exists(args.output):
    raise FileExistsError("Output file exists: " + args.output)

input_type = parent_parser.input_type(args.input_charge_type)

respin1, respin2 = _get_input_files(args.respin_location, respin1_fn="",
                                    respin2_fn="")
# Note: the molecule is taken from one of the respin files...
molecule = _read_respin(respin1)[-1]
# ... but this will throw an error if the molecule in log is not the same one
charges._get_charges(args.input_charge_type, args.input_charge_file,
                     input_type, molecule)

averaged_charges = equivalence(molecule, args.input_charge_type,
                               args.respin_location, respin1_fn="",
                               respin2_fn="")
charges._update_molecule_with_charges(molecule, averaged_charges, "avgd")

message = charges.compare_charges(args.input_charge_type, "avgd", molecule,
                                  thresh=args.thresh)
if message:
    print(message)
else:
    print("No charges vary by more than the threshold of " + str(args.thresh))

charges.dump_charges_to_qout(molecule, "avgd", args.output)
