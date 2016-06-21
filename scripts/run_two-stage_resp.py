#!/usr/bin/env python3

from repESP.resp import run_resp
from repESP import charges
from repESP import resp_helpers

import argparse
import resp_parser

import os
import shutil

help_description = """Evaluate the fit quality of the given charges on the
    provided mesh of fitting points."""

parser = argparse.ArgumentParser(
    parents=[resp_parser.parser],
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description=help_description)

parser.add_argument("esp_filename",
                    help=resp_parser.esp_file_help,
                    metavar="ESP_FILE_NAME")

parser.add_argument("--save_resp_to",
                    help="save the files of the RESP calculation to the given "
                    "directory. Note that RESP calculations are not performed "
                    "if all charges are to be scaled.",
                    metavar="SAVE_RESP_TO")

parser.add_argument("-o", "--output",
                    help="output file name",
                    default='resp_charges.txt')

args = parser.parse_args()

if args.save_resp_to is None:
    save_resp_to = "resp_temp_dir-dont_remove/"
else:
    save_resp_to = args.save_resp_to + "/"

if os.path.exists(args.output):
    raise FileExistsError("Output file exists: " + args.output)
if os.path.exists(save_resp_to):
    raise FileExistsError("Output directory exists: " + save_resp_to)

molecule = run_resp(args.respin_location, save_resp_to, resp_type='two_stage',
                    esp_fn=args.esp_filename)

if args.save_resp_to is None:
    shutil.rmtree(save_resp_to)

print("\nThe molecule with RESP charges is:")
for atom in molecule:
    atom.print_with_charge('resp')

charges.dump_charges_to_qout(molecule, "resp", args.output)
