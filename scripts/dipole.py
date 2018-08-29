#!/usr/bin/env python3

from repESP.resp_helpers import G09_esp
from repESP import charges, rep_esp

import argparse
import charges_parser
import resp_parser
import numpy as np

help_description = """Calculate dipole moment in Debye."""

parser = argparse.ArgumentParser(
    parents=[charges_parser.parser, resp_parser.parser],
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description=help_description)

# TODO: The ESP file is only necessary to create the molecule but that's not
# clear from the help message. Also, `respin_location` is only necessary to
# locate that ESP file.
parser.add_argument("esp_file",
                    help=resp_parser.esp_file_help,
                    metavar="ESP_FILE")

args = parser.parse_args()

input_type = charges_parser.input_type(args.input_charge_type)

# Note: the molecule is taken from the esp file
molecule = G09_esp(args.respin_location + '/' + args.esp_file).molecule
# ... but this will throw an error if the molecule in log is not the same one
charges._get_charges(args.input_charge_type, args.input_charge_file,
                     input_type, molecule)

dipole = rep_esp.calc_dipole(molecule, args.input_charge_type)

for coord, component in zip(['X', 'Y', 'Z'], dipole):
    print("{0}:   {1: .6f}".format(coord, component))

print("TOT: {0: .6f}".format(np.linalg.norm(dipole)))
