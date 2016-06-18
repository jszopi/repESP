#!/usr/bin/env python3

from repESP import cube_helpers, charges
from repESP.rep_esp import calc_grid_field

import argparse
import charges_parser
import os

help_description = """
    Create a Gaussian cube file (.cub) reproducing the ESP from given charges.

    If you do not want to process the charge in any way, raw charges can be
    extracted from Gaussian .log or AIMAll .sumviz files. However, it is
    probably sensible to average or equivalence the charges first, using the
    'average.py' script from this suite. The charges can then be read in from
    the output text file.
    """

parser = argparse.ArgumentParser(
    parents=[charges_parser.parser],
    description=help_description,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("template_cube",
                    help="template cube file (content doesn't matter, only the"
                    " grid and atomic coordinates are used)",
                    metavar="TEMPLATE_CUBE")

parser.add_argument("-o", "--output",
                    help="output file name",
                    default='rep_esp.cub')

args = parser.parse_args()

if os.path.exists(args.output):
    raise FileExistsError("Output file exists: " + args.output)

template_cube = cube_helpers.Cube(args.template_cube)
molecule = template_cube.molecule

input_type = charges_parser.input_type(args.input_charge_type)
charges._get_charges(args.input_charge_type, args.input_charge_file,
                     input_type, molecule)

rep_field = calc_grid_field(molecule, template_cube.field.grid, 'rep_esp',
                            [args.input_charge_type])[0]

rep_field.write_cube(args.output, molecule, args.input_charge_type)
