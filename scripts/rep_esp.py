#!/usr/bin/env python3

from repESP import cube_helpers, charges
from repESP.rep_esp import calc_grid_field

import argparse
import os

help_description = """
    Create a Gaussian cube file (.cub) reproducing the ESP from given charges.

    If you do not want to process the charge in any way, raw charges can be
    extracted from Gaussian .log or AIMAll .sumviz files. However, it is
    probably sensible to average or equivalence the charges first, using the
    appropriate script from this suite. The charges can then be read in from
    the output text file.
    """
# TODO: give the name of the appropriate script once it's been written

parser = argparse.ArgumentParser(
    description=help_description,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("template_cube",
                    help="template cube file (content doesn't matter, only the"
                    " grid and atomic coordinates are used)")

charge_choices = ['list', 'mulliken', 'nbo', 'aim', 'mk', 'chelp', 'chelpg',
                  'hly']
charge_choices_text = "'" + "', '".join(charge_choices[1:-1]) + "' and '" + \
                      charge_choices[-1] + "'."
parser.add_argument(
    "input_charge_type",
    help="type of input charges. To extract charges from a text file, use the "
    "'list' option. Other available options are: " + charge_choices_text,
    choices=charge_choices, metavar="input_charge_type")

parser.add_argument(
    "input_charge_file",
    help=""""input file to extract charges from. This can be a Gaussian .log,
    AIMAll .sumviz or a text file. The appropriate format for the text file is
    used throughout this suite of scripts and is very simple:

    1. Values should be space-separated
    2. Atom order should follow that in the template cube file
    3. Line breaks are not significant
    """)

parser.add_argument("-o", "--output",
                    help="output file name",
                    default='rep_esp.cub')

args = parser.parse_args()

if os.path.exists(args.output):
    raise FileExistsError("Output file exists: " + args.output)

template_cube = cube_helpers.Cube(args.template_cube)
molecule = template_cube.molecule

if args.input_charge_type == 'list':
    input_type = 'qout'
elif args.input_charge_type == 'aim':
    input_type = 'sumviz'
else:
    input_type = 'log'

charges._get_charges(args.input_charge_type, args.input_charge_file,
                     input_type, molecule)

rep_field = calc_grid_field(molecule, template_cube.field.grid, 'rep_esp',
                            [args.input_charge_type])[0]

rep_field.write_cube(args.output, molecule, args.input_charge_type)
