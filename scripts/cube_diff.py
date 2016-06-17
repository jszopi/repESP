#!/usr/bin/env python3

from repESP.cube_helpers import Cube
from repESP.field_comparison import difference, filter_by_dist

import argparse
import os

help_description = """Create a Gaussian cube file (.cub) with the difference
                   between two cube files. The difference is calculated as
                   FIELD1 - FIELD2."""

parser = argparse.ArgumentParser(description=help_description)

parser.add_argument("field1", metavar="FIELD1")
parser.add_argument("field2", metavar="FIELD2")

parser.add_argument("--relative",
                    help="At each point divide the difference by the value of "
                    "the first field at this point.",
                    action="store_true")

parser.add_argument("--absolute",
                    help="Return the absolute value of the difference",
                    action="store_true")

parser.add_argument("--exclude",
                    help="Exclude points lying closer from the electron "
                    "density (ED) isosurface of the given isovalue closer than"
                    " the given distance.",
                    nargs=3,
                    metavar=("ED_CUBE", "ED_ISOVALUE", "EXCLUSION_DISTANCE"))

parser.add_argument("--exclusion_as_zero",
                    help="""When excluding points, excluded points are by
                    default replaced with 'NAN'. If this flag is set, zeros are
                    written instead. This has a cosmetic effect on the
                    appearance of excluded edges of isosurface plots. If this
                    flag is set, the excluded areas are filled, making the
                    plots appear continuous.""",
                    action="store_true")

parser.add_argument("-o", "--output",
                    help="output file name",
                    default='difference.cub')

args = parser.parse_args()

if os.path.exists(args.output):
    raise FileExistsError("Output file exists: " + args.output)

if args.exclusion_as_zero and not args.exclude:
    raise ValueError("The --exclusion_as_zero flag should only be set when "
                     "--exclude options are set.")

field1 = Cube(args.field1)
molecule = field1.molecule
field1 = field1.field
field2 = Cube(args.field2).field

# Exclusion: run checks before the (possibly costly) difference is calculated
if args.exclude:
    ed_cube, isoval, excl_dist = args.exclude
    ed_cube = Cube(ed_cube)
    isoval = float(isoval)
    excl_dist = float(excl_dist)

diff = difference(field1, field2, relative=args.relative,
                  absolute=args.absolute)

if args.exclude:
    diff.check_nans = False
    dist = ed_cube.field.distance_transform(isoval)
    assign_val = 0 if args.exclusion_as_zero else None
    _dist, diff.values = filter_by_dist(excl_dist, dist, diff,
                                        assign_val=assign_val)

diff.write_cube(args.output, molecule)
