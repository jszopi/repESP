import sys, os
sys.path.append(f"{os.path.dirname(__file__)}")
sys.path.append(f"{os.path.dirname(__file__)}/../repESP_old")

from cube_helpers import Cube
from field_comparison import difference, filter_by_dist
import resp_helpers

import argparse

def main():

    help_description = """Calculate the difference between two fields FIELD1 -
                       FIELD2. The fields can be given either as Gaussian cube
                       files (.cub) or ESP fitting points (.esp)."""

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
                        " the given distance. This only works with cube files.",
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
                        help="output file path and name without extension",
                        default='field_diff')


    args = parser.parse_args()

    filetype = args.field1[-4:]
    filetype2 = args.field2[-4:]
    if filetype not in ['.cub', '.esp']:
        print("The fields must have the extension .cub or .esp")
        sys.exit(1)
    if filetype != filetype2:
        print("The fields given have different extensions")
        sys.exit(1)

    if os.path.exists(args.output + filetype):
        raise FileExistsError("Output file exists: " + args.output + filetype)

    if args.exclude is not None and filetype == '.esp':
        raise ValueError("The --exclude option can only be used with .cub files.")
    if args.exclusion_as_zero and args.exclude is None:
        raise ValueError("The --exclusion_as_zero flag should only be set when "
                         "--exclude options are set.")

    if filetype == '.cub':
        field1 = Cube(args.field1)
        molecule = field1.molecule
        field1 = field1.field
        field2 = Cube(args.field2).field

        # Exclusion: run checks before the (possibly costly) difference is
        # calculated
        if args.exclude:
            ed_cube, isoval, excl_dist = args.exclude
            ed_cube = Cube(ed_cube)
            isoval = float(isoval)
            excl_dist = float(excl_dist)
    elif filetype == '.esp':
        esp_file1 = resp_helpers.G09_esp(args.field1)
        esp_file2 = resp_helpers.G09_esp(args.field2)
        # TODO: Molecules not cross-checked. They should be because .esp files in
        # Antechamber format produced by repESP will have no identitity info. This
        # is fine for now for the difference (actual - reproduced), because actual
        # will be Gaussian .esp but not the other way round.
        molecule = esp_file1.molecule
        field1 = esp_file1.field
        field2 = esp_file2.field

    diff = difference(field1, field2, relative=args.relative,
                      absolute=args.absolute)

    if filetype == '.cub':
        if args.exclude:
            diff.check_nans = False
            dist = ed_cube.field.distance_transform(isoval)
            assign_val = 0 if args.exclusion_as_zero else None
            _dist, diff.values = filter_by_dist(excl_dist, dist, diff,
                                                assign_val=assign_val)

        diff.write_cube(args.output + '.cub', molecule)

    elif filetype == '.esp':
        diff.write_to_file(args.output + '.esp', molecule)
